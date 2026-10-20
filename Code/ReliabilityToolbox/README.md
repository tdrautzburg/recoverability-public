# ReliabilityToolbox

> **Note for this public repository.** This file is the toolbox's technical
> documentation and development log. Several sections refer to models, folders and
> scripts of the research project that are not part of this repository (HLW,
> Bauer–Rudebusch, `Miscellaneous`, importance-sampling drivers, tests on their
> workspaces). The files present here are those needed by the UC3 and
> Smets–Wouters workflows; see the top-level `README.md` for how to run them.

The single home of the reliability-diagnostics code. It started (2026-07) as the
active pipeline scripts of `../Miscellaneous/`, **renamed to meaningful names**
(table below); in 2026-08/09 two copies evolved in parallel — this folder (the
posterior-predictive / importance-sampling layer) and `HLW_new/Miscellaneous`
(true-states support, eq.-(3.17) smoother R², tail index, Riccati flag) — and on
**2026-09-15 the two were merged back here**. `HLW_new/Miscellaneous` is retired;
`HLW_new/*/main_run_*.m` and `model_assessment.m` now `addpath` this folder.

## Merge of 2026-09-15 (Miscellaneous ⇄ ReliabilityToolbox)

Common ancestor: the `HLW_new/Miscellaneous/attic/August 25 2026` snapshot.

| taken from `HLW_new/Miscellaneous` | taken from `ReliabilityToolbox` |
|---|---|
| both Kalman filters: relative Riccati stopping rule, `MaxIter` 1e5, non-finite guard, `ric_flag` output | engine switches `M_.fast_mode` / `M_.loglik_only` |
| engine: optional 6th input `true_states`; `R2_h` / `smoother_lag_m` / `smoother_base` via `fixed_lag_smoother_map`; `nu` / `nu_order` via `tail_unreliability_index`; `KF_V`, `KF_ric_flag` on `ModelMatrices`; one (TV) data pass with a common burn-in mask | engine: `obs_ss` computed in step 1 (available to `loglik_only` callers), `X0_pruned` / `B0_pruned` on `ModelMatrices`, TV data pass via `apply_tv_filter` |
| `reliability_draw_wrapper`: constant 19-cell output | posterior-predictive layer (see below), `augment_ssm.m`, tests |
| helpers `clip_series`, `fixed_lag_smoother_map`, `load_rrstar`, `report_corr`, `shade_band`, `tail_unreliability_index` | |

Both numerics are preserved: with the switches off the merged engine is the
Miscellaneous engine; with `fast_mode` / `loglik_only` on it is the
ReliabilityToolbox engine on every field both produced. Skipped objects are
returned as NaN (or `[]`), never 0.

**Interface after the merge.**

- `compute_state_reliability_diagnostics` returns **7 outputs**
  (`LatentStateDiagnostics, ModelMatrices, obs_indices, key_state_indices,
  TV_state_estimate, TV_surprises, TV_log_llhood`). The steady-state-gain
  data pass and its outputs `state_estimate`, `surprises`, `log_llhood`,
  `.steady_CorrY_H`, `.steady_RMSE` are gone everywhere; `uc3_level.mod`,
  `HLW_2023_US_TD_ME_small.mod`, `fast_mode_test.m`, `test_R2_Y*.m`,
  `script_model_test_plot.m`, `evaluate_at_reweighted_means.m` and
  `visualize_reweighting_results.m` were updated accordingly.
- `reliability_draw_wrapper` returns 19 cells in every branch (see its header).
- `compute_state_reliability_vs_truth.m` is superseded by the engine's
  `true_states` input and must be deleted from both folders (the file could
  not be removed by the merge tooling).

**Posterior-predictive / importance-sampling layer (only here):**
`simulate_artificial_data_wrapper`, `artificial_data_wrapper`,
`posterior_predictive_wrapper`, `apply_tv_filter`, `build_qscale`,
`reweight_draws`, `evaluate_at_reweighted_means`, `plot_reweighting`,
`visualize_reweighting_results`, `manual_posterior_functions.m`, `run.sh`
(unfinished SLURM skeleton). `posterior_predictive_wrapper` runs the engine
with `M_.loglik_only = true` and applies `apply_tv_filter` to each data set
against the one `Sequence`.

**Dynare dependency.** `durbin_koopman_diffuse_filter` calls Dynare's
`kalman_transition_matrix` and `compute_Pinf_Pstar` (not in this repo;
`diffuse_filter_initialization` is an in-repo reimplementation that the
production filter deliberately does not use). Running the engine therefore
needs Dynare's `matlab/` folder on the path even for population-only calls.

**Verification (2026-09-15, Octave 11.3 + Dynare from
`/opt/homebrew/opt/dynare/lib/dynare/matlab`, on the `HLW_2023_US_no_lags`
posterior-mode workspace `HLW_new/hlw_level/tmp.mat`, states of interest
`ygap`, `rrstar`, `drrstar`; scripts in `../../merge_check/`).**

1. Merged engine (defaults) vs retired Miscellaneous engine: every field of
   `LatentStateDiagnostics` (22), every common field of `ModelMatrices` (19)
   and all five remaining outputs **bit-identical**. Merged adds only
   `X0_pruned`, `B0_pruned`.
2. Merged `fast_mode` vs merged default: `CorrH_Observables`, `TV_CorrY_H`,
   `TV_RMSE`, `TV_KF_Sequence`, all model matrices and the filtered states
   identical; every skipped object is NaN / empty as documented.
3. Merged `fast_mode` vs pre-merge ReliabilityToolbox `fast_mode`: identical
   **except `TV_RMSE`** (max abs diff 40.8 pp). Expected and intended: the
   pre-merge RT engine averaged the squared forecast errors over the whole
   sample, burn-in included; the merged engine applies the Miscellaneous
   2026-08 rule — one `IDX_ok` mask (burn-in + heteroskedastic periods) for
   RMSEs and correlations alike (see "Burn-in" below). Any RMSE produced by
   the RT-line scripts before the merge is not comparable to the new numbers.
4. Merged `loglik_only` vs pre-merge RT `loglik_only`: `TV_KF_Sequence` and
   all model matrices identical; data pass correctly skipped.
5. Riccati flag on this workspace: converged in 1387 iterations, relative
   step 9.9e-11, rho(A-KC) = 0.9931.

## Lag augmentation in the engine and the smoother-revision regression (2026-09-15)

The states to diagnose are given as a J x 3 list, `M_.CellList_States`, and
`M_.states_of_interest` is derived from it so that every script indexing
results by name keeps working:

```matlab
M_.CellList_States    = { 'ygap',    [0 1 4 20], [];
                          'rrstar',  [0 1 4 20], [4 20];
                          'drrstar', [0 1],      [] };
M_.states_of_interest = expand_states_of_interest(M_.CellList_States);
%  -> {'ygap';'l1ygap';'l4ygap';'l20ygap';'rrstar';'l1rrstar';'l4rrstar';
%      'l20rrstar';'d4rrstar';'d20rrstar';'drrstar';'l1drrstar'}
```

column 2 = fixed-lag smoother lags m (0 = filter), column 3 = change horizons
h (diagnosed at m = 0). `reliability_draw_wrapper` passes `M_.CellList_States`
to the engine when the field exists (else `M_.states_of_interest`, the legacy
name list, which still works unchanged); the engine's second argument accepts
either form and returns the expanded names as
`LatentStateDiagnostics.state_labels`. The lag/change states are built by `augment_ssm`
inside the engine (step 2b) with the labels the `.mod` files used to
declare — `l<m><NAME>` (`l20rrstar`), `d<h><NAME>` (`d20rrstar`) — so those
declarations can be dropped from the `.mod`. The eq.-(3.17) smoother map is
then read off the spec (no structural detection). Diagnosed states are
listed in `LatentStateDiagnostics.state_labels`, one per row of `CorrH` /
`R2_h` / `nu`; the pruned state names are in `ModelMatrices.state_names`.
The legacy name-list input keeps its behaviour (same states, same pruning,
same outputs) up to the two numerical changes documented under "Exact match
with Dynare" below.

**Smoother-revision (proxy) regression** `LatentStateDiagnostics.fr` (also
wrapper cell 20): for every base state, h = 1..40 and m in
`M_.fr_smoother_lags` (default `[4 20]`, `[]` switches it off). Paper
notation (fixed 2026-09-16): its population R² is **P^{2,m}_{i,h}**
(`fr.R2_theory`), the OLS value **P̂^{2,m}_{i,h}** (`fr.R2_emp`) with
coefficients α̂, β̂ (`fr.alpha_emp`, `fr.beta_emp`); P is non-increasing in m,
equals 1 at m = 0 and bounds the infeasible R²_{i,h} from above. It is a
different object from the **hindsight reliability R^{2,m}_{i,h}** of
eq. (3.17) (`R2_h` on the `l<m>` rows, `CorrH` its square root), which is
non-decreasing in m and equals R²_{i,h} at m = 0.

```
f_{i,h,m} = alpha + beta r_{i,h},   r = xhat_{i,t+h|t+h} - xhat_{i,t+h|t},
                                    f = xhat_{i,t+h|t+h+m} - xhat_{i,t+h|t},
                                    xhat_{t+h|t} = A^h xhat_{t|t}
```

- population (steady-state filter, `fr_regression_theory`): `.R2_theory`,
  `.beta_theory` (+ `.cov_fr_theory`, `.var_r_theory`, `.var_f_theory`);
  alpha is identically 0. Derivation in the file header: both revisions are
  finite sums of innovations, so the moments are well defined for unit-root
  states;
- empirical (diffuse/TV filtered data, `fr_regression_empirical`): `.R2_emp`,
  `.alpha_emp`, `.beta_emp`, `.N_emp`, on the common burn-in /
  heteroskedastic mask applied to the whole window t..t+h+m.
- A state whose filter gain is zero (HLW's `drrstar`, a pure combination of
  date-t shocks) has `var_r = 0` and returns NaN, deliberately.

Verified (`merge_check/run_spec.m`, HLW no-lags workspace, 96-state augmented
system): base-state population objects unchanged by the augmentation to
Riccati tolerance, base filtered states unchanged at every t to 6e-13,
log-likelihood unchanged, `d4rrstar = rrstar − l4rrstar` exact in the
filtered states (5e-14) and in `KF_Omega`, `R2_h` weakly increasing in m, a
60,000-period Monte Carlo reproduces `R2_theory` / `beta_theory` / alpha = 0
within 0.008 / 2.6% / 0.02 sd, fast_mode and loglik_only behave as
documented. `test_augment_ssm.m` passes (it found and fixed one Octave shape
bug in `augment_ssm`).

### Exact match with Dynare (2026-09-15)

Reference: `oo_.UpdatedVariables` (x̂_{t|t}) from `calib_smoother` on the
HLW no-lags workspace (`kalman_algo = 4`, `lik_init = 3`,
`heteroskedastic_filter`); script `merge_check/run_dynare_match.m`.

| engine | max |x̂ − Dynare| over t, ygap / rrstar |
|---|---|
| before these changes (legacy names) | 2.0e-4 / 2.9e-4 |
| now, legacy names | 2.8e-13 / 5.0e-13 |
| now, 3-column spec (96-state augmented system) | 5.5e-13 / 7.9e-13 |

Two changes produced this; both alter legacy numbers slightly, toward Dynare:

1. **`A` is no longer rounded.** Step 2 used to do `A = round(A*1e6)/1e6`
   whenever a non-state variable of interest was appended (HLW: `drrstar`),
   to keep the exact-zero redundancy tests of step 3 from being fooled by
   round-off. That rounded the *dynamics* and was the entire 2.9e-4 gap. The
   rounding now applies to a pattern copy used by the pruning test only
   (same pruning decisions). Legacy population objects move by ~1e-6
   relative (`CorrH`, `R2_h`, `KF_Omega`, `K_gain`), filtered states by up
   to 2.9e-4, the log-likelihood in the 5th decimal.
2. **Augmented diffuse prior tied to Dynare's.** `durbin_koopman_diffuse_filter`
   step 4 now (a) copies the *full* rows of appended states into the stacked
   system (the old code copied their Dynare-state columns only, exact for
   non-state variables but dropping the `l_k → l_{k−1}` couplings of a lag
   chain) and (b) sets the diffuse part to `[I; E]·Pinf·[I; E]'`, with
   `Pinf` Dynare's own (un-augmented) matrix and `E = (I − A_aa)\A_ab` the
   exact map from the base states' diffuse component to the appended
   states' (a lag: `E = e_i'`; a change: `E v = 0` on unit-root directions;
   a non-state variable: its `ghx` row). Previously Dynare's
   `compute_Pinf_Pstar` was re-run on the extended system, whose Pinf is the
   orthogonal projector re-normalised by the appended components: with a
   lag chain that moved the t = 1 filtered states by ~1e2 and added 2.05 to
   the log-likelihood (diffuse −½·log F∞ constants). `Pstar` still comes
   from `compute_Pinf_Pstar` on the extended system; the exact diffuse
   filter depends on it only outside range(Pinf). For the legacy HLW input
   (`drrstar`, no unit-root component, `E = 0`) this change is numerically
   neutral. `Sequence.Pinf_ext_dynare` keeps the discarded matrix for
   auditing.

**Plotting scripts** (`HLW_plots*.m`, `Figure17_HLW_SPF.m`) index results by
`strcmp(M_.states_of_interest, name)`; with `M_.states_of_interest =
expand_states_of_interest(M_.CellList_States)` set in the `.mod` they need no
change, and the lag/change labels are the ones the old `.mod` declarations
produced (`l20rrstar`, `d20rrstar`).

**Known degenerate cells:** `TV_CorrY_H(strobs, h ≳ 12)` is `corr(error,
round-off)`: the `strobs` series is zero-filled for t ≤ 240 and its
informative stretch (t = 241..260) is inside the excluded COVID window for
h ≥ 12, so the random-walk forecast revision is identically 0 on the
admissible sample. The engine currently returns the round-off correlation
rather than NaN.

## Prior / posterior draws with `parfor` (`run_reliability_draws.m`, 2026-09-15)

`[results, draws, info] = run_reliability_draws('prior' | 'posterior', n, M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info[, opts])`
replaces Dynare's `prior_function` / `posterior_function`, whose
`execute_prior_posterior_function` evaluates the wrapper in a serial loop
(and needs the MATLAB-written MH history file). It draws the parameters
itself — priors with Dynare's `dprior`, posteriors straight from the
`*_mh*_blck*.mat` chain files with `mh_drop` per chain — and runs the wrapper
in a `parfor` over draws (each draw is independent). `results` has the layout
of `oo_.<type>_function_results`; assign it there and the plotting scripts
are unchanged. `HLW_2023_US_no_lags.mod` uses it.

**Engine speed (2026-09-15).** Two hot spots were vectorised, original code
left in place as comments: the error autocorrelation now accumulates `AKC^hh`
as an iterated product and reads all states off one diagonal per horizon, and
`theoretical_R2_wedge` accepts a vector `state_index` so the engine makes one
call instead of `n_states`. On the 77-state HLW system a full engine call went
from 3.7 s to 2.1 s; `Error_ACF` and `fh` change by ≤ 1.2e-15, everything
else is bit-identical. Where `parfor` pays is one level up — over draws
(`run_reliability_draws`) and over datasets (`model_assessment`'s `d`-loop,
a `parfor` since 2026-09-15) — not inside the engine, whose remaining cost is
chains of small dense products (the 250-horizon doubled-system
recursion, `eig`). The engine's `parfor` over shocks (step 5b) was replaced
by a doubling solve on 2026-09-16 (see the next section), which also removed
the Riccati iteration from the per-call cost (2.1 s → 0.8 s).

## Direct Riccati solve (`solve_steady_state_riccati.m`, 2026-09-16)

Every steady-state Riccati fixed point in the toolbox —
`kalman_filter_steady_and_diffuse` (the engine's filter), `steady_state_kalman_filter`
(used by `R2_wedge_factor` for `.fh_direct`) and `theoretical_R2_wedge`'s
`steady_state_lagged_measurement_filter` — now goes through one function,

```matlab
[Omega, K, V, Omega_tp1_t, ric_flag] = solve_steady_state_riccati(A, B, C, D[, opts]);
```

**Why.** The fixed-point iteration converges at rate `rho(A-KC)^2`, i.e. it
needs about `19/(1-rho^2)` steps for a 1e-10 relative step (1 400 at the HLW
posterior mean, 2 900 at the mode, 3e5 — beyond the 1e5 cap — for a pile-up
draw with `sigma_g`, `sigma_z` ≈ 0 and `rho ≈ 0.99997`). Three such solves
plus seven capped Lyapunov recursions per draw were about half of the per-draw
cost for prior draws (median `rho` 0.98) and three quarters for posterior draws
(median 0.996), see `ReadMe_Dynare_Run.md`; the capped pile-ups themselves are
rare (4 of 7 875 successful prior draws in the Windows log) but each cost
minutes and was returned with `converged = 0`. The equation is a DARE and is
solved directly, without the Control System Toolbox (`dare`/`idare` are not
available on the production machine):

1. `'dynare'` — Dynare's `kalman_steady_state` MEX (SLICOT), which ships with
   every Dynare build and solves the DARE without cross term. The ABCD form
   has the cross term `N = B D'`; with `R = D D'` invertible the standard
   substitution `Abar = A - N R^{-1} C`, `Qbar = B B' - N R^{-1} N'` reduces
   it to the MEX's form (same `Omega`). Used when the MEX is on the path and
   `rcond(R) > 1e-10`.
2. `'qz'` — own solver on the extended `(2n+ny)` pencil (Van Dooren 1981 /
   Arnold–Laub 1984) with base `qz`/`ordqz`; no invertibility of `R` or `A`
   needed. `Omega = U2 U1^{-1}` from the stable deflating subspace.
3. `'iterate'` — the original fixed-point iteration, kept as fallback and
   reference (the old code also remains as comments at the three call sites).

A candidate is accepted only if it is symmetric, PSD to tolerance and its
relative residual `r = ‖Omega − f(Omega)‖/‖Omega‖`, **scaled by
`1/(1 − rho(A−KC)^2)`**, is below 1e-8. The scaling is the point: `r` is one
step of a contraction with rate `rho^2`, so the distance to the fixed point is
about `r/(1−rho^2)`; a capped iteration on the pile-up draw shows `r = 2e-9`
while being 4e-5 off. `ric_flag` keeps the four fields the wrapper stores in
cell 19, with two changes of meaning: `.iterations` is 0 for a direct solve
(> 0 only when the fallback ran) and `.rel_step` is that scaled residual, i.e.
an **error estimate**, not the last iteration step. `.method` and
`.residuals` are added.

**Step 5b (per-shock error-variance decomposition)** had the same problem one
level down: seven Lyapunov recursions `Sigma = AKC Sigma AKC' + b_u b_u'`, each
capped at 1e4 steps and hitting the cap at `rho ≈ 1`. They are now solved by
the doubling (Smith) iteration `S ← S + M S M'`, `M ← M^2` (local
`dlyap_doubling` in the engine, ~25 doublings at `rho = 0.99995`); the
`parfor` loop is kept as a comment.

**Verification (Octave 11.3 + Dynare 7.1, HLW posterior-mean and posterior-mode
workspaces, `merge_check/check_dare.m`).** Both direct methods agree with each
other to ~1e-15 and with a 1e-13-converged iteration to ~1e-11. Against the
pre-change engine (iteration at 1e-10 step) every Riccati-dependent output —
`CorrH`, `R2_h`, `Error_ACF`, `FEV`, `weighted_HL`, `nu`, `fh`, `K_gain`,
`KF_Omega`, `fr.*` — differs by 1e-9 to 2e-8 relative, which is the *old*
iteration's own error (`step/(1−rho^2)` ≈ 1e-10/0.014); `shock_contributions`
(5b) agrees to 6e-15; `TV_*` and `Omega` are unaffected. The only larger
differences are 28 entries of `fh_direct` for Dynare `AUX_ENDO_LAG_*` states at
horizons where revision and forecast variance are both exactly zero (0/0): the
old iteration from `Omega = 0` produced roundoff garbage there (NaN, −0.7,
4.5), the direct solution gives ~1e-15; these are not states of interest.
Timing per engine call: 2.13 s → 0.83 s at the mean, 3.73 s → 0.93 s at the
mode, and the pile-up variant (`sigma_g`, `sigma_z` ÷ 100, `rho = 0.999969`,
`nu(rrstar) = 8359` quarters) 1.3 s where the old code ran three capped
Riccati loops plus seven capped Lyapunov loops (minutes). On 30 random stable
systems the methods agree to 4e-13; when `D D'` is singular the `'dynare'`
route is skipped and `'qz'` solves.

### Square, invertible systems: `'invertible'` shortcut (2026-10-16)

**Problem.** Smets–Wouters 2007 *without* measurement error (7 shocks, 7
noise-free observables, four of them first differences: `dy dc dinve dw`) has
transfer-function zeros at z = 1, so the DARE has no stabilising solution
(`rho(A-KC)` = 0.99999 in every draw). `'dynare'` and `'qz'` cannot return one,
and the `'iterate'` fallback approaches the limit only at rate ~1/t: it hit its
1e5 cap after ~20 s per engine pass with residual ≈ 0.5 and `converged = 0`.
Downstream, 0/500 posterior draws converged, half-lives were artefacts and the
Prop. 5.1 identity check failed (gap ≈ 1.1).

**Fix: guess and verify Omega = 0.** With D square and invertible, (*) at
Omega = 0 gives `B B' - B D' (D D')^{-1} D B' = 0`, so Omega = 0 is always a
fixed point. The gain is K = B D^{-1} and the closed loop is
`F = A - B D^{-1} C`, whose eigenvalues are the transmission zeros:
- rho(F) < 1: the stabilising solution (`ric_flag.method = 'invertible'`);
- rho(F) = 1: the limit of the Riccati recursion (maximal solution), approached
  at ~1/t (`'invertible_critical'`; SW: exactly four eigenvalues of F at 1);
- rho(F) > 1: not the limit; the shortcut is skipped and the other methods run.

`solve_steady_state_riccati` tries this first when `size(D,1) == size(D,2)` and
`rcond(D) > 1e-10`, and accepts it when `rho(F) <= 1 + opts.tol_unit` (1e-6) and
`||f(0)|| / ||B B'||` is below `tol_res_raw` (1e-10). It needs its own acceptance
test because the general one divides by ||Omega|| and 1 − rho². Every call site
passes its own ABCD for (*), so the shortcut applies wherever D is square.

A side finding: before this change, exactly invertible systems with rho(F) < 1
were also flagged as not converged. QZ returns Omega ≈ 0, but the relative
residual test rejects it.

**Effect on the SW baseline (posterior mean, one core).**
- Engine pass: 8.2 s → 0.44 s.
- Riccati flag: `converged = 1`, `rho_cl` = 1.
- `KF_Omega` = 0, `CorrH` = `R2_h` = 1 exactly.

Outputs that are 0/0 at Omega = 0 are NaN or `[]`, the engine's convention for
undefined values: `Error_ACF`, the FEV shares, `max_HL`/`weighted_HL`/`all_HL`,
`shock_contributions`, `nu`.

**Interpretation.** In the critical case the time-varying filter on the actual
sample (`TV_*`) still has a nonzero error at T; R² = 1 is the steady-state
limit.

**Test.** `test_invertible_riccati.m` covers four cases:
1. random systems with rho(F) < 1 (9–11 of 30 draws qualify), comparing the
   shortcut with QZ;
2. an MA(1) with its zero at 2, where the shortcut must not fire;
3. an MA(1) with a unit root (`'invertible_critical'`, and the iteration still
   at 1/(k+1));
4. a non-square D.

### QZ under MATLAB with complex eigenvalues (FIXED 2026-10-20)

`dare_qz` reads the pencil's eigenvalues with `ordeig`, which handles the 2×2 blocks of
the real generalised Schur form. It tested `exist('ordeig', 'file')`, which is 0 in
MATLAB: there `ordeig` is a built-in (`exist` returns 5; in Octave it is an m-file, 2).
So MATLAB always used the fallback `diag(AA)./diag(BB)`, which gives wrong moduli for
complex pairs. Whenever the closed loop had complex eigenvalues, the count of
eigenvalues inside the unit circle was off, `'qz'` was rejected, and the draw fell back
to the fixed-point iteration.

Example: the UC3-i2 posterior mean counted 115 instead of 114, although the stable and
unstable sets are cleanly separated (0.988 against 1.012). Every posterior draw iterated
(median 716 steps), while most prior draws, with real spectra, solved directly. The check
now accepts a file or a built-in.

**Effect** (short UC3-i2 run, one core):
- posterior stage 116 s → 94 s;
- error estimate 3e-9 → 3e-14;
- every draw solved directly;
- summary table unchanged.

On a Windows machine without the Parallel Computing Toolbox, posterior draws had taken
25 s each, about 4× the prior draws there, against 1.3× here. That is consistent with
the iteration being the expensive step on that machine.

Every model whose closed loop has complex eigenvalues was affected under MATLAB (UC3
levels too). Octave runs were not. Test: `test_riccati_qz_complex.m`.

## Implied reliability from the proxy regression (`fr_implied_R2.m`, 2026-09-17)

Proposition 5.1 of the paper ("Feasible smoother regression: bound and
identification") gives the population R² of the feasible regression of the
smoother revision on the filter revision as
`P^{2,m}_{i,h} = (V_h − Ω_ii)/(V_h − [Ω_m]_ii)`, `V_h = Var[x_{i,t+h}|y^t]`, with
α = 0, β = 1. Eliminating `V_h` against `R²_{i,h} = 1 − Ω_ii/V_h`,

```
R²_{i,h} = P (1 − ρ) / (1 − ρ P),      ρ = [Ω_m]_ii / Ω_ii .
```

`fr_implied_R2(P, omega_ratio)` evaluates this (P any array as in the `fr`
struct; `omega_ratio` scalar, same-size, or J × n_m broadcast over h). To make
it usable from a saved draw alone, `fr_regression_theory` now takes `KF_Omega`
as its 4th input and fills two more `fr` fields: `omega_ratio` (J × n_m, the ρ
above) and `R2_h_base` (J × H, the base state's eq.-(3.11) `R²_{i,h}` =
`Var(r_h)/(Var(r_h)+Ω_ii)`). Checks on the HLW posterior mean: `R2_h_base`
equals the engine's `R2_h` rows to 1e-14; `fr_implied_R2(fr.R2_theory,
fr.omega_ratio)` reproduces it to 1e-15 for m = 4 and 20 (draw by draw over
500 posterior draws: 6e-13). Feeding the OLS `fr.R2_emp` instead gives the
*data's* estimate of the filtered state's reliability under the model's ρ;
`HLW_plots.m` ("Implied R^2 from the empirical proxy regression") plots it
draw by draw against the model's `R²_{i,h}` and writes
`Tables_level/table_implied_R2.tex`; `model_assessment_plots.m` does the same
at the mode with the sampling bands of the assessment. Results saved before
this change (no `omega_ratio` in cell 20) are handled: `HLW_plots` rebuilds ρ
from cell 6, `(sd(l<m>NAME)/sd(NAME))²`.

`octave_compat/` holds stand-ins for `contains` and `histogram`; the plotting
scripts add it to the path only under Octave (`exist('OCTAVE_VERSION',
'builtin')`) — never `addpath` it in MATLAB.

## Posterior-predictive assessment with a re-estimated mode (`model_assessment.m`, 2026-09-15; in this folder since 2026-09-18)

`model_assessment.m` and `model_assessment_plots.m` live here (moved from
`HLW_new/hlw_level` on 2026-09-18; the stale copies in `hlw_spf/` and `hlw_ecm/`
were deleted, and `main_run_hlw_ecm.m` / `main_run_hlw_spf.m` call the function
API with `opts.model_dir` = their own folder). `model_assessment(D, opts)` reads the Dynare workspace from
`opts.model_dir` (default `Code/HLW_new/hlw_level`) and writes its result there;
`model_assessment_plots(Rfile, opts)` resolves relative `graph_dir`/`table_dir`
against the folder of `Rfile`.

`model_assessment(D, opts)` replaces the old 50,000-replication script. For each
of D datasets simulated at the posterior mode θ̂ it (i) computes the empirical
statistics at θ̂ on the data, (ii) at θ̂ on the simulated data, and (iii) at a
**re-estimated** θ̂^(d) obtained by damped Newton steps from θ̂ with a fixed
curvature matrix — the one-step (Le Cam) idea of `onestep_mode_note.tex`, made
robust after the plain one-step failed on this model. Findings that shaped it
(all verified numerically on `HLW_2023_US_no_lags`, see the function header):

- **Frozen heteroskedastic scales.** Dynare evaluates the `scales` expressions
  of `heteroskedastic_shocks` once at the calibrated values; the κ parameters
  never enter its likelihood (κ2021 = κ2022 = 1.038172 at the mode, the prior
  mode). θ̂ is the mode of the objective with the *stored* Qscale, so that is
  what the assessment uses (`opts.qscale = 'stored'`; the mode condition
  s(y) + ∇π(θ̂) = 0 then holds to 0.000 sd in every parameter). `build_qscale`
  (`'rebuild'`) assumes θ-dependent scales that the estimation never used.
- **Information-matrix failure.** Dynare's Σ is exactly the inverse observed
  Hessian on y (checked against our own numerical Hessian), but the Fisher
  information at θ̂ — covariance of the scores over datasets simulated from the
  model — disagrees: generalised eigenvalues 0.09–1.7 (`phi` sd 0.011 vs 0.029,
  `by` 0.024 vs 0.042). The note's step θ̂ + Σ[s(ỹ) − s(y)] therefore
  overshoots by up to 10× on model-generated data (log posterior on ỹ fell by
  100+). Default curvature is `Σ_F = [I_F − ∇²log π]⁻¹` from the simulated
  scores, steps are damped by backtracking on the log posterior of ỹ, up to 3
  iterations. The information comparison is printed as a diagnostic in its own
  right.
- Parameterisation is natural units (where Dynare's `csminwel` works and Σ
  lives), with the prior support enforced in the line search.

## Models solved outside Dynare: user solver and proper initialisation (2026-09-26)

Two options make the pipeline run on a model whose parameters are NOT mapped
to a state space by Dynare (first use: the Bauer–Rudebusch ESE term-structure
model, sampled in R, `Code/BauerRudebusch/`).

**User-supplied solver.** Every per-draw routine (`reliability_draw_wrapper`,
`simulate_artificial_data_wrapper`, `artificial_data_wrapper`,
`posterior_predictive_wrapper`, `model_assessment`'s local solvers,
`evaluate_at_reweighted_means`) now obtains the decision rule through one
function,

    [info, oo_, options_, M_] = solve_model_draw(M_, options_, oo_)

which runs `stoch_simul(M_, options_, oo_, M_.var_list_)` exactly as before
(verbose or under `evalc`) **unless** `M_.solve_function` is set, in which case
it calls the user map `S = f(M_.params, M_)` (name or handle; a name survives
save/load) and writes `S.A, S.B, S.C, S.D, S.obs_ss` into `oo_.dr` through
`dr_from_abcd` — the toolbox timing `x_t = A x_{t-1} + B e_t`,
`y_t = obs_ss + C x_{t-1} + D e_t`, states first in `M_.endo_names`,
observables after. `dr_from_abcd` also sets the smoother bookkeeping
(`oo_.bayestopt_*`, `M_.nstatic/npred`, `Sigma_e = I`, `H = 0`, …), so the
engine and the diffuse filter run unchanged (they still need Dynare's
`kalman_transition_matrix` / `compute_Pinf_Pstar` / `set_all_parameters` on
the path). `S.ok = false`, an empty `A` or an error mark a rejected draw
(`info = 1`, the callers' NaN path). Optional `S.X0`/`S.B0` override the
initial distribution per draw. `info` is now a SCALAR in both paths (Dynare's
first element), so the `info ~= 0` tests of the callers are robust.
`estim_params_from_names(param_names, M_)` builds the minimal
`estim_params_`/`bayestopt_` pair (and adds the `Sigma_e`, `Correlation_matrix`,
`Skew_e`, `H`, `Correlation_matrix_ME` fields Dynare 7's `set_all_parameters`
reads unconditionally) so that `set_all_parameters(x, …)` writes the draw into
`M_.params` one-to-one. `run_reliability_draws(..., opts)` takes the draws
themselves through `opts.draws` (k × n; all columns, or a random subset when
n_draws is smaller) instead of Dynare's prior sampler / MH files.
`build_qscale` now starts from the STORED `Qscale` (missing-data devices
survive) before overwriting the HLW kappa rows.

**Proper initialisation.** `M_.X_InitialDistribution.type` = `'diffuse'`
(default; exact diffuse unit roots, stationary block at its unconditional
covariance, `X0` the initial mean of the data pass, `B0` only seeds the
Riccati solve) or `'proper'`: the time-varying filter starts from
`x_0 ~ N(X0, B0 B0')` with no diffuse phase (`t_proper = 1`), the convention
of a sampler that initialises its Kalman filter with a wide finite variance
(ese.R: `a1 = 0.02` for i*, `P1 = I`). Implemented by reviving the legacy
`Omega_star_init`/`Omega_inf_init` arguments of `durbin_koopman_diffuse_filter`
(state block of the stacked `Pstar` set, `Pinf = 0`; a non-zero user `Pinf`
is refused) through the new 11th argument `init_type` of
`kalman_filter_steady_and_diffuse`. The population objects (`K`, `Omega`,
`CorrH`, half-lives, …) do not depend on the choice; `Sequence.init_type`
records it. Timing: the toolbox places the prior on `x_0`, the state that
`y_1` loads on with one lag; a sampler in the `y_t = Z α_t` convention places
it on `α_1 = x_1`. To match such a sampler exactly, back-shift its prior,
`X0 = A \ a1`, `P0 = A^{-1} (P1 − B B') A^{-T}` (a user solver can return
these per draw through `S.X0`/`S.B0`, as `ese_solve_draw` does); placing
`N(a1, P1)` on `x_0` directly is NOT immaterial even for a wide `P1` —
measured on the BR model against KFAS: 0.28 pp in the first filtered i* and
1.65 in the log-likelihood, versus 2e-7 pp and 6e-6 with the back-shift
(`Code/BauerRudebusch/check_kfas_vs_toolbox.{R,m}`, which also confirms the
Qscale missing-observation device against KFAS's NA handling to the same
precision).

**Model assessment for such a model.** `model_assessment` takes the objects a
Dynare workspace would have supplied through options: `opts.theta_hat` (mode,
free parameters), `opts.Sigma_H` (its covariance), `opts.lprior` (log prior
handle) and `opts.bounds` (support for the line search); it keeps a
workspace's own `M_.X_InitialDistribution` and reads the smoother bookkeeping
from `oo_` when `bayestopt_` lacks it. Parameters fixed by a normalisation must
be excluded from the free vector (zero posterior variance breaks the finite
differences) — `estim_params_from_names` with a subset of names does that, the
fixed values staying in `M_.params`. Also fixed on 2026-09-26, for every model:
the artificial datasets now set the cells of a missing-data device
(`opts.mask_ignore_shocks`, Qscale ≠ 1) to the placeholder `obs_ss` instead of
simulating them at the device's scale; the SPF assessment's `spf_obs` RMSE
band (median 199 against 0.86 in the data) was an artefact of the old
behaviour, every other observable was unaffected. BR driver:
`Code/BauerRudebusch/main_run_br_assessment.m` with `ese_log_prior.m`.
The mode itself comes from `find_posterior_mode` (2026-09-26): it maximises the
toolbox's own log-likelihood (solve → `loglik_only` engine → `apply_tv_filter`,
stored Qscale, the model's initialisation) plus the log prior over the free
parameters, in posterior-s.d. units, with Dynare's `csminwel1` (fallback
`fminsearch`), infeasible points (bounds, prior support, solver rejection)
handed to the optimiser through csminwel's penalty contract; a central-difference
Hessian at the mode gives `Sigma_H`. A mode found by another program (the R
sampler's optimiser) is only the starting point, since its objective differs in
constants, tolerances and initialisation; `model_assessment`'s mode-condition
line then checks the result (it should read ≈ 0 s.d. in every parameter).

**Tests.** `test_dynare_match` (regression against Dynare on the HLW fixture
`Code/HLW_new/hlw_level/tmp.mat`: the default diffuse/TV filter still equals
`oo_.UpdatedVariables` to 6e-13 for ygap, rrstar and drrstar after the diffuse
phase, `build_qscale` is unchanged on HLW, and both sampler paths of
`run_reliability_draws` run — re-run it after any change to the solve or
filter plumbing); `test_solve_model_draw` (Part A: `solve_model_draw` ≡ direct
`stoch_simul` on the HLW fixture, ghx/ghu/ys identical; Part B: toy
trend + cycle model through an anonymous map — proper init reproduces a
hand-written Kalman filter to 9e-16 in the states and 3e-13 in the
log-likelihood, diffuse and proper agree to 9e-16 after 100 periods,
population objects identical, rejected draws, wrapper and
`run_reliability_draws` with `opts.draws`);
`Code/BauerRudebusch/test_ese_solve_draw` (the MATLAB JSZ port `jsz_loadings` +
`ese_solve_draw` reproduce the R export's A, B, C, D, obs_ss for all 10 000
posterior and 10 000 prior draws to 3e-15).

## Name mapping (2026-07 renaming)

| This folder | Formerly (Miscellaneous/) |
|---|---|
| `reliability_draw_wrapper.m` | `f_wrapper_Extract_ABCD_par_while.m` |
| `compute_state_reliability_diagnostics.m` | `f_Extract_ABCD_par_while.m` |
| `steady_state_kalman_filter.m` | `my_KalmanFilter_while.m` |
| `kalman_filter_steady_and_diffuse.m` | `my_KalmanFilter_Heteroskedastic_jan27.m` |
| `durbin_koopman_diffuse_filter.m` | `f_kalman_filter_diffuse_hybrid.m` |
| `diffuse_filter_initialization.m` | `f_compute_Pinf_Pstar_dynare_style.m` |
| `robust_cholesky.m` | `f_robust_chol.m` |
| `forecast_errors_and_revisions.m` | `f_forecast_stats.m` |
| `estimation_error_half_life.m` | `half_life_by_shock.m` |
| `theoretical_R2_wedge.m` | `compute_newdeltaR2_wedge.m` |
| `make_positive_semidefinite.m` | `f_make_psd.m` |

### Simulation ("vs truth") support, added 2026-07-14, MERGED IN 2026-08

| This folder | Formerly |
|---|---|
| *(merged into `compute_state_reliability_diagnostics.m`)* | `f_Extract_ABCD_par_while_andrea.m` → `compute_state_reliability_vs_truth.m`, **superseded 2026-08** |
| `state_forecast_errors_and_revisions.m` | `f_forecast_stats_states.m` (recovered from Archive) |
| `R2_wedge_factor.m` | `compute_wedge_f.m` (recovered from Archive) |

`compute_state_reliability_vs_truth.m` was a separate engine that duplicated
the main one: 446 of its 486 code lines were verbatim copies, and every fix
had to be written twice. Two of its three differences are now part of
`compute_state_reliability_diagnostics`:

- an **optional 6th input `true_states`** (T x n_states, columns ordered like
  the pruned state vector). When supplied, the engine also returns
  state-level h-step forecast RMSEs and error/revision correlations against
  the truth, as `.TV_RMSE_states` / `.TV_CorrY_H_states` (`.steady_*_states`
  in the superseded file). NaN when it is absent, so nothing changes for
  existing callers;
- the **direct eq.-(26) wedge** from `R2_wedge_factor`, as `.fh_direct`.

⚠️ **The wedge names are not the same as in the superseded file.** Both
engines compute both objects, under opposite names:

| object | here | superseded `*_vs_truth` |
|---|---|---|
| `theoretical_R2_wedge` ("newdeltaR2") | `.fh` | `.fh_shift` |
| `R2_wedge_factor` (direct eq.-26) | `.fh_direct` | `.fh` |

`.fh` keeps the meaning it has always had in this engine, so no current
caller changes; scripts written against the superseded file must read
`.fh` → `.fh_direct`.

The third difference is **deliberately not carried over**: the superseded
engine forced the time-varying filter to be homoskedastic (`Qscale = 1`
regardless of `M_.heteroskedastic_shocks`). The merged engine honours
`M_.heteroskedastic_shocks` as it always has. For a model with no
`heteroskedastic_shocks` block (`uc3_level_i2`) the two are identical —
`Qscale` defaults to ones — so this only matters for models that actually
rescale shock variances (HLW/LW), where the old behaviour was arguably wrong
anyway.

### Tail-unreliability index, added 2026-08

| This folder | Formerly |
|---|---|
| `tail_unreliability_index.m` | *(new)* |

Implements the horizon-free reliability index of the note
*"A horizon-free reliability index for non-stationary states of arbitrary
integration order"* (`../Rate_of_convergence_R2.pdf`, Definition 2.1 and
Theorem 2.2, eq. (5)):

```
h^(2d+1) * (1 - R2_{h,i})  ->  nu_i^(d) = (2d+1) * Omega_ii / L_i^(d)
```

where `Omega` is the steady-state filtered error variance (`KF_Omega`),

```
d_i      = max{ k >= 0 : e_i' N^k B ~= 0 },        N^0 = Pi
L_i^(d)  = || e_i' N^d B ||^2 / (d!)^2
```

`Pi` is the Riesz projector on the *generalized* unit-root eigenspace of `A`
and `N = (A - I)*Pi` is the nilpotent part that records how the unit roots are
chained. `d_i` is state `i`'s integration order minus one (state `i` is
I(`d_i`+1)), and `L_i^(d)` is the innovation variance of the top rung of its
cumulation ladder. `nu` is scale-invariant and is the replacement for the
mechanically-converging long-horizon `R2_h` columns.

**Units.** `nu_i` is in quarters^(2`d_i`+1) — in quarters only where `d_i = 0`.
Indices of states with different `d` are **not comparable** and must not share
a table column (note, section 2.3).

**Revised 2026-09.** The original implementation computed the `d = 0` case
only, `h(1 - R2_{h,i}) -> Omega_ii / L_i` with `L_i = [Pi B B' Pi']_ii`, and
returned `NaN` for every state of higher order. Eq. (5) reduces to that at
`d = 0`, so all previously reported numbers are unchanged (verified to 4e-16
over 400 random systems, 1129 states, no state losing an index); what is new is
that I(2)+ states are indexed instead of dropped.

`compute_state_reliability_diagnostics` calls it in step 5e and returns the
result as `LatentStateDiagnostics.nu` (n_key_states x 1, ordered like
`states_of_interest`, i.e. like the rows of `CorrH`) together with
`LatentStateDiagnostics.nu_order`, the matching orders `d_i`. `nu` is `NaN` for
states without a stochastic trend at any rung (e.g. the output gap) — report
the paper's unconditional `R2` for those instead. Only `nu` and `nu_order` are
reported: the other objects of the note (`L_i^(d)`, `rho(A_s)`, the ladder of
rung norms, the onset horizon) stay in `info`. Skipped (NaN) in `fast_mode`.

### Filtering vs smoothing: eq. (3.17), added 2026-08

| This folder | Formerly |
|---|---|
| `fixed_lag_smoother_map.m` | *(new)* |

Table 3 of the draft puts the **filter** and the **5-yr fixed-lag smoother**
of `r*` in adjacent columns, and Section 3.5 says reliability must be weakly
increasing in the hindsight `m`. It was not: the smoothed `R2` came out
*below* the filtered one (0.02 vs 0.07 at h = 4, 0.20 vs 0.28 at h = 20).

The cause is the **denominator**, not the filter. The fixed-lag smoother is
implemented by augmenting the model with a lagged copy of the state,

```
l20rrstar(t) = rrstar(t-20)      // HLW_2023_US_TD_small.mod
```

so the *filtered* estimate of that coordinate is the m = 20 *smoothed*
estimate of `r*`, and `KF_Omega`'s diagonal entry for it is already the
smoothing error variance `[Omega_20]_ii`. But `CorrH.^2` is eq. (3.11)
applied coordinate by coordinate, so it normalised that numerator by
`Var[r*_{t+h-20} | y^t]` — the conditional variance of the **earlier** state,
short by 20 quarters of accumulated innovation variance for an I(1) trend.
This is precisely the case Section 3.5 warns about ("if each estimate is
allowed to generate its own denominator, the ratios need not inherit that
ordering"), and at h = m the reported number degenerates into
`1 - Omega_20/Omega_0`, a ratio of error variances rather than a reliability.

Step 5c of the engine now also returns **eq. (3.17)**,

```
R2^(m)_{i,h} = 1 - [Omega_m]_ii / Var[x_{i,t+h} | y^t],
```

which keeps the numerator and takes the denominator — eq. (3.10) — from the
**contemporaneous base state** at the **same horizon**, for every `m`. With a
common denominator the error-variance ordering `Omega_l <= Omega_m`
(eq. (3.16)) carries over one-for-one, so reliability is weakly increasing in
`m` by construction.

`fixed_lag_smoother_map` supplies the missing link — which states of interest
are lagged copies, of what, and at which lag. It does **not** parse names: in
the pruned system it compares impulse responses, `Psi_j(s) = e_j' A^(s-1) B`,
and accepts `j` as the m-period lag of `k` iff `Psi_j(s) = Psi_k(s-m)` for all
`s`. The number of leading all-zero response blocks pins `m` down uniquely and
the match at that `m` is the root of the lag chain. This works whether the
base is a genuine Dynare state or one of the augmented non-state variables of
interest (`rrstar = c*gg + zz` is the latter), and survives renaming. Set
`M_.fixed_lag_smoother_map = {'l20rrstar','rrstar',20; 'ldrrstar','drrstar',1}`
to pin the map by hand instead, or `cell(0,3)` to switch the correction off
and restore the eq.-(3.11) numbers.

**`CorrH` for smoothed states (changed 2026-09-15).** For a filtered state
`CorrH` is eq. (3.11) and `CorrH.^2 == R2_h` (verified to 4e-13). For a
fixed-lag smoothed state `l<m><NAME>` the doubled-system recursion gives
eq. (3.11) applied to the *lagged copy*, corr(x_{t+h−m}, x̂_{t+h−m|t+h} | y^t),
with the earlier state's variance in the denominator — the object Section 3.5
warns about. `CorrH` for those rows is now instead the conditional
correlation corr(x_{i,t+h}, x̂_{i,t+h|t+h+m} | y^t) = `sqrt(R2_h)`: since the
smoothed estimate is a conditional expectation under a finer information
set, Cov(x, x̂ | y^t) = Var(x̂ | y^t) and the squared conditional correlation
(the R² of regressing x_{t+h} on x̂ given y^t) is 1 − [Ω_m]_ii / Var[x_{t+h} | y^t]
= R^{2,m}_{i,h}. So `CorrH.^2 == R2_h` on every row, and the h-profile,
histogram and Bayes-factor figures built from cell 2 are comparable across
information sets. (The engine keeps the old eq.-(3.11) value only for its
verbose comparison print.)

New outputs:

| where | field / cell |
|---|---|
| `LatentStateDiagnostics` | `.R2_h` (n_key x MaxH), `.smoother_lag_m`, `.smoother_base`, `.smoother_map_info` |
| `reliability_draw_wrapper` | cell **10** `R2_h`, cell **11** `smoother_lag_m`, cell **12** `smoother_base` |

These three **reuse the placeholder slots 10-12**, which had been left empty
since the steady-state data pass was removed (see "One filter" below), so
every other index keeps its meaning and position. A failed draw writes a
scalar `NaN` in 10-11 and `{}` in 12. Like `nu`, the three are population
objects and are filled whether or not `M_.data` exists (NaN / empty in
`fast_mode`).

⚠️ Any script still reading cells 10-12 as the *pre-2026-08* steady-state
filter output was already getting empties; it will now get these instead. The
shapes differ, so it will fail loudly rather than quietly. `HLW_plots.m` /
`LW_plots.m` / `UC3_plots.m` are the ones to check.

`HLW_plots.m` reads cell 10 for rows 4-5 of the summary table, prints an
eq.-(3.17)-vs-eq.-(3.11) comparison, and checks the Section 3.5 ordering draw
by draw at every horizon. It tests whether cell 10 *contains* a matrix the
size of `CorrH`, and falls back to eq. (3.11) with a warning when it does not.
The engine runs the same two consistency checks internally and writes to
stderr if either fails.

Verified (2026-08, `test_fixed_lag_smoother_R2`, 9/9 pass, Octave 11.3; no
Dynare, no data) on an HLW-shaped test model — I(1) `r* = c*g + z`, a 20-link
lag chain, an IS/PC block, `r*` identified only through a noisy output gap,
`A` rounded to 1e-6 as step 2 does:

- eq. (3.17) reproduces `CorrH.^2` on every m = 0 state to 4e-13;
- `R2^(20) >= R2^(4) >= R2^(0)` at all 250 horizons (gaps +3.1e-2 and +1.3e-2
  at the tightest horizon), and `R2^(1) >= R2^(0)` for `drrstar`;
- the old measure reproduces the Table 3 pathology (smoother *minus* filter
  = -0.056 at its worst horizon; -0.041 at h = 4 and -0.006 at h = 20, the
  same sign and shape as the published 0.02/0.07 and 0.20/0.28).

(That test script is not in the repository.)

### I(2) states / defective unit roots (added 2026-08, resolved 2026-09)

HLW's unit eigenvalue is *not* semisimple: `ystar = ystar(-1) + gg(-1)/4 + shock`
with `gg` a random walk makes the trend block `[1 1/4; 0 1]`, a **Jordan block
at 1**, because `ystar` is I(2). Assumption 1 of the revised note admits that
case on purpose. Three stages of this file's history:

1. the first implementation detected the defect and returned `NaN` for *every*
   state;
2. the 2026-08 version switched to the **Riesz projector on the generalized
   eigenspace**, which exists either way, and rescued every state off the
   chain (`e_i'N = 0`) while still returning `NaN` on it;
3. the 2026-09 version indexes the chained states too, at their own order.

The order is now read off the ladder of row norms `||e_i' N^k Pi B||`,
`k = 0..(#unit roots - 1)`, stored in `info.rung_norms`: `d_i` is the highest
rung above the relative threshold, `NaN` if the whole row is zero.

In HLW, `gg`, `zz` and `rrstar = c*gg + zz` come out at `d = 0` — the last
being the note's own worked example, `L_{r*}^(0) = c^2 sigma_g^2 + sigma_z^2`,
reproduced exactly — and `ystar` at `d = 1` with `L^(1) = sigma_g^2/16`, where
before it had no index. When the unit eigenvalue is semisimple, `N = 0`, every
`d_i` is 0, `Pi` is the ordinary spectral projector and results are unchanged
(UC3 and the note's designs verified bit-for-bit).
`info.unit_root_defective`, `info.jordan_chain_states`, `info.d` and
`info.max_order` report the diagnosis; the projector is computed by ordered
real Schur + one Sylvester solve, with the eigenvector and `lim A^j` routes
kept as fallbacks and cross-checks.

`reliability_draw_wrapper` returns `nu` as **cell 17**, its orders as
**cell 18** and the Riccati flag as **cell 19** (the width grew 17 → 18 → 19
by appending, every existing index unchanged), so prior/posterior draws of the
index are collected the same way as the other diagnostics:

```matlab
nu_draws = cat(2, oo_.posterior_function_results{:, 17});   % n_key_states x n_draws
nu_prior = cat(2, oo_.prior_function_results{:,     17});
nu_ord   = cat(2, oo_.posterior_function_results{:, 18});   % same shape
ric      = cat(1, oo_.posterior_function_results{:, 19});   % n_draws x 4
```

Check that each row of `nu_ord` is constant across draws before taking medians:
rows with different `d` are in different units. `HLW_plots` does this and errors
out if any order is nonzero, since its summary tables print `nu` as quarters.
Drop failed draws with the usual `isscalar`/`isreal` masks on columns 2-3
first. The constant width is required — Dynare's
`execute_prior_posterior_function` sizes its storage from the first draw and
then assigns `results_cell(i,:)`, so successes and failures must agree — and
appending rather than inserting keeps the positional indexing used by
`HLW_plots.m` / `LW_plots.m` / `HLW_plots_TV.m` valid.

Verified against the three numerical designs of Section 1.8 of the note:
`Pi`, `L_1` and `nu_1` reproduce the published values (Design I
`nu_1 = 10.0154`, Design II `nu_1 = 0.3062` with `L_1 = 0.40` against an own
innovation variance of `0.04`, Design III `nu_1 = 6.8327`), Design I's full
`R2_h` / `h(1-R2_h)` table matches to four decimals, and `h(1-R2_{h,i})`
converges to `nu_i` in all three designs and in the UC3 model.

## One filter: steady-state data pass removed (2026-08)

The engine used to filter the data **twice** — once with the constant
steady-state gain, once with the exact diffuse / time-varying filter — and
return two parallel sets of empirical statistics. The steady-state **data
pass** is gone. What remains:

- the steady-state **Riccati solve** stays (since 2026-09-16 solved directly,
  see "Direct Riccati solve"). It is not a filter run
  over data, it is the definition of the paper's `K`, `Omega`, `A-KC`
  objects, and it still drives every population diagnostic: `CorrH` (R²_h),
  `Error_ACF`, `FEV`, the half-lives, `nu`, `CorrH_Observables`, the wedge
  `fh`. None of these is defined for a time-varying gain;
- the **data** are filtered once, by the diffuse/TV filter (`apply_tv_filter`).

Why the diffuse filter and not the constant-gain one:

- its empirical statistics are **exactly invariant** to the level of the
  non-stationary states; the constant-gain filter's are not (it starts at
  `X0 = 0` while the trends sit at the level of the data — `100*log(GDP)`
  ≈ 816 in `uc3_level_i2` — and that transient decays at `max|eig(A-KC)|`,
  i.e. at the reported estimation-error half-life). The invariance is what
  makes the simulation-based null of `model_assessment.m` — whose
  replications draw random trend levels — the right reference distribution;
- it is the filter Dynare itself used (`lik_init=3`), so these statistics
  are consistent with the reported states, likelihood and posterior;
- it is the only filter defined when shocks are heteroskedastic
  (`Qscale ~= 1`), as in the HLW/LW applications.

**Burn-in.** One `IDX_ok` mask is built once and used by *every* empirical
statistic — correlations and RMSEs alike. Previously the burn-in was applied
only to the TV correlations: the two filters' correlations were computed on
different samples (UC3: R² gaps up to 0.47 that vanish to 0.0000 once the
samples match), and neither RMSE excluded anything, so both measured the
initialization rather than forecast accuracy (yobs at h = 4: 133pp before,
2.8pp after). The length is `M_.reliability_burn_in`, default 50. **50 is not
enough when the estimation-error half-life is long** — it is 72–81 quarters in
`uc3_level_i2`, so an 816pp initial error is still ~500pp at t = 50. Set the
field explicitly; `ModelMatrices.AKC` gives the decay rate.

**Interface changes.**

| before | after |
|---|---|
| outputs 5–7 `state_estimate`, `surprises`, `log_llhood` | removed (engine returns 7 outputs, `TV_*` at 5–7) |
| `.steady_CorrY_H`, `.steady_RMSE` | removed; use `.TV_CorrY_H`, `.TV_RMSE` |
| `.steady_CorrY_H_states`, `.steady_RMSE_states` (vs_truth) | `.TV_CorrY_H_states`, `.TV_RMSE_states`, now on the main engine |
| — | `ModelMatrices.KF_V`, `.KF_ric_flag`, `.obs_ss`, `.X0_pruned`, `.B0_pruned` |
| wrapper cells 10–12 | hold `R2_h`, `smoother_lag_m`, `smoother_base`; 13–16 keep their positions |

`DeltaR2 = CorrH_Observables.^2 - TV_CorrY_H.^2`: the theoretical side stays
the steady-state population object, the empirical side is the diffuse
filter's. In finite samples the two filters differ, so DeltaR2 is not zero
under correct specification — but the test compares the empirical DeltaR2 to
the *simulated* distribution of DeltaR2, and the simulations run through the
identical code path, so the wedge is absorbed.

## Verification history (Octave 11.3 + Dynare 7.1, 2026-07 renaming)

- All 11 renamed files parse; `Tests/unit_tests_renamed.m`: 18/18 pass.
- `Tests/uc3_test_renamed.mod` ran the full reduced pipeline end-to-end.
- A/B against `Miscellaneous/`: **all `posterior_function` outputs and the
  direct posterior-mean run were bit-identical** (max abs diff = 0).
  `prior_function` outputs can only be compared within a run: the prior
  draws themselves are freshly randomized on each Dynare run (Dynare 7's
  posterior subsample is deterministic given the MH chain; its prior
  sampler is not), which is unrelated to the renaming.

(The `Tests/` folder referred to above is not in the repository.)

## Using this folder from a .mod file

1. `addpath(genpath('<path to>/ReliabilityToolbox'));`
2. the two estimation-function strings:
   `prior_function(function='reliability_draw_wrapper', ...)` and
   `posterior_function(function='reliability_draw_wrapper', ...)`
3. the direct call near the end:
   `[LatentStateDiagnostics, ModelMatrices, obs_indices, key_state_indices, TV_state_estimate, TV_surprises, TV_log_llhood] = compute_state_reliability_diagnostics(M_.observable_names, ...)`

Under Dynare 7, also copy `reliability_draw_wrapper.m` next to the `.mod`
file before running (Dynare 7 looks for prior/posterior functions in the
current directory only). Keep only ONE copy of these functions on the MATLAB
path at a time — do not `addpath` `Miscellaneous_Archive` together with this
folder.
