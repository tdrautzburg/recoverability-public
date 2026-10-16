# recoverability-public

Replication code for the reliability ("recoverability") of latent-state estimates from linear
state-space models. A model is estimated with [Dynare](https://www.dynare.org); a Kalman-filter
toolbox (`Code/ReliabilityToolbox`) then computes, draw by draw over prior and posterior samples,
how well the filter and fixed-lag smoothers recover the latent states:

- $R^2_h$ — reliability of the state estimate $h$ quarters ahead (eq. 3.17; filter and fixed-lag
  smoothers share the denominator, the conditional variance of the contemporaneous state),
- $\nu$ — tail-unreliability index (defined for non-stationary states only; "--" otherwise),
- half-life of the estimation errors (first lag at which their autocorrelation falls below 1/2),
- $\Omega_{ii}$ — steady-state variance of the estimation error (reported as its s.d.),
- empirical RMSEs, forecast-error and revision regressions, and a posterior-predictive
  ("Le Cam") model assessment with fan charts and p-value tables.

The code accompanies

> Amisano, Gianni, Thorsten Drautzburg, and Andrea Stella (2026). "How to catch a star?
> Reliability of filtering estimates in linear state space systems." RTE Conference,
> Philadelphia, PA.

Please cite this paper when using the code.

Two applications:

| Application | Folder | Model |
|---|---|---|
| UC3 | `Code/UC3` | Unobserved-components model of output, inflation and unemployment (random-walk trends, common AR(2) gap), US 1960Q1–2019Q4 |
| Smets–Wouters, baseline | `Code/SmetsWouters` | Smets and Wouters (2007) DSGE model, 7 observables, no measurement error, US 1947Q3–2004Q4 |
| Smets–Wouters, calibrated measurement error | `Code/SmetsWouters` | Same model with fixed measurement error on output, consumption, investment and wage growth and hours |

## Software requirements

| Software | Version used for the checks in this repository | Notes |
|---|---|---|
| MATLAB | R2026b | The checks ran with the Statistics and Machine Learning and Parallel Computing toolboxes installed (the reliability draws use `parfor`; without the Parallel Computing Toolbox MATLAB runs `parfor` loops serially). Running without the Statistics toolbox was not tested. No Control System Toolbox is needed. |
| Dynare | 7.2 (official macOS arm64 package) | Other Dynare versions were not checked here. |
| Internet access | — | Only once, to download the Smets–Wouters data (`get_sw_data.m`). |

GNU Octave (11.x) with Dynare is supported by the toolbox code (it is written to stay
Octave-compatible), but the workflows in this repository were run under MATLAB only.

### Configuring Dynare

The drivers call `Code/ReliabilityToolbox/setup_dynare_path.m`, which uses, in this order:

1. a Dynare already on the MATLAB path,
2. the environment variable `DYNARE_ROOT` (the Dynare folder that contains `matlab/`),
3. under MATLAB on macOS, the newest `/Applications/Dynare/<version>`,
4. the Homebrew installation `/opt/homebrew/opt/dynare/lib/dynare`.

On Windows or Linux set `DYNARE_ROOT`, e.g. `set DYNARE_ROOT=C:\dynare\7.2` or
`export DYNARE_ROOT=/opt/dynare-7.2`. If Dynare's compiled `kalman_steady_state` is not
available, the toolbox prints a warning and uses its own (slower) Riccati solver; results are
unaffected.

## Data

- **UC3**: `Code/UC3/UC3_data_yobs_1960_2019.csv` (included).
- **Smets–Wouters**: `usmodel_data.mat` (observables `dy dc dinve labobs pinfobs dw robs`, 230
  quarters, 1947Q3–2004Q4) and `usmodel_mode.mat` (the published posterior mode, used as the
  starting point of the mode search) from Johannes Pfeifer's replication files,
  <https://github.com/JohannesPfeifer/DSGE_mod/tree/master/Smets_Wouters_2007>. They are not
  included. `Code/SmetsWouters/get_sw_data.m` downloads them from
  `https://raw.githubusercontent.com/JohannesPfeifer/DSGE_mod/master/Smets_Wouters_2007/usmodel_data.mat`
  (and `.../usmodel_mode.mat`) into `Code/SmetsWouters/`; the SW drivers call it automatically
  and skip the download when the files are present. Without internet access, place the two files
  in `Code/SmetsWouters/` by hand.

## Running the workflows

All commands are run from the repository root. Each workflow runs in one MATLAB session; run
the workflows in separate sessions (the plotting stage keeps figures and global state).

```bash
# 1. UC3 (full run: 2 x 50,000 MH draws, 10,000 + 10,000 reliability draws, 300 datasets)
matlab -sd Code/UC3 -batch "main_run_uc3_levels"

# 2. Smets-Wouters baseline (short illustrative chain)
matlab -sd Code/SmetsWouters -batch "main_run_sw_baseline"

# 3. Smets-Wouters with calibrated measurement error (short illustrative chain)
matlab -sd Code/SmetsWouters -batch "main_run_sw_me"
```

(`-sd` sets MATLAB's start folder; equivalently, start MATLAB, `cd` into the folder and type the
script name.) The three workflows are independent and can be run in any order.

Measured run times on the machine used for the checks (18 cores, MATLAB R2026b, Dynare 7.2;
Dynare + draw stages / `model_assessment`, plotting adds a few minutes):

| Workflow | Settings | Run time |
|---|---|---|
| UC3, quick check | `UC3_MH_REPLIC=2000 UC3_NDRAWS=200 UC3_MA_D=20` | 2 min / 0.2 min |
| UC3, full | defaults | not re-run for this repository (hours) |
| SW with measurement error | defaults | 2 min / 0.5 min |
| SW baseline | defaults, one core (`-singleCompThread`, no parallel pool) | 15 min / 3 min |

Most of the SW baseline time is in the prior draws. A posterior draw takes about 0.1 s because
the steady-state Riccati equation is solved in closed form ($\Omega = 0$; see Limitations).
Prior draws whose transmission zeros lie outside the unit circle fall back to the general
solvers.

Each workflow has the same three stages:

1. **Dynare**: estimation (posterior mode, then random-walk Metropolis–Hastings with 2 chains and
   20% burn-in), smoother at the posterior mean, then `run_reliability_draws` computes the
   reliability diagnostics for prior draws and posterior (MH) draws;
2. **`model_assessment`**: posterior-predictive assessment with a damped-Newton re-estimated mode
   on artificial datasets simulated at the posterior mode;
3. **`HLW_plots`** (main plots only) and `model_assessment_plots`: figures and LaTeX tables.

### Run length

Set environment variables before starting MATLAB (unset = default):

| Variable | Meaning | UC3 default | SW default |
|---|---|---|---|
| `UC3_MH_REPLIC` / `SW_MH_REPLIC` | MH draws per chain (2 chains) | 50000 | 5000 |
| `UC3_NDRAWS` / `SW_NDRAWS` | prior and posterior reliability draws (each) | 10000 | 500 |
| `UC3_MA_D` / `SW_MA_D` | artificial datasets in `model_assessment` | 300 | 50 |

Example (a quick UC3 check): `UC3_MH_REPLIC=2000 UC3_NDRAWS=200 UC3_MA_D=20 matlab -sd Code/UC3 -batch "main_run_uc3_levels"`.
The drivers pass the first two to Dynare as macros (`-DMH_REPLIC=... -DN_DRAWS=...`); the
defaults are also set inside the `.mod` files, so `dynare uc3_levels.mod` alone reproduces the
full UC3 settings.

### Smets–Wouters settings

`Code/SmetsWouters/Smets_Wouters_2007_45.mod` is Pfeifer's implementation with two variants,
selected by the macro `SW_ME` (set by the drivers; the file is never edited):

- **Estimated** in both variants: the 36 parameters of the replication file (7 structural shock
  standard deviations and 29 structural parameters), with the priors of Smets and Wouters (2007).
  Mode: csminwel (`mode_compute=4`) started at the published mode. Sample and likelihood
  settings as in the replication file (`first_obs=1`, `presample=4`, `lik_init=2`), i.e. all 230
  quarters; Smets and Wouters (2007) estimate on 1966Q1–2004Q4 (`first_obs=71`).
- **Fixed** in both variants: `ctou=0.025`, `clandaw=1.5`, `cg=0.18`, `curvp=curvw=10`.
- **Measurement error** (`main_run_sw_me.m`, `SW_ME=1`): five extra shocks `me_dy`, `me_dc`,
  `me_dinve`, `me_dw`, `me_lab` added to the observation equations of `dy`, `dc`, `dinve`, `dw`,
  `labobs` (none on inflation and the interest rate). Their standard deviations are **calibrated,
  not estimated**:

  | Variable | Meaning | Default |
  |---|---|---|
  | `SW_ME_SHARE` | measurement-error variance as a share of each series' sample variance (full sample); s.d. = `sqrt(SW_ME_SHARE)` × sample s.d. | 0.10 |
  | `SW_ME_SD` | one absolute s.d. for all five series (overrides `SW_ME_SHARE`) | unset |

  With the default the standard deviations are 0.321 (`dy`), 0.266 (`dc`), 0.827 (`dinve`),
  0.197 (`dw`) and 0.917 (`labobs`); they are printed at the start of the run and stored in
  `M_.sw_me`. The 10% share is a calibration choice of this repository, not taken from Smets
  and Wouters (2007).
- **States of interest**: the flexible-price real rate `rrf` ("natural rate") and the output gap
  `ygap = y - yf` (defined as a linear combination in column 4 of `M_.CellList_States`), each with
  fixed-lag smoothers at lags 1, 4 and 20 (`l1rrf`, `l4rrf`, `l20rrf`, ...) and, for `rrf`, its
  4- and 20-quarter changes (`d4rrf`, `d20rrf`).

## Outputs

| | UC3 | SW baseline | SW with measurement error |
|---|---|---|---|
| Run folder | `Code/UC3/` | `Code/SmetsWouters/SW_baseline/` | `Code/SmetsWouters/SW_ME/` |
| Figures | `Graphs_level/` (`LW_KF_*.png`, `MA_*.png`) | `Graphs_level/` (`SW_*.png`, `MA_*.png`) | same |
| Tables (LaTeX) | `Tables_level/` | `Tables_level/` | same |
| Results | `UC3_results.mat` (+ `_prior`, `_posterior`) | `SW_results.mat` (+ `_prior`, `_posterior`) | same |
| Dynare workspace | `tmp.mat` | `tmp.mat` | same |
| Model assessment | `model_assessment_lecam.mat` | `model_assessment_lecam.mat` | same |
| Dynare output, MH chains, log | `uc3_levels/`, `uc3_levels.log` | `Smets_Wouters_2007_45/`, `Smets_Wouters_2007_45.log` | same |

The SW run folder is created by the driver and holds a copy of the `.mod` and the data; set
`SW_RUN_DIR` to use another folder name. Main tables: `table_summary.tex` (state-estimate s.d.,
$R^2_4$, $R^2_{20}$, $\nu$, half-life), `table_corrh_*.tex`, `table_halflife_*.tex`,
`table_implied_R2.tex`, `table_prior_posterior_params.tex`, `table_rmse_correlations.tex`,
`MA_table_*.tex`. Main figures: `*_R2_ptiles_*` (distribution of $R^2_h$ over horizons),
`*_R2H<h>_*` (histograms), `*_HalfLife_*`, `*_<state>.png` (filtered state with 95% band),
`MA_*` (model assessment). All outputs are git-ignored.

## Using another model instead of `uc3_levels.mod`

The pipeline is model-agnostic; the model-specific choices sit in the `.mod` file and in the
driver. To run another linear model:

1. **`.mod` file** — after `estimation`/`calib_smoother` and `stoch_simul`, copy the block of
   `uc3_levels.mod` that starts at `save tmp.mat`: it sets `M_.observable_names`, `M_.data`,
   `M_.CellList_States = {name, smoother lags, change horizons [, definition]}`,
   `M_.states_of_interest`, `M_.var_list_`, `M_.X_InitialDistribution`, the `oo_.bayestopt_*`
   fields and `oo_.dr.qz_crit`, then calls `save_results_split`, `run_reliability_draws('prior', ...)`
   and `run_reliability_draws('posterior', ...)`. Conventions the toolbox requires:
   - shocks: either declare `stderr 1` and put the standard deviations as parameters in the
     equations (UC3), or estimate `stderr` in `estimated_params` (SW); Dynare's measurement-error
     matrix (`varobs` with `stderr` on observables, `M_.H`) is **not** supported — declare
     measurement errors as shocks (as in the SW `.mod`);
   - a state of interest that is not a model variable can be defined as a linear combination in
     column 4, e.g. `{'ygap', [0 1 4 20], [], 'y - yf'}`.
2. **Driver** — copy `Code/UC3/main_run_uc3_levels.m` (or `Code/SmetsWouters/run_sw_variant.m`
   for a run in its own folder) and change: the `.mod` name in the `dynare` call, `CLS_uc3`
   (same list as `M_.CellList_States`), `HLW_RESULTS_FILE` (the name given to
   `save_results_split`), and `HLW_PLOT_CONFIG`: `PlotVar_Names` (states to plot, engine labels
   such as `l20ustar`), `time_start` (first quarter of the data, e.g. `1960` or `1947.5` for
   Q3), `common_scale_base`, `scatter_pairs`, `rmse_vars`, `tracking_vars`, `target_params`
   (parameter names of the new model) and `rstar_overlay = false`.
3. Run each model (or each variant of one `.mod`) in its own folder: Dynare writes a
   `+<modname>` package and an output folder next to the `.mod`.

## License

GNU General Public License v3.0 (`LICENSE`). `Code/SmetsWouters/Smets_Wouters_2007_45.mod` is
adapted from Johannes Pfeifer's DSGE_mod replication file (© 2007–2013 Frank Smets and Raf
Wouters, © 2013–15 Johannes Pfeifer, GPL v3 or later). The Smets–Wouters data are not
redistributed here; `get_sw_data.m` downloads them from the original source.

## Limitations

- **Short SW chains.** The default SW chains (2 × 5,000 draws, 20% burn-in) and 500 reliability
  draws illustrate the workflow; they are not a converged posterior and no convergence
  assessment is claimed. Increase `SW_MH_REPLIC` and `SW_NDRAWS` for substantive work.
- **SW baseline is a degenerate benchmark.** Without measurement error the SW model has seven
  shocks for seven noise-free observables, so $D$ is square and invertible and the shocks can
  be read off the data. The steady-state filter recovers the states exactly: $\Omega = 0$ and
  $R^2_h = 1$.
  - The toolbox detects this case directly (`solve_steady_state_riccati`, method
    `invertible_critical`; see "Square, invertible systems" in
    `Code/ReliabilityToolbox/README.md`). The four differenced observables put four
    eigenvalues of $A - BD^{-1}C$ on the unit circle, so this limit is not a stabilising
    solution: the time-varying filter reaches it only at rate $1/t$, and on the actual sample
    the filtered estimates still carry error.
  - Quantities that are 0/0 at $\Omega = 0$ are undefined: the error autocorrelation, the FEV
    shares, the half-life (shown as "--" in the tables) and the implied-reliability identity
    check (Prop. 5.1), which prints NaN.
- **SW states are stationary.** $\nu$ is not defined ("--" in the tables). $R^2_h$ rises with
  $h$ partly because its denominator, the conditional variance of the future state, rises
  towards the unconditional variance. The "s.d. of state estimate / innovation s.d." rows
  normalise by the impact standard deviation of the state, which is not comparable with the
  random-walk trends of UC3. `rrf` is in quarterly percent, as a deviation from steady state.
- **SW model assessment.** `model_assessment` re-estimates on artificial data with the toolbox's
  Kalman filter, whose likelihood is not exactly Dynare's under the replication settings
  (`presample=4`, `lik_init=2`): the mode condition it prints is 3.2 (measurement error) and 1.6
  (baseline) posterior s.d. rather than 0 (UC3: 0.000). The Newton steps therefore start off the
  toolbox's mode. Matching the initialisation (`lik_init=1`, `presample=0`) did not remove the
  gap and made the mode search fail.
