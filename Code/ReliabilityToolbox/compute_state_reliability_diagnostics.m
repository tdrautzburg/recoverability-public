function [LatentStateDiagnostics, ModelMatrices, obs_indices, key_state_indices, TV_state_estimate, TV_surprises, TV_log_llhood] = compute_state_reliability_diagnostics(observable_names, states_of_interest, oo_, M_, data, true_states)
%COMPUTE_STATE_RELIABILITY_DIAGNOSTICS  Build the ABCD state space from a solved Dynare
%   model and compute the latent-state reliability diagnostics of the paper.
%
%   RENAMED (2026-07): formerly f_Extract_ABCD_par_while.m in Miscellaneous/. All call
%   sites inside ReliabilityToolbox/ use the new names; the full mapping
%   is in README.md ("Suggested file names").
%
%   MERGED (2026-09-15): this file is the union of the two copies that had
%   diverged since 2026-08 - HLW_new/Miscellaneous (true_states input, eq.-(3.17)
%   R2_h via fixed_lag_smoother_map, tail index nu, Riccati flag, one data pass)
%   and ReliabilityToolbox (M_.fast_mode / M_.loglik_only switches, obs_ss and
%   X0_pruned/B0_pruned on ModelMatrices, the TV data pass factored out into
%   apply_tv_filter). Numerics of both are preserved: with the switches off the
%   result equals the Miscellaneous engine; with them on it equals the
%   ReliabilityToolbox engine on every field both produced. See README.md.
%
%   Pipeline:
%     1. Extract the ABCD representation from the Dynare decision rules:
%            x(t) = A x(t-1) + B e(t)
%            y(t) = C x(t-1) + D e(t),   e(t) ~ N(0, I)
%        (observables load on LAGGED states because Dynare's decision rule
%        writes every variable as ghx*s(t-1) + ghu*e(t));
%     2. augment the state vector with any "state of interest" that is not
%        a Dynare state variable, and prune states that neither feed the
%        observables nor any surviving state;
%     3. run the exact diffuse + time-varying Kalman filter
%        (kalman_filter_steady_and_diffuse), which also returns the
%        steady-state gain and covariances from its own Riccati iteration;
%        (until 2026-08 steady_state_kalman_filter was additionally called
%        here, but its outputs were overwritten unused - see section 4);
%     4. compute population reliability diagnostics:
%          - CorrH: corr(x(t), xhat(t)) as information accumulates over h
%            periods (the "reliability R" of the paper; CorrH.^2 is the
%            conditional R2 of eq. (3.11)). For a fixed-lag SMOOTHED state
%            'l<m><NAME>' it is, since 2026-09-15, the conditional
%            correlation corr(x_{i,t+h}, xhat_{i,t+h|t+h+m} | y^t) =
%            sqrt(R2_h) - the regression analogue of eq. (3.17), normalised
%            by the contemporaneous state's variance (step 5c);
%          - R2_h: the reliability under the eq.-(3.17) normalisation, which
%            is the one that is comparable ACROSS INFORMATION SETS and
%            therefore the one to report when a table puts filtered and
%            fixed-lag SMOOTHED estimates side by side (ADDED 2026-08; see
%            step 5c and fixed_lag_smoother_map). Equal to CorrH.^2 for
%            every row;
%          - Error_ACF: autocorrelation of the filtering error at lag h;
%          - FEV: fraction of the error variance resolved as h grows;
%          - half-lives of estimation errors, overall and by shock;
%          - nu: the horizon-free tail-unreliability index of the
%            non-stationary states (tail_unreliability_index), ADDED 2026-08,
%            REVISED 2026-09 to the general-integration-order formula, and
%            .nu_order, the integration order that fixes its units;
%          - theoretical R2 wedge objects (theoretical_R2_wedge);
%     5. if data is supplied: filter the data with the exact diffuse /
%        time-varying filter and compute empirical forecast-error RMSEs and
%        error/revision correlations by horizon
%        (forecast_errors_and_revisions), used for the DeltaR2 comparison
%        with the data.
%
%   FAST MODES (MERGED 2026-09 from ReliabilityToolbox):
%     M_.fast_mode   = true skips every population diagnostic that the
%                      empirical DeltaR2 comparison does not need (FEV, CorrH,
%                      Error_ACF, half-lives, R2_h, nu, the wedges fh/fh_direct)
%                      and runs the 5c covariance recursion only up to MaxH_Y,
%                      which is all CorrH_Observables needs. The data pass and
%                      CorrH_Observables are computed in full.
%     M_.loglik_only = true is strictly stronger: sections 5 and 6 are skipped
%                      entirely and only the state space, the Riccati solution
%                      and the diffuse/TV Sequence are built - what
%                      posterior_predictive_wrapper needs to evaluate
%                      log p(y|theta) on many data sets with apply_tv_filter.
%     Anything skipped is returned as NaN (or [] where the full engine returns
%     []), never as 0, so an uncomputed object cannot be mistaken for a result.
%
%   CHANGED (2026-08): the parallel STEADY-STATE-GAIN data pass (the old
%   step 6a) has been REMOVED. The steady-state Riccati solution itself is
%   untouched and still drives every population diagnostic of step 4 - it is
%   the definition of the paper's K, Omega, A-KC objects - but the data are
%   now filtered once, by the diffuse/TV filter only. Rationale:
%     - the diffuse filter's empirical statistics are EXACTLY invariant to
%       the level of the non-stationary states, the constant-gain filter's
%       are not (it starts at X0 = 0 while the trends sit at the level of
%       the data, e.g. 100*log(GDP) ~ 816 in uc3_level_i2, and that transient
%       decays at max|eig(A-KC)|, i.e. at the reported estimation-error
%       half-life). That invariance is what makes the simulation-based null
%       of model_assessment.m - whose replications draw random trend levels -
%       the right reference distribution for the empirical statistic;
%     - it is the filter Dynare itself used (lik_init=3), so these statistics
%       are consistent with the reported states, likelihood and posterior;
%     - it is the only filter that is defined when shocks are heteroskedastic
%       (Qscale ~= 1), as in the HLW/LW applications;
%     - having a single filter removes the sample/burn-in asymmetry that made
%       the two sets of statistics incomparable (see step 6).
%
%   INPUTS
%     observable_names   cell array, names of the observables (must be
%                        endogenous variables of the model)
%     states_of_interest EITHER (legacy) a cellstr of names of the latent
%                        states to diagnose, e.g. {'ygap';'rrstar'},
%                        OR (ADDED 2026-09-15) the J x 3 cell LIST kept in
%                        M_.CellList_States
%
%                          { 'ygap',    [0 1 4 20], []
%                            'rrstar',  [0 1 4 20], [4 20]
%                            'drrstar', [0 1],      []    }
%
%                        In the .mod file set both
%                          M_.CellList_States    = { ... };
%                          M_.states_of_interest = expand_states_of_interest(M_.CellList_States);
%                        so that M_.states_of_interest holds the names of the
%                        diagnosed states in the row order of every output
%                        (backward compatible with the plotting scripts);
%                        reliability_draw_wrapper hands M_.CellList_States to
%                        this function when the field exists.
%
%                        column 1: name of a model variable (state or not);
%                        column 2: fixed-lag SMOOTHER lags m to diagnose,
%                          0 = the filtered state itself, m >= 1 = the
%                          m-period smoothed estimate, realised as the
%                          filtered lag-augmented state 'l<m><NAME>'
%                          (l20rrstar, as in the old .mod files);
%                        column 3: change horizons h, diagnosed at m = 0
%                          only, as the augmented state 'd<h><NAME>' =
%                          NAME_t - NAME_{t-h} (d20rrstar).
%                        A J x 2 spec is read as lags only. The lag and
%                        change states are built by augment_ssm (step 2b)
%                        INSTEAD of being declared in the .mod file, and
%                        the eq.-(3.17) smoother map is then known by
%                        construction (no structural detection). The
%                        diagnosed states, in order, are
%                        LatentStateDiagnostics.state_labels; the legacy
%                        form gives the old behaviour exactly.
%     oo_                Dynare results structure; REQUIRED extra fields
%                        (set in the .mod file):
%                          oo_.bayestopt_smoother_var_list
%                          oo_.bayestopt_smoother_restrict_columns
%                          oo_.bayestopt_mf
%                          oo_.dr.qz_crit   (e.g. 1 + 1e-6)
%     M_                 Dynare model structure; REQUIRED extra fields:
%                          M_.X_InitialDistribution.X0/.B0 (and/or
%                          M_.heteroskedastic_shocks.Qscale)
%                        optional: M_.X_InitialDistribution.type (ADDED
%                        2026-09-26), 'diffuse' (default) or 'proper':
%                        'diffuse' initialises the unit roots exactly
%                        diffusely and the stationary block at its
%                        unconditional covariance (X0 is the initial mean of
%                        the data pass, B0 only seeds the Riccati solve);
%                        'proper' starts the time-varying filter from the
%                        finite prior x_0 ~ N(X0, B0*B0') with no diffuse
%                        phase, matching a sampler initialised with a wide
%                        finite variance (Bauer-Rudebusch ese.R). The
%                        population objects do not depend on this choice
%                        (see kalman_filter_steady_and_diffuse)
%                        optional: M_.solve_function - not read here, but by
%                        the per-draw routines through solve_model_draw
%                        optional: M_.verbose (default true)
%                        optional: M_.fast_mode, M_.loglik_only (default
%                        false), see FAST MODES above
%                        optional: M_.reliability_burn_in (default 50), the
%                        number of initial periods excluded from EVERY
%                        empirical statistic (see step 6)
%                        optional: M_.fixed_lag_smoother_map, a k x 3 cell
%                        array {state of interest, contemporaneous base
%                        state, m} that pins down which states of interest
%                        are fixed-lag SMOOTHED estimates instead of letting
%                        step 5c detect it from the impulse responses (or,
%                        with a spec, instead of reading it off the spec).
%                        Use cell(0,3) to switch the eq.-(3.17) correction off.
%                        optional: M_.fr_smoother_lags (default [4 20]), the
%                        smoother lags m of the smoother-revision regression
%                        (below); [] switches the regression off. Spec
%                        form only.
%     data               (optional) T x n_obs data matrix. When present,
%                        the outputs 5-7 are filled (unless M_.loglik_only).
%     true_states        (optional) T x n_states matrix of the TRUE state
%                        paths, columns ordered like the PRUNED state vector
%                        (i.e. like all_states / ModelMatrices.A). For
%                        simulation experiments where the states are known.
%                        When supplied, the engine additionally returns
%                        state-level h-step forecast RMSEs and error/revision
%                        correlations against the truth
%                        (.TV_RMSE_states, .TV_CorrY_H_states); those fields
%                        are NaN otherwise. Requires `data`.
%
%   OUTPUTS
%     LatentStateDiagnostics  struct of reliability statistics (see below).
%                             Includes .R2_h (ADDED 2026-08): the eq.-(3.17)
%                             reliability, n_key_states x MaxH, together with
%                             .smoother_lag_m / .smoother_base / .smoother_map_info.
%                             Includes .nu (ADDED 2026-08, REVISED 2026-09):
%                             the horizon-free tail-unreliability index, one
%                             entry per state of interest, NaN for states
%                             without a stochastic trend, and .nu_order (ADDED
%                             2026-09), the order d of each entry. nu is in
%                             quarters^(2d+1), so .nu_order is what makes the
%                             units readable and is required before entries
%                             with different d are compared. See step 5e.
%                             Includes .state_labels (ADDED 2026-09-15): the
%                             names of the diagnosed states (base names,
%                             'l<m><NAME>', 'd<h><NAME>'), one per row of
%                             CorrH / R2_h / nu, and .state_spec, the parsed
%                             spec. Includes .fr (ADDED 2026-09-15): the
%                             SMOOTHER-REVISION REGRESSION
%                                 f_{i,h,m} = alpha + beta r_{i,h},
%                                 r_{i,h}   = xhat_{i,t+h|t+h}   - xhat_{i,t+h|t},
%                                 f_{i,h,m} = xhat_{i,t+h|t+h+m} - xhat_{i,t+h|t},
%                             xhat_{t+h|t} = A^h xhat_{t|t}, for every base
%                             state of the spec, h = 1..MaxH_Y and m in
%                             M_.fr_smoother_lags: population R2 / beta
%                             under the steady-state filter (.R2_theory,
%                             .beta_theory; fr_regression_theory) and OLS
%                             R2 / alpha / beta on the diffuse-TV-filtered
%                             data (.R2_emp, .alpha_emp, .beta_emp, .N_emp;
%                             fr_regression_empirical). NaN when not
%                             computed (legacy input, fast/loglik modes as
%                             documented in those files, no data).
%     ModelMatrices           struct with A,B,C,D, gains and covariances, the
%                             observables' steady state obs_ss and the pruned
%                             initial distribution X0_pruned/B0_pruned (what
%                             simulate_artificial_data needs)
%     obs_indices             DR-order indices of the observables
%     key_state_indices       positions of states_of_interest in the pruned
%                             state vector
%     TV_state_estimate       T x n_key filtered states, diffuse/TV filter
%     TV_surprises            T x n_obs forecast errors, TV filter
%     TV_log_llhood           log-likelihood, TV filter (NaN if not computed)
%                             (the steady-state-gain counterparts of these
%                             three, outputs 5-7 until 2026-08, are gone)
%
%   MERGED (2026-08): compute_state_reliability_vs_truth.m, a 92%-identical
%   copy of this file, has been folded in here and DELETED. Two of its three
%   differences are now part of this engine:
%     - the optional `true_states` input and the state-level statistics it
%       enables (.TV_RMSE_states / .TV_CorrY_H_states, named .steady_*_states
%       in the deleted file);
%     - the eq.-(26) wedge from R2_wedge_factor, returned as .fh_direct. NOTE
%       THE NAMING: the deleted engine called that object .fh and put the
%       theoretical_R2_wedge ("newdeltaR2") version in .fh_shift; here .fh
%       keeps ITS existing meaning (theoretical_R2_wedge) and the direct
%       factor is the new field. Code written against the old file must be
%       re-read: .fh there == .fh_direct here.
%   The third difference is deliberately NOT carried over: the deleted file
%   forced the time-varying filter to be homoskedastic (Qscale = 1 regardless
%   of M_.heteroskedastic_shocks). This engine honours M_.heteroskedastic_shocks
%   as it always has. For a model without a heteroskedastic_shocks block
%   (uc3_level_i2) the two are identical, since Qscale defaults to ones; they
%   differ only for models that actually rescale shock variances (HLW/LW).
%
%   ADDED (2026-08): eq. (3.17), the filter-vs-smoother comparable R2.
%   Section 3.5 of the draft requires reliability to be weakly increasing in
%   the amount of hindsight m, but Table 3 reported CorrH.^2 - eq. (3.11) -
%   for the m = 20 smoother too, and eq. (3.11) applied to the augmented
%   lagged copy normalises by the conditional variance of the EARLIER state.
%   That is the "each estimate generates its own denominator" case Section 3.5
%   warns about, and it is why the smoothed R2 came out below the filtered
%   one. Step 5c now also returns .R2_h, which keeps the numerator (the
%   smoothing error variance [Omega_m]_ii, already on the diagonal of
%   KF_Omega) but takes the denominator from the contemporaneous base state.
%   See fixed_lag_smoother_map for how the base is identified.
%
%   See also: reliability_draw_wrapper, steady_state_kalman_filter,
%             kalman_filter_steady_and_diffuse, apply_tv_filter,
%             estimation_error_half_life, forecast_errors_and_revisions,
%             theoretical_R2_wedge, tail_unreliability_index,
%             fixed_lag_smoother_map, posterior_predictive_wrapper.

%% ------------------------------------------------------------------
%  0. Input handling
% -------------------------------------------------------------------
if ~isfield(M_, 'verbose')
    M_.verbose = true;
end

% Fast modes (MERGED 2026-09): fast_mode skips the population diagnostics
% that do not feed the empirical DeltaR2 comparison; loglik_only is strictly
% stronger and skips sections 5 and 6 in full, leaving only what is needed to
% evaluate the likelihood (the Riccati solution and the diffuse/TV Sequence).
% Skipped outputs are NaN or [], never 0.
if ~isfield(M_, 'fast_mode')
    M_.fast_mode = false;
end
if ~isfield(M_, 'loglik_only')
    M_.loglik_only = false;
end
LOGLIK_ONLY = M_.loglik_only;
FAST        = M_.fast_mode || LOGLIK_ONLY;

HAS_DATA          = exist('data', 'var') && ~isempty(data);
DO_ESTIMATE_STATE = HAS_DATA && ~LOGLIK_ONLY;

% MERGED (2026-08): the optional 6th input arrives from the former
% compute_state_reliability_vs_truth (see the header). When absent the engine
% behaves exactly as before.
HAVE_TRUE_STATES  = exist('true_states', 'var') && ~isempty(true_states);

MaxH   = 250;              % maximum horizon for the population correlation diagnostics
MaxH_Y = min(40, MaxH);    % horizon of the observable-level (data-comparable) objects

%% ------------------------------------------------------------------
%  1. Locate states and observables in the decision rules
% -------------------------------------------------------------------
% oo_.dr.state_var: indices (declaration order) of the state variables for
% the current parameterization. inv_order_var maps declaration order into
% DR (decision-rule) order, i.e. into the rows of ghx/ghu.
state_indices     = oo_.dr.inv_order_var(oo_.dr.state_var);
n_states_original = length(state_indices);
shock_names       = M_.exo_names;
n_shocks          = length(shock_names);
state_var_names   = M_.endo_names(oo_.dr.state_var);

% Full policy matrices: every variable = ghx*s(t-1) + ghu*e(t)
ghx_full = oo_.dr.ghx;
ghu_full = oo_.dr.ghu;

% SHOCK COVARIANCE (ADDED 2026-09-30). The toolbox timing has e_t ~ N(0, I), while Dynare's ghu loads
% on shocks with covariance M_.Sigma_e (the shocks block / estimated stderr). The toolbox's own .mod
% files declare unit-variance shocks and carry every s.d. as a parameter inside the equations, so
% Sigma_e = I there and nothing changes; a model written the Dynare way (SW07: 'stderr ea, ...' in
% estimated_params, measurement errors with their own stderr) had its s.d. silently ignored until
% now. ghu is scaled by a square root of Sigma_e: sqrt of the diagonal when Sigma_e is diagonal
% (keeps the per-shock attribution of the FEV), the lower Cholesky factor otherwise (symmetric
% square root if only PSD). Dynare's heteroskedastic scales (M_.heteroskedastic_shocks) multiply
% the s.d., so the Qscale treatment below is consistent with this. Dynare's measurement-error
% matrix M_.H (stderr on observed variables) is NOT supported: declare measurement errors as shocks.
if isfield(M_, 'Sigma_e') && ~isempty(M_.Sigma_e) && size(M_.Sigma_e, 1) == n_shocks
    Se = (M_.Sigma_e + M_.Sigma_e') / 2;
    if any(abs(Se(:) - reshape(eye(n_shocks), [], 1)) > 0)
        if isdiag(Se)
            S_e = diag(sqrt(max(diag(Se), 0)));
        else
            [L_e, flag_e] = chol(Se, 'lower');
            if flag_e == 0, S_e = L_e; else, [V_e, E_e] = eig(Se); S_e = V_e * diag(sqrt(max(diag(E_e), 0))) * V_e'; end
        end
        ghu_full = ghu_full * S_e;
    end
end
if isfield(M_, 'H') && ~isempty(M_.H) && any(M_.H(:) ~= 0)
    warning('compute_state_reliability_diagnostics:measurementError', ...
        'M_.H (Dynare measurement errors on observed variables) is ignored by the toolbox; declare measurement errors as exogenous shocks in the .mod file.');
end

% State transition: rows of ghx/ghu belonging to the states themselves
A_original = ghx_full(state_indices, :);   % square: n_states x n_states
B_original = ghu_full(state_indices, :);

% Observables
n_obs = length(observable_names);
obs_indices = zeros(n_obs, 1);
for i = 1:n_obs
    obs_indices(i) = oo_.dr.inv_order_var(find(strcmp(M_.endo_names, observable_names{i}))); %#ok<FNDSB>
end
C_original = ghx_full(obs_indices, :);
D_original = ghu_full(obs_indices, :);

% Steady state of the observables. Computed here rather than in section 6
% (where it used to live) so that it is available to every caller, including
% loglik_only ones and simulate_artificial_data, which adds it back to the
% simulated observables (MERGED 2026-09 from ReliabilityToolbox).
obs_ss = zeros(n_obs, 1);
for i = 1:n_obs
    var_idx   = find(strcmp(M_.endo_names, observable_names{i}), 1);
    obs_ss(i) = oo_.dr.ys(var_idx);
end

%% ------------------------------------------------------------------
%  2. Parse the states of interest; augment the state vector with
%     non-state variables of interest
% -------------------------------------------------------------------
% ADDED (2026-09-15): states_of_interest may be a J x 3 list (M_.CellList_States,
% see header). expand_states_of_interest - the same function the .mod file
% uses to fill M_.states_of_interest for the plotting scripts - returns the
% base variables, the smoother lags and the change horizons per base, and
% the expanded list of diagnosed states (key_names) with, for each, the
% smoother lag m and the base state whose conditional variance is the
% eq.-(3.17) denominator. The legacy cellstr input parses to lags {0}, no
% changes, no regression.
[key_names, sp] = expand_states_of_interest(states_of_interest);
spec_base = sp.base;  spec_lags = sp.lags;  spec_changes = sp.changes;  is_spec = sp.is_spec;
key_lag_m = sp.lag_m;  key_base = sp.base_of;
J_base    = numel(spec_base);
% DEFINED STATES (ADDED 2026-09-30): a base state may be a linear combination of endogenous
% variables given in column 4 of the list (e.g. {'ygap', [0], [], 'y - yf'}). The definitions
% come with the spec; a caller handing in the legacy name list (the cheap likelihood systems of
% make_loglik_ctx / model_assessment / find_posterior_mode pass M_.CellList_States(:, 1)) gets
% them from M_.CellList_States itself, so a defined state never has to exist in the .mod.
spec_defs = sp.defs;
if isfield(M_, 'CellList_States') && iscell(M_.CellList_States) && size(M_.CellList_States, 2) == 4
    for i = 1:J_base
        if isempty(spec_defs{i})
            k = find(strcmp(M_.CellList_States(:, 1), spec_base{i}), 1);
            if ~isempty(k) && ~isempty(M_.CellList_States{k, 4}), spec_defs{i} = M_.CellList_States{k, 4}; end
        end
    end
end

% Smoother lags of the smoother-revision regression (spec form only)
if is_spec
    if isfield(M_, 'fr_smoother_lags')
        fr_m = M_.fr_smoother_lags;
    else
        fr_m = [4 20];
    end
    fr_m = unique(double(fr_m(:)'));
    if any(~isfinite(fr_m)) || any(mod(fr_m, 1) ~= 0) || any(fr_m < 1)
        error('compute_state_reliability_diagnostics:badFrLags', ...
            'M_.fr_smoother_lags must contain finite integers >= 1 (or be empty).');
    end
else
    fr_m = zeros(1, 0);
end

n_key_states = numel(key_names);

% 2a. A base variable that is not a Dynare state (e.g. a static
% transformation of states, rrstar = c*gg + zz) is appended to the state
% vector with dynamics
%   z(t) = ghx_z * s(t-1) + ghu_z * e(t),
% i.e. a zero column block in the augmented A (nothing feeds back from z).
%   A DEFINED state (column 4 of the list, ADDED 2026-09-30) is appended the same way with
%   the rows coef' * ghx_full / coef' * ghu_full of the variables in its definition: exact
%   for any parameter draw, and its lags / changes follow through augment_ssm as usual.
%   Its name must not be a state variable of the model (ambiguous); a non-state
%   endogenous variable of the same name is shadowed by the definition.
non_state_vars_of_interest = {};
ghx_nonstate = zeros(0, size(ghx_full, 2));
ghu_nonstate = zeros(0, size(ghu_full, 2));

for i = 1:J_base
    if ~isempty(spec_defs{i})
        if any(strcmp(state_var_names, spec_base{i}))
            error('compute_state_reliability_diagnostics:definedState', ...
                'state of interest "%s" has a definition ("%s") but is also a state variable of the model.', spec_base{i}, spec_defs{i});
        end
        [cf, vars] = parse_state_definition(spec_defs{i}, M_.endo_names);
        rows = zeros(numel(vars), 1);
        for k = 1:numel(vars), rows(k) = oo_.dr.inv_order_var(find(strcmp(M_.endo_names, vars{k}), 1)); end
        non_state_vars_of_interest{end+1} = spec_base{i};                       %#ok<AGROW>
        ghx_nonstate(end+1, :) = cf' * ghx_full(rows, :);                       %#ok<AGROW>
        ghu_nonstate(end+1, :) = cf' * ghu_full(rows, :);                       %#ok<AGROW>
    elseif ~any(strcmp(state_var_names, spec_base{i}))
        endo_pos = find(strcmp(M_.endo_names, spec_base{i}), 1);
        if isempty(endo_pos)
            error('compute_state_reliability_diagnostics:unknownState', ...
                ['state of interest "%s" is not an endogenous variable of the model ' ...
                 '(define it in column 4 of M_.CellList_States if it is a linear combination of variables).'], spec_base{i});
        end
        non_state_vars_of_interest{end+1} = spec_base{i};                       %#ok<AGROW>
        ghx_nonstate(end+1, :) = ghx_full(oo_.dr.inv_order_var(endo_pos), :);   %#ok<AGROW>
        ghu_nonstate(end+1, :) = ghu_full(oo_.dr.inv_order_var(endo_pos), :);   %#ok<AGROW>
    end
end

n_non_state = length(non_state_vars_of_interest);
if n_non_state > 0
    if M_.verbose
        fprintf('Found %d variables of interest that are not state variables:\n', n_non_state);
        fprintf('  %s\n', non_state_vars_of_interest{:});
    end

    A_augmented = [A_original,                          zeros(n_states_original, n_non_state);
                   ghx_nonstate(:, 1:n_states_original), zeros(n_non_state, n_non_state)];
    B_augmented = [B_original; ghu_nonstate];

    % Extend the initial distribution (if provided) to the augmented state.
    % NOTE: the lower-left block reuses the loadings ghx_nonstate as the
    % covariance FACTOR; this is only an approximation of the implied
    % initial cross-covariance (kept from the original implementation; the
    % augmented branch is not exercised by uc3_level.mod, whose states of
    % interest are all genuine state variables).
    if isfield(M_, 'X_InitialDistribution')
        if isfield(M_.X_InitialDistribution, 'B0')
            M_.X_InitialDistribution.B0 = [M_.X_InitialDistribution.B0,            zeros(n_states_original, n_non_state);
                                           ghx_nonstate(:, 1:n_states_original),    eye(n_non_state)];
        end
        if isfield(M_.X_InitialDistribution, 'X0')
            M_.X_InitialDistribution.X0 = [M_.X_InitialDistribution.X0; zeros(n_non_state, 1)];
        end
    end

    state_var_names = [state_var_names; non_state_vars_of_interest'];

    C_augmented = [ghx_full(obs_indices, :), zeros(n_obs, n_non_state)];
    D_augmented = ghu_full(obs_indices, :);

    % CHANGED (2026-09-15): the line
    %     A_augmented = round(A_augmented*1e6)/1e6;
    % that used to sit here (to kill numerical dust before the exact-zero
    % redundancy tests of step 3) rounded the DYNAMICS, not just the zero
    % pattern. Measured on HLW it moved the filtered states by up to 2.9e-4
    % relative to Dynare's own updated states (oo_.UpdatedVariables) and was
    % the whole of that discrepancy: with A exact the two agree to 3e-13.
    % The rounding now applies to a pattern copy used by step 3 only (see
    % A_pattern below); A_augmented stays exactly Dynare's decision rule.
    n_states_augmented = n_states_original + n_non_state;
else
    A_augmented = A_original;
    B_augmented = B_original;
    C_augmented = ghx_full(obs_indices, :);
    D_augmented = ghu_full(obs_indices, :);
    n_states_augmented = n_states_original;
end

% 2b. ADDED (2026-09-15): lag and change states via augment_ssm, replacing
% the lagged copies that used to be declared in the .mod file - the labels
% are the same (l1rrstar..l20rrstar, d4rrstar, ...). Every base state gets
% the lags the spec asks for AND the lags the smoother-revision regression
% needs (fr_m); augment_ssm fills in the intermediate lags itself. The
% appended states are exact linear identities in the history of x, so:
%   - the initial-distribution factor B0 is extended by copying the base
%     state's row for a lag (l<k>NAME at t = 0 is a lagged value of the same
%     variable, so it starts with the base's spread and is perfectly
%     correlated with it) and by a zero row for a change. B0 only seeds the
%     Riccati iteration of the steady-state filter, whose fixed point does
%     not depend on it; the diffuse/TV filter builds its own exact
%     initialization from the augmented system inside
%     durbin_koopman_diffuse_filter (step 4 there);
%   - X0 is extended likewise (base value for lags, 0 for changes); it is
%     the initial mean of the data pass only.
if is_spec
    aug_spec = cell(J_base, 3);
    for j = 1:J_base
        aug_spec{j, 1} = spec_base{j};
        aug_spec{j, 2} = reshape(unique([spec_lags{j}(spec_lags{j} > 0), fr_m]), 1, []);
        aug_spec{j, 3} = spec_changes{j};
    end
    n_before = n_states_augmented;
    [A_augmented, B_augmented, C_augmented, D_augmented, state_var_names] = ...
        augment_ssm(A_augmented, B_augmented, C_augmented, D_augmented, state_var_names, aug_spec);
    state_var_names    = state_var_names(:);
    n_states_augmented = size(A_augmented, 1);
    q_app = n_states_augmented - n_before;

    if q_app > 0 && isfield(M_, 'X_InitialDistribution')
        Rsel = zeros(q_app, n_before);       % row selector: base state of each lag, 0 for changes
        for a = 1:q_app
            lab = state_var_names{n_before + a};
            tok = regexp(lab, '^l(\d+)([A-Za-z_].*)$', 'tokens', 'once');   % l<k><NAME>; names never start with a digit
            if ~isempty(tok)
                Rsel(a, find(strcmp(state_var_names(1:n_before), tok{2}), 1)) = 1;
            end
        end
        if isfield(M_.X_InitialDistribution, 'B0')
            B0 = M_.X_InitialDistribution.B0;
            M_.X_InitialDistribution.B0 = [B0,        zeros(n_before, q_app);
                                           Rsel * B0, zeros(q_app,    q_app)];
        end
        if isfield(M_.X_InitialDistribution, 'X0')
            X0 = M_.X_InitialDistribution.X0(:);
            M_.X_InitialDistribution.X0 = [X0; Rsel * X0];
        end
    end
    if M_.verbose
        fprintf('augment_ssm appended %d lag/change state(s)\n', q_app);
    end
end

% Positions of the diagnosed states in the augmented (pre-pruning) vector
key_state_indices_augmented = zeros(n_key_states, 1);
for i = 1:n_key_states
    pos = find(strcmp(state_var_names, key_names{i}), 1);
    if isempty(pos)
        error('compute_state_reliability_diagnostics:missingKey', ...
            'diagnosed state "%s" was not found in the augmented state vector.', key_names{i});
    end
    key_state_indices_augmented(i) = pos;
end

% States the smoother-revision regression needs (base state and its Lag_m
% for every m in fr_m); protected from pruning like the diagnosed states.
fr_names_needed = {};
if is_spec && ~isempty(fr_m)
    for j = 1:J_base
        fr_names_needed{end+1} = spec_base{j}; %#ok<AGROW>
        for m = fr_m
            fr_names_needed{end+1} = sprintf('l%d%s', m, spec_base{j}); %#ok<AGROW>
        end
    end
end
protected = false(n_states_augmented, 1);
protected(key_state_indices_augmented) = true;
for i = 1:numel(fr_names_needed)
    protected(find(strcmp(state_var_names, fr_names_needed{i}), 1)) = true;
end

%% ------------------------------------------------------------------
%  3. Prune redundant states
% -------------------------------------------------------------------
% A state is redundant iff (i) it is not a state of interest (nor needed by
% the smoother-revision regression), (ii) it does not enter any observable
% (zero column of C), and (iii) no other state depends on it (zero
% off-diagonal column of A). Removing such states keeps the filter exactly
% equivalent while shrinking the system.
keep_states = true(n_states_augmented, 1);
redundant_states = [];

% Zero pattern for the redundancy tests. When a non-state variable of
% interest was appended in step 2a its ghx row carries numerical dust, so the
% pattern is rounded to 1e-6 - exactly the rounding that used to be applied
% to A_augmented itself (same pruning decisions as before, exact dynamics).
if n_non_state > 0
    A_pattern = round(A_augmented*1e6)/1e6;
else
    A_pattern = A_augmented;
end

for i = 1:n_states_augmented
    if ~protected(i)
        influences_observables = any(C_augmented(:, i) ~= 0);
        other_states_depend = false;
        for j = 1:n_states_augmented
            if j ~= i && A_pattern(j, i) ~= 0
                other_states_depend = true;
                break;
            end
        end
        if ~influences_observables && ~other_states_depend
            keep_states(i) = false;
            redundant_states = [redundant_states, i]; %#ok<AGROW>
            if M_.verbose
                fprintf('Found redundant state: %s (index %d)\n', state_var_names{i}, i);
            end
        end
    end
end

A = A_augmented(keep_states, keep_states);
B = B_augmented(keep_states, :);
C = C_augmented(:, keep_states);
D = D_augmented;

% Prune the initial distribution accordingly
if isfield(M_, 'X_InitialDistribution')
    if isfield(M_.X_InitialDistribution, 'B0')
        M_.X_InitialDistribution.B0 = M_.X_InitialDistribution.B0(keep_states, keep_states);
    end
    if isfield(M_.X_InitialDistribution, 'X0')
        M_.X_InitialDistribution.X0 = M_.X_InitialDistribution.X0(keep_states, 1);
    end
end

n_states = sum(keep_states);
state_var_names = state_var_names(keep_states);

% Re-index the states of interest within the pruned state vector
key_state_indices = zeros(n_key_states, 1);
current_idx = 0;
for i = 1:n_states_augmented
    if keep_states(i)
        current_idx = current_idx + 1;
        key_state_indices(key_state_indices_augmented == i) = current_idx;
    end
end

% Rows of the pruned state vector used by the smoother-revision regression
fr_base_idx = zeros(J_base, 1);
fr_lag_idx  = zeros(J_base, numel(fr_m));
if is_spec && ~isempty(fr_m)
    for j = 1:J_base
        fr_base_idx(j) = find(strcmp(state_var_names, spec_base{j}), 1);
        for im = 1:numel(fr_m)
            fr_lag_idx(j, im) = find(strcmp(state_var_names, sprintf('l%d%s', fr_m(im), spec_base{j})), 1);
        end
    end
end

if M_.verbose
    fprintf('Removed %d redundant state(s)\n', length(redundant_states));
    disp('Remaining states:');
    for i = 1:n_states
        disp(['Row ' num2str(i) ': ' state_var_names{i}]);
    end
    disp('Variables of interest and their indices in the state vector:');
    for i = 1:n_key_states
        disp(['  ' key_names{i} ': ' num2str(key_state_indices(i))]);
    end
    disp(['A matrix dimensions: ' num2str(size(A,1)) 'x' num2str(size(A,2))]);
    disp(['B matrix dimensions: ' num2str(size(B,1)) 'x' num2str(size(B,2))]);
    disp(['C matrix dimensions: ' num2str(size(C,1)) 'x' num2str(size(C,2))]);
    disp(['D matrix dimensions: ' num2str(size(D,1)) 'x' num2str(size(D,2))]);
end

%% ------------------------------------------------------------------
%  4. Kalman filters
% -------------------------------------------------------------------
% REMOVED (2026-08): a preceding standalone call to steady_state_kalman_filter
% (the old "step 4a"). It was dead code: kalman_filter_steady_and_diffuse below
% runs its own Riccati iteration and overwrote all four of its outputs
% (K_gain, KF_Omega, KF_V, KF_Omega_tp1_t) before any of them was read, and the
% call below is unconditional. Consequence to be aware of: the steady-state
% covariance used from here on is the one computed inside
% kalman_filter_steady_and_diffuse, which iterates from B0*B0' (since 2026-09
% to a RELATIVE tolerance of 1e-10, see that file), not from 0 to 1e-12. On the
% UC3 calibration the two solutions agree to ~2e-6 in relative terms; since
% nu = (2d+1)*Omega_ii / L_i^(d) is linear in Omega_ii at every order d, that
% is also the accuracy of nu (nu = 3.116279 vs 3.116289 for ustar).
%
% Exact diffuse + time-varying filter. Requires the initial distribution
% and/or heteroskedasticity information on M_ plus the smoother bookkeeping
% fields on oo_ (see header).
if ~isfield(M_, 'X_InitialDistribution') && ~isfield(M_, 'heteroskedastic_shocks')
    error('compute_state_reliability_diagnostics:missingInit', ...
        ['M_ must carry M_.X_InitialDistribution (X0/B0) and/or ' ...
         'M_.heteroskedastic_shocks.Qscale: the time-varying diffuse filter ' ...
         'needs them (set them in the .mod file, as in uc3_level.mod).']);
end

if ~isfield(M_, 'heteroskedastic_shocks')
    % Homoskedastic default: unit variance scales. The Sequence must cover
    % the whole sample even when loglik_only skips the data pass here, since
    % the caller then applies it to the data itself (HAS_DATA, not
    % DO_ESTIMATE_STATE).
    if HAS_DATA
        T = size(data, 1);
    else
        T = 250;   % population-only call: any horizon long enough to pass the diffuse phase
    end
    M_.heteroskedastic_shocks.Qscale = ones(n_shocks, T);
    Qscale = M_.heteroskedastic_shocks.Qscale;
else
    try
        Qscale = M_.heteroskedastic_shocks.Qscale;
    catch
        Qscale = M_.heteroskedastic_shocks.Qscale_orig;
    end
    if isempty(Qscale)
        Qscale = ones(n_shocks, size(data, 1));
    end
    T = size(Qscale, 2);
end

% Hand the smoother bookkeeping to the diffuse filter through oo_.dr
oo_.dr.restrict_var_list = oo_.bayestopt_smoother_var_list;
oo_.dr.restrict_columns  = oo_.bayestopt_smoother_restrict_columns;
oo_.dr.mf                = oo_.bayestopt_mf;
oo_.dr.nstatic           = M_.nstatic;
oo_.dr.npred             = M_.npred;
oo_.dr.keep_states       = keep_states;

if DO_ESTIMATE_STATE
    tv_data = data';
else
    tv_data = [];
end
% Initialization of the time-varying filter (ADDED 2026-09-26): 'diffuse'
% (default, as always) or 'proper' (x_0 ~ N(X0, B0*B0'), no diffuse phase); see
% the header and kalman_filter_steady_and_diffuse.
init_type = 'diffuse';
if isfield(M_.X_InitialDistribution, 'type') && ~isempty(M_.X_InitialDistribution.type)
    init_type = lower(char(M_.X_InitialDistribution.type));
end
if ~any(strcmp(init_type, {'diffuse', 'proper'}))
    error('compute_state_reliability_diagnostics:initType', ...
        'M_.X_InitialDistribution.type must be ''diffuse'' or ''proper'', got ''%s''.', init_type);
end
if strcmp(init_type, 'proper') && (~isfield(M_.X_InitialDistribution, 'B0') || isempty(M_.X_InitialDistribution.B0) ...
                                   || ~isfield(M_.X_InitialDistribution, 'X0') || isempty(M_.X_InitialDistribution.X0))
    error('compute_state_reliability_diagnostics:initType', ...
        'M_.X_InitialDistribution.type = ''proper'' needs both X0 and B0 (x_0 ~ N(X0, B0*B0'')).');
end
[K_gain, ~, ~, KF_Omega, KF_V, KF_Omega_tp1_t, ~, Sequence, KF_ric_flag] = kalman_filter_steady_and_diffuse(A, B, C, D, tv_data, ...
    M_.X_InitialDistribution.X0, M_.X_InitialDistribution.B0, Qscale, T, oo_.dr, init_type);

%% ------------------------------------------------------------------
%  5. Population diagnostics
% -------------------------------------------------------------------
% Filtering-error recursion: with xtilde = x - xhat,
%   xtilde(t) = (A - K C) xtilde(t-1) + (B - K D) e(t),
% so AKC and BKD govern how estimation errors persist and are fed by shocks.
AKC = A - K_gain*C;
BKD = B - K_gain*D;

% Doubled system tracking (x, xhat) jointly:
%   [x(t); xhat(t)] = A_tilde [x(t-1); xhat(t-1)] + B_tilde e(t)
A_tilde = [A,        zeros(n_states, n_states);
           K_gain*C, AKC];
B_tilde = [B; K_gain*D];

% Numerical cleanup of BKD*BKD' (clip tiny entries, symmetrize, and project
% on the PSD cone) so the Lyapunov iterations below stay well behaved.
BKD(abs(BKD(:)) < 1e-7) = 0;
[V_eig, D_eig] = eig(BKD * BKD');
D_eig = diag(max(diag(D_eig), 0));
BKD_BKDt = V_eig*D_eig*V_eig';
BKD_BKDt = (BKD_BKDt + BKD_BKDt')/2;

% Zero out rows/columns corresponding to error components that receive no
% shocks at all
for i = 1:size(BKD, 1)
    if all(abs(BKD(i,:)) < 1e-10)
        BKD_BKDt(i,:) = 0;
        BKD_BKDt(:,i) = 0;
    end
end
BKD_BKDt(abs(BKD_BKDt(:)) < 1e-10) = 0;
BKD_BKDt = (BKD_BKDt + BKD_BKDt')/2;

n_BKD = size(BKD, 1);
if min(eig(BKD_BKDt)) < 0
    FixBKD = 1e-8;   % tiny diagonal loading if projection left a negative eig
else
    FixBKD = 0;
end

% Grid and defaults shared by all of section 5. Everything below is
% pre-set to "not computed" (NaN / []) and overwritten only by the blocks
% that actually run under the current fast-mode setting.
MaxIter  = 1e4;
conv_tol = 1e-12;
h_grid   = [1:100, 150:50:MaxIter, MaxIter];
All_STD  = diag(KF_Omega);   % filtered state variances (steady state)

Omega          = nan(n_states, n_states);
FEV            = nan(n_key_states, numel(h_grid));
CorrH          = nan(n_key_states, MaxH);
CorrH_knownX0  = nan(n_key_states, MaxH);
Error_AR_h     = nan(n_key_states, MaxH);
CorrY_H        = nan(n_obs, MaxH_Y);
R2_h           = nan(n_key_states, MaxH);
smoother_lag_m = nan(n_key_states, 1);
smoother_base  = cell(n_key_states, 1);
smoother_map_info = struct('source', 'not computed (fast_mode / loglik_only)');
nu             = nan(n_key_states, 1);
nu_order       = nan(n_key_states, 1);
max_HL = []; all_HL = []; weighted_HL = []; weights = [];
Omega_Breakdown = zeros(n_states, n_states, size(BKD, 2));

% Smoother-revision regression container (ADDED 2026-09-15); NaN unless a
% spec with a non-empty M_.fr_smoother_lags was given. The population side
% is filled here (it is a few matrix products, so it is not subject to
% fast_mode), the empirical side in step 6c.
fr = fr_regression_empty(spec_base, fr_m, MaxH_Y);
DO_FR = is_spec && ~isempty(fr_m) && ~LOGLIK_ONLY;

% Nothing in section 5 is used to evaluate the likelihood.
if ~LOGLIK_ONLY

if DO_FR
    fr = fr_regression_theory(A, K_gain, KF_V, KF_Omega, fr_base_idx, fr_lag_idx, fr_m, MaxH_Y, fr);
end

% --- 5a. FEV: error variance of xhat(t) built from information windows ---
% Iterating Var(err) from 0 gives, at iteration h, the variance of the
% estimation error when the filter has only run for h periods. FEV is then
% 1 - errvar(h)/errvar(inf) (normalized by the larger of the filter
% uncertainty and the limit), i.e. the share of the eventual information
% already accumulated after h observations.
if ~FAST
    Omega = zeros(n_states, n_states);   % start at zero: important for the FEV interpretation
    tmp_FEV = zeros(n_states, MaxIter);
    for ii = 1:MaxIter
        NewOmega = AKC*Omega*AKC' + BKD_BKDt + eye(n_BKD)*FixBKD;
        NewOmega(abs(NewOmega(:)) < 1e-10) = 0;
        NewOmega = (NewOmega + NewOmega')/2;
        tmp_FEV(:, ii) = diag(NewOmega);
        if norm(NewOmega - Omega, 'fro') < conv_tol
            % Converged: the remaining path is flat, fill and stop early
            tmp_FEV(:, ii+1:MaxIter) = repmat(diag(NewOmega), 1, MaxIter-ii);
            Omega = NewOmega;
            break;
        end
        Omega = NewOmega;
    end

    FEV = 1 - tmp_FEV(key_state_indices, h_grid) ./ ...
          repmat(max(All_STD(key_state_indices), tmp_FEV(key_state_indices, end)), [1, numel(h_grid)]);

    % Guard against round-off leaving Omega slightly indefinite
    if min(eig(Omega)) < 0
        disp(eig(Omega));
        [V_eig, D_eig] = eig(Omega);
        D_eig = diag(max(diag(D_eig), 0));
        NewOmega = real(V_eig*D_eig*V_eig');
        Omega = (NewOmega + NewOmega')/2;
        disp(eig(Omega));
    end
end

% --- 5b. Error-variance decomposition by shock ---
% Same Lyapunov equation, but feeding one shock at a time. Used to weight the
% half-lives.
% CHANGED (2026-09-16): each fixed point Sigma = AKC Sigma AKC' + b_u b_u' is
% obtained by the doubling (Smith) iteration - S <- S + M S M', M <- M^2 -
% which converges in log2 of the number of steps the plain recursion needs
% (~25 doublings at rho(AKC) = 0.99995, where the recursion hit its 1e4 cap
% on every pile-up draw, 7 times over). The old recursion also clipped
% entries below 1e-10 to zero at every step; the doubling result is the exact
% fixed point, so the two differ at that level for converged draws.
nshocks = size(BKD, 2);
if ~FAST
    for uu = 1:nshocks
        Omega_Breakdown(:,:,uu) = dlyap_doubling(AKC, BKD(:,uu)*BKD(:,uu)');
    end
end

% ORIGINAL CODE (per-shock Lyapunov recursion from KF_Omega), kept for comparison:
% if ~FAST
%     parfor uu = 1:nshocks
%         TmpOmega = KF_Omega;
%         for ii = 1:MaxIter
%             NewOmega = AKC*TmpOmega*AKC' + BKD(:,uu)*BKD(:,uu)';
%             NewOmega(abs(NewOmega(:)) < 1e-10) = 0;
%             NewOmega = (NewOmega + NewOmega')/2;
%             if norm(NewOmega - TmpOmega, 'fro') < conv_tol
%                 TmpOmega = NewOmega;
%                 break;
%             end
%             TmpOmega = NewOmega;
%         end
%         Omega_Breakdown(:,:,uu) = TmpOmega;
%     end
% end

% --- 5c. Correlation of states (and observables) with their estimates ---
% Propagate the joint covariance of (x, xhat) forward: after h periods of
% data the correlation corr(x_i, xhat_i) measures how much of the state the
% filter has recovered ("reliability"); squared it is the R2 of the paper.
%
% ADDED (2026-08): eq. (3.17), the reliability measure that is comparable
% ACROSS INFORMATION SETS, i.e. between the filter (m = 0) and the fixed-lag
% smoothers (m > 0) of Section 3.5.
%
%   R2^(m)_{i,h} = 1 - [Omega_m]_ii / Var[x_{i,t+h} | y^t],     eq. (3.17)
%
% Both estimates are residualised against the same baseline information set
% y^t and normalised by the SAME conditional total variance (3.10) - that of
% the CONTEMPORANEOUS state i, at horizon h. Since nested information sets
% order the error variances (Omega_l <= Omega_m for l >= m, eq. (3.16)), a
% common denominator makes reliability weakly increasing in the amount of
% hindsight, as Section 3.5 requires.
%
% WHAT WAS WRONG BEFORE. CorrH.^2 below is eq. (3.11) applied coordinate by
% coordinate, which for a fixed-lag smoother means applying it to the
% AUGMENTED LAGGED COPY (l20rrstar(t) = rrstar(t-20)). Its own denominator is
% then Var[r*_{t+h-m} | y^t], the conditional variance of the EARLIER state -
% exactly the failure mode flagged in Section 3.5 ("the resulting comparison
% is not the one covered by (3.16)"). For an I(1) state that denominator is
% short by m periods' worth of accumulated innovation variance, which is why
% the smoothed R2 in Table 3 of the draft came out BELOW the filtered one even
% though the smoothing error variance is uniformly smaller. At h = m the
% reported number even degenerates to 1 - Omega_m/Omega_0, a RATIO OF ERROR
% VARIANCES rather than a reliability.
%
% Note that KF_Omega already carries the right NUMERATOR for both cases: it is
% the steady-state filtered variance of the whole augmented state vector, so
% its diagonal entry for the lagged copy IS [Omega_m]_ii = Var[r*_{t-m} | y^t].
% Only the denominator has to be redirected, which is what
% fixed_lag_smoother_map resolves, just below.
CorrY_H = zeros(n_obs, MaxH_Y);
if FAST
    % Only CorrH_Observables is needed, and nothing beyond MaxH_Y is used.
    H_loop = MaxH_Y;
else
    H_loop        = MaxH;
    CorrH         = zeros(n_key_states, MaxH);
    CorrH_knownX0 = zeros(n_key_states, MaxH);
    Error_AR_h    = zeros(n_key_states, MaxH);
    R2_h          = zeros(n_key_states, MaxH);

    % Which states of interest are m-period lagged copies, and of what?
    %   - explicit M_.fixed_lag_smoother_map: always wins;
    %   - spec input (2026-09-15): known by construction - every 'l<m><NAME>'
    %     key is the m-period smoother of NAME, everything else is m = 0 - so
    %     the map is written from the spec;
    %   - legacy cellstr input: detected structurally from the impulse
    %     responses of the pruned system (fixed_lag_smoother_map), so no
    %     naming convention is assumed.
    if isfield(M_, 'fixed_lag_smoother_map')
        [smoother_lag_m, smoother_base_idx, smoother_map_info] = ...
            parse_fixed_lag_smoother_override(M_.fixed_lag_smoother_map, key_names, ...
                                              key_state_indices, state_var_names, M_.verbose);
    elseif is_spec
        sel = key_lag_m > 0;
        spec_map = [key_names(sel), key_base(sel), num2cell(key_lag_m(sel))];
        [smoother_lag_m, smoother_base_idx, smoother_map_info] = ...
            parse_fixed_lag_smoother_override(spec_map, key_names, ...
                                              key_state_indices, state_var_names, M_.verbose);
        smoother_map_info.source = 'states_of_interest spec (augment_ssm labels)';
    else
        [smoother_lag_m, smoother_base_idx, smoother_map_info] = ...
            fixed_lag_smoother_map(A, B, key_state_indices, state_var_names, [], [], M_.verbose);
    end
    smoother_base = smoother_map_info.base_name;
end

% Initial condition: x(0) uncertain with the steady-state filter variance,
% xhat(0) known (the doubled system starts at the filter's information set)
Cov0 = A_tilde*[KF_Omega, zeros(n_states, n_states); zeros(n_states, 2*n_states)]*A_tilde' + B_tilde*B_tilde';
Cov  = Cov0;
C_aug = [C, zeros(n_obs, n_states); zeros(n_obs, n_states), C];
CovYYhat = C_aug*Cov*C_aug' + [D*D', zeros(n_obs, n_obs); zeros(n_obs, 2*n_obs)];
AutoCov0 = AKC*KF_Omega*AKC' + BKD_BKDt + eye(n_BKD)*FixBKD;
Cov_knownX0 = B_tilde*B_tilde';   % alternative: initial state known exactly

% CHANGED (2026-09-15): the error autocorrelation used AKC^hh inside the
% (hh, nn) double loop - a fresh matrix power per horizon and per state
% (18,838 powers of a 77x77 matrix on the augmented HLW system, ~0.5 s of a
% 3.7 s call). AKC^hh is now accumulated as an iterated product once per hh
% and all states of interest are read off the diagonal of AKC^hh*AutoCov0 at
% once. Same quantity; the two differ only in floating-point rounding
% (matrix power uses repeated squaring), ~1e-15 relative.
AKCh     = eye(n_states);
dAutoCov = diag(AutoCov0);

for hh = 1:H_loop
    if ~FAST
        AKCh = AKCh * AKC;                       % AKC^hh
        dAA  = diag(AKCh * AutoCov0);
        Error_AR_h(:, hh) = dAA(key_state_indices) ./ dAutoCov(key_state_indices);

        for nn = 1:n_key_states
            kk = key_state_indices(nn);

            % Autocorrelation of the filtering error at lag h - ORIGINAL CODE,
            % kept for comparison with the vectorised version above:
            % sel_vec = zeros(1, n_states);
            % sel_vec(1, kk) = 1;
            % Error_AR_h(nn, hh) = sel_vec*AKC^hh*AutoCov0*sel_vec' / (sel_vec*AutoCov0*sel_vec');

            % corr(x, xhat) after h periods (block 1 = x, block 2 = xhat)
            if abs(Cov(kk + n_states, kk + n_states)) > 1e-12
                CorrH(nn, hh) = Cov(kk, kk + n_states) / ...
                    sqrt(Cov(kk, kk) * Cov(kk + n_states, kk + n_states));
            else
                CorrH(nn, hh) = nan(1);
            end
            if abs(Cov_knownX0(kk + n_states, kk + n_states)) > 1e-12
                CorrH_knownX0(nn, hh) = Cov_knownX0(kk, kk + n_states) / ...
                    sqrt(Cov_knownX0(kk, kk) * Cov_knownX0(kk + n_states, kk + n_states));
            else
                CorrH_knownX0(nn, hh) = nan(1);
            end

            % Eq. (3.17). Numerator: the estimation error variance of the
            % coordinate that is actually being estimated - [Omega_m]_ii for a
            % fixed-lag smoother, [Omega_0]_ii = Omega_ii for the filter.
            % Denominator: eq. (3.10), Var[x_{t+h}|y^t], for the CONTEMPORANEOUS
            % base state, which is Cov's (1,1) block at row/column bb. When the
            % state of interest is not a lagged copy, bb == kk and this is
            % algebraically identical to CorrH(nn,hh)^2, i.e. eq. (3.11).
            bb = smoother_base_idx(nn);
            if Cov(bb, bb) > 1e-12
                R2_h(nn, hh) = 1 - KF_Omega(kk, kk) / Cov(bb, bb);
            else
                R2_h(nn, hh) = nan(1);
            end
        end
    end

    % Same object for the observables (used to benchmark against data)
    if hh <= MaxH_Y
        for nn = 1:n_obs
            if abs(CovYYhat(nn, nn)) > 1e-12
                CorrY_H(nn, hh) = CovYYhat(nn, nn + n_obs) / sqrt(CovYYhat(nn, nn)*CovYYhat(nn + n_obs, nn + n_obs));
            else
                CorrY_H(nn, hh) = nan(1);
            end
        end
    end

    % Advance the joint covariances one period.
    % (The recursion adds B_tilde*B_tilde' each step; an earlier version
    % incorrectly re-added Cov0, which embeds the initial-state uncertainty.)
    Cov = A_tilde*Cov*A_tilde' + B_tilde*B_tilde';
    if ~FAST
        Cov_knownX0 = A_tilde*Cov_knownX0*A_tilde' + B_tilde*B_tilde';
    end
    if hh <= MaxH_Y
        CovYYhat = C_aug*Cov*C_aug' + [D*D', zeros(n_obs, n_obs); zeros(n_obs, 2*n_obs)];
    end
end

if M_.verbose
    disp('KF hidden state standard deviation: t|t');
    disp(sqrt(diag(KF_Omega(key_state_indices, key_state_indices))));
    disp('KF hidden state standard deviation: t+1|t');
    disp(sqrt(diag(KF_Omega_tp1_t(key_state_indices, key_state_indices))));
    disp('Filtering error standard deviation');
    disp(sqrt(diag(Omega(key_state_indices, key_state_indices))));
end

if ~FAST
    % CHANGED (2026-09-15): CorrH of a fixed-lag SMOOTHED state (m > 0) is now
    % the conditional correlation that is the regression analogue of eq.
    % (3.17): conditional on y^t, xhat_{i,t+h|t+h+m} is the conditional
    % expectation of x_{i,t+h} under a finer information set, so
    %     Cov(x_{i,t+h}, xhat | y^t) = Var(xhat | y^t)
    % and the squared conditional correlation - the R2 of regressing
    % x_{i,t+h} on the smoothed estimate given y^t - equals
    %     Var(xhat | y^t) / Var(x_{i,t+h} | y^t) = 1 - [Omega_m]_ii / Var(x_{i,t+h} | y^t) = R2_h.
    % Hence CorrH = sqrt(R2_h) for those rows, with the CONTEMPORANEOUS base
    % state's variance in the denominator. The value the doubled-system
    % recursion produced above for such a row - eq. (3.11) applied to the
    % lagged copy, i.e. corr(x_{i,t+h-m}, xhat_{i,t+h-m|t+h} | y^t) with the
    % EARLIER state's variance - is kept in CorrH_eq311 for the verbose
    % comparison below only. Rows with m = 0 are untouched (CorrH.^2 == R2_h
    % there already, checked below). max(.,0) guards the sqrt against a
    % round-off-negative R2_h at short horizons.
    CorrH_eq311 = CorrH;
    for nn = find(smoother_lag_m(:)' > 0)
        CorrH(nn, :) = sqrt(max(R2_h(nn, :), 0));
    end

    % Consistency checks on eq. (3.17). Both are properties the measure must have
    % by construction, so a violation means a bug (or a numerically degenerate
    % draw), not an economic finding.
    %   (i)  for m = 0 it must reproduce eq. (3.11), i.e. CorrH.^2;
    %   (ii) at every horizon it must be weakly higher for a smoother than for the
    %        filter on the same base state - the point of Section 3.5.
    % Tolerance 1e-6 (was 1e-8 until 2026-09-15): near-unit-root draws put
    % round-off of 1e-8..2e-8 into the 250-step covariance recursion, which
    % is noise for this purpose, not a bug.
    R2_h_eq311_gap = max(abs(R2_h(smoother_lag_m == 0, :) - CorrH(smoother_lag_m == 0, :).^2), [], 2);
    if ~isempty(R2_h_eq311_gap) && any(R2_h_eq311_gap > 1e-6)
        fprintf(2, ['compute_state_reliability_diagnostics: eq. (3.17) and eq. (3.11) disagree by ' ...
                    '%.2e on a FILTERED state; they are algebraically identical there.\n'], max(R2_h_eq311_gap));
    end
    R2_h_hindsight_gap = 0;
    for nn = 1:n_key_states
        if smoother_lag_m(nn) > 0
            mm = find(key_state_indices == smoother_base_idx(nn), 1);   % is the base itself reported?
            if ~isempty(mm)
                R2_h_hindsight_gap = min(R2_h_hindsight_gap, min(R2_h(nn, :) - R2_h(mm, :)));
            end
        end
    end
    if R2_h_hindsight_gap < -1e-8
        fprintf(2, ['compute_state_reliability_diagnostics: eq. (3.17) is NOT weakly increasing in the ' ...
                    'amount of hindsight (worst gap %.2e). Section 3.5 says it must be.\n'], R2_h_hindsight_gap);
    elseif M_.verbose
        for nn = 1:n_key_states
            if smoother_lag_m(nn) > 0
                fprintf(['Eq. (3.17) R2 for %s (m = %d hindsight on %s): %.4f at h = 4, %.4f at h = 20 ' ...
                         '(eq. (3.11) on the lagged copy gave %.4f and %.4f)\n'], ...
                    key_names{nn}, smoother_lag_m(nn), smoother_base{nn}, ...
                    R2_h(nn, min(4, MaxH)), R2_h(nn, min(20, MaxH)), ...
                    CorrH_eq311(nn, min(4, MaxH))^2, CorrH_eq311(nn, min(20, MaxH))^2);
            end
        end
    end

    % --- 5d. Half-lives of estimation errors ---
    % Skip when the error variance is identically zero (states observed
    % exactly, no estimation error to decay).
    if any(diag(Omega) ~= 0)
        [max_HL, all_HL, shock_idx] = estimation_error_half_life(AKC, BKD, key_state_indices, 1e3);
        if ~any(isnan(max_HL))
            if M_.verbose
                disp('half life');
                disp(max_HL);
            end
            % Weight each shock's half-life by its share of the error variance
            weighted_HL = zeros(1, n_key_states);
            for nn = 1:n_key_states
                weights = squeeze(Omega_Breakdown(key_state_indices(nn), key_state_indices(nn), :)) / ...
                          sum(squeeze(Omega_Breakdown(key_state_indices(nn), key_state_indices(nn), :)));
                weighted_HL(nn) = all_HL(:, nn)' * weights(shock_idx);
            end
        else
            max_HL = []; all_HL = []; weighted_HL = []; weights = [];
        end
    else
        max_HL = []; all_HL = []; weighted_HL = []; weights = [];
    end

    % --- 5e. Tail-unreliability index nu (ADDED 2026-08, REVISED 2026-09) ----
    % Horizon-free reliability index for the non-stationary states of interest,
    % from the note "A horizon-free reliability index for non-stationary states of
    % arbitrary integration order" (Rate_of_convergence_R2.pdf, Theorem 2.2).
    %
    % MOTIVATION. The conditional reliability computed in 5c satisfies
    %   R2_{h,i} = 1 - Omega_ii / Var[x_{i,t+h} | y^t]     (paper eq. (3.11)),
    % whose numerator does NOT depend on h. For a unit-root state the denominator
    % diverges in h, so R2_{h,i} -> 1 for every model, state and parameterization
    % (paper's Theorem 3.4): the long-horizon columns of the headline tables are
    % pinned against 1 and cannot discriminate between models. All of the
    % discriminating content lies in the RATE of that convergence, which Theorem
    % 2.2 of the note isolates as a single horizon-free number:
    %
    %   h^(2d+1) * (1 - R2_{h,i})  -->  nu_i = (2d+1) * Omega_ii / L_i^(d),
    %
    % as h -> Inf, where d = d_i is the integration ORDER of state i (state i is
    % I(d_i + 1)) and L_i^(d) = ||e_i' N^d B||^2/(d!)^2 is its generalised long-run
    % innovation variance, built from the unit-root projector Pi and the nilpotent
    % N = (A - I)*Pi. Note that L_i^(d) is NOT the state's own innovation variance
    % [B*B']_ii - see tail_unreliability_index.m for the derivation and the
    % counterexamples. Larger nu_i means LESS reliable.
    %
    % CHANGED (2026-09). This step used to call the d = 0 specialisation
    % h*(1 - R2_{h,i}) -> Omega_ii/L_i, which is valid only when every
    % non-stationary state is I(1) and returned NaN otherwise. HLW's ystar is I(2)
    % (ystar = ystar(-1) + gg(-1)/4 + shock with gg a random walk) and used to come
    % back empty for that reason; it now receives an index. Every state with d_i = 0
    % - which is all six of HLW's states_of_interest, ygap through d20rrstar - keeps
    % the number it had before, since eq. (5) reduces to the old formula at d = 0.
    %
    % UNITS. nu_i is in QUARTERS^(2*d_i+1), so it is in quarters for the I(1) states
    % and only for those. `nu_order` below carries d_i alongside nu, and entries
    % with different d_i must NOT be pooled into one table column or compared to
    % each other (note, section 2.3). At d_i = 0, and when row i of A equals e_i',
    % nu_i is exactly the half-reliability horizon of the R2_{h,i} curve.
    %
    % WHICH Omega. KF_Omega is the steady-state FILTERED error variance
    % Var(x_t | y^t), i.e. exactly the numerator of 1 - R2_{h,i}, and it is the
    % same matrix that seeds Cov0 in step 5c. Using it (rather than the Lyapunov
    % iterate `Omega` of step 5a, which has been clipped and PSD-projected) makes
    % nu the exact tail coefficient of the R2_{h,i} path reported above - i.e. of
    % CorrH.^2, and equally of the variable `R2_h` added in step 5c, since the two
    % differ only by a shift of the horizon in the denominator and that shift
    % washes out in the tail. nu is therefore already comparable across
    % information sets: a fixed-lag smoother inherits the same order d_i and the
    % same L_i^(d) as its base state, so nu_smoother/nu_filter = Omega_m/Omega_0,
    % which is the eq.-(3.16) ordering.
    %
    % NaN ENTRIES are expected and meaningful: nu is undefined for a state with no
    % stochastic trend at any order (the output gap ygap, for which
    % h*(1 - R2_{h,i}) diverges instead of converging). For those states the
    % paper's UNCONDITIONAL R2 is already the right horizon-free summary, so nu and
    % R2 must not be pooled into a single table column. nu is also NaN if
    % Assumption 1 of the note fails for this draw (an explosive root, or a root of
    % unit modulus away from +1), which keeps a pathological draw from crashing the
    % sampler. A defective unit root is NO LONGER such a case.
    %
    % NOTE. Do not try to read nu off h^(2d+1)*(1 - R2_{h,i}) at a moderate
    % horizon: the convergence is slow. The expression above is exact and
    % closed-form, so in the Bayesian implementation the posterior of nu comes from
    % evaluating it at each draw.
    %
    % Only nu and nu_order are packed into LatentStateDiagnostics below. The
    % underlying pieces are available from the same call should the other objects
    % of the note be wanted later (generalised long-run innovation variances L, the
    % unit-root projector Pi, and an `info` struct carrying rho(A_s), the ladder of
    % row norms ||e_i' N^k Pi B|| from which the orders are read, and the numerical
    % checks): [nu, L, Pi, nu_info, nu_order] = tail_unreliability_index(...).
    [nu, ~, ~, nu_info, nu_order] = tail_unreliability_index(A, B, KF_Omega, ...
                                                            key_state_indices, M_.verbose);

    if M_.verbose
        disp('Tail-unreliability index nu (NaN = no stochastic trend at any order)');
        disp('  [state index, order d, nu in quarters^(2d+1)]:');
        disp([key_state_indices(:), nu_order(:), nu(:)]);
        if nu_info.max_order > 0
            disp(['NOTE: orders differ across states, so nu is NOT in a common ' ...
                  'unit - do not pool these in one column.']);
        end
    end
end % ~FAST

end % ~LOGLIK_ONLY

%% ------------------------------------------------------------------
%  6. Data pass: filter the observed data (only when data is supplied)
% -------------------------------------------------------------------
% Defaults so that the output structs are always well defined
TV_state_estimate = [];
TV_surprises      = [];
TV_log_llhood     = NaN;
TV_CorrY_H     = nan(n_obs, MaxH_Y);
TV_RMSE        = nan(n_obs, MaxH_Y);
fh             = nan(n_obs, MaxH_Y, n_states);
fh_direct      = nan(n_obs, MaxH_Y, n_states);
TV_CorrY_H_states = nan(n_states, MaxH_Y);
TV_RMSE_states    = nan(n_states, MaxH_Y);

if DO_ESTIMATE_STATE
    tol = 1e-6;
    [T, n_data] = size(data);

    if isequal(n_obs, n_data)

        % Demean the data with the model's steady state of the observables
        % (obs_ss is computed at the end of section 1)
        data = data - obs_ss';   % T x n_obs

        %% ============================================================
        %  6a. Sample selection, shared by EVERY empirical statistic
        % =============================================================
        % Two kinds of period are excluded, and the SAME mask is then used
        % for the correlations and for the RMSEs.
        %
        %   (i) BURN-IN. The filter starts at X0 (zero unless set on M_)
        %       while the non-stationary states sit at the level of the data
        %       (yobs ~ 816 in uc3_level_i2), so the early forecast errors
        %       are dominated by the initialization transient. That transient
        %       decays at max|eig(A - K C)|, i.e. at exactly the estimation-
        %       error half-life computed in step 5d, and over the
        %       contaminated window both the forecast error and the forecast
        %       revision are driven by the same decaying path, which pushes
        %       corr(error, revision) towards one and inflates the RMSE by
        %       orders of magnitude.
        %
        %   (ii) HETEROSKEDASTIC PERIODS. Qscale ~= 1 breaks the constant-
        %        variance interpretation of these statistics.
        %
        % CHANGED (2026-08): the burn-in used to be applied ONLY to the TV
        % correlations - not to the (now removed) steady-state pass, and to
        % neither RMSE. The two filters' correlations were therefore computed
        % on different samples (on UC3, R2 gaps of up to 0.47 that vanish to
        % 0.0000 once the samples match), and both RMSEs measured the
        % initialization rather than forecast accuracy (yobs, h = 4: 133pp
        % vs 2.8pp once the burn-in is applied).
        %
        % CAVEAT ON THE DEFAULT: 50 quarters is NOT enough when the
        % estimation-error half-life is long. In uc3_level_i2 it is 72-81
        % quarters (Table 1 of the draft), so a 816pp initial error is still
        % ~500pp at t = 50. Raise M_.reliability_burn_in accordingly, or
        % lengthen the sample; ModelMatrices.AKC gives the decay rate.
        if isfield(M_, 'reliability_burn_in') && ~isempty(M_.reliability_burn_in)
            burn_in = M_.reliability_burn_in;
        else
            burn_in = 50;
        end
        burn_in = max(0, min(burn_in, T));

        if ~isfield(M_.heteroskedastic_shocks, 'Qscale')
            M_.heteroskedastic_shocks.Qscale = ones(n_shocks, T);
        end
        Q = M_.heteroskedastic_shocks.Qscale;
        Q(:, T+1:end) = [];
        Q(isnan(Q)) = 1;
        % Shocks whose rescaling is a DATA DEVICE rather than a variance break
        % are exempted from the exclusion through M_.reliability_mask_ignore_shocks
        % (cell of exo names; default none, i.e. unchanged behaviour). ADDED
        % 2026-09-19 for HLW_2023_US_TD_SPF.mod, which scales the SPF measurement
        % error by 10000 before 1992 and by 4/8/12 in non-survey quarters: with
        % the plain rule only 34 of 264 periods survived and EVERY empirical
        % statistic (TV_RMSE, TV_CorrY_H, fr.*_emp) was NaN. The exempted
        % shock still enters the time-varying filter with its stored scales;
        % only the sample-selection mask ignores it.
        Qmask = Q;
        if isfield(M_, 'reliability_mask_ignore_shocks') && ~isempty(M_.reliability_mask_ignore_shocks)
            ign = ismember(cellstr(M_.exo_names), cellstr(M_.reliability_mask_ignore_shocks));
            if ~any(ign)
                warning('compute_state_reliability_diagnostics:maskIgnore', ...
                        'M_.reliability_mask_ignore_shocks names no shock of the model; mask unchanged.');
            end
            Qmask(ign, :) = 1;
        end
        TV_var = any(abs(Qmask - 1) > 1e-4);   % 1 x T flag: any (non-exempted) shock rescaled at t
        IDX_ok = (TV_var == 0);
        IDX_ok(1:burn_in) = false;

        if HAVE_TRUE_STATES && ~isequal(size(true_states), [T, n_states])
            error('compute_state_reliability_diagnostics:badTrueStates', ...
                ['true_states must be T x n_states = %d x %d (columns ' ...
                 'ordered like the pruned state vector); got %d x %d.'], ...
                T, n_states, size(true_states, 1), size(true_states, 2));
        end

        %% ========================================================
        %  6b. Diffuse / time-varying filter: apply the pre-computed
        %      gains to the data. This is the ONLY data pass.
        %      - diffuse phase: exact univariate DK update on the
        %        stacked state [y-only; x], using Sequence.K_stacked;
        %      - post-diffuse: standard multivariate update on x with
        %        Sequence.K_gain.
        %      The recursion lives in apply_tv_filter (extracted verbatim
        %      2026-08 so that posterior_predictive_wrapper can reuse it on
        %      many data sets against one Sequence).
        % =========================================================
        if ~isfield(M_, 'X_InitialDistribution') || ...
                ~isfield(M_.X_InitialDistribution, 'X0') || ...
                isempty(M_.X_InitialDistribution.X0)
            X0_use = [];                  % apply_tv_filter reads [] as zeros(n_states,1)
        else
            X0_use = M_.X_InitialDistribution.X0;
        end

        [TV_log_llhood, all_states, TV_surprises, filter_ok] = ...
            apply_tv_filter(A, C, Sequence, data, X0_use, tol);

        if filter_ok
            TV_state_estimate = all_states(:, key_state_indices);
        else
            % Same contract as the inline recursion had: empty on failure
            TV_state_estimate = [];
            TV_surprises      = [];
        end

        if filter_ok
            %% ====================================================
            %  6c. Empirical forecast RMSE and error/revision
            %      correlations by horizon (the DeltaR2 building
            %      blocks; eq. (5.2)-(5.3) of the paper)
            % =====================================================
            for hh = 1:MaxH_Y
                [TV_yerr, TV_rev] = forecast_errors_and_revisions(A, C, all_states, data, hh);

                % Exclude t if ANY period in [t-1-hh, ..., t] is excluded:
                % forecast origin and all intermediate periods must be clean
                IDX_ok_hh = ~backward_window_max(double(~IDX_ok), hh+1);

                for nn = 1:n_obs
                    X1 = TV_yerr(IDX_ok_hh, nn);
                    X2 = TV_rev(IDX_ok_hh, nn);

                    % Joint NaN mask keeps the pairs aligned (the original
                    % code masked X1 and X2 separately, which errors or
                    % misaligns when their NaN patterns differ)
                    ok = ~isnan(X1) & ~isnan(X2);
                    TV_CorrY_H(nn, hh) = corr(X1(ok), X2(ok));

                    % Same sample for the RMSE (it used to run over the whole
                    % series, burn-in included - see the note in step 6a)
                    TV_RMSE(nn, hh) = sqrt(mean(X1(~isnan(X1)).^2));
                end
            end

            % Smoother-revision regression on the filtered data (ADDED
            % 2026-09-15): OLS of f = xhat_{t+h|t+h+m} - xhat_{t+h|t} on
            % r = xhat_{t+h|t+h} - xhat_{t+h|t}, per base state, h and m,
            % on the same admissible sample (window t..t+h+m clean).
            if DO_FR
                fr = fr_regression_empirical(all_states, A, fr_base_idx, fr_lag_idx, fr_m, MaxH_Y, IDX_ok, fr);
            end

            % State-level forecast statistics against the TRUE states
            % (simulation experiments only; MERGED 2026-08 from the deleted
            % compute_state_reliability_vs_truth). Same objects as above but
            % for the states instead of the observables, with the truth in
            % place of the data, and on the same sample.
            if HAVE_TRUE_STATES
                for hh = 1:MaxH_Y
                    [TV_xerr, TV_xrev] = state_forecast_errors_and_revisions(A, all_states, true_states, hh);
                    IDX_ok_hh = ~backward_window_max(double(~IDX_ok), hh+1);
                    for nn = 1:n_states
                        X1 = TV_xerr(IDX_ok_hh, nn);
                        X2 = TV_xrev(IDX_ok_hh, nn);
                        ok = ~isnan(X1) & ~isnan(X2);   % joint NaN mask
                        TV_CorrY_H_states(nn, hh) = corr(X1(ok), X2(ok));
                        TV_RMSE_states(nn, hh)    = sqrt(mean(X1(~isnan(X1)).^2));
                    end
                end
            end

            % Theoretical R2 wedge between each state and the observables
            % (eq. (26) objects). The steady-state filter inside
            % theoretical_R2_wedge is identical across states, so it is
            % computed once and reused via 'filter_cache'. Both wedges are
            % population objects and are skipped in fast mode.
            if ~FAST
                % CHANGED (2026-09-15): one call for all states. theoretical_R2_wedge
                % now accepts a vector 'state_index' and returns the wedge as
                % n_obs x n_h x n_states; the observable-side objects and the
                % forward covariance path, which do not depend on the state, are
                % computed once instead of n_states times (0.75 s -> ~0.05 s on the
                % augmented HLW system). Same numbers up to rounding (~1e-15).
                out_wedge = theoretical_R2_wedge(A, B, C, D, 1:MaxH_Y, 'state_index', 1:n_states);
                fh = out_wedge.wedge;

                % ORIGINAL CODE, kept for comparison with the vectorised call above:
                % filter_cache = [];
                % for i_states = 1:n_states
                %     if isempty(filter_cache)
                %         out_wedge = theoretical_R2_wedge(A, B, C, D, 1:MaxH_Y, 'state_index', i_states);
                %         filter_cache = out_wedge.steady_state;
                %     else
                %         out_wedge = theoretical_R2_wedge(A, B, C, D, 1:MaxH_Y, 'state_index', i_states, ...
                %             'filter_cache', filter_cache);
                %     end
                %     fh(:, :, i_states) = out_wedge.wedge;
                % end

                % The OTHER version of the same eq.-(26) object: the direct
                % adjustment factor of R2_wedge_factor, returned as .fh_direct
                % (MERGED 2026-08; it was .fh in the deleted vs_truth engine,
                % which in turn called the block above .fh_shift - see the header).
                % Its steady-state filter does not depend on the horizon, so it is
                % computed once and passed as a cache. Kept as a separate Riccati
                % call, exactly as in the deleted file, so the numbers reproduce
                % bit-for-bit (steady_state_kalman_filter iterates from 0, the
                % engine's own solution from B0*B0').
                [wf_K, ~, ~, wf_Sigma, wf_V] = steady_state_kalman_filter(A, B, C, D);
                wf_cache = struct('K_gain', wf_K, 'Sigma', wf_Sigma, 'Sigma_v', wf_V);
                for hh = 1:MaxH_Y
                    wedge_temp = R2_wedge_factor(A, B, C, D, hh, wf_cache);   % n_states x n_obs
                    fh_direct(:, hh, :) = wedge_temp';
                end
            end
        else
            fprintf(2, 'Cholesky Decomp of V failed \n');
        end
    else
        fprintf(2, 'wrong number of data series \n');
    end
end

%% ------------------------------------------------------------------
%  7. Pack the outputs
% -------------------------------------------------------------------
% `nu` (ADDED 2026-08) is n_key_states x 1, ordered like states_of_interest,
% i.e. like the ROWS of CorrH; it is NaN for states without a stochastic trend
% (see step 5e). `nu_order` (ADDED 2026-09) is the matching vector of
% integration orders d_i: nu(i) is in quarters^(2*nu_order(i)+1), so it is in
% quarters exactly where nu_order is 0 - which is every state of interest in
% HLW and UC3 - and entries with different nu_order are in different units.
% REMOVED (2026-08): the fields 'steady_CorrY_H' and 'steady_RMSE', together
% with the steady-state-gain data pass that produced them (see step 6).
% 'TV_CorrY_H' / 'TV_RMSE' keep their names and are now the only empirical
% statistics; 'CorrH_Observables' remains their THEORETICAL (steady-state,
% population) counterpart, so DeltaR2 = CorrH_Observables.^2 - TV_CorrY_H.^2.
%
% ADDED (2026-08, merged from the deleted compute_state_reliability_vs_truth):
%   'fh_direct'          n_obs x MaxH_Y x n_states, the eq.-(26) wedge from
%                        R2_wedge_factor. 'fh' is UNCHANGED (theoretical_R2_wedge).
%                        In the deleted engine these two were called 'fh' and
%                        'fh_shift' RESPECTIVELY SWAPPED - see the header.
%   'TV_CorrY_H_states'  n_states x MaxH_Y, state-vs-truth error/revision
%   'TV_RMSE_states'     n_states x MaxH_Y, state-vs-truth forecast RMSE
%                        Both are NaN unless `true_states` was supplied.
%
% ADDED (2026-08, eq. (3.17) / Section 3.5):
%   'R2_h'            n_key_states x MaxH, the reliability measure that is
%                     comparable across information sets,
%                     R2^(m)_{i,h} = 1 - [Omega_m]_ii / Var[x_{i,t+h}|y^t].
%                     THIS, not CorrH.^2, is what belongs in a table that puts
%                     filtered and smoothed estimates side by side (Table 3).
%                     Rows whose state of interest is not a lagged copy carry
%                     m = 0 and are identical to CorrH.^2 (eq. (3.11)), so
%                     nothing changes for the filtered columns.
%   'smoother_lag_m'  n_key_states x 1, the m of eq. (3.17); 0 = filter.
%   'smoother_base'   n_key_states x 1 cell, name of the contemporaneous base
%                     state supplying the denominator (itself when m = 0).
%   'smoother_map_info'  the diagnostics of fixed_lag_smoother_map (accepted
%                     match quality and all candidate matches), for auditing.
%
% ADDED (2026-09-15):
%   'state_labels'    n_key_states x 1 cell, the diagnosed states in row
%                     order (base names, 'l<m><NAME>', 'd<h><NAME>').
%   'state_spec'      the parsed states_of_interest: .base, .lags, .changes,
%                     .fr_smoother_lags, .is_spec.
%   'fr'              the smoother-revision regression (fr_regression_empty
%                     for the field list); population side from step 5,
%                     empirical side from step 6c.
state_spec = struct('base', {spec_base}, 'lags', {spec_lags}, 'changes', {spec_changes}, 'defs', {spec_defs}, ...
                    'fr_smoother_lags', fr_m, 'is_spec', is_spec);
fieldNames = {'CorrH', 'Error_ACF', 'max_HL', 'weighted_HL', 'all_HL', 'shock_contributions', 'FEV', 'CorrH_knownX0', 'TV_KF_Sequence', 'CorrH_Observables', 'TV_CorrY_H', 'TV_RMSE', 'fh', 'fh_direct', 'TV_CorrY_H_states', 'TV_RMSE_states', 'nu', 'nu_order', 'R2_h', 'smoother_lag_m', 'smoother_base', 'smoother_map_info', ...
    'state_labels', 'state_spec', 'fr'};
variables  = {CorrH, Error_AR_h, max_HL, weighted_HL, all_HL, weights, FEV, CorrH_knownX0, Sequence, CorrY_H, TV_CorrY_H, TV_RMSE, fh, fh_direct, TV_CorrY_H_states, TV_RMSE_states, nu, nu_order, R2_h, smoother_lag_m, smoother_base, smoother_map_info, ...
    key_names, state_spec, fr};
LatentStateDiagnostics = cell2struct(variables, fieldNames, 2);

% KF_V (ADDED 2026-08): the steady-state innovation covariance C*Omega*C'+D*D'.
% It used to be consumed only by the removed steady-state data pass; it is a
% population object of the paper (the variance term of eq. (5.1)), so it is
% returned here instead of being discarded.
% KF_ric_flag (ADDED 2026-09): the steady-state Riccati convergence flag
% (.converged/.iterations/.rel_step/.rho_cl). Carried out of here so the
% per-draw wrapper can record it: KF_Omega, KF_Omega_tp1_t, KF_V, AKC, BKD and
% every half-life and FEV derived from them are only meaningful when
% .converged is true, and a pile-up draw can exhaust the iteration cap.
% X0_pruned / B0_pruned (MERGED 2026-09): the initial distribution after the
% augmentation of step 2 and the pruning of step 3, i.e. exactly what the
% filter above was started from. simulate_artificial_data needs them to draw
% the initial state of the DGP consistently with the filter.
if isfield(M_, 'X_InitialDistribution') && isfield(M_.X_InitialDistribution, 'X0')
    X0_pruned = M_.X_InitialDistribution.X0;
else
    X0_pruned = [];   % apply_tv_filter treats [] as zeros(n_states,1)
end
if isfield(M_, 'X_InitialDistribution') && isfield(M_.X_InitialDistribution, 'B0')
    B0_pruned = M_.X_InitialDistribution.B0;
else
    B0_pruned = [];
end

% state_names (ADDED 2026-09-15): labels of the rows of the pruned augmented
% state vector, i.e. of A, B, K_gain, KF_Omega and of the columns of
% TV_KF_Sequence.Omega. Needed to address the augment_ssm states by name.
fieldNames = {'AKC', 'BKD', 'A', 'B', 'C', 'D', 'K_gain', 'A_tilde', 'B_tilde', 'KF_Omega', 'KF_Omega_tp1_t', ...
    'KF_V', 'Omega', 'A_original', 'B_original', 'C_original', 'D_original', 'obs_ss', 'KF_ric_flag', ...
    'X0_pruned', 'B0_pruned', 'state_names'};
variables  = {AKC, BKD, A, B, C, D, K_gain, A_tilde, B_tilde, KF_Omega, KF_Omega_tp1_t, KF_V, Omega, ...
    A_original, B_original, C_original, D_original, obs_ss, KF_ric_flag, X0_pruned, B0_pruned, state_var_names};
ModelMatrices = cell2struct(variables, fieldNames, 2);

end


%% ========================================================================
function [lag_m, base_idx, info] = parse_fixed_lag_smoother_override(map, states_of_interest, key_state_indices, state_var_names, verbose)
%PARSE_FIXED_LAG_SMOOTHER_OVERRIDE  Read M_.fixed_lag_smoother_map.
%   Bypasses the structural detection of fixed_lag_smoother_map when the
%   lag structure is known and should be pinned down in the .mod file:
%
%       M_.fixed_lag_smoother_map = { 'l20rrstar', 'rrstar',  20
%                                     'ldrrstar',  'drrstar',  1 };
%
%   k x 3 cell array of {state of interest, contemporaneous base state, m}.
%   States of interest not listed are treated as filtered estimates (m = 0),
%   so M_.fixed_lag_smoother_map = cell(0,3) disables the correction entirely
%   and restores the pre-2026-08 eq.-(3.11) numbers.
%
%   The base must be a state of the PRUNED state vector, which in practice
%   means it must itself be one of M_.states_of_interest (otherwise step 3
%   may prune it); an unknown name is an error rather than a silent m = 0.
    n_key    = numel(key_state_indices);
    lag_m    = zeros(n_key, 1);
    base_idx = key_state_indices(:);

    info = struct();
    info.base_name  = state_var_names(base_idx);
    info.mismatch   = nan(n_key, 1);
    info.candidates = cell(n_key, 1);
    info.max_lag    = 0;
    info.tol_rel    = NaN;
    info.source     = 'M_.fixed_lag_smoother_map';

    if isempty(map)
        return
    end
    if ~iscell(map) || size(map, 2) ~= 3
        error('compute_state_reliability_diagnostics:badSmootherMap', ...
            'M_.fixed_lag_smoother_map must be a k x 3 cell array {state, base, lag}.');
    end

    for rr = 1:size(map, 1)
        nn = find(strcmp(states_of_interest, map{rr, 1}), 1);
        if isempty(nn)
            error('compute_state_reliability_diagnostics:badSmootherMap', ...
                'M_.fixed_lag_smoother_map row %d names "%s", which is not in M_.states_of_interest.', ...
                rr, map{rr, 1});
        end
        bb = find(strcmp(state_var_names, map{rr, 2}), 1);
        if isempty(bb)
            error('compute_state_reliability_diagnostics:badSmootherMap', ...
                ['M_.fixed_lag_smoother_map row %d names base "%s", which is not in the pruned ' ...
                 'state vector. Add it to M_.states_of_interest so that step 3 keeps it.'], ...
                rr, map{rr, 2});
        end
        lag_m(nn)          = map{rr, 3};
        base_idx(nn)       = bb;
        info.base_name{nn} = state_var_names{bb};
    end

    if verbose
        disp('Fixed-lag smoother map (eq. (3.17) denominators), from M_.fixed_lag_smoother_map:');
        for nn = 1:n_key
            if lag_m(nn) > 0
                fprintf('  %-14s = %s(-%d)  ->  denominator: Var[%s_{t+h}|y^t]\n', ...
                    states_of_interest{nn}, info.base_name{nn}, lag_m(nn), info.base_name{nn});
            else
                fprintf('  %-14s  filtered estimate (m = 0), eq. (3.17) = eq. (3.11)\n', ...
                    states_of_interest{nn});
            end
        end
    end
end


%% ========================================================================
function out = backward_window_max(x, k)
%BACKWARD_WINDOW_MAX  out(t) = max(x(max(1,t-k):t)).
%   Portable replacement for movmax(x, [k, 0]) (movmax needs MATLAB R2016a+
%   and does not exist in Octave). Used to flag any period whose backward
%   window of length k+1 touches an excluded observation.
    T = numel(x);
    out = x;
    for t = 1:T
        out(t) = max(x(max(1, t-k):t));
    end
end


%% ========================================================================
function S = dlyap_doubling(M, W)
%DLYAP_DOUBLING  Solve S = M S M' + W (rho(M) < 1) by the doubling / Smith iteration.
%   S_0 = W, S_{k+1} = S_k + M_k S_k M_k', M_{k+1} = M_k^2: after k steps S_k is
%   the sum of the first 2^k terms of sum_j M^j W M'^j, so the error decays like
%   rho(M)^(2^(k+1)) - quadratic convergence, base MATLAB/Octave only. Stops
%   when the increment is negligible relative to S or M_k has underflowed.
    S  = (W + W') / 2;
    Mk = M;
    for k = 1:200
        dS = Mk * S * Mk';
        S  = S + dS;
        S  = (S + S') / 2;
        if norm(dS, 'fro') <= 1e-14 * max(norm(S, 'fro'), realmin), break; end
        Mk = Mk * Mk;
        if ~all(isfinite(Mk(:))) || norm(Mk, 'fro') < 1e-300, break; end
    end
end
