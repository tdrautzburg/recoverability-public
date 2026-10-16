function output_cell = reliability_draw_wrapper(xparam1, original_M_, original_options_, original_oo_, estim_params_, bayestopt_, dataset_, dataset_info) %#ok<INUSD>
%RELIABILITY_DRAW_WRAPPER  Per-draw reliability diagnostics for prior_function/posterior_function.
%
%   RENAMED (2026-07): formerly f_wrapper_Extract_ABCD_par_while.m in Miscellaneous/. All call
%   sites inside ReliabilityToolbox/ use the new names; the full mapping
%   is in README.md ("Suggested file names").
%
%   This function is the entry point that Dynare's
%       prior_function(function='reliability_draw_wrapper', ...)
%       posterior_function(function='reliability_draw_wrapper', ...)
%   commands call once for every parameter draw. For each draw it:
%
%     1. writes the drawn parameter vector into the model,
%     2. re-solves the model to get the decision rules - with stoch_simul, or
%        (CHANGED 2026-09-26) with the user-supplied map M_.solve_function when
%        the model is not a Dynare model; both through solve_model_draw,
%     3. calls compute_state_reliability_diagnostics, which builds the ABCD state-space
%        representation, runs the (diffuse/heteroskedastic) Kalman filter and
%        computes the reliability diagnostics used in the paper
%        (correlation with the truth by horizon, half-lives of estimation
%        errors, forecast-error variance shares, state uncertainty, ...).
%
%   REQUIRED EXTRA FIELDS ON M_ (set in the .mod file before the
%   prior_function/posterior_function call):
%       M_.observable_names   e.g. {'yobs';'piobs';'uobs'}
%       M_.states_of_interest e.g. {'ygap';'ustar';'ystar';'pistar'}
%       M_.CellList_States    optional (ADDED 2026-09-15): J x 3 list
%                             {name, smoother lags, change horizons}; when
%                             present it is what the engine diagnoses, and
%                             M_.states_of_interest should be set to
%                             expand_states_of_interest(M_.CellList_States)
%       M_.var_list_          the var_list_ from the .mod file
%       M_.data               T x n_obs data matrix (its presence switches on
%                             the data-based/filtering outputs, see below)
%       M_.X_InitialDistribution.X0 / .B0   initial state mean / std factor
%                             (.type 'diffuse'|'proper', see the engine)
%       M_.solve_function     optional (ADDED 2026-09-26): user map from deep
%                             parameters to (A,B,C,D,obs_ss) replacing
%                             stoch_simul, see solve_model_draw
%       M_.verbose            optional, default true; false suppresses output
%
%   INPUTS (fixed by the Dynare prior/posterior_function interface):
%       xparam1        vector of estimated parameters for this draw, ordered
%                      as in the estimated_params block (= bayestopt_.name)
%       original_M_    model structure (a private copy for this draw)
%       original_options_, original_oo_, estim_params_, bayestopt_,
%       dataset_, dataset_info  : usual Dynare structures
%
%   OUTPUT: output_cell, a 1x20 cell array. Cells 13-16 are filled only when
%   M_.data exists and are left empty ([]) otherwise; the WIDTH is 20 either
%   way (see "CONSTANT WIDTH" below). Cell {20} (ADDED 2026-09-15) is
%   LatentStateDiagnostics.fr, the smoother-revision regression struct
%   (fr_regression_empty lists its fields); M_.states_of_interest may now be
%   the J x 3 spec described in compute_state_reliability_diagnostics, in
%   which case the rows of cells 1-8, 10-12, 17-18 follow
%   LatentStateDiagnostics.state_labels rather than M_.states_of_interest.
%       { 1} max_HL        : worst-case half-life of estimation errors, per state
%       { 2} CorrH         : corr(state, filtered estimate) by horizon
%       { 3} Error_ACF     : autocorrelation of estimation errors by horizon
%       { 4} weighted_HL   : half-life weighted by shock contributions
%       { 5} all_HL        : half-life by shock (shocks x states)
%       { 6} sqrt diag KF_Omega        : filtered state std, t|t
%       { 7} sqrt diag KF_Omega_tp1_t  : predicted state std, t+1|t
%       { 8} FEV           : 1 - relative forecast-error variance path (or NaN on failure)
%       { 9} params        : parameter vector actually used for this draw
%       -- eq. (3.17) objects, filled with or without data --
%       {10} R2_h          : n_key_states x MaxH, the eq.-(3.17) reliability
%                            R2^(m)_{i,h} = 1 - [Omega_m]_ii / Var[x_{i,t+h}|y^t]
%                            (ADDED 2026-08). This is the R2 to report when a
%                            table compares FILTERED and fixed-lag SMOOTHED
%                            estimates, because both are normalised by the
%                            conditional variance of the same contemporaneous
%                            state at the same horizon, which is what makes
%                            reliability weakly increasing in the hindsight m
%                            (Section 3.5). Since 2026-09-15 every row equals
%                            cell {2} squared: CorrH of a smoothed state
%                            'l<m><NAME>' is sqrt(R2_h), the conditional
%                            correlation corr(x_{t+h}, xhat_{t+h|t+h+m} | y^t).
%       {11} smoother_lag_m: n_key_states x 1, the m of eq. (3.17), 0 for a
%                            filtered estimate. Says which rows of {10} the
%                            correction actually bites on (ADDED 2026-08).
%       {12} smoother_base : n_key_states x 1 cell of names, the
%                            contemporaneous state whose conditional variance
%                            (3.10) is the denominator of that row of {10}
%                            (the state's own name when m = 0). Lets a
%                            plotting script pair each smoother column with
%                            its filter column without parsing names
%                            (ADDED 2026-08).
%            REUSED SLOTS.  Until 2026-08 cells 10-12 held the steady-state-
%            gain filter's states, surprises and log-likelihood; that data
%            pass was removed and the cells were left EMPTY as placeholders so
%            that 13-17 would not shift. The eq.-(3.17) objects now occupy
%            them, which keeps the width at 17 and every other index where it
%            was. A script written against the pre-2026-08 wrapper that still
%            reads 10-12 will now get these instead of empties - it was
%            already broken, but it will now fail loudly rather than quietly.
%       -- the following only when M_.data exists --
%       {13} TV_state_estimates : filtered states (diffuse + time-varying filter)
%       {14} TV_surprises       : forecast errors, time-varying filter
%       {15} TV_log_llhood      : log-likelihood, time-varying filter
%       {16} TV state std paths : sqrt of predictive variances P(t|t-1) of the
%                                 states of interest (n_key_states x T)
%       -- population diagnostic again, filled with or without data --
%       {17} nu            : tail-unreliability index, n_key_states x 1,
%                            ordered like M_.states_of_interest; NaN for states
%                            without a stochastic trend (ADDED 2026-08, see
%                            tail_unreliability_index.m). REVISED 2026-09 to the
%                            general-integration-order formula, so the UNIT is
%                            quarters^(2d+1) with d the matching entry of {18} -
%                            quarters exactly where that d is 0, which is every
%                            state of interest in HLW and UC3
%       {19} riccati_flag  : 1x4 numeric row [converged, iterations, rel_step,
%                            rho_cl] for the steady-state Riccati solve
%                            (ADDED 2026-09), or nan(1,4) for a failed draw.
%                            converged is 0/1. Everything derived from the
%                            steady-state covariance - cells {1}-{8} and {10}
%                            - is only trustworthy when it is 1, so screen on
%                            it before taking medians across draws. rho_cl is
%                            rho(A-K*C); at 1 - 1e-5 or closer the draw is a
%                            sigma_g/sigma_z pile-up. CHANGED 2026-09-16: the
%                            equation is solved directly (solve_steady_state_
%                            riccati), so iterations is 0 for a direct solve
%                            (> 0 only for the fixed-point fallback) and
%                            rel_step is the relative ERROR estimate
%                            r/(1-rho_cl^2) of the accepted solution (residual
%                            check, tolerance 1e-8), no longer the last step.
%       {18} nu_order      : the integration orders d of {17}, n_key_states x 1,
%                            NaN wherever {17} is (ADDED 2026-09). Carried per
%                            draw because it is what fixes the units of {17}: a
%                            median over draws is only meaningful within a
%                            constant d, and the structural order can in
%                            principle move if a draw sends a chain coefficient
%                            to zero
%
%   POSTERIOR/PRIOR DRAWS OF nu.  Cell 17 is what makes the posterior of the
%   tail index available: after the posterior_function call, collect it across
%   draws exactly as the plotting scripts do for the other diagnostics,
%
%       nu_draws = cat(2, oo_.posterior_function_results{:, 17});   % n_key x n_draws
%       nu_prior = cat(2, oo_.prior_function_results{:,     17});
%       nu_ord   = cat(2, oo_.posterior_function_results{:, 18});   % same shape
%
%   and then take medians/quantiles across the second dimension - after
%   checking that each row of nu_ord is constant across draws, since rows with
%   different d are in different units and a median across them is meaningless. Failed draws
%   contribute a scalar NaN in every cell and are dropped by the usual
%   isscalar/isreal masks on columns 2-3 before this concatenation; note that
%   cat(2, ...) requires the surviving cells to have equal row counts, which
%   is why the failure path below writes a scalar NaN rather than a vector.
%   Genuine NaNs INSIDE a successful draw's nu (stationary states such as the
%   output gap) are expected and should be reported as the paper's
%   unconditional R2 instead, not imputed.
%
%   CONSTANT WIDTH.  Dynare's execute_prior_posterior_function sizes its
%   storage from the FIRST draw, results_cell = cell(n_draws, columns(junk)),
%   and then assigns results_cell(i,:) = wrapper(...). Every draw must
%   therefore return the SAME number of cells, successes and failures alike,
%   or the sampling loop errors out part-way through. That is why nu sits at a
%   fixed index 17, nu_order at 18, and the eq.-(3.17) objects at fixed
%   indices 10-12, in both
%   the data and the no-data branch (and why the no-data branch pre-allocates
%   the full width, leaving 13-16 empty) rather than being packed at the first
%   free slot of each branch. Keeping every OTHER index unchanged matters
%   because the downstream plotting scripts address these results positionally
%   (e.g. cat(3, oo_.posterior_function_results{:, 8}) for the FEV).
%
%   On solver failure (indeterminacy/BK violation for a draw) the cells are
%   filled with NaN so that the draw can be discarded downstream.
%
%   See also: compute_state_reliability_diagnostics, execute_prior_posterior_function.

% -------------------------------------------------------------------------
% 1. Apply the parameter draw to THIS DRAW'S copy of the model structure.
% -------------------------------------------------------------------------
% NOTE (2026-07, corrected assessment): Dynare's own
% execute_prior_posterior_function already applies each draw via
%     M_ = set_all_parameters(parameter_mat(:,i), estim_params_, M_)
% before calling this wrapper (verified in Dynare 7.1 and in the 2015-era
% version used for the paper), so original_M_ arrives here ALREADY carrying
% the draw and the line below is, under stock Dynare, a numerically
% neutral re-application. The original code instead looped over
%     set_param_value(est_param_names{ii}, xparam1(ii))
% which writes into the *global* M_ (never used by this function) - a
% redundant no-op that did NOT corrupt results, but that (a) depended on
% oo_.posterior_median existing (crashes for prior_function without a
% prior MH run), (b) silently mutated global state, and (c) would break
% under parallel execution or under any Dynare version that stops
% pre-applying the draw. The explicit set_all_parameters call below keeps
% the wrapper correct on its own, independent of Dynare internals.
original_M_ = set_all_parameters(xparam1, estim_params_, original_M_);

% Which states to diagnose (ADDED 2026-09-15): the J x 3 list
% M_.CellList_States = {name, smoother lags, change horizons} when the .mod
% file defines it, else the legacy M_.states_of_interest name list. With the
% list, M_.states_of_interest is expected to hold
% expand_states_of_interest(M_.CellList_States) - the diagnosed states in the
% row order of every cell below - which is what the plotting scripts index.
if isfield(original_M_, 'CellList_States')
    state_list = original_M_.CellList_States;
else
    state_list = original_M_.states_of_interest;
end

% The presence of a data matrix decides whether we only compute
% population/theoretical diagnostics (8 outputs) or also run the filters on
% the data (16 outputs).
DO_ESTIMATE = isfield(original_M_, 'data');
if ~isfield(original_M_, 'verbose')
    original_M_.verbose = true;
end

% -------------------------------------------------------------------------
% 2. Re-solve the model for this draw.
% -------------------------------------------------------------------------
% info ~= 0 signals a solver failure (e.g. Blanchard-Kahn violation, or a
% draw the user map rejects) for this particular draw; we then return NaNs
% instead of crashing the whole prior/posterior sampling loop.
% CHANGED 2026-09-26: stoch_simul (verbose or under evalc) or the user map
% M_.solve_function, both inside solve_model_draw, which returns a scalar info.
[info, new_oo_, new_options_, new_M_] = solve_model_draw(original_M_, original_options_, original_oo_); %#ok<ASGLU>

% -------------------------------------------------------------------------
% 3. Pre-allocate the output cell.
% -------------------------------------------------------------------------
% CHANGED (2026-08): the width is 17 in BOTH branches. It used to be 8 without
% data (grown to 9 by the params assignment below) and 16 with data; nu is
% appended at index 17 and the eq.-(3.17) objects fill the placeholder slots
% 10-12, in both cases so that a single, fixed column holds each whichever
% branch ran, and so that Dynare's fixed-width storage of the draws stays
% valid (see "CONSTANT WIDTH" in the header). Cells 13-16 stay empty in the
% no-data branch, exactly as they already did on the failure path.
% CHANGED (2026-09): the width is 18. nu_order is APPENDED at 18 rather than
% inserted next to nu, so that every existing index - 17 included - keeps its
% position and the positional reads in the plotting scripts are unaffected.
% CHANGED (2026-09): the width is 19. The Riccati convergence flag is APPENDED
% at 19, again leaving every existing index untouched.
% CHANGED (2026-09-15): the width is 20. The smoother-revision regression
% struct LatentStateDiagnostics.fr is APPENDED at 20 (scalar NaN on failure).
N_OUTPUT_CELLS = 20;
output_cell = cell(1, N_OUTPUT_CELLS);

if info == 0
    % ---------------------------------------------------------------------
    % 4. Successful solve: run the diagnostics engine.
    % ---------------------------------------------------------------------
    output_cell{1, 9} = new_M_.params;   % record the params actually used

    if DO_ESTIMATE
        % qz_crit is used downstream to classify unit roots when building
        % the diffuse initialization (Pinf/Pstar).
        new_oo_.dr.qz_crit = 1 + 1e-6;

        [LatentStateDiagnostics, ModelMatrices, obs_indices, key_state_indices, ...
            TV_state_estimates, TV_surprises, TV_log_llhood] = ...
            compute_state_reliability_diagnostics(original_M_.observable_names, state_list, new_oo_, new_M_, original_M_.data); %#ok<ASGLU>

        % Cells 10-12 held the steady-state-gain filter's state estimates,
        % surprises and log-likelihood until that data pass was removed in
        % 2026-08 (see compute_state_reliability_diagnostics, step 6). They are
        % now the eq.-(3.17) objects, assigned below with the other population
        % diagnostics because they do not depend on the data.
        output_cell{1, 13} = TV_state_estimates;
        output_cell{1, 14} = TV_surprises;
        output_cell{1, 15} = TV_log_llhood;

        % Time paths of the predictive standard deviations of the states of
        % interest: Sequence.Omega(:,:,t) holds P(t|t-1), so we extract the
        % diagonal of the (key_state x key_state) block for t = 1..T using
        % linear indexing (avoids a loop over T).
        [T, m] = size(TV_state_estimates);
        sub = LatentStateDiagnostics.TV_KF_Sequence.Omega(key_state_indices, key_state_indices, :);
        diagIdx   = (1:m) + (0:m-1)*m;            % diagonal positions within one m x m slice
        sliceSize = m * m;
        allIdx    = diagIdx' + (0:T-1) * sliceSize; % m x T linear indices across slices
        output_cell{1, 16} = sqrt(reshape(sub(allIdx), m, T));
    else
        new_oo_.dr.qz_crit = 1 + 1e-6;
        [LatentStateDiagnostics, ModelMatrices, obs_indices, key_state_indices] = ...
            compute_state_reliability_diagnostics(original_M_.observable_names, state_list, new_oo_, new_M_); %#ok<ASGLU>
    end

    % Population diagnostics (available with or without data)
    output_cell{1, 1} = LatentStateDiagnostics.max_HL;
    output_cell{1, 2} = LatentStateDiagnostics.CorrH;
    output_cell{1, 3} = LatentStateDiagnostics.Error_ACF;
    output_cell{1, 4} = LatentStateDiagnostics.weighted_HL;
    output_cell{1, 5} = LatentStateDiagnostics.all_HL;
    % Kalman-filter state uncertainty, filtered (t|t) and predicted (t+1|t)
    output_cell{1, 6} = sqrt(diag(ModelMatrices.KF_Omega(key_state_indices, key_state_indices)));
    output_cell{1, 7} = sqrt(diag(ModelMatrices.KF_Omega_tp1_t(key_state_indices, key_state_indices)));
    output_cell{1, 8} = LatentStateDiagnostics.FEV;

    % Tail-unreliability index nu (ADDED 2026-08, REVISED 2026-09): n_key_states
    % x 1, ordered like M_.states_of_interest, in quarters^(2d+1) with d the
    % integration order in cell 18. It is a population object (a closed-form
    % function of the steady-state filtered variance and of the generalised
    % long-run innovation variance of the state's top cumulation rung), so it is
    % returned for every successful draw whether or not data was supplied.
    % Collecting this cell across draws is what delivers the posterior of nu -
    % see "POSTERIOR/PRIOR DRAWS OF nu" in the header. Cell 18 travels with it
    % because a median over draws is only in a single unit within a constant d.
    output_cell{1, 17} = LatentStateDiagnostics.nu;
    output_cell{1, 18} = LatentStateDiagnostics.nu_order;

    % Riccati convergence flag (ADDED 2026-09). Stored as a fixed-width 1x4
    % NUMERIC row rather than the struct returned by the filter, so that the
    % failure path can write nan(1,4) and keep the constant-width contract,
    % and so the whole posterior collapses to a matrix in one call:
    %     ric = cat(1, oo_.posterior_function_results{:, 19});   % n_draws x 4
    %     mean(ric(:,1) == 0)     % share of draws that did NOT converge
    %     max(ric(ric(:,1)==1,4)) % worst rho(A-KC) among converged draws
    rf = ModelMatrices.KF_ric_flag;
    output_cell{1, 19} = [double(rf.converged), rf.iterations, rf.rel_step, rf.rho_cl];

    % Smoother-revision regression (ADDED 2026-09-15): the whole struct
    % LatentStateDiagnostics.fr (population R2/beta under the steady-state
    % filter, OLS R2/alpha/beta on the TV-filtered data, per base state x
    % horizon x smoother lag; see fr_regression_empty). One cell rather than
    % one per array so the width contract is a single new index. Collect it as
    %     FR = [oo_.posterior_function_results{:, 20}];   % struct array over draws
    %     R2 = cat(4, FR.R2_emp);                          % J x H x n_m x n_draws
    % after dropping failed draws (scalar NaN in this cell).
    output_cell{1, 20} = LatentStateDiagnostics.fr;

    % Eq. (3.17) (ADDED 2026-08): the reliability measure that may be compared
    % across information sets, so that the fixed-lag SMOOTHER columns of the
    % summary table are weakly more reliable than the FILTER column, as
    % Section 3.5 requires. Like nu, these are population objects, returned for
    % every successful draw whether or not data was supplied. Cell 2 (CorrH) is
    % deliberately left alone: it is eq. (3.11), it is what the h-profile,
    % histogram and Bayes-factor figures of the filtered states are built from,
    % and for those states the two measures coincide anyway. Cell 11 records
    % which rows differ.
    output_cell{1, 10} = LatentStateDiagnostics.R2_h;
    output_cell{1, 11} = LatentStateDiagnostics.smoother_lag_m;
    output_cell{1, 12} = LatentStateDiagnostics.smoother_base;

    fprintf('. \n')   % progress marker: one dot per successful draw
else
    % ---------------------------------------------------------------------
    % 5. Failed solve: return NaNs so the draw can be filtered out later.
    % ---------------------------------------------------------------------
    output_cell{1, 1} = nan(1);
    output_cell{1, 2} = nan(1);
    output_cell{1, 3} = nan(1);
    output_cell{1, 4} = nan(1);
    output_cell{1, 5} = nan(1);
    output_cell{1, 6} = nan(1);
    output_cell{1, 7} = nan(1);
    output_cell{1, 8} = nan(1);
    output_cell{1, 9} = xparam1;   % keep the rejected draw for reference
    % Eq. (3.17) cells for a rejected draw: a SCALAR NaN in 10-11 and an empty
    % cell in 12, so the draw stays detectable by the same isscalar test the
    % plotting scripts apply to columns 2-3 (the draw failed before the state
    % space was built, so n_key_states is not known here).
    output_cell{1, 10} = nan(1);
    output_cell{1, 11} = nan(1);
    output_cell{1, 12} = {};
    % cells 13-16 are left empty ([]), as in the original implementation
    % nu for a rejected draw: a SCALAR NaN, matching cells 1-8. Deliberately
    % not nan(n_key_states,1) - the draw failed before the state space was
    % built, so n_key_states is not known here, and the scalar keeps the
    % failure detectable by the same isscalar test the plotting scripts apply
    % to columns 2-3.
    output_cell{1, 17} = nan(1);
    output_cell{1, 18} = nan(1);   % nu_order, same scalar-NaN convention as 17
    % Riccati flag: nan(1,4), NOT a scalar, so that the cat(1, ...) idiom in
    % the header still returns an n_draws x 4 matrix with the failed draws as
    % all-NaN rows. A failed draw never reached the Riccati solve.
    output_cell{1, 19} = nan(1, 4);
    output_cell{1, 20} = nan(1);   % fr struct: scalar NaN, same convention as 1-8
    fprintf(2, 'x \n');  % progress marker: one red x per failed draw
end

end
