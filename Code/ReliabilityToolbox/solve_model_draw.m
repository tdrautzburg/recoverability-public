function [info, oo_, options_, M_] = solve_model_draw(M_, options_, oo_)
%SOLVE_MODEL_DRAW  Decision rule for one parameter draw: stoch_simul, or a user-supplied map.
%
%   [info, oo_, options_, M_] = SOLVE_MODEL_DRAW(M_, options_, oo_)
%
%   ADDED 2026-09-26. Every per-draw routine of the toolbox (reliability_draw_wrapper,
%   simulate_artificial_data_wrapper, artificial_data_wrapper,
%   posterior_predictive_wrapper, model_assessment, evaluate_at_reweighted_means)
%   used to call Dynare's stoch_simul directly to turn the draw written into M_ by
%   set_all_parameters into the decision rule oo_.dr from which the engine reads its
%   A, B, C, D matrices. They now all call this function, which
%
%     - by DEFAULT runs  stoch_simul(M_, options_, oo_, M_.var_list_)  exactly as
%       before (console output swallowed by evalc unless M_.verbose is true);
%
%     - when M_.solve_method = 'backward' (ADDED 2026-09-28), solves a purely
%       backward-looking model in closed form through backward_looking_dr (Dynare's
%       generated Jacobian, one linear solve, exact; ~60x faster than stoch_simul on
%       HLW / UC3) and falls back to stoch_simul if the model has leads or the solve fails;
%
%     - when M_.solve_function is set, calls a USER-SUPPLIED map from the deep
%       parameters to the state space and writes its result into oo_.dr through
%       dr_from_abcd, so that a model estimated outside Dynare (the Bauer-Rudebusch
%       ESE model, sampled in R) runs through the same pipeline with no other
%       change.
%
%   USER MAP.  M_.solve_function is the NAME (char; preferred, it survives
%   save/load and parfor) or a handle of a function
%
%       S = f(params, M_)
%
%   params is M_.params AFTER set_all_parameters has written the draw. For a model
%   outside Dynare whose whole parameter vector is estimated, params IS the draw,
%   in M_.param_names = bayestopt_.name order (estim_params_from_names builds the
%   estim_params_/bayestopt_ pair that makes set_all_parameters do that). Shock
%   standard deviations that Dynare's stderr blocks estimate live in M_.Sigma_e,
%   not in params; the M_ argument gives the map access to them and to any fixed
%   inputs the driver stored on M_ (e.g. M_.ese: maturities, portfolio weights).
%
%   S is a struct with fields
%       A (n x n), B (n x q), C (p x n), D (p x q), obs_ss (p x 1)
%   in the toolbox timing
%       x_t = A x_{t-1} + B e_t,   y_t = obs_ss + C x_{t-1} + D e_t,   e_t ~ N(0, I),
%   y in the units of M_.data. The n states are the FIRST n names of M_.endo_names
%   and the p observables (M_.observable_names) the remaining ones - see
%   dr_from_abcd, which also sets the smoother bookkeeping the engine needs.
%   Optional fields:
%       X0, B0   per-draw initial state mean (n x 1) and covariance factor (n x n,
%                P0 = B0*B0'), written into M_.X_InitialDistribution for this draw;
%                the model-constant X_InitialDistribution set by the driver is used
%                otherwise (its .type field is kept either way)
%       ok       false marks a draw the map rejects (non-stationary block, covariance
%                not positive definite, ...)
%   A rejected draw may also be signalled by returning [] (or an empty A) or by
%   throwing an error; all three give info = 1, which sends the callers down their
%   failure path (NaN cells, draw dropped downstream).
%
%   OUTPUTS mirror stoch_simul's: info (SCALAR: 0 = solved, the first element of
%   Dynare's info vector otherwise), oo_ with oo_.dr (plus, in the user path, the
%   oo_.bayestopt_* bookkeeping of dr_from_abcd), options_ (unchanged in the user
%   path), M_ (with the dimension fields of dr_from_abcd in the user path). On
%   failure the input structures are returned unchanged. oo_.dr.qz_crit is set to
%   1 + 1e-6 after every successful solve in both paths (stoch_simul rebuilds
%   oo_.dr without it; the callers also set it themselves, as before).
%
%   See also: dr_from_abcd, estim_params_from_names, reliability_draw_wrapper,
%             run_reliability_draws, compute_state_reliability_diagnostics.

info = 0;

if isfield(M_, 'solve_function') && ~isempty(M_.solve_function)
    % ---------------------------------------------------------------- user map
    f = M_.solve_function;
    if ischar(f) || isa(f, 'string')
        f = str2func(char(f));
    elseif ~isa(f, 'function_handle')
        error('solve_model_draw:solve_function', 'M_.solve_function must be a function name or handle.');
    end
    oo_in = oo_;  M_in = M_;
    try
        S  = f(M_.params, M_);
        ok = isstruct(S) && isfield(S, 'A') && ~isempty(S.A);
        if ok && isfield(S, 'ok') && ~isempty(S.ok), ok = logical(S.ok(1)); end
        if ~ok
            info = 1;
            oo_ = oo_in;  M_ = M_in;
            return
        end
        [oo_, M_] = dr_from_abcd(oo_, M_, S.A, S.B, S.C, S.D, S.obs_ss);
        if isfield(S, 'X0') && ~isempty(S.X0), M_.X_InitialDistribution.X0 = S.X0(:); end
        if isfield(S, 'B0') && ~isempty(S.B0), M_.X_InitialDistribution.B0 = S.B0;    end
    catch ME
        if isfield(M_in, 'verbose') && M_in.verbose
            fprintf(2, 'solve_model_draw: user map failed: %s\n', ME.message);
        end
        info = 1;
        oo_ = oo_in;  M_ = M_in;
    end
    return
end

% ------------------------------------------- closed form for backward-looking models (opt-in)
% ADDED 2026-09-28: M_.solve_method = 'backward' replaces stoch_simul by backward_looking_dr
% (Dynare's generated Jacobian + one linear solve; exact for models without leads, tested
% bit-for-bit against stoch_simul on HLW and UC3 by test_backward_looking_dr; 0.2-0.3 ms
% instead of 12-18 ms). Models with forward-looking variables, or any failure, fall
% through to stoch_simul below. Default (field absent or 'stoch_simul'): unchanged.
if isfield(M_, 'solve_method') && ~isempty(M_.solve_method) && strcmpi(char(M_.solve_method), 'backward')
    [dr_bw, ok_bw, why_bw] = backward_looking_dr(M_, oo_);
    if ok_bw
        oo_.dr = dr_bw;
        oo_.dr.qz_crit = 1 + 1e-6;
        oo_.steady_state = dr_bw.ys;
        return
    elseif isfield(M_, 'verbose') && M_.verbose
        fprintf('solve_model_draw: backward-looking solver not used (%s); stoch_simul instead\n', why_bw);
    end
end

% ------------------------------------------------------------------ Dynare path
if ~isfield(M_, 'var_list_')
    error('solve_model_draw:var_list', 'M_.var_list_ is missing (set it in the .mod epilogue) and M_.solve_function is not set.');
end
oo_in = oo_;  options_in = options_;  M_in = M_;
try
    if isfield(M_, 'verbose') && M_.verbose
        [info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, M_.var_list_);
    else
        % evalc swallows stoch_simul's console output
        [~, info, oo_, options_, M_] = evalc('stoch_simul(M_, options_, oo_, M_.var_list_);');
    end
catch
    info = 1;
end
if isempty(info), info = 0; end
info = info(1);
if info ~= 0
    oo_ = oo_in;  options_ = options_in;  M_ = M_in;
else
    % stoch_simul rebuilds oo_.dr and drops the unit-root threshold the diffuse
    % filter reads; the callers re-set it after every solve (same value), and
    % since 2026-09-26 it is set here as well so that no caller can forget it.
    oo_.dr.qz_crit = 1 + 1e-6;
end
end
