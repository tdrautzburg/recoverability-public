function [oo_, M_] = dr_from_abcd(oo_, M_, A, B, C, D, obs_ss)
%DR_FROM_ABCD  Write a state space (A, B, C, D, obs_ss) into Dynare-style oo_.dr / M_.
%
%   [oo_, M_] = DR_FROM_ABCD(oo_, M_, A, B, C, D, obs_ss)
%
%   ADDED 2026-09-26 (generalised from Code/BauerRudebusch/ese_dr_from_abcd.m).
%   The reliability engine reads the model from Dynare's decision rule,
%       every variable = ys + ghx * s(t-1) + ghu * e(t),
%   in DR order. This function fills oo_.dr (ghx, ghu, order_var, inv_order_var,
%   state_var, ys, eigval, qz_crit) and the smoother bookkeeping the engine and the
%   diffuse filter need (oo_.bayestopt_smoother_var_list, _restrict_columns,
%   oo_.bayestopt_mf, the mirror fields on oo_.dr, M_.nstatic/npred/nspred/...,
%   M_.Sigma_e = I, M_.H = 0) from a state space in the toolbox timing
%       x_t = A x_{t-1} + B e_t,   y_t = obs_ss + C x_{t-1} + D e_t,   e_t ~ N(0, I).
%   solve_model_draw calls it on the output of a user-supplied solver, so the
%   function is safe to call once per draw on a copy of oo_/M_.
%
%   CONVENTIONS the caller must respect
%     - M_.endo_names lists the n states FIRST, in the row order of A, then the
%       p observables in the row order of C (n = size(A,1), p = size(C,1),
%       numel(M_.endo_names) = n + p): the observables are separate variables,
%       declared after the states (an observed state must be given its own
%       observable row, as in Dynare's measurement equations);
%     - M_.exo_names has q = size(B,2) entries, the columns of B and D;
%     - the DR order is [observables (static); states (purely backward)], so
%       oo_.dr.order_var = [n+1:n+p, 1:n] and oo_.dr.ys = [zeros(n,1); obs_ss];
%     - states are in deviation from their (zero) steady state; obs_ss carries
%       the level of the observables (added back by the engine and the
%       simulators).
%   M_.X_InitialDistribution is NOT touched here: the driver sets the initial
%   state distribution (X0, B0, and since 2026-09-26 the optional .type) once, and
%   a user solver may override X0/B0 per draw through solve_model_draw.
%
%   See also: solve_model_draw, estim_params_from_names, compute_state_reliability_diagnostics.

[n_states, n_shocks] = size(B);
n_obs  = size(C, 1);
n_endo = n_states + n_obs;
obs_ss = obs_ss(:);

if ~isequal(size(A), [n_states, n_states]), error('dr_from_abcd:dim', 'A must be %d x %d.', n_states, n_states); end
if ~isequal(size(C), [n_obs, n_states]),    error('dr_from_abcd:dim', 'C must be %d x %d.', n_obs, n_states); end
if ~isequal(size(D), [n_obs, n_shocks]),    error('dr_from_abcd:dim', 'D must be %d x %d.', n_obs, n_shocks); end
if numel(obs_ss) ~= n_obs,                  error('dr_from_abcd:dim', 'obs_ss must have %d entries.', n_obs); end
if ~isfield(M_, 'endo_names') || numel(M_.endo_names) ~= n_endo
    error('dr_from_abcd:names', 'M_.endo_names must list the %d states followed by the %d observables.', n_states, n_obs);
end
if ~isfield(M_, 'exo_names') || numel(M_.exo_names) ~= n_shocks
    error('dr_from_abcd:names', 'M_.exo_names must have %d entries (columns of B and D).', n_shocks);
end
if ~all(all(isfinite(A))) || ~all(all(isfinite(B))) || ~all(all(isfinite(C))) || ~all(all(isfinite(D))) || ~all(isfinite(obs_ss))
    error('dr_from_abcd:nonfinite', 'the state-space matrices contain non-finite entries.');
end

order_var     = [n_states+1:n_endo, 1:n_states]';
inv_order_var = zeros(n_endo, 1); inv_order_var(order_var) = 1:n_endo;

ghx = zeros(n_endo, n_states); ghu = zeros(n_endo, n_shocks);
ghx(1:n_obs, :)     = C;  ghu(1:n_obs, :)     = D;
ghx(n_obs+1:end, :) = A;  ghu(n_obs+1:end, :) = B;

oo_.dr.ghx           = ghx;
oo_.dr.ghu           = ghu;
oo_.dr.order_var     = order_var;
oo_.dr.inv_order_var = inv_order_var;
oo_.dr.state_var     = 1:n_states;
oo_.dr.ys            = [zeros(n_states, 1); obs_ss];
oo_.dr.eigval        = eig(A);
oo_.dr.qz_crit       = 1 + 1e-6;
oo_.steady_state     = oo_.dr.ys;
oo_.exo_steady_state = zeros(n_shocks, 1);

M_.endo_nbr = n_endo;  M_.exo_nbr = n_shocks;
M_.nstatic = n_obs;    M_.npred = n_states;  M_.nboth = 0;  M_.nfwrd = 0;
M_.nspred  = n_states; M_.nsfwd = 0;
M_.state_var = 1:n_states;
M_.Sigma_e = eye(n_shocks);
M_.H       = zeros(n_obs);
M_.Correlation_matrix    = eye(n_shocks);   % read unconditionally by Dynare 7's set_all_parameters
M_.Correlation_matrix_ME = eye(n_obs);
if ~isfield(M_, 'Skew_e'), M_.Skew_e = zeros(0, 4); end   % no skewed shocks (csn_update_specification)

% smoother bookkeeping (DR-order indices); every variable is a state or an observable
k_states = (M_.nstatic + (1:M_.nspred))';
k_obs    = inv_order_var(n_states+1:n_endo);
svl      = unique([k_states; k_obs]);
[~, restrict_columns] = ismember(k_states, svl);
[~, mf]               = ismember(k_obs, svl);
oo_.bayestopt_smoother_var_list         = svl;
oo_.bayestopt_smoother_saved_var_list   = svl;
oo_.bayestopt_smoother_restrict_columns = restrict_columns;
oo_.bayestopt_smoother_mf               = mf;
oo_.bayestopt_mf  = k_obs;  oo_.bayestopt_mf1 = k_obs;  oo_.bayestopt_mf2 = k_obs;
oo_.bayestopt_mfys = (n_states+1:n_endo)';
oo_.dr.restrict_var_list = svl;  oo_.dr.restrict_columns = restrict_columns;  oo_.dr.mf = mf;
oo_.dr.nstatic = M_.nstatic;  oo_.dr.npred = M_.npred;  oo_.dr.nboth = 0;  oo_.dr.nfwrd = 0;
oo_.dr.nspred  = M_.nspred;   oo_.dr.nsfwd = 0;
end
