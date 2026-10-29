function [dr, ok, info] = backward_looking_dr(M_, oo_)
%BACKWARD_LOOKING_DR  First-order decision rule of a purely backward-looking Dynare model without QZ.
%
%   [dr, ok, info] = BACKWARD_LOOKING_DR(M_, oo_)
%
%   ADDED 2026-09-28. For a model with no forward-looking variables (M_.nfwrd = 0,
%   M_.nboth = 0; HLW and its ECM/SPF variants, UC3) the linearised equations read
%       f_0 y_t + f_- y_{t-1} + f_u u_t = 0,
%   so the solution is y_t = G_y y_{t-1} + G_u u_t with G_y = -f_0^{-1} f_-, G_u = -f_0^{-1} f_u:
%   one evaluation of Dynare's generated Jacobian (+<fname>/dynamic_g1.m, sparse
%   layout: columns [y_{t-1}, y_t, y_{t+1}, u_t]) and one linear solve, instead of
%   stoch_simul's generalised Schur decomposition (about 0.2 ms against 30 ms for
%   HLW). The steady state is recomputed from the model's steadystate.m block and
%   set_auxiliary_variables.m at the current M_.params. The result is written into
%   a copy of oo_.dr in Dynare's own layout (ghx over dr.state_var, rows in
%   dr.order_var), so everything downstream (engine, augment, plots) is unchanged.
%   ok = false (dr untouched) when the model has leads, f_0 is singular, or the
%   generated files are missing - the caller then falls back to stoch_simul.
%   Exactness: test_backward_looking_dr compares ghx/ghu with stoch_simul at the
%   posterior mode and at posterior draws (agreement to ~1e-13 expected).
%
%   See also: solve_model_draw (opts / M_.solve_method = 'backward'), dr_from_abcd.
ok = false;  dr = oo_.dr;  info = '';
if ~isfield(M_, 'nfwrd') || M_.nfwrd > 0 || M_.nboth > 0, info = 'model has forward-looking variables'; return; end
n = M_.endo_nbr;  nu = M_.exo_nbr;
% (package functions +<fname>/... are not visible to exist() under Octave: every call is wrapped in try/catch)
fn_g1 = [M_.fname '.dynamic_g1'];  fn_ss = [M_.fname '.steadystate'];  fn_aux = [M_.fname '.set_auxiliary_variables'];
if ~isfield(M_, 'dynamic_g1_sparse_rowval'), info = 'M_ has no sparse Jacobian indices (Dynare < 6?)'; return; end
% steady state at the current parameters
ys = zeros(n, 1);  exo_ss = zeros(nu, 1);  params = M_.params;
try
    [ys, params, chk] = feval(fn_ss, ys, exo_ss, params);
    if chk ~= 0 || any(~isfinite(ys)), info = 'steady state failed'; return; end
catch err
    if isempty(regexp(err.message, 'undefined|Undefined|Unrecognized|not found', 'once')), info = ['steadystate failed: ' err.message]; return; end
end
try, ys = feval(fn_aux, ys, exo_ss, params); catch, end
% Jacobian at the steady state (Dynare 7 sparse representation)
try
    g1 = feval(fn_g1, [ys; ys; ys], exo_ss, params, ys, M_.dynamic_g1_sparse_rowval, M_.dynamic_g1_sparse_colval, M_.dynamic_g1_sparse_colptr);
catch err
    info = ['dynamic_g1 failed: ' err.message];  return;
end
g1 = full(g1);
if size(g1, 2) ~= 3 * n + nu, info = sprintf('unexpected Jacobian width %d (3n + nu = %d)', size(g1, 2), 3 * n + nu); return; end
f_m = g1(:, 1:n);  f_0 = g1(:, n + 1:2 * n);  f_u = g1(:, 3 * n + 1:end);
if rcond(f_0) < 1e-12, info = 'f_0 (contemporaneous Jacobian) singular'; return; end
G_y = -(f_0 \ f_m);  G_u = -(f_0 \ f_u);
if any(~isfinite(G_y(:))) || any(~isfinite(G_u(:))), info = 'non-finite solution'; return; end
% Dynare layout: rows in DR order, ghx columns = the state variables (declaration indices)
dr.ys  = ys;
dr.ghx = G_y(dr.order_var, dr.state_var);
dr.ghu = G_u(dr.order_var, :);
if isfield(dr, 'eigval'), dr.eigval = eig(G_y(dr.state_var, dr.state_var)); end
ok = true;
end
