function fr = fr_regression_empty(states, m_grid, H)
%FR_REGRESSION_EMPTY  NaN-filled container for the smoother-revision regression.
%
%   fr = FR_REGRESSION_EMPTY(states, m_grid, H)
%
%   states  J x 1 cellstr, base-state names (rows of every array below)
%   m_grid  1 x n_m smoother lags (third dimension)
%   H       forecast horizons h = 1..H (second dimension)
%
%   Fields (J x H x n_m unless noted; NaN = not computed):
%     .states, .h (1 x H), .m (1 x n_m)
%     .R2_theory, .beta_theory, .cov_fr_theory, .var_r_theory, .var_f_theory
%           population objects, steady-state filter (fr_regression_theory);
%           population alpha is identically 0 and is not stored
%     .omega_ratio  J x n_m, [Omega_m]_ii / Omega_ii: smoothing over filtering
%           error variance of the base state (ADDED 2026-09-17). With it the
%           empirical P^2 identifies R^2 through fr_implied_R2 (Prop. 5.1 (ii)):
%               R^2_{i,h} = P (1 - rho) / (1 - rho P),   rho = [Omega_m]_ii / Omega_ii
%     .R2_h_base    J x H, the model's eq.-(3.11) reliability R^2_{i,h} of the
%           base state = Var(r_h) / (Var(r_h) + Omega_ii) (ADDED 2026-09-17), so
%           that fr alone carries both sides of that identity
%     .R2_emp, .alpha_emp, .beta_emp, .N_emp
%           OLS on the filtered data (fr_regression_empirical); N_emp is the
%           number of usable observations
%
%   NOTATION (paper, fixed 2026-09-16). This is the PROXY regression: the
%   smoother revision f = xhat_{i,t+h|t+h+m} - xhat_{i,t+h|t} stands in for
%   the infeasible forecast error x_{i,t+h} - xhat_{i,t+h|t} on the LEFT, the
%   filter revision r = xhat_{i,t+h|t+h} - xhat_{i,t+h|t} is the regressor.
%       P^{2,m}_{i,h}      population R2         -> .R2_theory
%       Phat^{2,m}_{i,h}   OLS R2 on the data    -> .R2_emp
%       alphahat, betahat  OLS coefficients      -> .alpha_emp, .beta_emp
%                          (population values 0 and 1; .beta_theory stores 1)
%   Closed form (draft.tex Prop. 5.1 "Feasible smoother regression", eq. R2m,
%   corrected 2026-09-16), with V_h = Var[x_{i,t+h}|y^t] and U_h = V_h - Omega_ii:
%       P^{2,m}_{i,h} = (V_h - Omega_ii) / (V_h - [Omega_m]_ii)
%                     = U_h / (U_h + Omega_ii - [Omega_m]_ii),
%   population alpha = 0, beta = 1, residual variance Omega_ii - [Omega_m]_ii;
%   equals 1 at m = 0, is NON-increasing in m, and P^{2,infty} >= R^2_{i,h}; it
%   identifies the infeasible reliability through (Prop. 5.1 (ii))
%       1 - R^2_{i,h} = phi^m (1 - P^{2,m}_{i,h}) (V_h - [Omega_m]_ii) / V_h,
%       phi^m = Omega_ii / (Omega_ii - [Omega_m]_ii).
%   fr_regression_theory does NOT use this closed form: it builds Var(r),
%   Var(f), Cov(f,r) from the innovations representation and forms
%   Cov^2/(Var(f) Var(r)); on the HLW posterior mean the two agree to 2e-14
%   and .beta_theory = Cov/Var(r) is 1 to machine precision (check of
%   2026-09-16). A base state the data never reveal contemporaneously (K row
%   zero, e.g. HLW's drrstar) has Var(r) = 0 and is left NaN: the regressor is
%   degenerate, which is the h >= 1 non-degeneracy condition of the proposition.
%   Do NOT confuse it with the HINDSIGHT reliability R^{2,m}_{i,h} of eq.
%   (3.17), 1 - [Omega_m]_ii / Var(x_{i,t+h}|y^t), which has the infeasible
%   error on the left and the smoother revision on the RIGHT, equals R^2_{i,h}
%   at m = 0 and is non-DEcreasing in m: that one is LatentStateDiagnostics.R2_h
%   on the 'l<m><NAME>' rows (and CorrH there is its square root).
%
%   See also: fr_regression_theory, fr_regression_empirical,
%             compute_state_reliability_diagnostics.

J   = numel(states);
n_m = numel(m_grid);
z   = nan(J, H, n_m);
fr = struct( ...
    'states',        {reshape(states, [], 1)}, ...
    'h',             1:H, ...
    'm',             reshape(m_grid, 1, []), ...
    'R2_theory',     z, 'beta_theory', z, 'cov_fr_theory', z, ...
    'var_r_theory',  z, 'var_f_theory', z, ...
    'omega_ratio',   nan(J, n_m), 'R2_h_base', nan(J, H), ...
    'R2_emp',        z, 'alpha_emp', z, 'beta_emp', z, 'N_emp', z);
end
