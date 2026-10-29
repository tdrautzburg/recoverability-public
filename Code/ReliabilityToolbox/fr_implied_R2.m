function R2 = fr_implied_R2(P, omega_ratio)
%FR_IMPLIED_R2  Reliability R^2_{i,h} implied by a proxy-regression R^2 (Prop. 5.1).
%
%   R2 = FR_IMPLIED_R2(P, omega_ratio)
%
%   Proposition 5.1 of the paper ("Feasible smoother regression: bound and
%   identification", draft.tex) gives, for the regression of the smoother
%   revision f = xhat_{t+h|t+h+m} - xhat_{t+h|t} on the filter revision
%   r = xhat_{t+h|t+h} - xhat_{t+h|t},
%
%       P^{2,m}_{i,h} = (V_h - Omega_ii) / (V_h - [Omega_m]_ii),   V_h = Var[x_{i,t+h}|y^t],
%
%   and R^2_{i,h} = 1 - Omega_ii / V_h. Eliminating V_h,
%
%       R^2_{i,h} = P * (1 - rho) / (1 - rho * P),     rho = [Omega_m]_ii / Omega_ii,
%
%   which maps the FEASIBLE regression R^2 (population P, fr.R2_theory, or the
%   OLS estimate Phat, fr.R2_emp) into the INFEASIBLE reliability using only
%   the model-implied error-variance ratio rho (fr.omega_ratio). With P =
%   fr.R2_theory the result reproduces fr.R2_h_base exactly (identity check);
%   with P = fr.R2_emp it is the data's estimate of R^2_{i,h} under the model's
%   rho, the object plotted against the model's own R^2_{i,h}.
%
%   SHAPES. P is J x H x n_m (x D ...) as in the fr struct, or any array;
%   omega_ratio is either
%     - the same size as P (elementwise), or a scalar, or
%     - J x n_m (as fr.omega_ratio), broadcast over the horizon dimension 2
%       and any trailing dimensions of P.
%   NaN propagates. rho = 1 (no smoothing gain, m = 0) gives R2 = 0 for P < 1
%   and NaN for P = 1, both correct limits of the identity.
%
%   See also: fr_regression_theory, fr_regression_empty, fr_regression_empirical.

if isscalar(omega_ratio) || isequal(size(omega_ratio), size(P))
    rho = omega_ratio;
else
    [J, n_m] = size(omega_ratio);
    szP = size(P);
    if szP(1) ~= J || numel(szP) < 3 || szP(3) ~= n_m
        error('fr_implied_R2:shape', ...
              'omega_ratio (%d x %d) must be scalar, the size of P, or J x n_m with P of size J x H x n_m (x ...).', J, n_m);
    end
    rho = reshape(omega_ratio, [J, 1, n_m]);          % broadcast over h (dim 2) and trailing dims
end

R2 = P .* (1 - rho) ./ (1 - rho .* P);
end
