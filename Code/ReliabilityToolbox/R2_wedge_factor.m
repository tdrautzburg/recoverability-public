function [f_wedge, details] = R2_wedge_factor(A, B, C, D, h, filter_cache)
%R2_WEDGE_FACTOR  Adjustment factor f_{i,j}(h) of eq. (26), states x observables.
%
%   f_wedge = R2_WEDGE_FACTOR(A, B, C, D, h) computes, at horizon h
%   (default 1), the matrix of adjustment factors linking the theoretical
%   R^2 of state i to the theoretical R^2 of observable j:
%
%       f_{i,j}(h) = [RevVar_x(i) / RevVar_y(j)] * [FEV_y(j) / FEV_x(i)],
%
%   where RevVar is the h-period forecast-revision variance
%   Omega_{h-1} = sum_{k=0}^{h-1} A^k (K V K') A^k'  (state block; C-loaded
%   for the observables) and FEV the corresponding h-step forecast-error
%   variance. Rows index states, columns index observables.
%
%   [f_wedge, details] = R2_WEDGE_FACTOR(...) also returns the intermediate
%   objects (revision/forecast variances, per-variable R^2's, K, Sigma, V).
%
%   OPTIONAL 6th ARGUMENT filter_cache: a struct with fields K_gain, Sigma,
%   Sigma_v (the outputs of steady_state_kalman_filter(A,B,C,D): gain,
%   filtered state covariance, innovation covariance). Since the filter
%   does not depend on h, callers that loop over horizons should compute it
%   once and pass it here to skip the Riccati iteration on every call.
%
%   RENAMED (2026-07): formerly compute_wedge_f.m (recovered from
%   Archive/Miscellaneous for the *_vs_truth variant of the diagnostics
%   engine, merged into compute_state_reliability_diagnostics and deleted in
%   2026-08; this file now produces that engine's .fh_direct. Original note:
%   engine). Logic unchanged apart from the optional filter_cache.
%
%   See also: theoretical_R2_wedge (the "newdeltaR2.pdf" variant of the
%   same object). NAMING: compute_state_reliability_diagnostics stores THIS
%   function's output as .fh_direct and theoretical_R2_wedge's as .fh. The
%   deleted *_vs_truth engine used the opposite convention (.fh here,
%   .fh_shift there), so old scripts must be re-read.
%   See also: steady_state_kalman_filter.

narginchk(4, 6);

if nargin < 5 || isempty(h)
    h = 1;
end

validateattributes(A, {'numeric'}, {'2d', 'square'}, mfilename, 'A', 1);
validateattributes(B, {'numeric'}, {'2d'}, mfilename, 'B', 2);
validateattributes(C, {'numeric'}, {'2d'}, mfilename, 'C', 3);
validateattributes(D, {'numeric'}, {'2d'}, mfilename, 'D', 4);
validateattributes(h, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, mfilename, 'h', 5);

nstates = size(A, 1);
nobs    = size(C, 1);
nshocks = size(B, 2);

if size(B, 1) ~= nstates
    error('R2_wedge_factor:DimensionMismatch', ...
        'B must have the same number of rows as A.');
end
if size(C, 2) ~= nstates
    error('R2_wedge_factor:DimensionMismatch', ...
        'C must have the same number of columns as A has states.');
end
if size(D, 1) ~= nobs || size(D, 2) ~= nshocks
    error('R2_wedge_factor:DimensionMismatch', ...
        'D must be size(C,1)-by-size(B,2).');
end

% Steady-state Kalman filter objects (compute, or reuse the caller's cache)
if nargin >= 6 && ~isempty(filter_cache)
    K_gain  = filter_cache.K_gain;
    Sigma   = filter_cache.Sigma;
    Sigma_v = filter_cache.Sigma_v;
else
    [K_gain, ~, ~, Sigma, Sigma_v] = steady_state_kalman_filter(A, B, C, D);
end

% Revision covariance accumulated over h periods:
%   Omega_{h-1} = sum_{k=0}^{h-1} A^k (K V K') A^k'
Omega_hminus1      = zeros(nstates, nstates);
A_power            = eye(nstates);
shock_revision_cov = K_gain * Sigma_v * K_gain';

terms_remaining = h;
while terms_remaining > 0
    Omega_hminus1 = Omega_hminus1 + A_power * shock_revision_cov * A_power';
    A_power = A * A_power;
    terms_remaining = terms_remaining - 1;
end

state_revision_var = diag(Omega_hminus1);
obs_revision_var   = diag(C * Omega_hminus1 * C');

state_forecast_var = diag(Sigma + Omega_hminus1);
obs_forecast_var   = diag(Sigma_v + C * Omega_hminus1 * C');

delta   = bsxfun(@rdivide, state_revision_var, obs_revision_var.');
nu      = bsxfun(@rdivide, obs_forecast_var.', state_forecast_var);
f_wedge = delta .* nu;

if nargout > 1
    details = struct();
    details.h = h;
    details.K_gain = K_gain;
    details.Sigma = Sigma;
    details.Sigma_v = Sigma_v;
    details.Omega_hminus1 = Omega_hminus1;
    details.delta = delta;
    details.nu = nu;
    details.state_revision_var = state_revision_var;
    details.obs_revision_var = obs_revision_var;
    details.state_forecast_var = state_forecast_var;
    details.obs_forecast_var = obs_forecast_var;
    details.R2_state = state_revision_var ./ state_forecast_var;
    details.R2_obs = obs_revision_var ./ obs_forecast_var;
end
end
