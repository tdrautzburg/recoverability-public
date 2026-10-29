function [K_gain, filtered_X, log_llhood, Omega, V, Omega_tp1_t, surprises, ric_flag] = steady_state_kalman_filter(A, B, C, D, data)
%STEADY_STATE_KALMAN_FILTER  Steady-state Kalman filter for the ABCD representation.
%
%   RENAMED (2026-07): formerly my_KalmanFilter_while.m in Miscellaneous/. All call
%   sites inside ReliabilityToolbox/ use the new names; the full mapping
%   is in README.md ("Suggested file names").
%
%   State-space convention (Fernandez-Villaverde, Rubio-Ramirez, Sargent,
%   Watson, 2007 "ABCD" form, as produced from Dynare decision rules):
%
%       x(t) = A x(t-1) + B e(t)
%       y(t) = C x(t-1) + D e(t),        e(t) ~ N(0, I)
%
%   Note that the observables load on the LAGGED state; this matches
%   Dynare's ghx/ghu decision-rule representation used in
%   compute_state_reliability_diagnostics.
%
%   The function solves the filter Riccati equation (directly, via
%   solve_steady_state_riccati, since 2026-09-16; it used to iterate from 0)
%   and returns the steady-state PREDICTIVE-form gain, i.e. the K in
%
%       xhat(t|t) = A xhat(t-1|t-1) + K v(t),
%       v(t)      = y(t) - C xhat(t-1|t-1),
%
%   with K = (A Omega C' + B D') / V and V = C Omega C' + D D'.
%
%   OUTPUTS
%     K_gain       steady-state Kalman gain (predictive form, see above)
%     filtered_X   n_states x T filtered states (only if data is provided)
%     log_llhood   Gaussian log-likelihood of the data, up to the
%                  -0.5*T*n_obs*log(2*pi) constant (only if data provided)
%     Omega        steady-state filtered state covariance, P(t|t)
%     V            steady-state innovation covariance
%     Omega_tp1_t  steady-state predicted state covariance, P(t+1|t)
%     surprises    n_obs x T one-step-ahead forecast errors (if data provided)
%     ric_flag     struct describing the Riccati solve (ADDED 2026-09, fields
%                  as returned by solve_steady_state_riccati since 2026-09-16):
%                    .converged   true iff the accepted solution passed the
%                                 scaled residual check
%                    .iterations  fixed-point iterations run (0 = direct solve)
%                    .rel_step    error estimate r/(1 - rho^2) of the solution
%                    .rho_cl      rho(A - K*C)
%                    .method      'dynare' | 'qz' | 'iterate'
%
%   INPUT data (optional): n_obs x T matrix (time in COLUMNS).
%
%   Notes on changes vs. the original version:
%     - the Riccati loop now stops on convergence (norm of the update below
%       tol) instead of always running a fixed 1000 iterations; a maximum
%       iteration count is kept as a safeguard for unit-root models;
%     - fixed nobs = size(C,1) (the original used size(C,2), i.e. the
%       number of states, which mis-sized the "surprises" array and made
%       the data branch error out whenever n_states ~= n_obs);
%     - the Cholesky factor of V is computed once outside the time loop
%       (V is constant in the steady-state filter);
%     - all optional outputs are initialized, so the function no longer
%       errors when called with an empty data argument.

nstates = size(A, 1);
nobs    = size(C, 1);   % FIX: number of observables = rows of C (was size(C,2))

% Default outputs when no data is supplied
filtered_X = [];
log_llhood = [];
surprises  = [];

% -------------------------------------------------------------------------
% 1. Steady-state covariance: Riccati fixed point.
% -------------------------------------------------------------------------
% CHANGED (2026-09-16): solved directly by solve_steady_state_riccati (DARE
% via Dynare's MEX or the toolbox's QZ solver, fixed-point iteration as the
% fallback and reference) - see kalman_filter_steady_and_diffuse for the
% reason. The iteration used to start from zero here; the fallback now starts
% from B*B' (same fixed point, no singular V when D*D' is singular).
[Omega, K_gain, V, Omega_tp1_t, ric_flag] = solve_steady_state_riccati(A, B, C, D);

% ORIGINAL CODE (fixed-point iteration from Omega = 0, 2026-09 version), kept
% for comparison:
% MaxIter     = 1e5;
% tol_riccati = 1e-10;                        % RELATIVE step, not absolute
% Omega   = zeros(nstates, nstates);
% change_norm = Inf;
% nonfinite   = false;
% ii = 0;
% while change_norm > tol_riccati && ii < MaxIter
%     ii = ii + 1;
%     V        = C*Omega*C' + D*D';           % innovation covariance
%     K_gain   = (A*Omega*C' + B*D') / V;     % predictive-form gain
%     NewOmega = A*Omega*A' + B*B' - K_gain*V*K_gain';
%     NewOmega = (NewOmega + NewOmega')/2;    % drop asymmetry accumulated by roundoff
%     if ~all(isfinite(NewOmega(:)))
%         nonfinite = true;                   % keep the last finite iterate
%         break
%     end
%     change_norm = norm(NewOmega - Omega, 'fro') / ...
%         max([norm(Omega, 'fro'), norm(NewOmega, 'fro'), realmin]);
%     Omega    = NewOmega;
% end
% Omega_tp1_t = A*Omega*A' + B*B';
% V      = C*Omega*C' + D*D';
% K_gain = (A*Omega*C' + B*D') / V;
% nonfinite = nonfinite || ~all(isfinite(Omega(:))) || ~all(isfinite(K_gain(:)));
% if nonfinite, rho_cl = NaN; else, rho_cl = max(abs(eig(A - K_gain*C))); end
% ric_flag = struct('converged', ~nonfinite && (change_norm <= tol_riccati), ...
%                   'iterations', ii, 'rel_step', change_norm, 'rho_cl', rho_cl);
% if ~ric_flag.converged
%     warning('riccati:noconv', ['%s: Riccati did not converge after %d iterations ' ...
%             '(relative step %.3e, rho(A-KC) = %.6f).%s'], 'steady_state_kalman_filter', ...
%             ii, change_norm, ric_flag.rho_cl, repmat(' Iterate became non-finite.', 1, double(nonfinite)));
% end

% -------------------------------------------------------------------------
% 2. Optional: run the filter on data with the steady-state gain.
% -------------------------------------------------------------------------
if exist('data', 'var') && ~isempty(data)
    TT = size(data, 2);
    filtered_X = nan(nstates, TT);
    surprises  = nan(nobs, TT);
    current_filtered_X = zeros(nstates, 1);

    L         = chol(V, 'lower');           % V constant => factor once
    log_det_V = 2 * sum(log(diag(L)));
    log_llhood = 0;
    for tt = 1:TT
        surprises(:, tt)   = data(:, tt) - C*current_filtered_X;
        z                  = L \ surprises(:, tt);
        log_llhood         = log_llhood - 0.5*(z'*z) - 0.5*log_det_V;
        current_filtered_X = A*current_filtered_X + K_gain*surprises(:, tt);
        filtered_X(:, tt)  = current_filtered_X;
    end
end

end
