function [xerr, rev] = state_forecast_errors_and_revisions(A, state_estimate, true_states, h)
%STATE_FORECAST_ERRORS_AND_REVISIONS  h-step forecast errors/revisions for the STATES.
%
%   State-level analogue of forecast_errors_and_revisions: compares state
%   forecasts against the TRUE state paths (available in simulation
%   experiments, where the states are known).
%
%       xerr(t,:) = x(t)         - xhat(t|t-1-h)    [h-step forecast error]
%       rev(t,:)  = xhat(t|t-1)  - xhat(t|t-1-h)    [forecast revision]
%
%   Forecast identities (states are dated t, so relative to the observation
%   equation y(t) = C x(t-1) + D e(t) one extra power of A is needed):
%       xhat(t|t-1-h) = A^(h+1) * xhat(t-1-h|t-1-h)
%       xhat(t|t-1)   = A       * xhat(t-1|t-1)
%
%   INPUTS
%     A              : nx x nx  state transition matrix
%     state_estimate : T  x nx  filtered states, row t = xhat(t|t)
%     true_states    : T  x nx  true (simulated) states, row t = x(t)'
%     h              : integer >= 1, forecast horizon
%
%   OUTPUT: T x nx matrices, NaN for t < h+2.
%
%   RENAMED (2026-07): formerly f_forecast_stats_states.m (recovered from
%   Archive/Miscellaneous for the *_vs_truth variant of the diagnostics
%   engine, which was merged into compute_state_reliability_diagnostics and
%   deleted in 2026-08; this file is now called from that engine's step 6c
%   whenever its optional `true_states` argument is supplied. Original note:
%   engine). Logic unchanged; header comments corrected (the original
%   carried copy-paste artifacts from the observable-level version).

T       = size(true_states, 1);
nx      = size(A, 1);
Ahplus1 = A^(h+1);          % nx x nx, precomputed
tval    = (h+2):T;          % valid observation indices

% Baseline forecast:  xhat(t|t-1-h) = A^(h+1) * xhat(t-1-h|t-1-h)
xhat_base  = state_estimate(tval-1-h, :) * Ahplus1';   % (T-h-1) x nx

% One-step forecast:  xhat(t|t-1) = A * xhat(t-1|t-1)
xhat_1step = state_estimate(tval-1, :) * A';           % (T-h-1) x nx

% Allocate with NaN for out-of-range periods
xerr = NaN(T, nx);
rev  = NaN(T, nx);

xerr(tval, :) = true_states(tval, :) - xhat_base;
rev(tval, :)  = xhat_1step - xhat_base;
