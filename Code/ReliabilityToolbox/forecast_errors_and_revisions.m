function [yerr, rev] = forecast_errors_and_revisions(A, C, state_estimate, data, h)
% Computes h-step forecast error and forecast revision for all t.
%
% RENAMED (2026-07): formerly f_forecast_stats.m in Miscellaneous/. All call
% sites inside ReliabilityToolbox/ use the new names; the full mapping
% is in README.md ("Suggested file names").
%
%   yerr(t,:) = y(t) - yhat(t|t-1-h)           [h-step forecast error]
%   rev(t,:)  = yhat(t|t-1) - yhat(t|t-1-h)    [forecast revision]
%
% Forecast identities used:
%   yhat(t|t-1-h) = C * A^h * xhat(t-1-h|t-1-h)
%   yhat(t|t-1)   = C       * xhat(t-1|t-1)
%
% Inputs:
%   A              : nx x nx   state transition matrix
%   C              : ny x nx   measurement matrix  [B,D unused: don't enter means]
%   state_estimate : T  x nx   filtered states,  row t = xhat(t|t)
%   data           : T  x ny   observations,     row t = y(t)'
%   h              : integer >= 1, forecast horizon
%
% Output: T x ny matrices, NaN for t < h+2

T    = size(data, 1);
ny   = size(C, 1);
CAh  = C * (A^h);          % precompute  ny x nx
tval = (h+2):T;            % valid observation indices

% Baseline forecast:  yhat(t|t-1-h) = CAh * xhat(t-1-h|t-1-h)
yhat_base  = state_estimate(tval-1-h, :) * CAh';   % (T-h-1) x ny

% One-step forecast:  yhat(t|t-1)   = C   * xhat(t-1|t-1)
yhat_1step = state_estimate(tval-1,   :) * C';     % (T-h-1) x ny

% Allocate with NaN for out-of-range periods
yerr = NaN(T, ny);
rev  = NaN(T, ny);

yerr(tval, :) = data(tval, :) - yhat_base;
rev(tval,  :) = yhat_1step    - yhat_base;