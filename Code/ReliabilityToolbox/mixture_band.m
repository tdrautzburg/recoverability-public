function [q_lo, q_hi, sd_eq] = mixture_band(mu, sig, p_lo, p_hi)
% MIXTURE_BAND  Percentile band of a Gaussian mixture over posterior draws, per period.
%
%   [q_lo, q_hi, sd_eq] = mixture_band(mu, sig)              2.5% and 97.5%
%   [q_lo, q_hi, sd_eq] = mixture_band(mu, sig, p_lo, p_hi)
%
%   mu, sig   T x n: per-draw filtered mean and filter standard deviation of one
%             state (draw d contributes N(mu(t,d), sig(t,d)^2) at period t)
%   q_lo,q_hi T x 1: quantiles of the equal-weight mixture, found by bisection
%             on F(x) = mean_d Phi((x - mu_d)/sig_d) (Phi via erfc, no toolbox)
%   sd_eq     T x 1: (q_hi - q_lo) / (2 * 1.96), the standard deviation of the
%             Gaussian with the same central band width
%
%   ADDED 2026-09-20. Replaces mean +/- 1.96 * sqrt(mean_d sig_d^2 + var_d mu_d)
%   throughout HLW_plots: that Gaussian "overall" band is driven by the few
%   pile-up draws whose steady-state filter variance is enormous (on the HLW
%   baseline the top 5% of posterior draws carried 80% of the mean variance,
%   sqrt(mean var) 3.8pp against a median of 1.5pp), so the mixture percentile is
%   the honest 95% band. Draws with NaN mean or sd are dropped period by period.
if nargin < 3 || isempty(p_lo), p_lo = 0.025; end
if nargin < 4 || isempty(p_hi), p_hi = 0.975; end
[T, n] = size(mu);
% squeeze() of a T x 1 x n slice and of a 1 x T x n slice differ when n is 0 or 1
% (Octave: 1 x T x 0 -> 1 x T x 0 stays 3-D, 1 x T x 1 -> 1 x T): accept the
% transposed / degenerate shapes (2026-09-28, proper-prior BR tracking panels)
if ~isequal(size(sig), [T, n])
    if isequal(size(sig'), [T, n]),                sig = sig';
    elseif numel(sig) == 0 && n == 0,              sig = zeros(T, 0);
    elseif numel(sig) == T * n,                     sig = reshape(sig, T, n);
    else, error('mixture_band: mu and sig must both be T x n.');
    end
end
if n == 0                                        % empty draw set (e.g. no draw meets a tracking-plot R^2 threshold): NaN band, as mean([]) would give
    q_lo = nan(T, 1);  q_hi = q_lo;  sd_eq = q_lo;
    return
end
ok  = isfinite(mu) & isfinite(sig) & sig >= 0;
sig = max(sig, 1e-12);
mu(~ok) = 0;  sig(~ok) = 1;                      % neutralised below through the weights
w   = double(ok);  nw = sum(w, 2);  nw(nw == 0) = NaN;
Phi = @(z) 0.5 * erfc(-z / sqrt(2));
F   = @(x) sum(w .* Phi((x - mu) ./ sig), 2) ./ nw;       % T x 1 mixture cdf at x (T x 1)
lo0 = min(mu - 8 * sig .* w - 1e6 * (1 - w), [], 2);      % excluded draws pushed out of the way
hi0 = max(mu + 8 * sig .* w + 1e6 * (1 - w), [], 2);
lo0 = min(lo0, min(mu, [], 2)); hi0 = max(hi0, max(mu, [], 2));
q_lo = bisect(F, p_lo, lo0, hi0);
q_hi = bisect(F, p_hi, lo0, hi0);
sd_eq = (q_hi - q_lo) / (2 * 1.959963984540054);
end

function q = bisect(F, p, a, b)
% vectorised bisection over periods: F is monotone in x for each row
for it = 1:70
    c  = 0.5 * (a + b);
    Fc = F(c);
    up = Fc < p;                                   % root lies above c
    a(up)  = c(up);
    b(~up) = c(~up);
end
q = 0.5 * (a + b);
end
