function tau = first_crossing_below(acf, level)
%FIRST_CROSSING_BELOW  Half-life of estimation errors from their autocorrelation function.
%
%   tau = FIRST_CROSSING_BELOW(acf, level)
%
%   ADDED 2026-09-30. acf is n_states x n_lags (x n_draws): the autocorrelation of the filtering
%   error at lags 1..n_lags (the engine's Error_ACF / wrapper cell 3). For every state (and draw)
%       tau = min{ h >= 1 : acf(h) < level },   tau = Inf if no such h
%   with level = 1/2 the half-life of the paper. Returns n_states x n_draws (n_states x 1 for a
%   matrix input). NaN autocorrelations never count as a crossing.
%
%   See also: estimation_error_half_life (impulse-response half-lives, shock by shock), HLW_plots.

if nargin < 2 || isempty(level), level = 0.5; end
[n, H, D] = size(acf);
tau = inf(n, D);
for d = 1:D
    for i = 1:n
        h = find(acf(i, :, d) < level, 1, 'first');
        if ~isempty(h), tau(i, d) = h; end
    end
end
if H == 0, tau(:) = Inf; end
end
