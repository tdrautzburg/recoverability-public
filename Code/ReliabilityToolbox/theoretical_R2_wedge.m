function out = theoretical_R2_wedge(A, B, C, D, horizons, varargin)
%THEORETICAL_R2_WEDGE Compute the wedge objects from newdeltaR2.pdf.
%
% RENAMED (2026-07): formerly compute_newdeltaR2_wedge.m in Miscellaneous/. All call
% sites inside ReliabilityToolbox/ use the new names; the full mapping
% is in README.md ("Suggested file names").
%
% This file matches the notation in newdeltaR2.pdf, which uses
%
%   x_t = A x_{t-1} + B epsilon_t
%   y_t = C x_{t-1} + D epsilon_t
%   e_t ~ N(0, I).
%
% For each horizon h, the function computes:
%   - SNR_ai(h) and align2_ai(h) for R^2(x_i, h)
%   - SNR_cj(h) and align2_cj(h) for R^2(y_j, h)
%   - R^2(x_i, h), R^2(y_j, h), and the wedge
%
% with
%   a_1' = \iota_i' A
%   b_1' = iota_i' B
%   c_j' = iota_j' C
%   d_j' = iota_j' D.
%
% The note defines
%   R^2(x_i,h) = SNR_ai(h) * align2_ai(h)
%   R^2(y_j,h) = SNR_cj(h) * align2_cj(h)
% and therefore
%   f_th(h,i,j) = [SNR_ai(h) / SNR_cj(h)] * [align2_ai(h) / align2_cj(h)].
%
% Under Kalman orthogonality, Sigma_er(h) = Sigma_r(h). This implementation
% uses that identity when building the second moments from the steady-state
% Kalman filter.
%
% Example:
%   out = theoretical_R2_wedge(A, B, C, D, [1 2 4 8]);
%   out.state.snr
%   out.state.align2
%   out.observables.snr
%   out.observables.align2
%   out.wedge
%
% CHANGED (2026-09-15): 'state_index' may be a VECTOR. With n_i > 1 indices
% the state-side objects (out.state.*, out.rows.a1/b1) are n_i x n_h and
% out.wedge / out.implied_state_R2_from_yj are n_obs x n_h x n_i, third
% dimension in the order of state_index; with a scalar index every output
% keeps its old shape (1 x n_h rows, n_obs x n_h wedge). The observable-side
% objects and the forward covariance path do not depend on the state and are
% computed once, which is why compute_state_reliability_diagnostics now makes
% one call with state_index = 1:n_states instead of n_states calls. The
% quadratic forms a_i' Sigma a_i are evaluated for all i at once as the
% diagonal of A*Sigma*A', numerically the same to ~1e-15.

parser = inputParser();
parser.addParameter('state_index', 1, @(x) validateattributes(x, {'numeric'}, ...
    {'vector', 'integer', 'positive'}));
parser.addParameter('obs_indices', [], @(x) validateattributes(x, {'numeric'}, ...
    {'vector', 'integer', 'positive'}));
parser.addParameter('tol', 1e-12, @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'positive'}));
parser.addParameter('max_iter', 100000, @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}));
parser.addParameter('P0', [], @(x) isnumeric(x));
% filter_cache: pass out.steady_state from a previous call with the SAME
% A,B,C,D to skip the steady-state Riccati iteration (the filter does not
% depend on state_index, so it can be reused when looping over states).
parser.addParameter('filter_cache', [], @(x) isempty(x) || isstruct(x));
parser.parse(varargin{:});
opts = parser.Results;

[nx, nx2] = size(A);
if nx ~= nx2
    error('A must be square.');
end

[nxB, ne] = size(B);
if nxB ~= nx
    error('B must have the same number of rows as A.');
end

[ny, nxC] = size(C);
if nxC ~= nx
    error('C must have the same number of columns as A.');
end

if ~isequal(size(D), [ny, ne])
    error('D must be ny-by-ne, where ne = size(B, 2).');
end

horizons = horizons(:).';
validateattributes(horizons, {'numeric'}, {'vector', 'integer', 'nonnegative'});

stateIdx = opts.state_index(:).';
if any(stateIdx > nx)
    error('state_index exceeds the number of states.');
end
numI = numel(stateIdx);

if isempty(opts.obs_indices)
    obsIdx = 1:ny;
else
    obsIdx = opts.obs_indices(:).';
end

if any(obsIdx > ny)
    error('obs_indices contains an index larger than the number of observables.');
end

if isempty(opts.P0)
    P0 = 1e6 * eye(nx);
else
    P0 = opts.P0;
    if ~isequal(size(P0), [nx, nx])
        error('P0 must be nx-by-nx.');
    end
end

if ~isempty(opts.filter_cache)
    % Reuse a previously computed steady-state filter (see header note)
    Pf    = opts.filter_cache.filtered_cov;
    Ppred = opts.filter_cache.predicted_cov;
    K     = opts.filter_cache.kalman_gain;
    V     = opts.filter_cache.innovation_cov;
    info  = opts.filter_cache.info;
else
    [Pf, Ppred, K, V, info] = steady_state_lagged_measurement_filter(A, B, C, D, P0, opts.tol, opts.max_iter);
end

aRows = A(stateIdx, :);          % numI x nx, rows a_i'
bRows = B(stateIdx, :);          % numI x ne, rows b_i'
cRows = C(obsIdx, :);
dRows = D(obsIdx, :);

stateNoise = row_squared_norms(bRows);   % numI x 1
obsNoise = row_squared_norms(dRows);     % numJ x 1

numH = numel(horizons);
numJ = numel(obsIdx);
maxH = max(horizons);

SigmaEStore = cell(1, maxH + 1);
SigmaEStore{1} = symmetrize(Pf); % h = 0
P = Pf;
for step = 1:maxH
    P = A * P * A' + B * B';
    SigmaEStore{step + 1} = symmetrize(P);
end

SigmaE = zeros(nx, nx, numH);
SigmaR = zeros(nx, nx, numH);
SigmaER = zeros(nx, nx, numH);

% VECTORISED VERSION (2026-09-15): all states in state_index and all
% observables at once per horizon. quadform(v, M) = max(real(v'Mv), 0) for
% every row v of a matrix R is max(real(diag(R*M*R')), 0); safe_divide is
% applied elementwise with the same threshold.
stateSignal = zeros(numI, numH);
stateRevision = zeros(numI, numH);
stateCovER = zeros(numI, numH);
obsSignal = zeros(numJ, numH);
obsRevision = zeros(numJ, numH);
obsCovER = zeros(numJ, numH);

for hIdx = 1:numH
    h = horizons(hIdx);

    SigmaE(:, :, hIdx) = SigmaEStore{h + 1};
    SigmaR(:, :, hIdx) = symmetrize(SigmaE(:, :, hIdx) - Pf);
    SigmaER(:, :, hIdx) = SigmaR(:, :, hIdx);

    stateSignal(:, hIdx)   = quadform_rows(aRows, SigmaE(:, :, hIdx));
    stateRevision(:, hIdx) = quadform_rows(aRows, SigmaR(:, :, hIdx));
    stateCovER(:, hIdx)    = quadform_rows(aRows, SigmaER(:, :, hIdx));

    obsSignal(:, hIdx)   = quadform_rows(cRows, SigmaE(:, :, hIdx));
    obsRevision(:, hIdx) = quadform_rows(cRows, SigmaR(:, :, hIdx));
    obsCovER(:, hIdx)    = quadform_rows(cRows, SigmaER(:, :, hIdx));
end

stateForecastErrorVar = stateSignal + stateNoise;                       % numI x numH
stateSNR    = safe_divide(stateSignal, stateForecastErrorVar);
stateAlign2 = safe_divide(stateCovER.^2, stateSignal .* stateRevision);
stateR2     = stateSNR .* stateAlign2;

obsForecastErrorVar = obsSignal + obsNoise;                             % numJ x numH
obsSNR    = safe_divide(obsSignal, obsForecastErrorVar);
obsAlign2 = safe_divide(obsCovER.^2, obsSignal .* obsRevision);
obsR2     = obsSNR .* obsAlign2;

% wedge f_th(j, h, i) = [SNR_ai(h) / SNR_cj(h)] * [align2_ai(h) / align2_cj(h)]
fth = zeros(numJ, numH, numI);
impliedStateR2 = zeros(numJ, numH, numI);
for iIdx = 1:numI
    fth(:, :, iIdx) = safe_divide(repmat(stateSNR(iIdx, :), numJ, 1), obsSNR) ...
        .* safe_divide(repmat(stateAlign2(iIdx, :), numJ, 1), obsAlign2);
    impliedStateR2(:, :, iIdx) = fth(:, :, iIdx) .* obsR2;
end
if numI == 1                     % scalar state_index: the old 2-D shapes
    fth = fth(:, :, 1);
    impliedStateR2 = impliedStateR2(:, :, 1);
end

% ORIGINAL CODE (scalar state_index), kept for comparison with the vectorised
% block above:
% a1 = A(stateIdx, :).';  b1 = B(stateIdx, :).';  stateNoise = squared_norm(b1);
% for hIdx = 1:numH
%     h = horizons(hIdx);
%     SigmaE(:, :, hIdx) = SigmaEStore{h + 1};
%     SigmaR(:, :, hIdx) = symmetrize(SigmaE(:, :, hIdx) - Pf);
%     SigmaER(:, :, hIdx) = SigmaR(:, :, hIdx);
%     stateSignal(hIdx) = quadform(a1, SigmaE(:, :, hIdx));
%     stateRevision(hIdx) = quadform(a1, SigmaR(:, :, hIdx));
%     stateCovER(hIdx) = quadform(a1, SigmaER(:, :, hIdx));
%     stateForecastErrorVar(hIdx) = stateSignal(hIdx) + stateNoise;
%     stateSNR(hIdx) = safe_divide(stateSignal(hIdx), stateForecastErrorVar(hIdx));
%     stateAlign2(hIdx) = safe_divide(stateCovER(hIdx)^2, stateSignal(hIdx) * stateRevision(hIdx));
%     stateR2(hIdx) = stateSNR(hIdx) * stateAlign2(hIdx);
%     for jIdx = 1:numJ
%         cj = cRows(jIdx, :).';
%         djNoise = obsNoise(jIdx);
%         obsSignal(jIdx, hIdx) = quadform(cj, SigmaE(:, :, hIdx));
%         obsRevision(jIdx, hIdx) = quadform(cj, SigmaR(:, :, hIdx));
%         obsCovER(jIdx, hIdx) = quadform(cj, SigmaER(:, :, hIdx));
%         obsForecastErrorVar(jIdx, hIdx) = obsSignal(jIdx, hIdx) + djNoise;
%         obsSNR(jIdx, hIdx) = safe_divide(obsSignal(jIdx, hIdx), obsForecastErrorVar(jIdx, hIdx));
%         obsAlign2(jIdx, hIdx) = safe_divide(obsCovER(jIdx, hIdx)^2, obsSignal(jIdx, hIdx) * obsRevision(jIdx, hIdx));
%         obsR2(jIdx, hIdx) = obsSNR(jIdx, hIdx) * obsAlign2(jIdx, hIdx);
%         fth(jIdx, hIdx) = safe_divide(stateSNR(hIdx), obsSNR(jIdx, hIdx)) ...
%             * safe_divide(stateAlign2(hIdx), obsAlign2(jIdx, hIdx));
%         impliedStateR2(jIdx, hIdx) = fth(jIdx, hIdx) * obsR2(jIdx, hIdx);
%     end
% end

out = struct();
out.horizons = horizons;
out.state_index = stateIdx;
out.obs_indices = obsIdx;

out.steady_state = struct();
out.steady_state.filtered_cov = Pf;
out.steady_state.predicted_cov = Ppred;
out.steady_state.kalman_gain = K;
out.steady_state.innovation_cov = V;
out.steady_state.info = info;

out.rows = struct();
out.rows.a1 = aRows;             % numI x nx (1 x nx for a scalar state_index, as before)
out.rows.b1 = bRows;
out.rows.cj = cRows;
out.rows.dj = dRows;

out.second_moments = struct();
out.second_moments.Sigma_e = SigmaE;
out.second_moments.Sigma_r = SigmaR;
out.second_moments.Sigma_er = SigmaER;

out.state = struct();
out.state.signal = stateSignal;
out.state.noise = stateNoise;
out.state.revision_var = stateRevision;
out.state.cov_er = stateCovER;
out.state.forecast_error_var = stateForecastErrorVar;
out.state.snr = stateSNR;
out.state.align2 = stateAlign2;
out.state.R2 = stateR2;

out.observables = struct();
out.observables.signal = obsSignal;
out.observables.noise = obsNoise;
out.observables.revision_var = obsRevision;
out.observables.cov_er = obsCovER;
out.observables.forecast_error_var = obsForecastErrorVar;
out.observables.snr = obsSNR;
out.observables.align2 = obsAlign2;
out.observables.R2 = obsR2;

out.wedge = fth;
out.implied_state_R2_from_yj = impliedStateR2;

if nx == 1
    out.scalar_case = struct();
    out.scalar_case.wedge_snr_only = bsxfun(@rdivide, stateSNR, obsSNR);
end

end

function x = quadform_rows(R, M)
% quadform for every row of R at once: x(i) = max(real(R(i,:) * M * R(i,:)'), 0)
x = real(sum((R * M) .* R, 2));
x = max(x, 0);
end

function [Pf, Ppred, K, V, info] = steady_state_lagged_measurement_filter(A, B, C, D, P0, tol, maxIter)
% Steady-state Kalman covariance for y_t = C x_{t-1} + D e_t.
% CHANGED (2026-09-16): the same fixed point as the engine's, obtained from
% solve_steady_state_riccati (direct DARE solve, iteration fallback started at
% P0) instead of the third independent fixed-point iteration this file used to
% run - which hit its cap alongside the other two on pile-up draws.
[Pf, K, V, Ppred, rf] = solve_steady_state_riccati(A, B, C, D, ...
    struct('Omega0', symmetrize(P0), 'tol_riccati', tol, 'MaxIter', maxIter));
Pf = symmetrize(Pf);  V = symmetrize(V);  Ppred = symmetrize(Ppred);

info = struct();
info.converged  = rf.converged;
info.iterations = rf.iterations;
info.tol        = tol;
info.max_iter   = maxIter;
info.method     = rf.method;
info.residual   = rf.rel_step;

% ORIGINAL CODE (own fixed-point iteration), kept for comparison:
% P = symmetrize(P0);
% Q = B * B';
% R = D * D';
% N = B * D';
% converged = false;
% for iter = 1:maxIter
%     V = symmetrize(C * P * C' + R);
%     S = A * P * C' + N;
%     K = right_solve(S, V);
%     Ppred = A * P * A' + Q;
%     Pnew = symmetrize(Ppred - K * V * K');
%     if norm(Pnew - P, 'fro') <= tol * (1 + norm(P, 'fro'))
%         P = Pnew;
%         converged = true;
%         break;
%     end
%     P = Pnew;
% end
% Pf = P;
% V = symmetrize(C * Pf * C' + R);
% S = A * Pf * C' + N;
% K = right_solve(S, V);
% Ppred = symmetrize(A * Pf * A' + Q);
% info = struct('converged', converged, 'iterations', iter, 'tol', tol, 'max_iter', maxIter);
% if ~converged
%     warning('theoretical_R2_wedge:noConvergence', ...
%         'Steady-state covariance recursion did not converge within max_iter.');
% end

end

function x = quadform(v, M)
x = real(v' * M * v);
x = max(x, 0);
end

function n2 = squared_norm(v)
n2 = real(v' * v);
n2 = max(n2, 0);
end

function n2 = row_squared_norms(X)
n2 = sum(X.^2, 2);
n2 = max(n2, 0);
end

function y = safe_divide(num, den)
% Elementwise (num and den of equal size, or one of them scalar); NaN where
% the denominator is negligible relative to the numerator - the same
% threshold as the original scalar version:
%   if abs(den) <= eps(max(1, abs(num))), y = NaN; else, y = num / den; end
y = num ./ den;
small = abs(den) <= eps(max(1, abs(num)));
y(small) = NaN;
end

function M = symmetrize(M)
M = 0.5 * (M + M');
end

function X = right_solve(A, B)
if rcond(B) < 1e-12
    X = A * pinv(B);
else
    X = A / B;
end
end
