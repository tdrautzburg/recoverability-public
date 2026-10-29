function [max_half_life, all_half_lives, shock_idx] = estimation_error_half_life(AKC, BKD, which_i, max_H)
%ESTIMATION_ERROR_HALF_LIFE  Half-life of Kalman filtering errors, shock by shock.
%
%   RENAMED (2026-07): formerly half_life_by_shock.m in Miscellaneous/. All call
%   sites inside ReliabilityToolbox/ use the new names; the full mapping
%   is in README.md ("Suggested file names").
%
%   The filtering error z(t) = x(t) - xhat(t) evolves as
%
%       z(t) = (A - K C) z(t-1) + (B - K D) e(t)  =  AKC z(t-1) + BKD e(t),
%
%   so the impulse response of z to shock u tells how long an estimation
%   error caused by that shock persists. For each target state (rows
%   which_i of z) and each shock, the "half-life" is measured as the last
%   horizon at which |IRF| still exceeds 50% of its peak absolute value.
%
%   INPUTS
%     AKC     n x n   error transition matrix, A - K*C
%     BKD     n x k   error shock loading, B - K*D
%     which_i         indices of the target states within z
%     max_H           horizon over which IRFs are computed
%
%   OUTPUTS
%     max_half_life   1 x n_targets  max over shocks (worst case), NaN if
%                     the error process is explosive
%     all_half_lives  n_shocks_active x n_targets  half-life per shock
%     shock_idx       indices of the shocks that actually load on the
%                     error (all-zero columns of BKD are dropped); use it
%                     to align all_half_lives with the original shock list
%
%   CONVENTIONS / CAVEATS (kept from the original implementation so that
%   published results are reproducible):
%     - the half-life is expressed in IRF-slice units where slice 1 is the
%       impact period (h = 0); the "+hh+1" indexing below adds one period
%       relative to the raw slice position of the last crossing;
%     - if the IRF never falls below 50% of its peak within max_H, the
%       half-life is right-censored at ~max_H;
%     - if the IRF is below 50% of its peak everywhere after the peak, the
%       half-life is set to 1;
%     - eigenvalues of AKC up to 1.0001 are tolerated (allows for unit
%       roots surviving in the error process due to round-off).

ntargets = numel(which_i);

% Drop shocks that do not load on the error at all (zero columns of BKD)
shock_idx = 1:size(BKD, 2);
shock_idx(:, prod(BKD == 0) == 1) = [];
BKD = BKD(:, shock_idx);
nshocks = size(BKD, 2);

if max(abs(eig(AKC))) > 1.0001
    % Explosive error process: half-lives are meaningless
    fprintf(2, 'system is not stable: Max eigenvalue of A-KC is explosive \n');
    all_half_lives = nan(nshocks, ntargets);
else
    all_half_lives = nan(nshocks, ntargets);

    % IRFs of the error process: slice 1 = impact (h = 0)
    nvars = size(AKC, 1);
    IRFs = zeros(nvars, nshocks, max_H + 1);
    IRFs(:, :, 1) = BKD;
    for hh = 2:max_H + 1
        IRFs(:, :, hh) = AKC * IRFs(:, :, hh - 1);
    end

    for nn = 1:ntargets
        % reshape (not squeeze): squeeze flips the orientation to a column
        % when nshocks == 1, which made the function silently return 1
        squeezed_IRFs = reshape(IRFs(which_i(nn), :, :), nshocks, []);   % nshocks x (max_H+1)
        [max_IRF, hh] = max(abs(squeezed_IRFs), [], 2);     % peak and its slice, per shock
        for uu = 1:nshocks
            % Last slice after the peak where |IRF| is still above half the peak
            tmp = find(abs(squeezed_IRFs(uu, hh(uu)+1:end)) / max_IRF(uu) > 0.5, 1, 'last') + hh(uu) + 1;
            if isempty(tmp)
                all_half_lives(uu, nn) = 1;
            else
                all_half_lives(uu, nn) = tmp;
            end
        end
    end
end

% Worst case across shocks (max ignores NaN only if some entries are finite;
% an all-NaN column propagates NaN, which callers test with isnan)
max_half_life = max(all_half_lives, [], 1);
