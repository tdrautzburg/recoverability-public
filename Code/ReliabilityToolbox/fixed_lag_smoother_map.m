function [lag_m, base_idx, info] = fixed_lag_smoother_map(A, B, key_state_indices, state_var_names, max_lag, tol_rel, verbose)
%FIXED_LAG_SMOOTHER_MAP  Which states of interest are fixed-lag SMOOTHED estimates?
%
%   ADDED (2026-08). Support for the eq.-(3.17) reliability measure of the
%   paper (Section 3.5, "Filtering versus smoothing: the value of hindsight").
%
%   WHY THIS EXISTS.  Section 2.2 implements the m-period fixed-lag smoother by
%   augmenting the model with a LAGGED COPY of the state,
%
%       l20rrstar(t) = rrstar(t-20)      (HLW_2023_US_TD_small.mod)
%
%   so that the Kalman FILTER estimate of the augmented state at t,
%   E[l20rrstar(t) | y^t] = E[rrstar(t-20) | y^t], is exactly the m = 20
%   fixed-lag SMOOTHED estimate of r*, and the filtered variance of that
%   coordinate, [KF_Omega]_jj, is exactly the smoothing error variance
%   [Omega_m]_ii of eq. (3.16)-(3.17).
%
%   The catch is the DENOMINATOR.  Applying eq. (3.11) mechanically to the
%   augmented coordinate normalises by Var[x_{j,t+h} | y^t] = Var[r*_{t+h-m} | y^t],
%   the conditional variance of the EARLIER state - which is precisely the
%   violation the paper warns about ("augmenting the state vector with lagged
%   copies of the states produces a smoothed estimate of an earlier state, so
%   the conditional variance in the denominator is that of the earlier state,
%   and the resulting comparison is not the one covered by (3.16)").  For an
%   I(1) state that denominator is smaller than the filter's by roughly m
%   periods' worth of innovation variance, so the smoothed R2 can - and in
%   Table 3 of the draft does - come out BELOW the filtered R2, even though
%   the smoother's error variance is uniformly smaller.
%
%   Eq. (3.17) fixes the normalisation by keeping the FILTER's denominator,
%
%       R2^(m)_{i,h} = 1 - [Omega_m]_ii / Var[x_{i,t+h} | y^t],
%
%   i.e. the conditional variance (3.10) of the CONTEMPORANEOUS base state i
%   (rrstar), at the same horizon h, for both m = 0 and m > 0.  With a common
%   denominator, Omega_l <= Omega_m for l >= m carries over one-for-one and
%   reliability is weakly increasing in the amount of hindsight.
%
%   This function supplies the missing link: for every state of interest it
%   reports the lag m and the index of the base state i whose conditional
%   variance belongs in the denominator.
%
%   HOW THE MAPPING IS FOUND.  Structurally, not by name.  In the pruned
%   state space x(t) = A x(t-1) + B e(t) the impulse responses of state j are
%   Psi_j(s) = e_j' A^(s-1) B, and state j is the m-period lag of state k iff
%
%       Psi_j(s) = Psi_k(s-m)  for all s >= 1,   with Psi_k(s) = 0 for s <= 0,
%
%   i.e. iff row j reproduces row k shifted right by m.  The number of leading
%   all-zero response blocks of row j pins m down uniquely, and the base is
%   then the (unique) row that matches at that m; it is by construction a root
%   of the lag chain, since a root responds on impact.  Name conventions
%   ('l20rrstar') are therefore never parsed - the test works whether the base
%   is a genuine Dynare state or one of the augmented non-state variables of
%   interest (rrstar = c*gg + zz is the latter), and it is immune to renaming.
%
%   INPUTS
%     A, B               pruned state-space matrices (n_states x n_states,
%                        n_states x n_shocks), as built in step 3 of
%                        compute_state_reliability_diagnostics
%     key_state_indices  n_key x 1, positions of the states of interest in the
%                        pruned state vector
%     state_var_names    n_states x 1 cell of names (only used for reporting)
%     max_lag            (optional) largest lag considered. Default
%                        n_states - 1: no longer chain is representable, so
%                        the default never truncates a real chain. Setting it
%                        BELOW the true lag does not disable the match, it
%                        makes it land on an intermediate link of the chain
%                        instead of the root - a denominator that is still
%                        dated wrong, just less so. Leave it empty.
%     tol_rel            (optional) relative tolerance of the response match.
%                        Default 1e-5. NOTE this cannot be tightened much:
%                        step 2 of the engine rounds A to 1e-6 absolute, so
%                        the base row and the lagged row pick up slightly
%                        different rounding and agree only to ~1e-6 in
%                        relative terms. A structurally DIFFERENT state misses
%                        by O(1) over the whole response window, so there is
%                        no practical risk of a false positive at 1e-5.
%     verbose            (optional) print the detected map. Default false.
%
%   OUTPUTS
%     lag_m     n_key x 1, the m of eq. (3.17). 0 for a state that is not a
%               lagged copy of anything, i.e. for a filtered estimate.
%     base_idx  n_key x 1, index IN THE PRUNED STATE VECTOR of the base state
%               whose conditional variance (3.10) is the denominator. Equal to
%               key_state_indices when lag_m = 0, so callers can use it
%               unconditionally.
%     info      struct with .base_name (n_key x 1 cell), .mismatch (n_key x 1,
%               the achieved relative mismatch of the accepted match, NaN when
%               lag_m = 0), .candidates (n_key x 1 cell of [k m err] rows, the
%               matches found AT THE ACCEPTED LAG - more than one row there
%               means two states of the model are observationally identical,
%               and the closest match was taken), .max_lag, .tol_rel.
%
%   A state of interest that LOOKS lagged (leading zero responses) but whose
%   base was pruned away in step 3 gets lag_m = 0 and a warning under verbose:
%   the caller then falls back to the eq.-(3.11) normalisation for that row,
%   which is the pre-2026-08 behaviour. Keep the base in states_of_interest -
%   as HLW does with rrstar - and this cannot happen.
%
%   See also: compute_state_reliability_diagnostics (step 5c),
%             reliability_draw_wrapper.

if nargin < 5 || isempty(max_lag), max_lag = size(A, 1) - 1; end
if nargin < 6 || isempty(tol_rel), tol_rel = 1e-5;           end
if nargin < 7 || isempty(verbose), verbose = false;          end

n_states = size(A, 1);
n_shocks = size(B, 2);
n_key    = numel(key_state_indices);

max_lag  = max(0, min(max_lag, n_states - 1));

lag_m    = zeros(n_key, 1);
base_idx = key_state_indices(:);

info = struct();
info.base_name  = state_var_names(base_idx);
info.mismatch   = nan(n_key, 1);
info.candidates = cell(n_key, 1);
info.max_lag    = max_lag;
info.tol_rel    = tol_rel;

if max_lag == 0 || n_shocks == 0
    return
end

%% ------------------------------------------------------------------
%  Impulse responses Psi(:,:,s) = A^(s-1) * B, s = 1 .. S
% -------------------------------------------------------------------
% Stored HORIZON-MAJOR and flattened, one row per state:
%   Psi2(j, (s-1)*n_shocks + (1:n_shocks)) = e_j' A^(s-1) B.
% With that layout "row j equals row k shifted right by m" is a plain
% comparison of two column ranges, so the search below needs no reshaping and
% runs one vectorised comparison per lag - it is executed once per posterior
% draw, so the constant matters.
%
% S must leave enough horizons AFTER the shift for the comparison to be
% informative at the largest lag considered: S - max_lag = n_states + 10.
S    = max_lag + n_states + 10;
Psi2 = zeros(n_states, S*n_shocks);
Mpow = B;
for s = 1:S
    Psi2(:, (s-1)*n_shocks + (1:n_shocks)) = Mpow;
    Mpow = A*Mpow;
end

% Scale of each row, used to normalise the comparison and to skip rows that
% carry no shock at all (those would match anything).
row_scale = max(abs(Psi2), [], 2);
zero_row  = ~(row_scale > 0);
self      = (1:n_states)';

%% ------------------------------------------------------------------
%  Match every state of interest against every state
% -------------------------------------------------------------------
% The lag is pinned down BEFORE any comparison: if state j is an m-period lag
% of anything, it cannot respond to a shock before period m+1, so m equals the
% number of leading all-zero response blocks of row j. The search therefore
% starts at that m - where the match is against the ROOT of the lag chain, the
% contemporaneous base state we want - and only walks down if the root is not
% in the pruned state vector, in which case the deepest surviving link of the
% chain is the best available denominator. In the intended setup (the base is
% itself one of M_.states_of_interest, as rrstar and drrstar are in HLW) the
% first test succeeds and the loop runs once.
for nn = 1:n_key
    jj = key_state_indices(nn);
    if zero_row(jj)
        continue
    end

    % First horizon at which row jj responds; m_top = leading zero blocks.
    first_nz = find(max(abs(reshape(Psi2(jj, :), n_shocks, S)), [], 1) > tol_rel*row_scale(jj), 1);
    m_top    = min(first_nz - 1, max_lag);

    for m = m_top:-1:1
        % Compare row jj from horizon m+1 onwards against EVERY row from
        % horizon 1 onwards, in one shot.
        tail    = Psi2(jj, m*n_shocks+1:end);
        err     = max(abs(bsxfun(@minus, Psi2(:, 1:(S-m)*n_shocks), tail)), [], 2) ./ ...
                  max(row_scale, row_scale(jj));
        matched = find(err <= tol_rel & ~zero_row & self ~= jj);

        if ~isempty(matched)
            [best_err, which]   = min(err(matched));
            lag_m(nn)           = m;
            base_idx(nn)        = matched(which);
            info.base_name{nn}  = state_var_names{matched(which)};
            info.mismatch(nn)   = best_err;
            info.candidates{nn} = [matched(:), repmat(m, numel(matched), 1), err(matched)];
            break
        end
    end

    if lag_m(nn) == 0 && m_top > 0 && verbose
        % Leading zeros but no base found: the base was pruned in step 3.
        fprintf(2, ['fixed_lag_smoother_map: %s does not respond before period %d but no ' ...
                    'base state was found - it is treated as a filtered estimate ' ...
                    '(eq. (3.11)). Add its contemporaneous counterpart to ' ...
                    'M_.states_of_interest to get eq. (3.17).\n'], state_var_names{jj}, m_top+1);
    end
end

if verbose
    disp('Fixed-lag smoother map (eq. (3.17) denominators):');
    for nn = 1:n_key
        if lag_m(nn) > 0
            fprintf('  %-14s = %s(-%d)  ->  denominator: Var[%s_{t+h}|y^t]  (rel. mismatch %.1e)\n', ...
                state_var_names{key_state_indices(nn)}, info.base_name{nn}, lag_m(nn), ...
                info.base_name{nn}, info.mismatch(nn));
        else
            fprintf('  %-14s  filtered estimate (m = 0), eq. (3.17) = eq. (3.11)\n', ...
                state_var_names{key_state_indices(nn)});
        end
    end
end

end
