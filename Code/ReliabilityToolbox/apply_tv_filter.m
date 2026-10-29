function [loglik, all_states, surprises, filter_ok] = apply_tv_filter(A, C, Sequence, data, X0, tol)
%APPLY_TV_FILTER  Apply pre-computed diffuse/time-varying Kalman gains to data.
%
%   Extracted verbatim (2026-08) from section 6b of
%   compute_state_reliability_diagnostics so that the same recursion can be
%   reused by the posterior-predictive reweighting code, which evaluates
%   log p(y|theta) for several data sets against ONE Sequence object.
%
%   The gains in Sequence depend only on (A,B,C,D,X0,B0,Qscale,T) and NOT on
%   the data, so a caller may build Sequence once and call this repeatedly.
%
%   INPUTS
%     A, C       pruned state transition / observation loading matrices
%     Sequence   output of kalman_filter_steady_and_diffuse
%     data       T x n_obs, ALREADY DEMEANED by the observables' steady state
%     X0         initial state mean; [] means zeros(n_states,1)
%     tol        (optional, default 1e-6) threshold on F below which a
%                univariate diffuse-phase observation contributes no
%                quadratic likelihood term
%     chol_cache (optional) 1 x T cell of Cholesky factors of Sequence.V(:,:,t);
%                pass the third output of build_chol_cache to skip refactorizing
%
%   OUTPUTS
%     loglik     scalar log-likelihood (NaN if a Cholesky failed)
%     all_states T x n_states filtered states xhat(t|t); rows from the
%                failure point onward are NaN if filter_ok is false
%     surprises  T x n_obs forecast errors; rows in the diffuse phase are
%                NaN by construction (see note below)
%     filter_ok  logical, false if any Cholesky failed

if nargin < 6 || isempty(tol),        tol = 1e-6;        end

[T, n_obs] = size(data);
n_states   = size(A, 1);

all_states = nan(T, n_states);
surprises  = nan(T, n_obs);
filter_ok  = true;

if isempty(X0)
    current_filtered_X = zeros(n_states, 1);
else
    current_filtered_X = X0(:);
end

% Stacked initial state: [0 for the y-only slots; x0]
alpha_filtered = [zeros(Sequence.n_obs_only, 1); current_filtered_X];
initial_state  = current_filtered_X';   % 1 x n_states, used after hand-off

% Start from the diffuse-phase log-likelihood contribution
% (-0.5 log F terms; see durbin_koopman_diffuse_filter)
loglik = Sequence.loglik_contribut;
tt     = 0;

while filter_ok && tt < T
    tt = tt + 1;

    if tt < Sequence.t_proper
        % ---------- DIFFUSE PHASE: univariate sequential ----------
        alpha_hat = Sequence.T_stacked * alpha_filtered;

        for i = 1:n_obs
            v_i       = data(tt, i) - Sequence.Z(i, :) * alpha_hat;
            alpha_hat = alpha_hat + Sequence.K_stacked(:, i, tt) * v_i;

            % Quadratic term only for NON-diffuse observations (DK 2012, eq. 7.14)
            if ~Sequence.obs_is_diffuse(i, tt) ...
                    && isfinite(Sequence.F_used(i, tt)) ...
                    && Sequence.F_used(i, tt) > tol
                loglik = loglik - 0.5 * (v_i^2 / Sequence.F_used(i, tt));
            end
        end

        alpha_filtered   = alpha_hat;
        all_states(tt,:) = alpha_hat(Sequence.n_obs_only + 1 : end);
        initial_state    = all_states(tt, :);

    else
        % ---------- STANDARD PHASE: multivariate batch ----------
        [L, chol_success] = robust_cholesky(Sequence.V(:,:,tt));
       

        if ~chol_success
            loglik    = NaN;
            filter_ok = false;
            break;
        end

        surprises(tt, :) = data(tt, :) - (C * initial_state')';
        new_state        = A * initial_state' + Sequence.K_gain(:,:,tt) * surprises(tt, :)';
        all_states(tt,:) = new_state;

        try
            log_det_V    = 2 * sum(log(diag(L)));
            inv_surprise = L \ surprises(tt, :)';
            loglik       = loglik - 0.5*(inv_surprise'*inv_surprise) - 0.5*log_det_V;
        catch
            loglik    = NaN;
            filter_ok = false;
            break;
        end

        initial_state = new_state';
    end
end

end