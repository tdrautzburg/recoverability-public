function [K_gain, filtered_X, log_llhood, Omega, V, Omega_tp1_t, surprises, Sequence, ric_flag] = kalman_filter_steady_and_diffuse(A, B, C, D, data, X0, B0, Qscale, T, dr, init_type)
%KALMAN_FILTER_STEADY_AND_DIFFUSE  Steady-state + time-varying (diffuse) Kalman filter.
%
%   RENAMED (2026-07): formerly my_KalmanFilter_Heteroskedastic_jan27.m in Miscellaneous/. All call
%   sites inside ReliabilityToolbox/ use the new names; the full mapping
%   is in README.md ("Suggested file names").
%
%   State-space convention (ABCD form, observables load on lagged states):
%
%       x(t) = A x(t-1) + B Q(t)^{1/2} e(t)
%       y(t) = C x(t-1) + D Q(t)^{1/2} e(t),    e(t) ~ N(0, I)
%
%   where Q(t) = diag(Qscale(:,t)) allows for heteroskedastic shock
%   variances (Qscale scales VARIANCES; Qscale = 1 everywhere recovers the
%   homoskedastic model).
%
%   The function returns two filters:
%     (a) the STEADY-STATE filter: gain and covariances from the Riccati
%         fixed point, solved directly by solve_steady_state_riccati since
%         2026-09-16 (outputs K_gain, Omega, V, Omega_tp1_t, and - if data is
%         supplied - filtered_X, log_llhood, surprises);
%     (b) the TIME-VARYING filter object `Sequence`, produced by
%         durbin_koopman_diffuse_filter: an exact Durbin-Koopman diffuse
%         initialization for the unit roots, per-period gains/variances,
%         and the diffuse-phase log-likelihood contribution.
%
%   INPUTS
%     A,B,C,D  ABCD matrices (see above)
%     data     n_obs x T data matrix (time in COLUMNS); may be empty
%     X0       initial state mean (used by the data pass below)
%     B0       factor of the initial state covariance: P0 = B0*B0'
%              (with init_type 'diffuse' used only to initialize the
%              steady-state Riccati iteration - the time-varying filter uses
%              the exact diffuse initialization; with 'proper' it IS the
%              time-varying filter's initial covariance, see init_type)
%     Qscale   n_shocks x T matrix of shock VARIANCE scales (NaN treated as 1)
%     T        sample length (defaults to size(Qscale,2))
%     dr       Dynare decision-rule structure; must carry the fields set by
%              the caller: restrict_var_list, restrict_columns, mf, nstatic,
%              npred, keep_states, qz_crit (see compute_state_reliability_diagnostics)
%     init_type  ADDED 2026-09-26, optional, default 'diffuse':
%              'diffuse'  exact Durbin-Koopman diffuse initialization of the
%                         unit roots, stationary block at its unconditional
%                         covariance (Dynare's compute_Pinf_Pstar) - as before;
%              'proper'   the time-varying filter starts from the PROPER prior
%                         x_0 ~ N(X0, B0*B0') with NO diffuse phase (t_proper
%                         = 1), the convention of a sampler that initialises
%                         its Kalman filter with a wide finite variance (the
%                         Bauer-Rudebusch ese.R: a1 = 2% for i*, P1 = I). The
%                         steady-state objects are unaffected either way.
%
%   OUTPUTS
%     ric_flag   : struct describing the Riccati solve (ADDED 2026-09), with
%                  fields .converged, .iterations, .rel_step, .rho_cl, and
%                  since 2026-09-16 .method ('dynare'|'qz'|'iterate') and
%                  .residuals; see solve_steady_state_riccati. .iterations is
%                  0 for a direct solve; .rel_step is the error estimate
%                  r/(1-rho^2) of the accepted solution. .converged is false
%                  if no method passed the residual check; Omega must not be
%                  treated as a fixed point in that case.
%     K_gain, Omega, V, Omega_tp1_t : steady-state filter objects
%     filtered_X, log_llhood, surprises : data pass results (see caveat below)
%     Sequence : time-varying filter object (see durbin_koopman_diffuse_filter)
%
%   CAVEAT on the data pass in this function: during the diffuse phase
%   (t < Sequence.t_proper) Sequence.K_gain(:,:,t) is zero, so the state is
%   updated by pure prediction; the exact diffuse mean update on the stacked
%   state is implemented by the caller (compute_state_reliability_diagnostics), which
%   uses Sequence.K_stacked. The filtered_X/log_llhood computed here are
%   therefore only used for standalone/diagnostic calls.
%
%   Changes vs. the original version:
%     - fixed `isempty('X0')` (always false: it tested the literal string
%       'X0') -> `isempty(X0)`;
%     - removed a dead Sequence pre-allocation that was overwritten by the
%       hybrid filter output, and a dead call to
%       diffuse_filter_initialization whose outputs were ignored by
%       durbin_koopman_diffuse_filter (it computes its own initialization
%       via Dynare's compute_Pinf_Pstar);
%     - with the log-likelihood fix in durbin_koopman_diffuse_filter
%       (Sequence.loglik_contribut now contains ONLY diffuse-phase terms),
%       the -0.5*log|V| added below for t >= t_proper is no longer double
%       counted.

nstates = size(A, 1);
nshocks = size(B, 2);
nobs    = size(C, 1);

% ------------------------------------------------------------------------
% 0. Input handling / defaults
% ------------------------------------------------------------------------
if ~exist('Qscale', 'var') || isempty(Qscale)
    if ~exist('T', 'var') || isempty(T), T = 250; end
    Qscale = ones(nshocks, T);
end
Qscale(isnan(Qscale)) = 1;          % missing scales = homoskedastic
if ~exist('T', 'var') || isempty(T)
    T = size(Qscale, 2);
end

% ------------------------------------------------------------------------
% 1. Steady-state filter: Riccati fixed point.
% ------------------------------------------------------------------------
% CHANGED (2026-09-16): solved directly by solve_steady_state_riccati (Dynare's
% DARE MEX or the toolbox's own QZ solver; the fixed-point iteration below is
% kept there as the fallback and reference). The iteration converges at rate
% rho(A-KC)^2 (~19/(1-rho^2) steps): thousands of steps for typical HLW draws,
% past the 1e5 cap for pile-up draws (sigma_g, sigma_z ~ 0, rho ~ 0.99997);
% with three such solves per draw it was half to three quarters of the
% per-draw cost. The direct solve takes milliseconds at any rho and is
% accepted only after a residual check (see that function). ric_flag keeps its
% four fields; .rel_step is now an ERROR estimate r/(1-rho^2) rather than the
% last iteration step, .iterations is 0 when a direct solver was used.
% B0*B0' is passed as the start of the fallback iteration, as before.
if ~exist('B0', 'var') || isempty(B0)
    Omega0 = [];
else
    Omega0 = B0*B0';
end
[Omega, K_gain, V, Omega_tp1_t, ric_flag] = solve_steady_state_riccati(A, B, C, D, struct('Omega0', Omega0));

% ORIGINAL CODE (fixed-point iteration, 2026-09 version), kept for comparison:
% if ~exist('B0', 'var') || isempty(B0)
%     Omega = zeros(nstates, nstates);
% else
%     Omega = B0*B0';
% end
% MaxIter     = 1e5;
% tol_riccati = 1e-10;                        % RELATIVE step, not absolute
% change_norm = Inf;
% nonfinite   = false;
% ii = 0;
% while change_norm > tol_riccati && ii < MaxIter
%     ii = ii + 1;
%     V        = C*Omega*C' + D*D';
%     K_gain   = (A*Omega*C' + B*D') / V;
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
%             '(relative step %.3e, rho(A-KC) = %.6f).%s'], 'kalman_filter_steady_and_diffuse', ...
%             ii, change_norm, ric_flag.rho_cl, repmat(' Iterate became non-finite.', 1, double(nonfinite)));
% end

% ------------------------------------------------------------------------
% 2. Time-varying filter. Default: exact diffuse initialization (Durbin-Koopman
%    univariate treatment of the unit roots, Dynare-style Pinf/Pstar), the
%    hybrid filter computing its own initialization internally from `dr` (the
%    6th/7th arguments empty). ADDED 2026-09-26: init_type 'proper' passes
%    P0 = B0*B0' as the finite initial covariance and a zero diffuse part, so
%    the recursion starts in the standard phase at t = 1 (see header).
% ------------------------------------------------------------------------
if ~exist('init_type', 'var') || isempty(init_type), init_type = 'diffuse'; end
switch lower(init_type)
    case 'diffuse'
        Sequence = durbin_koopman_diffuse_filter(A, B, C, D, Qscale, [], [], T, dr);
    case 'proper'
        if ~exist('B0', 'var') || isempty(B0)
            error('kalman_filter_steady_and_diffuse:init', 'init_type ''proper'' needs B0 (P0 = B0*B0'').');
        end
        Sequence = durbin_koopman_diffuse_filter(A, B, C, D, Qscale, B0*B0', zeros(nstates), T, dr);
    otherwise
        error('kalman_filter_steady_and_diffuse:init', 'init_type must be ''diffuse'' or ''proper'', got ''%s''.', init_type);
end
Sequence.init_type = lower(init_type);
Sequence.t_proper  = Sequence.d_diffuse + 1;   % first post-diffuse period (1 under 'proper')

% ------------------------------------------------------------------------
% 3. Optional data pass (steady-state mean update with time-varying gains;
%    see CAVEAT in the header).
% ------------------------------------------------------------------------
filtered_X = [];
log_llhood = [];
surprises  = [];

if exist('data', 'var') && ~isempty(data)
    TT = size(data, 2);
    filtered_X = nan(nstates, TT);
    surprises  = nan(nobs, TT);

    if ~exist('X0', 'var') || isempty(X0)      % FIX: was isempty('X0')
        current_filtered_X = zeros(nstates, 1);
    else
        current_filtered_X = X0;
    end

    % Start from the diffuse-phase contribution (log F terms and diffuse
    % quadratic terms are handled inside the hybrid filter / by the caller).
    log_llhood = Sequence.loglik_contribut;

    [L, chol_success] = robust_cholesky(Sequence.V(:,:,1));
    tt = 0;
    while chol_success && tt < TT
        tt = tt + 1;
        surprises(:, tt) = data(:, tt) - C*current_filtered_X;

        % Gaussian log-likelihood only over the post-diffuse sample;
        % L is refactored each period because V is time-varying.
        if tt >= Sequence.t_proper
            try
                log_det_V    = 2 * sum(log(diag(L)));
                inv_surprise = L \ surprises(:, tt);
                log_llhood   = log_llhood - 0.5*(inv_surprise'*inv_surprise) - 0.5*log_det_V;
            catch
                log_llhood   = NaN;
                surprises    = [];
                chol_success = false;
            end
        end

        current_filtered_X = A*current_filtered_X + Sequence.K_gain(:,:,tt)*surprises(:, tt);
        filtered_X(:, tt)  = current_filtered_X;
        [L, chol_success]  = robust_cholesky(Sequence.V(:,:,min(T, tt+1)));
    end
end

end
