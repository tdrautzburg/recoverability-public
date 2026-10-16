function Sequence = durbin_koopman_diffuse_filter(A, B, C, D, Qscale, Omega_star_init, Omega_inf_init, T, dr, tol)
%DURBIN_KOOPMAN_DIFFUSE_FILTER  Hybrid diffuse Kalman filter (Dynare-style).
%
%   RENAMED (2026-07): formerly f_kalman_filter_diffuse_hybrid.m in Miscellaneous/. All call
%   sites inside ReliabilityToolbox/ use the new names; the full mapping
%   is in README.md ("Suggested file names").
%
%   Computes the covariance/gain sequence of the exact diffuse Kalman filter
%   for a model with unit roots, then hands over to the standard filter once
%   the diffuse phase ends. Only covariances and gains are computed here (no
%   data pass): callers apply the stored gains to the data.
%
%   Model (ABCD form, observables load on lagged states):
%       x(t) = A x(t-1) + B Q(t)^{1/2} e(t)
%       y(t) = C x(t-1) + D Q(t)^{1/2} e(t),   e(t) ~ N(0, I),
%   with Q(t) = diag(Qscale(:,t)) a VARIANCE scale.
%
%   METHOD
%   During the diffuse phase the filter works on the Dynare-style "stacked"
%   state alpha(t) = [y-only variables; x(t)] so that observables become
%   exact selections (Z picks elements of alpha), which permits the
%   univariate (sequential) Durbin-Koopman diffuse update (DK 2012, ch. 5/7).
%   The stacked transition/selection matrices and the diffuse initialization
%   Pinf/Pstar are taken from Dynare's own kalman_transition_matrix and
%   compute_Pinf_Pstar applied to dr, then pruned to the states kept by
%   compute_state_reliability_diagnostics (dr.keep_states) and permuted to match the
%   ordering of the A,B,C,D matrices. Once rank(Z Pinf Z') = 0 the filter
%   switches to the standard multivariate recursion on x alone.
%
%   INPUTS
%     A,B,C,D          ABCD matrices AFTER pruning (as built in
%                      compute_state_reliability_diagnostics)
%     Qscale           n_shocks x T variance scales
%     Omega_star_init, Omega_inf_init
%                      BOTH EMPTY (default): the initialization is computed
%                      internally with Dynare's compute_Pinf_Pstar (exact
%                      diffuse treatment of the unit roots).
%                      REVIVED 2026-09-26 for the PROPER initialization: with
%                      Omega_star_init an m x m covariance (m = size(A,1),
%                      pruned ABCD order) and Omega_inf_init empty or all
%                      zeros, the state block of the stacked initial
%                      covariance is set to Omega_star_init and the diffuse
%                      part to zero, so the recursion starts in the standard
%                      phase at t = 1 (d_diffuse = 0) from P(0|0) =
%                      Omega_star_init. A non-zero Omega_inf_init (user-given
%                      partially diffuse prior) is not supported and errors.
%     T                sample length
%     dr               Dynare decision-rule structure carrying (set by caller):
%                        restrict_var_list  (= bayestopt_.smoother_var_list)
%                        restrict_columns   (= bayestopt_.smoother_restrict_columns)
%                        mf                 (= bayestopt_.mf, obs positions)
%                        keep_states        (pruning mask from f_Extract...)
%                        qz_crit            (unit-root threshold, e.g. 1+1e-6)
%     tol              numerical tolerance (default 1e-10)
%
%   OUTPUT: Sequence, a struct with
%     .K_stacked(:,:,t)   diffuse-phase FILTERED update gains on the stacked
%                         state, one column per observable (univariate DK)
%     .F_used(i,t)        innovation variance used for obs i at t
%     .obs_is_diffuse(i,t) true if obs i at t was updated with F_inf
%     .K_gain(:,:,t)      standard-phase PREDICTIVE gains on x
%                         (zero during the diffuse phase)
%     .V(:,:,t)           innovation covariance at t
%     .Omega(:,:,t)       predictive state covariance P(t|t-1) (slice 1 = P(1|0))
%     .Omega_star/.Omega_inf  finite/diffuse parts of the covariance
%     .loglik_contribut   DIFFUSE-PHASE log-likelihood contribution only:
%                         sum of -0.5*log(F) terms (callers must add the
%                         quadratic terms for non-diffuse observations, and
%                         the full Gaussian term for post-diffuse periods).
%     .d_diffuse          last diffuse period (0 if none)
%     .filter_type{t}     'diffuse' or 'standard'
%     .T_stacked, .Z, .n_obs_only, .our_idx, .Pinf_stacked, .Pstar_stacked
%
%   Changes vs. the original version:
%     - LOG-LIKELIHOOD FIX: the standard-phase -0.5*log|V| terms are no
%       longer accumulated into loglik_contribut. Both call sites
%       (compute_state_reliability_diagnostics and kalman_filter_steady_and_diffuse)
%       add the full Gaussian term -0.5*v'V^{-1}v - 0.5*log|V| for
%       t >= t_proper themselves, so the old behavior double counted
%       -0.5*log|V| over the post-diffuse sample.
%     - removed dead code: a hand-built stacked system that was immediately
%       overwritten, an obsolete state-index computation, and a debug block
%       referencing an undefined variable.
%     - the -0.5*log(2*pi) constants are omitted throughout (as before).
%
%   See also: compute_Pinf_Pstar, kalman_transition_matrix,
%             compute_state_reliability_diagnostics.

    if nargin < 10 || isempty(tol), tol = 1e-10; end

    [m_abcd, ~] = size(A);       % full ABCD state dimension
    n = size(C, 1);              % number of observables

    %% ------------------------------------------------------------------
    %  1. Stacked system from Dynare's own machinery
    % -------------------------------------------------------------------
    iv = dr.restrict_var_list;   % variables kept by the smoother
    ic = dr.restrict_columns;    % positions of the states within iv

    [dyn_T, dyn_R] = kalman_transition_matrix(dr, iv, ic);

    % Positions (within dyn_T) of the states and of the "pure" observables
    % (observables that are not themselves states).
    state_idx = dr.restrict_columns(:)';
    m_full    = length(state_idx);
    obs_idx   = dr.mf(:)';
    obs_only  = obs_idx(~ismember(obs_idx, state_idx));
    our_idx   = [obs_only(:)', state_idx(:)'];

    % Diffuse initialization, computed exactly as Dynare does it.
    [Pstar_dynare, Pinf_dynare] = compute_Pinf_Pstar(dr.mf, dyn_T, dyn_R, diag(Qscale(:,1)), dr.qz_crit);

    T_stacked     = dyn_T(our_idx, our_idx);
    R_stacked     = dyn_R(our_idx, :);
    Pinf_stacked  = Pinf_dynare(our_idx, our_idx);
    Pstar_stacked = Pstar_dynare(our_idx, our_idx);

    %% ------------------------------------------------------------------
    %  2. Observation selection matrix Z on the stacked state
    % -------------------------------------------------------------------
    n_obs      = length(obs_idx);
    n_obs_only = sum(~ismember(obs_idx, state_idx));
    n_stacked  = n_obs_only + m_full;

    Z = zeros(n_obs, n_stacked);
    for i = 1:n_obs
        pos = find(our_idx == obs_idx(i));
        Z(i, pos) = 1;
    end

    %% ------------------------------------------------------------------
    %  3. Align the stacked system with the pruned ABCD system:
    %     (a) drop the states removed by compute_state_reliability_diagnostics,
    %     (b) permute DR ordering -> declaration ordering (that of A,B,C,D).
    % -------------------------------------------------------------------
    % Map each stacked state position to its declaration-order index
    state_var_dr = dr.inv_order_var(dr.state_var);
    state_idx_dr = iv(state_idx);           % convert to DR positions
    stacking_to_decl = zeros(m_full, 1);
    for k = 1:m_full
        stacking_to_decl(k) = find(state_var_dr == state_idx_dr(k));
    end

    % (a) prune to the same states as the ABCD system
    keep_in_stacking = dr.keep_states(stacking_to_decl);
    kept_positions   = find(keep_in_stacking);

    prune_sel     = [1:n_obs_only, n_obs_only + kept_positions(:)'];
    T_stacked     = T_stacked(prune_sel, prune_sel);
    R_stacked     = R_stacked(prune_sel, :);
    Pinf_stacked  = Pinf_stacked(prune_sel, prune_sel);
    Pstar_stacked = Pstar_stacked(prune_sel, prune_sel);
    Z             = Z(:, prune_sel);

    % (b) permute the state block into declaration order
    kept_decl_idx = stacking_to_decl(kept_positions);
    [~, perm]     = sort(kept_decl_idx);
    full_perm     = [1:n_obs_only, n_obs_only + perm(:)'];

    T_stacked     = T_stacked(full_perm, full_perm);
    R_stacked     = R_stacked(full_perm, :);
    Pinf_stacked  = Pinf_stacked(full_perm, full_perm);
    Pstar_stacked = Pstar_stacked(full_perm, full_perm);
    Z             = Z(:, full_perm);
    n_stacked     = size(T_stacked, 1);

    %% ------------------------------------------------------------------
    %  4. Extend the stacking with augmented states (if the ABCD system was
    %     augmented with variables of interest that are not Dynare states,
    %     and/or with the lag and change states of augment_ssm; they live
    %     in the bottom rows of A and B, in that order).
    % -------------------------------------------------------------------
    m_stacking = n_stacked - n_obs_only;    % Dynare states in the stacking
    n_aug      = m_abcd - m_stacking;       % augmented state count

    if n_aug > 0
        % CHANGED (2026-09-15): the appended rows used to be copied with their
        % Dynare-state columns only, A(m_stacking+1:end, 1:m_stacking). That
        % was exact for the non-state variables of interest (their rows load
        % on Dynare states alone) but drops the couplings among appended
        % states that augment_ssm creates - Lag_k loads on Lag_{k-1} and
        % d_h on Lag_{h-1}. Copying the full appended rows is exact in both
        % cases; for the old augmentation the extra block is identically
        % zero, so nothing changes there.
        A_aug = A(m_stacking+1:end, :);              % n_aug x m_abcd
        B_aug = B(m_stacking+1:end, :);              % n_aug x n_shocks

        n_stacked_base   = n_stacked;
        Pinf_stacked_base = Pinf_stacked;

        T_ext = zeros(n_stacked + n_aug);
        T_ext(1:n_stacked, 1:n_stacked) = T_stacked;
        T_ext(n_stacked+1:end, n_obs_only+1:end) = A_aug;
        T_stacked = T_ext;

        R_stacked = [R_stacked; B_aug];
        Z         = [Z, zeros(size(Z,1), n_aug)];
        n_stacked = size(T_stacked, 1);

        % Stationary part of the initialization: recomputed by Dynare's
        % routine on the extended system (the appended states need their own
        % Pstar rows). Its Pinf is NOT used - see below.
        mf_ext = zeros(size(Z,1), 1);
        for i = 1:size(Z,1)
            mf_ext(i) = find(Z(i,:));
        end
        [Pstar_stacked, Pinf_ext_dynare] = compute_Pinf_Pstar(...
            mf_ext, T_stacked, R_stacked, diag(Qscale(:,1)), dr.qz_crit);

        % CHANGED (2026-09-15): diffuse part tied to Dynare's own Pinf.
        % compute_Pinf_Pstar sets Pinf to the orthogonal projector on the
        % unit-root subspace of the system it is given; on the EXTENDED
        % system that projector is re-normalised by the appended components,
        % so its base block no longer equals Dynare's Pinf for the model.
        % The filtered states at t = 1 and the diffuse-phase -0.5*log(F_inf)
        % terms of the likelihood then differ from Dynare's (measured on
        % HLW: 1e2 in the t = 1 states, 2.05 in the log-likelihood), while
        % everything from t = 2 on is unaffected. The appended states are
        % exact linear identities in the history of x, so their diffuse
        % component is a fixed linear map of the base states' diffuse
        % component: at the unit root z = A_ab x + A_aa z, i.e.
        %     z_diffuse = E x_diffuse,   E = (I - A_aa) \ A_ab,
        % (A_aa is nilpotent - lag chains, changes, static non-state
        % variables - so I - A_aa is invertible; a lag of x has E = e_i',
        % a change has E v = 0 on every unit-root direction v, a non-state
        % variable has E = its ghx row). The exact diffuse prior of the
        % extended system is therefore  [I; E] Pinf [I; E]', whose base block
        % IS Dynare's Pinf: Z Pinf Z' and the base rows of Pinf Z' are
        % unchanged, hence the base filter and the diffuse log-likelihood
        % constants coincide with Dynare's exactly. The stationary part
        % Pstar from compute_Pinf_Pstar is kept: the exact diffuse filter
        % depends on Pstar only through its component outside range(Pinf)
        % (the post-diffuse quantities are the GLS treating the diffuse
        % directions as fixed unknowns), which is what makes the two pieces
        % combinable. Verified against oo_.UpdatedVariables (merge_check).
        A_ab = A(m_stacking+1:end, 1:m_stacking);
        A_aa = A(m_stacking+1:end, m_stacking+1:end);
        E    = (eye(n_aug) - A_aa) \ A_ab;                          % n_aug x m_stacking
        Ext  = [eye(n_stacked_base); zeros(n_aug, n_obs_only), E];  % n_stacked x n_stacked_base
        Pinf_stacked = Ext * Pinf_stacked_base * Ext';
        Pinf_stacked = 0.5 * (Pinf_stacked + Pinf_stacked');
        Sequence.Pinf_ext_dynare = Pinf_ext_dynare;   % what compute_Pinf_Pstar would have used, for auditing
    end

    m = m_abcd;   % from here on, n_stacked = n_obs_only + m

    %% ------------------------------------------------------------------
    %  4b. ADDED 2026-09-26: proper (finite, non-diffuse) initialization.
    %      Replaces the Dynare-built Pstar/Pinf of steps 1-4 by the caller's
    %      P(0|0) on the state block and a zero diffuse part. The y-only
    %      slots and their cross terms are left at zero: with Pinf = 0 the
    %      diffuse phase never opens (rank(Z Pinf Z') = 0 at t = 1), and the
    %      standard phase reads only the state block Omega_star_0 below.
    % -------------------------------------------------------------------
    if ~isempty(Omega_star_init) || ~isempty(Omega_inf_init)
        if isempty(Omega_star_init)
            error('durbin_koopman_diffuse_filter:init', 'Omega_inf_init without Omega_star_init: a purely diffuse user prior is not supported.');
        end
        if ~isequal(size(Omega_star_init), [m, m])
            error('durbin_koopman_diffuse_filter:init', 'Omega_star_init must be %d x %d (pruned ABCD state order).', m, m);
        end
        if ~isempty(Omega_inf_init) && any(abs(Omega_inf_init(:)) > 0)
            error('durbin_koopman_diffuse_filter:init', 'a non-zero Omega_inf_init (partially diffuse user prior) is not supported; pass [] or zeros.');
        end
        Pstar_stacked = zeros(n_stacked);
        Pstar_stacked(n_obs_only+1:end, n_obs_only+1:end) = 0.5 * (Omega_star_init + Omega_star_init');
        Pinf_stacked  = zeros(n_stacked);
        Sequence.init_type = 'proper';
    else
        Sequence.init_type = 'diffuse';
    end

    %% ------------------------------------------------------------------
    %  5. Pre-allocate outputs
    % -------------------------------------------------------------------
    Sequence.our_idx        = our_idx;
    Sequence.K_stacked      = zeros(n_stacked, n, T);
    Sequence.F_used         = zeros(n, T);
    Sequence.obs_is_diffuse = false(n, T);
    Sequence.K_gain         = zeros(m, n, T);
    Sequence.V              = zeros(n, n, T);
    Sequence.Omega          = zeros(m, m, T+1);
    Sequence.Omega_star     = zeros(m, m, T+1);
    Sequence.Omega_inf      = zeros(m, m, T+1);
    Sequence.loglik_contribut = 0;
    Sequence.filter_type    = cell(T, 1);
    Sequence.Pinf_stacked   = Pinf_stacked;
    Sequence.Pstar_stacked  = Pstar_stacked;
    Sequence.Z              = Z;

    %% ------------------------------------------------------------------
    %  6. Initialize (state block of the stacked initialization)
    % -------------------------------------------------------------------
    Omega_star_0 = Pstar_stacked(n_obs_only+1:end, n_obs_only+1:end);  % m x m
    Omega_inf_0  = Pinf_stacked(n_obs_only+1:end, n_obs_only+1:end);   % m x m

    Sequence.Omega(:,:,1)      = Omega_star_0 + Omega_inf_0;
    Sequence.Omega_star(:,:,1) = Omega_star_0;
    Sequence.Omega_inf(:,:,1)  = Omega_inf_0;

    Pstar = Pstar_stacked;      % stacked covariances, updated in the loop
    Pinf  = Pinf_stacked;

    d_diffuse        = 0;
    Omega_star_state = Omega_star_0;   % filtered P_xx handed to the standard phase

    %% ------------------------------------------------------------------
    %  7. Main covariance/gain recursion
    % -------------------------------------------------------------------
    for tt = 1:T

        % Time-varying shock covariance for this period
        Q_t        = diag(Qscale(:, tt));
        QQ_stacked = R_stacked * Q_t * R_stacked';
        BD = B * Q_t * D';
        BB = B * Q_t * B';
        DD = D * Q_t * D';

        % Still in the diffuse phase while the diffuse part of the
        % innovation variance has positive rank.
        newRank    = rank(Z * Pinf * Z', tol);
        is_diffuse = (newRank > 0);

        if is_diffuse
            d_diffuse = tt;
            Sequence.filter_type{tt} = 'diffuse';

            % --- univariate (sequential) DK diffuse update on [y; x] ---
            [K_stacked_current, Pstar, Pinf, loglik_t, F_used, obs_is_diffuse] = ...
                diffuse_step_stacked(Pstar, Pinf, Z, tol);

            Sequence.K_stacked(:,:,tt)    = K_stacked_current;  % filtered gains
            Sequence.F_used(:,tt)         = F_used;
            Sequence.obs_is_diffuse(:,tt) = obs_is_diffuse;

            % Hand-over covariance: FILTERED P_xx(tt|tt), extracted before
            % the prediction step overwrites Pstar.
            Omega_star_state = Pstar(n_obs_only+1:end, n_obs_only+1:end);

            % Prediction step: advance stacked covariances to tt+1
            Pstar = T_stacked * Pstar * T_stacked' + QQ_stacked;
            Pinf  = T_stacked * Pinf  * T_stacked';

            % Store predicted state covariances (prior at tt+1)
            Omega_star = Pstar(n_obs_only+1:end, n_obs_only+1:end);
            Omega_inf  = Pinf(n_obs_only+1:end, n_obs_only+1:end);

            Sequence.V(:,:,tt) = Z * (Pstar + Pinf) * Z';
            Sequence.Omega(:,:,tt+1)      = f_make_psd_local(Omega_star + Omega_inf);
            Sequence.Omega_star(:,:,tt+1) = f_make_psd_local(Omega_star);
            Sequence.Omega_inf(:,:,tt+1)  = f_make_psd_local(Omega_inf);

            % Diffuse-phase log-likelihood contribution (-0.5*log F terms;
            % quadratic terms for non-diffuse obs are added by the caller).
            Sequence.loglik_contribut = Sequence.loglik_contribut + loglik_t;

        else
            Sequence.filter_type{tt} = 'standard';

            % --- standard multivariate recursion on x alone ---
            Omega_star = Omega_star_state;

            V = C * Omega_star * C' + DD;
            V = 0.5 * (V + V');

            K_state = (A * Omega_star * C' + BD) / V;   % predictive gain
            Sequence.K_gain(:,:,tt) = K_state;

            Omega_star_state = A * Omega_star * A' + BB - K_state * V * K_state';

            Sequence.V(:,:,tt) = V;
            Sequence.Omega(:,:,tt+1)      = f_make_psd_local(Omega_star_state);
            Sequence.Omega_star(:,:,tt+1) = f_make_psd_local(Omega_star_state);
            Sequence.Omega_inf(:,:,tt+1)  = zeros(m);

            % NOTE (log-likelihood fix): no contribution is accumulated here.
            % Callers add the full Gaussian term -0.5*v'inv(V)v - 0.5*log|V|
            % for every post-diffuse period; the old code also added
            % -0.5*log|V| here, double counting it.
        end
    end

    Sequence.d_diffuse  = d_diffuse;
    Sequence.n_obs_only = n_obs_only;
    Sequence.T_stacked  = T_stacked;
end


%% ========================================================================
function [K_current, Pstar, Pinf, loglik, F_used, is_diffuse] = ...
    diffuse_step_stacked(Pstar, Pinf, Z, tol)
%DIFFUSE_STEP_STACKED  Univariate Durbin-Koopman diffuse update, stacked state.
%
%   Processes the observations one at a time (valid because Z is a pure
%   selection matrix on the stacked state, so measurement errors are
%   uncorrelated across the sequential updates by construction).
%
%   Outputs:
%     K_current  - FILTERED update gains: alpha(t|t,i) = alpha(t|t,i-1) + K_i v_i.
%                  For diffuse obs K_i = Pinf Z_i'/F_inf; otherwise
%                  Pstar Z_i'/F_star.  (m_stacked x n_obs)
%     F_used     - innovation variance used per observation (n_obs x 1)
%     is_diffuse - true where the diffuse variance F_inf was used
%     loglik     - ONLY the -0.5*log(F) terms. The caller adds the quadratic
%                  term -0.5*v^2/F for NON-diffuse observations; no quadratic
%                  term enters for diffuse observations (DK 2012, eq. 7.14).

    n         = size(Z, 1);
    m_stacked = size(Pstar, 1);

    K_current  = zeros(m_stacked, n);
    F_used     = zeros(n, 1);
    is_diffuse = false(n, 1);
    loglik     = 0;

    for i = 1:n
        Zi = Z(i, :);

        f_inf  = Zi * Pinf * Zi';
        f_star = Zi * Pstar * Zi';
        k_inf  = Pinf * Zi';
        k_star = Pstar * Zi';

        if f_inf > tol
            % DIFFUSE observation: update both Pstar and Pinf (DK 2012, 5.7)
            K_current(:, i) = k_inf / f_inf;
            F_used(i)       = f_inf;
            is_diffuse(i)   = true;

            Pstar = Pstar ...
                + k_inf * (k_inf' / f_inf) * (f_star / f_inf) ...
                - k_star * (k_inf' / f_inf) ...
                - (k_inf / f_inf) * k_star';
            Pinf = Pinf - k_inf * k_inf' / f_inf;

            loglik = loglik - 0.5 * log(f_inf);

        elseif f_star > tol
            % STANDARD observation encountered within the diffuse phase
            K_current(:, i) = k_star / f_star;
            F_used(i)       = f_star;
            is_diffuse(i)   = false;

            Pstar = Pstar - k_star * k_star' / f_star;

            loglik = loglik - 0.5 * log(f_star);

        else
            % SINGULAR observation: carries no information, skip
            K_current(:, i) = zeros(m_stacked, 1);
            F_used(i)       = Inf;
            is_diffuse(i)   = false;
        end

        % Enforce symmetry after every rank-one update
        Pstar = 0.5 * (Pstar + Pstar');
        Pinf  = 0.5 * (Pinf + Pinf');
    end
end


%% ========================================================================
function M = f_make_psd_local(M)
%F_MAKE_PSD_LOCAL  Project a symmetric matrix onto the PSD cone.
%   Clips negative eigenvalues at zero (they arise from round-off in the
%   covariance recursions) and re-symmetrizes.
    M = 0.5 * (M + M');
    [V, D] = eig(M);
    d = max(real(diag(D)), 0);
    M = real(V * diag(d) * V');
    M = 0.5 * (M + M');
end
