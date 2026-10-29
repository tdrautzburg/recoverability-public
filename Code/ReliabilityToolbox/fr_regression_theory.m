function fr = fr_regression_theory(A, K, V, Omega, base_idx, lag_idx, m_grid, H, fr)
%FR_REGRESSION_THEORY  Population moments of the smoother-revision regression.
%
%   fr = FR_REGRESSION_THEORY(A, K, V, Omega, base_idx, lag_idx, m_grid, H, fr)
%
%   CHANGED 2026-09-17: 4th input Omega (steady-state filtered covariance,
%   ModelMatrices.KF_Omega) added; it fills fr.omega_ratio and fr.R2_h_base
%   (see fr_regression_empty), the two objects that turn an empirical P^2 into
%   an implied R^2 (fr_implied_R2). Nothing else changed.
%
%   For each base state i (row base_idx(j) of the PRUNED, lag-augmented state
%   vector), forecast horizon h = 1..H and smoother lag m = m_grid(im), the
%   population R2 and slope of the OLS regression
%
%       f_{i,h,m,t} = alpha + beta * r_{i,h,t},
%
%       r_{i,h,t}   = xhat_{i,t+h|t+h}   - xhat_{i,t+h|t}     (filter revision)
%       f_{i,h,m,t} = xhat_{i,t+h|t+h+m} - xhat_{i,t+h|t}     (smoother revision)
%
%   with xhat_{t+h|t} = A^h xhat_{t|t}, under the STEADY-STATE filter.
%
%   DERIVATION. With innovations nu_t = y_t - C xhat_{t-1|t-1} (white, with
%   covariance V) the steady-state filter is xhat_{t|t} = A xhat_{t-1|t-1} + K nu_t,
%   so for the augmented system
%
%       xhat_{t+s|t+s} = A^s xhat_{t|t} + sum_{j=1}^{s} A^{s-j} K nu_{t+j}.
%
%   Taking row i at s = h gives r; taking the l<m>NAME row at s = h+m gives f,
%   because l<m>NAME_{t+h+m} = NAME_{t+h} and the deterministic part
%   e_Lm' A^{h+m} xhat_{t|t} equals e_i' A^h xhat_{t|t} (the lag chain is an
%   exact identity, and the base block of the augmented filter coincides
%   with the un-augmented one since the appended states carry no
%   information). Hence, with G_s = A^s K,
%
%       r_{h}   = sum_{s=0}^{h-1}   e_i'  G_s nu_{t+h-s}
%       f_{h,m} = sum_{s=0}^{h+m-1} e_Lm' G_s nu_{t+h+m-s}
%
%   and, nu being white, both are stationary for ANY A (unit roots included):
%
%       Var(r)    = sum_{s=0}^{h-1}   e_i'  G_s V G_s'     e_i
%       Var(f)    = sum_{s=0}^{h+m-1} e_Lm' G_s V G_s'     e_Lm
%       Cov(f, r) = sum_{s=0}^{h-1}   e_Lm' G_{s+m} V G_s' e_i
%
%   (the innovation nu_{t+h-s} enters r through G_s and f through G_{s+m}).
%   Population alpha is 0 (both revisions have mean zero); beta = Cov/Var(r);
%   R2 = Cov^2 / (Var(f) Var(r)).
%
%   Var(r) = V_h - Omega_ii (the h-step forecast-error variance less the
%   filtering-error variance), so the base state's eq.-(3.11) reliability is
%   R^2_{i,h} = 1 - Omega_ii / V_h = Var(r) / (Var(r) + Omega_ii), stored as
%   fr.R2_h_base; and [Omega_m]_ii is the diagonal entry of Omega at the
%   l<m>NAME row, stored relative to Omega_ii as fr.omega_ratio.
%
%   INPUTS
%     A, K, V    pruned augmented transition matrix, steady-state gain and
%                innovation covariance (ModelMatrices.A, .K_gain, .KF_V)
%     Omega      steady-state filtered state covariance (ModelMatrices.KF_Omega)
%     base_idx   J x 1 row indices of the base states
%     lag_idx    J x n_m row indices of l<m>NAME of each base state, one
%                column per entry of m_grid
%     m_grid     1 x n_m smoother lags m >= 1
%     H          largest forecast horizon (h = 1..H)
%     fr         struct to fill (see fr_regression_empty); the fields
%                R2_theory, beta_theory, cov_fr_theory, var_r_theory,
%                var_f_theory (J x H x n_m) are written.
%
%   See also: fr_regression_empirical, fr_regression_empty,
%             compute_state_reliability_diagnostics.

J   = numel(base_idx);
n_m = numel(m_grid);
n_obs = size(K, 2);
n   = size(A, 1);

% G_s = A^s K and G_s V, for s = 0 .. H + max(m) - 1 (slice s+1 holds s)
S  = H + max(m_grid);
G  = zeros(n, n_obs, S);
GV = zeros(n, n_obs, S);
G(:, :, 1)  = K;
GV(:, :, 1) = K * V;
for s = 2:S
    G(:, :, s)  = A * G(:, :, s-1);
    GV(:, :, s) = G(:, :, s) * V;
end

for j = 1:J
    kk = base_idx(j);
    % q_r(s) = e_i' G_s V G_s' e_i,  s = 0..H-1  -> Var(r_h) = sum_{s<h} q_r(s)
    q_r = zeros(1, H);
    for s = 0:H-1
        q_r(s+1) = GV(kk, :, s+1) * G(kk, :, s+1)';
    end
    var_r = cumsum(q_r);                       % 1 x H, entry h
    fr.R2_h_base(j, :) = var_r ./ (var_r + Omega(kk, kk));   % eq. (3.11) for the base state

    for im = 1:n_m
        m  = m_grid(im);
        lm = lag_idx(j, im);
        fr.omega_ratio(j, im) = Omega(lm, lm) / Omega(kk, kk);

        q_f = zeros(1, H + m);                 % Var(f_{h,m}) = sum_{s < h+m} q_f(s)
        for s = 0:H+m-1
            q_f(s+1) = GV(lm, :, s+1) * G(lm, :, s+1)';
        end
        var_f_all = cumsum(q_f);               % entry h+m

        q_c = zeros(1, H);                     % Cov = sum_{s<h} e_Lm' G_{s+m} V G_s' e_i
        for s = 0:H-1
            q_c(s+1) = GV(lm, :, s+m+1) * G(kk, :, s+1)';
        end
        cov_fr = cumsum(q_c);

        for h = 1:H
            vr = var_r(h);  vf = var_f_all(h+m);  cv = cov_fr(h);
            fr.var_r_theory(j, h, im)  = vr;
            fr.var_f_theory(j, h, im)  = vf;
            fr.cov_fr_theory(j, h, im) = cv;
            if vr > eps && vf > eps
                fr.beta_theory(j, h, im) = cv / vr;
                fr.R2_theory(j, h, im)   = cv^2 / (vf * vr);
            end                                % else: left NaN (degenerate)
        end
    end
end

end
