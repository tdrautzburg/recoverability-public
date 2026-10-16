function fr = fr_regression_empirical(all_states, A, base_idx, lag_idx, m_grid, H, IDX_ok, fr)
%FR_REGRESSION_EMPIRICAL  OLS of the smoother revision on the filter revision, on data.
%
%   fr = FR_REGRESSION_EMPIRICAL(all_states, A, base_idx, lag_idx, m_grid, H, IDX_ok, fr)
%
%   Sample counterpart of fr_regression_theory, on the filtered states of
%   ONE data pass (the diffuse / time-varying filter in
%   compute_state_reliability_diagnostics). For every base state i, horizon
%   h = 1..H and smoother lag m in m_grid it runs
%
%       f_{i,h,m,t} = alpha + beta * r_{i,h,t} + error,     t = 1..T-h-m
%
%       r_{i,h,t}   = xhat_{i,t+h|t+h}   - xhat_{i,t+h|t}
%       f_{i,h,m,t} = xhat_{i,t+h|t+h+m} - xhat_{i,t+h|t}
%
%   where xhat_{t+h|t} = A^h xhat_{t|t} (row t of all_states*(A^h)'),
%   xhat_{i,t+h|t+h} is row t+h of all_states in column i, and
%   xhat_{i,t+h|t+h+m} is row t+h+m of all_states in the column of the
%   lag-augmented state l<m>NAME (the m-period fixed-lag smoothed estimate
%   of NAME_{t+h}).
%
%   SAMPLE. An observation t is used only if EVERY period in the window
%   t .. t+h+m is admissible under IDX_ok (the engine's common mask: burn-in
%   and heteroskedastic periods excluded), so that the forecast origin, the
%   filter update at t+h and the whole smoothing window are clean. At least
%   3 usable observations and a non-degenerate regressor are required;
%   otherwise the cell is left NaN.
%
%   INPUTS
%     all_states  T x n_states filtered states xhat(t|t), pruned augmented
%                 state vector (as built in the engine's step 6b)
%     A           pruned augmented transition matrix
%     base_idx    J x 1 column indices of the base states
%     lag_idx     J x n_m column indices of Lag_m of each base state
%     m_grid      1 x n_m smoother lags m >= 1
%     H           largest forecast horizon (h = 1..H)
%     IDX_ok      1 x T (or T x 1) logical, admissible periods
%     fr          struct to fill (see fr_regression_empty); the fields
%                 R2_emp, alpha_emp, beta_emp, N_emp (J x H x n_m) are written.
%
%   See also: fr_regression_theory, fr_regression_empty,
%             compute_state_reliability_diagnostics.

[T, n] = size(all_states);
J   = numel(base_idx);
n_m = numel(m_grid);

bad = double(~IDX_ok(:)');                 % 1 x T
cs  = [0, cumsum(bad)];                    % cs(t+1) = number of bad periods in 1..t

Ah = eye(n);
for h = 1:H
    Ah = A * Ah;                           % A^h
    Xh = all_states * Ah';                 % row t: xhat_{t+h|t}'

    for j = 1:J
        kk = base_idx(j);
        fc = Xh(:, kk);                    % forecast of state i made at t

        for im = 1:n_m
            m  = m_grid(im);
            lm = lag_idx(j, im);
            Tn = T - h - m;                % last usable origin
            if Tn < 3, continue; end

            t  = (1:Tn)';
            r  = all_states(t + h,     kk) - fc(t);
            f  = all_states(t + h + m, lm) - fc(t);

            % window t..t+h+m must be free of excluded periods
            nbad  = cs(t + h + m + 1) - cs(t);
            valid = (nbad(:) == 0) & isfinite(r) & isfinite(f);
            N = sum(valid);
            fr.N_emp(j, h, im) = N;
            if N < 3, continue; end

            rv = r(valid);  fv = f(valid);
            rbar = mean(rv);  fbar = mean(fv);
            Srr = sum((rv - rbar).^2);
            Sff = sum((fv - fbar).^2);
            Srf = sum((rv - rbar) .* (fv - fbar));
            if Srr <= eps * max(1, sum(rv.^2)) || Sff <= eps * max(1, sum(fv.^2))
                continue;                  % degenerate regressor / regressand: NaN
            end
            beta  = Srf / Srr;
            fr.beta_emp(j, h, im)  = beta;
            fr.alpha_emp(j, h, im) = fbar - beta * rbar;
            fr.R2_emp(j, h, im)    = Srf^2 / (Srr * Sff);
        end
    end
end

end
