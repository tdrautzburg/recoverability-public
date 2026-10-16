function [Omega, K, V, Omega_tp1_t, ric_flag] = solve_steady_state_riccati(A, B, C, D, opts)
%SOLVE_STEADY_STATE_RICCATI  Steady-state Kalman filter covariance (one solver for the toolbox).
%
%   [Omega, K, V, Omega_tp1_t, ric_flag] = SOLVE_STEADY_STATE_RICCATI(A, B, C, D)
%   [...] = SOLVE_STEADY_STATE_RICCATI(A, B, C, D, opts)
%
%   ABCD form with observables loading on LAGGED states,
%       x(t) = A x(t-1) + B e(t),   y(t) = C x(t-1) + D e(t),   e ~ N(0, I).
%   Omega = Var(x_t | y^t) is the stabilising solution of
%
%       Omega = A Omega A' + B B' - (A Omega C' + B D') (C Omega C' + D D')^{-1} (A Omega C' + B D')',   (*)
%
%   K = (A Omega C' + B D') V^{-1} with V = C Omega C' + D D', and
%   Omega_tp1_t = A Omega A' + B B' (= Var(x_{t+1} | y^t)). Every Riccati fixed
%   point in the toolbox (kalman_filter_steady_and_diffuse, steady_state_kalman_filter,
%   theoretical_R2_wedge) goes through this function since 2026-09-16.
%
%   WHY (2026-09-16). The fixed-point iteration used until now converges at rate
%   rho(A-KC)^2, about 19/(1-rho^2) steps for a 1e-10 relative step: 1.4e3 at
%   the HLW posterior mean, 3e3 at the mode, 3e5 (past the 1e5 cap, minutes)
%   for pile-up draws with sigma_g, sigma_z near zero. Three such solves per
%   draw were about half of the per-draw cost on prior draws and three
%   quarters on posterior draws. (*) is a discrete algebraic Riccati equation
%   and is solved directly here. MATLAB's dare/idare need the Control System
%   Toolbox, which the production machine does not have, so:
%
%   METHODS, tried in this order (opts.method forces one):
%     'invertible'  (ADDED 2026-10-16) square, invertible D (as many shocks as
%                observables, no noise-free direction missing): at Omega = 0, (*)
%                gives B B' - B D' (D D')^{-1} D B' = 0, so Omega = 0 IS a fixed
%                point, with K = B D^{-1} and closed loop A - K C = F = A - B D^{-1} C
%                (the eigenvalues of F are the transmission zeros). Hence
%                  rho(F) < 1   Omega = 0 is the stabilising solution ('invertible';
%                               the state is recovered exactly in the limit);
%                  rho(F) = 1   Omega = 0 is the limit of the Riccati recursion (the
%                               maximal solution) but not stabilising: the filter
%                               approaches it only at rate ~1/t and no solver below
%                               can return it - QZ finds no stable deflating
%                               subspace and the iteration hits its cap after ~20 s
%                               with residual ~0.5 (Smets-Wouters 2007 without
%                               measurement error: four eigenvalues of F at 1 from
%                               the differenced observables) ('invertible_critical');
%                  rho(F) > 1   Omega = 0 is a fixed point but not the limit; the
%                               methods below find the stabilising solution.
%                Accepted when rho(F) <= 1 + opts.tol_unit and the absolute residual
%                ||f(0)|| / ||B B'|| is below opts.tol_res_raw. ric_flag.converged is
%                true (the limit is exact); rho_cl = rho(F) flags the critical case.
%     'dynare'   Dynare's kalman_steady_state MEX (SLICOT SB02OD, generalised
%                Schur), which ships with every Dynare build and solves the
%                uncorrelated-noise DARE  P = T P T' + Q - T P Z'(Z P Z' + H)^{-1} Z P T'.
%                (*) has the cross term N = B D'; with R = D D' invertible the
%                standard substitution
%                    Abar = A - N R^{-1} C,   Qbar = B B' - N R^{-1} N'
%                turns (*) into that form with T = Abar, Z = C, H = R (the same
%                Omega solves both). Used when the MEX is on the path and
%                rcond(R) > 1e-10.
%     'qz'       own solver on the extended (2n+ny) pencil of the DARE
%                (Van Dooren 1981 / Arnold-Laub 1984), base MATLAB / Octave qz +
%                ordqz, no invertibility of R or A required. The stable deflating
%                subspace [U1; U2; U3] gives Omega = U2 U1^{-1}.
%     'iterate'  the original fixed-point iteration from opts.Omega0 (or B B'),
%                relative tolerance 1e-10, cap 1e5, kept as the fallback and as
%                the reference implementation.
%   Whatever method produced Omega, the result is ACCEPTED only if it is
%   symmetric, positive semidefinite (to tolerance) and satisfies (*): the
%   relative residual r = ||Omega - f(Omega)|| / ||Omega|| of the fixed-point
%   map f, SCALED by 1/(1 - rho(A-KC)^2), must be below opts.tol_res (1e-8).
%   The scaling matters: the residual is one step of a contraction with rate
%   rho^2, so the distance to the fixed point is about r/(1 - rho^2) - at
%   rho = 0.99997 a capped iteration shows r = 2e-9 while being 4e-5 off.
%   A candidate whose RAW residual is at solver precision (r <= opts.tol_res_raw,
%   1e-10) is accepted as well: at rho = 0.999995 (a prior draw) 1 - rho^2 is
%   1e-5 and an exact solution with r = 1e-12 has an estimate of 1e-7, which
%   no method can beat; rel_step reports that estimate. Otherwise the next
%   method is tried. If all fail, the candidate with the smallest residual is
%   returned with ric_flag.converged = false.
%
%   ric_flag (same four fields the wrapper stores in cell 19, plus two):
%     .converged   true iff the accepted solution passed the scaled residual check
%     .iterations  fixed-point iterations run (0 when a direct solver was used)
%     .rel_step    scaled relative residual r/(1 - rho^2) of the accepted
%                  solution, i.e. an estimate of the relative ERROR (was: last
%                  relative step of the iteration, a smaller number)
%     .rho_cl      rho(A - K C)
%     .method      'invertible' | 'invertible_critical' | 'dynare' | 'qz' | 'iterate'
%     .residuals   raw relative residual of every method that was tried
%                  (for 'invertible': ||f(0)|| / ||B B'||)
%
%   opts: .method ('auto'), .Omega0 (start of the iteration, [] = B*B'),
%         .tol_res (1e-8), .tol_res_raw (1e-10), .MaxIter (1e5),
%         .tol_riccati (1e-10), .tol_unit (1e-6, how close to 1 counts as the
%         unit circle in the 'invertible' check), .verbose (false)
%
%   See also: kalman_filter_steady_and_diffuse, steady_state_kalman_filter,
%             theoretical_R2_wedge.

if nargin < 5, opts = struct(); end
def = struct('method', 'auto', 'Omega0', [], 'tol_res', 1e-8, 'tol_res_raw', 1e-10, 'MaxIter', 1e5, 'tol_riccati', 1e-10, 'tol_unit', 1e-6, 'verbose', false);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = def.(fn{i}); end
end

n  = size(A, 1);
Q  = B * B';   Q = (Q + Q') / 2;
R  = D * D';   R = (R + R') / 2;
N  = B * D';

if strcmp(opts.method, 'auto')
    order = {'dynare', 'qz', 'iterate'};
else
    order = {opts.method};
    if ~strcmp(opts.method, 'iterate'), order{end+1} = 'iterate'; end   % always keep the fallback
end

ric_flag = struct('converged', false, 'iterations', 0, 'rel_step', NaN, 'rho_cl', NaN, 'method', '', 'residuals', struct());
Omega = []; K = []; V = [];

% 'invertible' shortcut (ADDED 2026-10-16, see header): exact Omega = 0 when D is square and
% invertible and rho(A - B D^{-1} C) <= 1. Its own acceptance test: the residual check below
% divides by ||Omega|| and by 1 - rho^2, both zero here.
if any(strcmp(opts.method, {'auto', 'invertible'})) && size(D, 1) == size(D, 2) && rcond(D) > 1e-10
    F_inv  = A - (B / D) * C;
    rho_F  = max(abs(eig(F_inv)));
    f0     = Q - N * (R \ N');                         % f(0) of (*)
    res0   = norm(f0, 'fro') / max(norm(Q, 'fro'), realmin);
    ric_flag.residuals.invertible = res0;
    if opts.verbose, fprintf('solve_steady_state_riccati: invertible check rho(F) = %.8f, residual at Omega = 0 %.2e\n', rho_F, res0); end
    if rho_F <= 1 + opts.tol_unit && res0 <= opts.tol_res_raw
        Omega = zeros(n);
        if rho_F < 1 - opts.tol_unit, ric_flag.method = 'invertible'; else, ric_flag.method = 'invertible_critical'; end
        ric_flag.converged = true;
        ric_flag.rel_step  = res0;
        ric_flag.rho_cl    = rho_F;
        V = R;
        K = N / R;                                     % = B D^{-1}
        Omega_tp1_t = Q;
        return
    end
end
if strcmp(opts.method, 'invertible'), order = {'iterate'}; end   % shortcut not applicable: keep the fallback
last_iter = struct('Omega', [], 'iterations', 0, 'rel_step', NaN);
% (the local function dynare_mex_usable, at the end of this file, remembers whether the MEX loads)
best      = struct('Omega', [], 'res', Inf, 'rho', NaN, 'method', '');

for m = 1:numel(order)
    method = order{m};
    Om = [];
    try
        switch method
            case 'dynare'
                if ~(exist('kalman_steady_state', 'file') == 3 || exist('kalman_steady_state', 'file') == 2), continue; end
                if ~dynare_mex_usable(), continue; end   % a MEX that cannot be loaded costs ~50 ms per failed attempt (2026-09-28)
                if rcond(R) < 1e-10, continue; end
                RiN  = R \ N';                        % R^{-1} N'
                Abar = A - N * (R \ C);               % A - N R^{-1} C
                Qbar = Q - N * RiN;   Qbar = (Qbar + Qbar') / 2;
                Qbar = psd_project(Qbar);             % round-off can leave it indefinite
                Om = kalman_steady_state(Abar', Qbar, C', R);   % call pattern of Dynare's get_pstar.m
            case 'qz'
                Om = dare_qz(A, C, Q, R, N);
            case 'iterate'
                [Om, it, step] = riccati_iterate(A, B, C, D, opts.Omega0, opts.tol_riccati, opts.MaxIter);
                last_iter = struct('Omega', Om, 'iterations', it, 'rel_step', step);
        end
    catch err
        if opts.verbose, fprintf('solve_steady_state_riccati: %s failed (%s)\n', method, err.message); end
        Om = [];
    end
    if isempty(Om) || ~all(isfinite(Om(:))), ric_flag.residuals.(method) = NaN; continue; end
    Om = (Om + Om') / 2;
    [res, rho] = riccati_residual(Om, A, C, Q, R, N);
    ric_flag.residuals.(method) = res;
    err_est = res / max(1 - rho^2, 1e-12);        % distance to the fixed point, not one step of it
    if opts.verbose, fprintf('solve_steady_state_riccati: %s residual %.2e, rho %.6f, error estimate %.2e\n', method, res, rho, err_est); end
    psd_ok = min(eig(Om)) >= -opts.tol_res * max(1, norm(Om, 'fro'));
    % Two ways to pass (ADDED 2026-09-16, after a prior draw with rho = 0.999995):
    % (a) the error estimate is below tol_res, or (b) the raw residual is at
    % solver precision (tol_res_raw, 1e-10). Near a unit closed-loop root
    % 1 - rho^2 is ~1e-5, so even a machine-precision residual of 1e-12 gives
    % an error estimate of 1e-7 > tol_res; rejecting it would discard an exact
    % solution and run the fallback iteration for 1e5 steps to a far worse one
    % (residual 3e-6). rel_step still reports the estimate, so the caller sees
    % the attainable precision.
    if psd_ok && (err_est <= opts.tol_res || res <= opts.tol_res_raw)
        Omega = Om;
        ric_flag.converged = true;
        ric_flag.method    = method;
        ric_flag.rel_step  = err_est;
        ric_flag.rho_cl    = rho;
        if strcmp(method, 'iterate'), ric_flag.iterations = last_iter.iterations; end
        break;
    end
    if psd_ok && (isempty(best.Omega) || res < best.res)   % remember the best rejected candidate
        best = struct('Omega', Om, 'res', res, 'rho', rho, 'method', method);
    end
end

if isempty(Omega)
    % every method failed the check: return the candidate with the smallest
    % residual (CHANGED 2026-09-16: was always the iteration's last iterate,
    % which can be the worst of the three), flagged
    if isempty(best.Omega)
        if isempty(last_iter.Omega)
            [Om, it, step] = riccati_iterate(A, B, C, D, opts.Omega0, opts.tol_riccati, opts.MaxIter);
            last_iter = struct('Omega', Om, 'iterations', it, 'rel_step', step);
        end
        Om = (last_iter.Omega + last_iter.Omega') / 2;
        [res, rho] = riccati_residual(Om, A, C, Q, R, N);
        best = struct('Omega', Om, 'res', res, 'rho', rho, 'method', 'iterate');
    end
    Omega = best.Omega;
    ric_flag.method     = best.method;
    ric_flag.iterations = last_iter.iterations * strcmp(best.method, 'iterate');
    ric_flag.rel_step   = best.res / max(1 - best.rho^2, 1e-12);
    ric_flag.rho_cl     = best.rho;
    ric_flag.converged  = false;
end

V = C * Omega * C' + R;   V = (V + V') / 2;
K = (A * Omega * C' + N) / V;
Omega_tp1_t = A * Omega * A' + Q;
if ~ric_flag.converged
    warning('riccati:noconv', ...
        'solve_steady_state_riccati: no method met the residual tolerance (best %s, residual %.3e, rho(A-KC) = %.6f).', ...
        ric_flag.method, ric_flag.rel_step, ric_flag.rho_cl);
end
end


%% ========================================================================
function [res, rho] = riccati_residual(Om, A, C, Q, R, N)
% relative residual of (*) and rho(A - K C) at Om
V  = C * Om * C' + R;
S  = A * Om * C' + N;
K  = S / V;
F  = A * Om * A' + Q - K * S';
res = norm(Om - F, 'fro') / max(norm(Om, 'fro'), realmin);
if all(isfinite(K(:)))
    rho = max(abs(eig(A - K * C)));
else
    rho = NaN;
end
end

function [Omega, ii, change_norm] = riccati_iterate(A, B, C, D, Omega0, tol_riccati, MaxIter)
% ORIGINAL fixed-point iteration (kalman_filter_steady_and_diffuse, 2026-09 version):
% relative stopping rule, non-finite guard. Default start B*B' (PSD, keeps V
% non-singular when D*D' is); callers that want the old start pass Omega0.
if isempty(Omega0), Omega = B * B'; else, Omega = Omega0; end
change_norm = Inf;  ii = 0;
while change_norm > tol_riccati && ii < MaxIter
    ii = ii + 1;
    V        = C*Omega*C' + D*D';
    K_gain   = (A*Omega*C' + B*D') / V;
    NewOmega = A*Omega*A' + B*B' - K_gain*V*K_gain';
    NewOmega = (NewOmega + NewOmega')/2;
    if ~all(isfinite(NewOmega(:)))
        break                                   % keep the last finite iterate
    end
    change_norm = norm(NewOmega - Omega, 'fro') / max([norm(Omega, 'fro'), norm(NewOmega, 'fro'), realmin]);
    Omega = NewOmega;
end
end

function X = dare_qz(A, C, Q, R, N)
% Stabilising solution of  X = A X A' + Q - (A X C' + N)(C X C' + R)^{-1}(A X C' + N)'
% via the extended pencil (Van Dooren 1981). Write it in control form with
% F = A', G = C', S = N, i.e. X = F'XF + Q - (F'XG + S)(R + G'XG)^{-1}(F'XG + S)'.
% The (2n+p) x (2n+p) pencil  lambda*E - M  with
%     M = [F  0  G ;  -Q  I  -S ;  S'  0  R],   E = [I 0 0 ;  0 F' 0 ;  0 -G' 0]
% has n finite eigenvalues inside the unit circle (those of the closed loop),
% n outside (their reciprocals) and p at infinity; the deflating subspace of
% the inside ones, spanned by the first n columns Z(:,1:n) = [U1; U2; U3] of
% the ordered generalised Schur form, gives X = U2 * inv(U1).
n = size(A, 1);  p = size(C, 1);
F = A';  G = C';  S = N;
M = [F,            zeros(n),  G;
     -Q,           eye(n),   -S;
     S',           zeros(p, n), R];
E = [eye(n),       zeros(n),  zeros(n, p);
     zeros(n),     F',        zeros(n, p);
     zeros(p, n),  -G',       zeros(p)];
[AA, BB, ~, Z] = qz(M, E, 'real');
% eigenvalues of the pencil: alpha./beta; "inside" = |alpha| < |beta| (finite, in the unit disc)
alpha = diag(AA);  beta = diag(BB);
% ordqz needs a logical selector for the eigenvalues to move to the top; use the
% 2x2-block-safe eigenvalue extraction via ordeig where available
if exist('ordeig', 'file')
    ev = ordeig(AA, BB);
else
    ev = alpha ./ beta;
end
inside = abs(ev) < 1;
if sum(inside) ~= n
    error('dare_qz:count', 'pencil has %d eigenvalues inside the unit circle, expected %d.', sum(inside), n);
end
[~, ~, ~, Z] = ordqz(AA, BB, eye(size(AA)), Z, inside);
U1 = Z(1:n, 1:n);
U2 = Z(n+1:2*n, 1:n);
if rcond(U1) < 1e-14
    error('dare_qz:singularU1', 'deflating subspace is not a graph subspace (rcond(U1) = %.2e).', rcond(U1));
end
X = U2 / U1;
X = (X + X') / 2;
end

function M = psd_project(M)
[Vv, Dd] = eig(M);
M = Vv * diag(max(real(diag(Dd)), 0)) * Vv';
M = (M + M') / 2;
end

function ok = dynare_mex_usable()
%DYNARE_MEX_USABLE  Does Dynare's kalman_steady_state MEX load and solve a scalar DARE? Checked once per session.
%   ADDED 2026-09-28: the MATLAB MEX of the official Dynare 7.2 macOS (arm64) package is
%   linked without SLICOT (dlopen: symbol '_sb02od_' not found in flat namespace), and
%   every failed load attempt inside the try/catch above cost ~50 ms - the BR likelihood
%   went from 7 ms to 59 ms. The Octave (Homebrew) MEX works. The check solves
%   P = T P T + Q - T P Z (Z P Z + H)^{-1} Z P T with T = 0.5, Q = Z = H = 1.
persistent usable
if isempty(usable)
    try
        P = kalman_steady_state(0.5, 1, 1, 1);
        usable = isscalar(P) && isfinite(P) && abs(P - (0.25 * P + 1 - 0.25 * P^2 / (P + 1))) < 1e-8;
    catch
        usable = false;
    end
    if ~usable
        warning('solve_steady_state_riccati:mex', ['Dynare''s kalman_steady_state MEX does not load or gives a wrong result in this ' ...
                'session; the toolbox QZ solver is used for every steady-state Riccati equation.']);
    end
end
ok = usable;
end
