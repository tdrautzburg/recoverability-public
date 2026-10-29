% TEST_INVERTIBLE_RICCATI  The 'invertible' shortcut of solve_steady_state_riccati (ADDED 2026-10-16).
%
%   ABCD timing x(t) = A x(t-1) + B e(t), y(t) = C x(t-1) + D e(t). With D square and invertible,
%   Omega = 0 solves the Riccati equation and the closed loop is F = A - B D^{-1} C. Checks:
%     1. rho(F) < 1 (random systems): method 'invertible', Omega = 0, agrees with the QZ solver.
%     2. rho(F) > 1 (zeros moved outside): the shortcut must NOT fire; the result equals QZ's.
%     3. rho(F) = 1, MA(1) with a unit root  y(t) = e(t) - e(t-1): method 'invertible_critical',
%        Omega = 0, while the old fixed-point iteration is still far from 0 at its cap
%        (convergence at rate 1/t: Omega after k steps of the scalar recursion is 1/(k+1)).
%     4. D not square (more shocks than observables): the shortcut is skipped.
%     5. Outputs K = B D^{-1}, V = D D', Omega_tp1_t = B B' in the shortcut.
%   Self-contained (no Dynare, no fixtures); prints a pass message or errors.

%   Without the shortcut the existing methods cannot ACCEPT Omega = 0 even when rho(F) < 1: QZ
%   returns Omega ~ 0, but the acceptance test measures the residual relative to ||Omega||, so the
%   solution is flagged as not converged. The forced 'qz' / 'iterate' calls below therefore warn;
%   those warnings are expected and switched off here.

here = fileparts(mfilename('fullpath'));
addpath(here);
warn_state = warning('off', 'riccati:noconv');
cleanup_warn = onCleanup(@() warning(warn_state));
if exist('rng', 'file'), rng(20261016); else, rand('state', 20261016); randn('state', 20261016); end %#ok<RAND>

% ---- 1. invertible, rho(F) < 1 ----
n_ok = 0;
for trial = 1:30
    n = 3 + mod(trial, 4);  p = 2 + mod(trial, 2);
    A = randn(n);  A = 0.9 * A / max(abs(eig(A)));
    B = randn(n, p);  D = randn(p) + 2 * eye(p);
    C = 0.3 * randn(p, n);                       % keep the draw only if rho(F) < 1

    F = A - (B / D) * C;
    if max(abs(eig(F))) >= 0.99, continue; end
    [Om, K, V, Otp1, fl] = solve_steady_state_riccati(A, B, C, D);
    assert(strcmp(fl.method, 'invertible'), 'case 1: method %s, expected invertible', fl.method);
    assert(fl.converged && all(Om(:) == 0), 'case 1: Omega must be exactly 0');
    Oq = solve_steady_state_riccati(A, B, C, D, struct('method', 'qz'));
    assert(norm(Oq, 'fro') <= 1e-9 * max(1, norm(B * B', 'fro')), 'case 1: QZ disagrees (||Omega_qz|| = %.2e)', norm(Oq, 'fro'));
    assert(norm(K - B / D, 'fro') <= 1e-12 * max(1, norm(K, 'fro')) && norm(V - D * D', 'fro') <= 1e-12 * norm(V, 'fro') ...
           && norm(Otp1 - B * B', 'fro') <= 1e-12 * norm(Otp1, 'fro'), 'case 5: K, V, Omega_tp1_t at Omega = 0');
    n_ok = n_ok + 1;
end
assert(n_ok >= 5, 'case 1: only %d usable random systems', n_ok);

% ---- 2. non-invertible, rho(F) > 1: scalar MA(1) y = e - 2 e(-1) (zero at 2) ----
A = 0; B = 1; C = -2; D = 1;                      % x(t) = e(t), y(t) = -2 x(t-1) + e(t), F = 2
[Om, ~, ~, ~, fl] = solve_steady_state_riccati(A, B, C, D);
Oq = solve_steady_state_riccati(A, B, C, D, struct('method', 'qz'));
assert(~any(strcmp(fl.method, {'invertible', 'invertible_critical'})), 'case 2: shortcut fired with rho(F) = 2');
assert(fl.converged && abs(Om - Oq) <= 1e-10 && Om > 0.5, 'case 2: expected the nonzero stabilising solution (Omega = %.4f, QZ %.4f)', Om, Oq);

% ---- 3. critical, rho(F) = 1: y = e - e(-1) ----
A = 0; B = 1; C = -1; D = 1;                      % F = 1
[Om, K, ~, ~, fl] = solve_steady_state_riccati(A, B, C, D);
assert(strcmp(fl.method, 'invertible_critical') && fl.converged && Om == 0 && K == 1, 'case 3: method %s, Omega %.3g', fl.method, Om);
% the old iteration from B B' = 1: Omega_k = 1/(k+1), i.e. still 1e-3 after 1000 steps
Oi = solve_steady_state_riccati(A, B, C, D, struct('method', 'iterate', 'MaxIter', 1000));
assert(abs(Oi - 1/1001) < 1e-6, 'case 3: iteration after 1000 steps should be 1/1001, is %.6g', Oi);

% ---- 4. D not square ----
A = 0.5 * eye(2); B = [1 0 0.3; 0 1 0.2]; C = eye(2); D = [0.5 0 0; 0 0.5 0];
[~, ~, ~, ~, fl] = solve_steady_state_riccati(A, B, C, D);
assert(~any(strcmp(fl.method, {'invertible', 'invertible_critical'})), 'case 4: shortcut fired with non-square D');

fprintf('test_invertible_riccati: all checks passed (%d random invertible systems).\n', n_ok);
