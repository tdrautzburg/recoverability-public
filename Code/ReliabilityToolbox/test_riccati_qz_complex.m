% TEST_RICCATI_QZ_COMPLEX  The QZ route of solve_steady_state_riccati with complex closed-loop eigenvalues.
%
%   ADDED 2026-10-20. Until then dare_qz tested exist('ordeig', 'file'), which is 0 in MATLAB
%   (ordeig is a built-in there), and read the eigenvalues off the diagonal of the real Schur form
%   instead. Complex pairs sit in 2x2 blocks, so their moduli were wrong, the count of eigenvalues
%   inside the unit circle was off, and 'qz' was rejected in favour of the slow iteration (the UC3-i2
%   posterior draws). This test forces method 'qz' on random stable systems whose closed loop has
%   complex eigenvalues, with more shocks than observables (the 'invertible' shortcut does not
%   apply), and checks that 'qz' is accepted and agrees with a tightly converged iteration.
%   Self-contained (no Dynare, no fixtures); prints a pass message or errors.

here = fileparts(mfilename('fullpath'));
addpath(here);
if exist('rng', 'file'), rng(20261020); else, rand('state', 20261020); randn('state', 20261020); end %#ok<RAND>
warn_state = warning('off', 'riccati:noconv');
cleanup_warn = onCleanup(@() warning(warn_state));

n_tested = 0;
for trial = 1:20
    n = 6 + mod(trial, 5);  p = 2;  ne = 4;
    % rotation blocks give complex eigenvalues of A (and, generically, of the closed loop)
    A = zeros(n);
    for k = 1:2:n-1
        r = 0.6 + 0.35 * rand;  th = pi * rand;
        A(k:k+1, k:k+1) = r * [cos(th), -sin(th); sin(th), cos(th)];
    end
    if mod(n, 2) == 1, A(n, n) = 0.5; end
    T = randn(n);  A = T * A / T;                       % hide the block structure
    B = randn(n, ne);  C = randn(p, n);  D = randn(p, ne);
    [Om, K, ~, ~, fl] = solve_steady_state_riccati(A, B, C, D, struct('method', 'qz'));
    cl = eig(A - K * C);
    if ~any(abs(imag(cl)) > 1e-8), continue; end        % keep systems with complex closed-loop roots
    assert(strcmp(fl.method, 'qz') && fl.converged, 'trial %d: QZ not accepted (method %s, converged %d)', trial, fl.method, fl.converged);
    Oi = solve_steady_state_riccati(A, B, C, D, struct('method', 'iterate', 'tol_riccati', 1e-14, 'MaxIter', 1e6));
    assert(norm(Om - Oi, 'fro') <= 1e-8 * max(1, norm(Oi, 'fro')), 'trial %d: QZ and iteration disagree (%.2e)', trial, norm(Om - Oi, 'fro'));
    n_tested = n_tested + 1;
end
assert(n_tested >= 5, 'only %d systems with complex closed-loop eigenvalues', n_tested);
fprintf('test_riccati_qz_complex: all checks passed (%d systems with complex closed-loop eigenvalues).\n', n_tested);
