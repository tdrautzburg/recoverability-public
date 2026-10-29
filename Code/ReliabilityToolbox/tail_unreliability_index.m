function [nu, L, Pi, info, d] = tail_unreliability_index(A, B, Omega, state_indices, verbose)
%TAIL_UNRELIABILITY_INDEX  Horizon-free reliability index nu for non-stationary
%   states of ARBITRARY integration order.
%
%   ADDED (2026-08). REVISED (2026-09) to the general-order formula.
%   Implements the tail-unreliability index of the note "A horizon-free
%   reliability index for non-stationary states of arbitrary integration
%   order" (Rate_of_convergence_R2.pdf, Definition 2.1 and Theorem 2.2), whose
%   eq. (5) supersedes the I(1)-only limit this function used to compute. It is
%   the non-stationary counterpart of the paper's unconditional R2 and replaces
%   the long-horizon conditional R2 columns of the headline tables.
%
%   WHY IT EXISTS.  For a unit-root state the paper's conditional reliability
%
%       R2_{h,i} = 1 - Omega_ii / D_i(h),                                (note eq. (1))
%       D_i(h)   = [A^h Omega (A^h)']_ii
%                  + sum_{j=1..h} [A^{j-1} B B' (A^{j-1})']_ii,          (note eq. (2))
%
%   converges to 1 for EVERY model, state and parameterization (paper's
%   Theorem 3.4): the numerator does not depend on h, while the denominator
%   diverges. A statistic whose limit is the same for all objects being
%   compared cannot discriminate between them, so all of the model- and
%   state-specific content of {R2_{h,i}} lies in the RATE at which it
%   approaches 1. Theorem 2.2 of the note extracts that rate:
%
%       h^(2d+1) * (1 - R2_{h,i})  -->  nu_i^(d) = (2d+1) * Omega_ii / L_i^(d)
%
%   as h --> Inf, with d = d_i the ORDER of state i (state i is I(d_i + 1)) and
%   L_i^(d) its generalised long-run innovation variance of order d. This is
%   eq. (5) of the note, and it is the unique non-degenerate limit obtainable
%   from {R2_{h,i}} by rescaling in the horizon: lim_h h^alpha (1 - R2_{h,i})
%   is 0 for alpha < 2d+1, is nu_i^(d) for alpha = 2d+1, and diverges for
%   alpha > 2d+1. Larger nu_i means LESS reliable.
%
%   WHAT CHANGED (2026-09).  The previous implementation computed the d = 0
%   specialisation only,
%
%       h * (1 - R2_{h,i})  -->  nu_i = Omega_ii / L_i,      L_i = [Pi B B' Pi']_ii,
%
%   and returned NaN for every state of higher order, because for those
%   h*(1 - R2_{h,i}) -> 0 and there is no index at that rescaling. HLW's ystar
%   is exactly such a state and used to come back empty. The general formula
%   below reduces to the old one when d_i = 0, so every number previously
%   reported is unchanged; what is new is that I(2) and higher states now get
%   an index instead of a NaN.
%
%   THE INGREDIENTS.  Let Pi be the Riesz (spectral) projector onto the
%   GENERALISED eigenspace of the eigenvalue 1 and
%
%       N = (A - I) * Pi,        A = Pi + N + A_s,   N nilpotent, rho(A_s) < 1
%
%   (note eq. (3)). N records how the unit roots are CHAINED: N = 0 says the
%   unit eigenvalue is semisimple, so every non-stationary state is I(1);
%   N ~= 0 says some state carries a unit-root state as a drift and is
%   therefore of higher order. Then (note, Definition 2.1, eq. (4)):
%
%     d_i      : the ORDER of state i,
%                    d_i = max{ k >= 0 : e_i' N^k B ~= 0 },     N^0 = Pi,
%                i.e. the highest rung of state i's cumulation ladder that
%                still carries shocks. State i is I(d_i + 1). d_i = 0 is the
%                I(1) case in which N plays no role for state i.
%     L_i^(d)  : the generalised long-run innovation variance of order d,
%                    L_i^(d) = || e_i' N^d B ||^2 / (d!)^2
%                            = [N^d B B' (N^d)']_ii / (d!)^2.
%                It is the innovation variance of the TOP rung of state i's
%                ladder; equivalently the sum of squared normalised long-run
%                impulse responses e_i' N^d B / d! = lim_k k^{-d} dx_{i,t+k}/de_t';
%                equivalently, up to (d!)^2, the long-run variance of the
%                (d_i+1)-differenced state. At d = 0 it is the old
%                L_i = [Pi B B' Pi']_ii = ||e_i' Pi B||^2, the per-period
%                variance of the permanent (random-walk) component of state i.
%     Omega_ii : the steady-state filtered error variance Var(x_{i,t} | y^t),
%                i.e. the diagonal of the object the engines call KF_Omega.
%                It is the numerator of 1 - R2_{h,i} at every horizon.
%
%   The single reason the order matters is visible in the growth rate of the
%   denominator (note eq. (6)),
%
%       D_i(h) = L_i^(d) / (2d+1) * h^(2d+1) + O(h^(2d)),
%
%   cumulation along the ladder makes the shock-accumulation term of eq. (2)
%   grow like h^(2d+1) rather than like h, while Omega_ii is horizon-invariant
%   regardless.
%
%   NOT [B*B']_ii.  The natural but incorrect normalization is the state's own
%   one-step innovation variance: the two coincide only when state i is a pure
%   random walk driven by its own shock. In general (i) part of a shock's
%   impact on state i is transitory and does not enter L_i^(d), and (ii) shocks
%   to OTHER states can have permanent effects on state i through the
%   off-diagonal dynamics and then enter L_i^(d) without appearing in [B*B']_ii
%   at all. This is exactly the case for a composite trend such as HLW's
%   r* = c*g + z, where L_{r*}^(0) = c^2 sigma_g^2 + sigma_z^2, and even more so
%   for HLW's ystar, whose top rung is carried entirely by the shock to g.
%
%   UNITS, AND WHY nu IS NOT COMPARABLE ACROSS ORDERS.  nu_i^(d) is measured in
%   periods^(2d+1) (quarters^(2d+1) in the applications), so indices of states
%   with different d must NOT be put side by side or pooled in a single table
%   column: an r* with nu^(0) = 40 and a ystar with nu^(1) = 2449 are not
%   "sixty times worse" - they are not in the same units (note, section 2.3).
%   The order is returned in `d` (and in info.d) precisely so that the caller
%   can label, group or convert them; a statistic measured in periods for every
%   d is constructed in Appendix A.13 of the note and is not computed here.
%
%   READING nu.  At d = 0 the decomposition D_i(h) = Omega_ii + h*L_i is a
%   fixed STOCK of standing uncertainty plus a FLOW of genuine movement
%   accumulating at a constant rate, so nu_i = Omega_ii / L_i is a stock
%   divided by a flow, hence a horizon; when row i of A equals e_i' (state i's
%   own equation is a driftless random walk) the whole reliability profile
%   collapses to the exact one-parameter hyperbola R2_{h,i} = h / (h + nu_i)
%   and nu_i is precisely the half-reliability horizon.
%
%   THAT READING DOES NOT SURVIVE d >= 1 (note, section 2.3). With d_i >= 1 the
%   state's variance accumulates at a rate that is itself growing in the
%   horizon, so there is no constant flow to divide by; what replaces it is
%
%       D_i(h) = Omega_ii  +  (D_i(h) - Omega_ii),
%                \_______/     \_________________/
%                 standing      accumulated movement, growing like h^(2d+1)
%                uncertainty
%
%   with both terms nonnegative and the second no longer proportional to h. The
%   correct verbal statement is that forecasting an I(d+1) trend far out is hard
%   not because one becomes worse at locating it, but because it has more time
%   to move AND because the rate at which it moves is itself accumulating.
%
%   DO NOT estimate nu_i by evaluating h^(2d+1)*(1 - R2_{h,i}) at a moderate
%   horizon: the convergence is slow (the statistic is ~33% below its limit at
%   h = 20 and ~9% below at h = 100 in the note's Design I, and slower still at
%   higher orders, where the neglected term is O(h^(2d)) against h^(2d+1)).
%   nu_i is available exactly and in closed form from Omega, N and Pi, which is
%   what this function computes. In a Bayesian implementation nu_i is a smooth
%   function of the structural parameters, so its posterior follows by
%   evaluating this function at each draw.
%
%   STATE-SPACE CONVENTION (identical to the rest of the toolbox, ABCD form
%   with observables loading on the LAGGED state):
%
%       x(t) = A x(t-1) + B e(t)
%       y(t) = C x(t-1) + D e(t),     e(t) ~ N(0, I).
%
%   INPUTS
%     A              n_states x n_states transition matrix (PRUNED state
%                    vector, i.e. ModelMatrices.A)
%     B              n_states x n_shocks shock loading matrix
%     Omega          n_states x n_states steady-state FILTERED error variance
%                    Var(x_t | y^t) (i.e. KF_Omega, not the predicted variance
%                    KF_Omega_tp1_t)
%     state_indices  (optional) indices of the states for which nu is wanted,
%                    in the pruned state vector; default: all states
%     verbose        (optional, default false) print/raise warnings when
%                    Assumption 1 of the note fails. Diagnostics are recorded
%                    in `info` regardless, so per-draw loops can stay silent.
%
%   OUTPUTS
%     nu     numel(state_indices) x 1 tail-unreliability indices, in
%            periods^(2*d+1) STATE BY STATE (see UNITS above), in the order of
%            state_indices. NaN where the index is undefined, see "WHEN nu IS
%            NaN" below.
%     L      n_states x 1 generalised long-run innovation variances
%            L_i^(d_i) = ||e_i' N^{d_i} B||^2 / (d_i!)^2, for ALL states, each
%            at ITS OWN order d_i (0 where the state carries no trend).
%            CHANGED 2026-09: this used to be [Pi B B' Pi']_ii for every state,
%            i.e. the d = 0 rung only; the two agree wherever d_i = 0.
%     Pi     n_states x n_states Riesz projector on the unit-root GENERALIZED
%            eigenspace. It equals the ordinary spectral projector whenever the
%            unit eigenvalue is semisimple, and still exists when it is not
%     info   struct of numerical diagnostics:
%              .n_unit_roots       number of eigenvalues of A at 1 (algebraic
%                                  multiplicity)
%              .rho_stable         max modulus of the remaining eigenvalues,
%                                  i.e. rho(A_s); a by-product of the spectrum
%                                  check, kept because it governs how far out
%                                  the asymptotic regime starts
%              .d                  n_states x 1 order d_i of EVERY state, NaN
%                                  where the state carries no stochastic trend.
%                                  d = info.d(state_indices) (ADDED 2026-09)
%              .rung_norms         n_states x (n_unit_roots) matrix of the row
%                                  norms ||e_i' N^k Pi B||, k = 0,1,2,...; the
%                                  ladder from which d_i is read off. Column
%                                  k+1 holds rung k (ADDED 2026-09)
%              .max_order          max(info.d), the highest order present; 0
%                                  means every non-stationary state is I(1)
%                                  (ADDED 2026-09)
%              .assumption1_ok     true if A has semisimple unit eigenvalues and
%                                  all other eigenvalues strictly inside the
%                                  unit circle. FALSE IS NOT A FAILURE: it means
%                                  some state has d_i >= 1, which eq. (5)
%                                  handles. Assumption 1 of the note itself only
%                                  asks that every eigenvalue be 1 or inside the
%                                  unit circle, which is the .spectrum_ok test
%              .spectrum_ok        true if every eigenvalue of A is either 1 or
%                                  strictly inside the unit circle, i.e.
%                                  Assumption 1 of the note. FALSE means nu is
%                                  NaN throughout (ADDED 2026-09)
%              .assumption2_ok     DEPRECATED alias of .assumption1_ok, kept so
%                                  that existing callers do not break
%              .unit_root_defective true if N = (A - I)*Pi is nonzero, i.e.
%                                  there is a Jordan block at 1
%              .jordan_chain_states indices of the states with e_i'N ~= 0. These
%                                  are the states of order d_i >= 1. CHANGED
%                                  2026-09: they no longer get nu = NaN
%              .projector_source   'schur', 'eig', 'power', 'trivial' (no unit
%                                  root, Pi = 0) or 'none' (failure)
%              .projector_residual max residual of Pi*Pi = Pi, A*Pi = Pi*A and
%                                  trace(Pi) = .n_unit_roots
%              .power_check_gap    ||Pi - A^J|| for large J (cross-check against
%                                  Pi = lim A^j), NaN when that limit does not
%                                  exist, which is the case precisely when the
%                                  unit eigenvalue is defective
%              .messages           cell array of diagnostic messages
%     d      numel(state_indices) x 1 order d_i of each requested state, in the
%            order of state_indices, NaN wherever nu is NaN. ADDED 2026-09:
%            nu(kk) is in periods^(2*d(kk)+1), so this output is what makes the
%            returned indices interpretable. Appended LAST so that every
%            existing [nu, L, Pi, info] call site keeps working unchanged.
%
%   WHEN nu IS NaN
%     - L_i^(d) = 0 at every rung, i.e. e_i' N^k B = 0 for all k >= 0. The
%       direction carries no stochastic trend at any order: either state i is
%       stationary (the HLW/UC3 output gap, for which D_i(h) converges to a
%       finite Sigma_ii and h*(1 - R2_{h,i}) diverges) or e_i is a
%       cointegrating direction. There the mechanical-convergence critique does
%       not arise and the paper's UNCONDITIONAL R2_i is already the right
%       horizon-free summary, so nu must not be reported and must not be pooled
%       in a single column with the trends' nu.
%     - The spectrum itself is unusable (an explosive eigenvalue, or an
%       eigenvalue of unit modulus away from +1, e.g. -1 or a complex pair on
%       the circle), i.e. Assumption 1 of the note fails. Then D_i(h) has no
%       power-law growth and the limit does not exist for ANY coordinate, so nu
%       is NaN throughout.
%
%   NO LONGER A NaN CASE: DEFECTIVE UNIT ROOTS (I(2) AND HIGHER STATES).
%   Assumption 1 of the note deliberately does NOT require the unit eigenvalue
%   to be semisimple, and estimated models routinely make it defective. HLW is
%   the case in point: ystar = ystar(-1) + gg(-1)/4 + shock with gg a random
%   walk gives the trend block [1 1/4; 0 1], a Jordan block at 1, because ystar
%   is I(2). This function works with the RIESZ projector on the generalized
%   eigenspace, which exists regardless, and then reads each state's order off
%   the ladder N^k Pi B:
%     - e_i'N = 0  : e_i'A^j -> e_i'Pi geometrically, d_i = 0, and nu_i is the
%                    familiar Omega_ii / L_i;
%     - e_i'N ~= 0 : state i accumulates a random walk. D_i(h) grows like
%                    h^(2d+1) with d = d_i >= 1, and nu_i = (2d+1)Omega_ii/L_i^(d)
%                    is its index - in periods^(2d+1), NOT in periods.
%   In HLW this gives gg, zz and rrstar = c*gg + zz an index at d = 0 (the
%   note's worked example, L_{r*}^(0) = c^2 sigma_g^2 + sigma_z^2) and ystar an
%   index at d = 1 with L_{ystar}^(1) = sigma_g^2/16, where before 2026-09 it
%   had none. When the unit eigenvalue is semisimple N = 0, every d_i is 0, Pi
%   is the ordinary spectral projector and nothing changes.
%
%   LINEAR FUNCTIONALS.  For a composite object a'x_t (e.g. HLW's r* = c*g + z
%   with a = (c,1)' on the (g,z) block) the note's Definition 2.1 applies to the
%   direction a verbatim: d(a) = max{k >= 0 : a'N^k B ~= 0} and
%   nu(a) = (2d+1) * a'*Omega*a / (||a' N^d B||^2 / (d!)^2), provided the
%   denominator is positive. Given the Pi returned here that is a few lines; the
%   coordinate case is a = e_i. In this toolbox composite states of interest are
%   appended to the state vector in compute_state_reliability_diagnostics
%   (section 2), so they already appear as coordinates and the coordinate
%   formula covers them.
%
%   SCALE INVARIANCE.  nu is invariant to a diagonal rescaling of the states
%   (percent vs. decimals vs. log points): under x -> S x we have A -> S A S^-1,
%   hence Pi -> S Pi S^-1 and N -> S N S^-1, so N^d B -> S N^d B and both
%   Omega_ii and L_i^(d) scale by s_i^2 while d_i is unchanged. This is what
%   makes nu comparable across states whose raw uncertainties are not - within
%   a common order d, see UNITS above.
%
%   See also: compute_state_reliability_diagnostics,
%             steady_state_kalman_filter.

%% ------------------------------------------------------------------
%  0. Input handling, tolerances and default (failure) outputs
% -------------------------------------------------------------------
n = size(A, 1);

if nargin < 4 || isempty(state_indices)
    state_indices = (1:n)';
end
state_indices = state_indices(:);

if nargin < 5 || isempty(verbose)
    verbose = false;
end

if size(A, 2) ~= n || size(B, 1) ~= n || any(size(Omega) ~= n)
    error('tail_unreliability_index:badInput', ...
        'A must be square and B, Omega must have %d rows/columns.', n);
end

% Tolerances -------------------------------------------------------------
% UNIT_TOL: an eigenvalue counts as a unit root when |lambda - 1| <= UNIT_TOL.
%   It is set to the same 1e-6 that the .mod files use for Dynare's
%   oo_.dr.qz_crit = 1 + 1e-6, so that a root treated as a unit root when the
%   diffuse filter is initialized is also treated as one here.
% RESID_TOL: relative tolerance for the projector identities Pi^2 = Pi,
%   A*Pi = Pi*A and trace(Pi) = (number of unit roots), which is how we verify
%   that the computed Pi really is the Riesz projector. Note that A*Pi = Pi is
%   NOT among them: it holds only in the semisimple case, and is exactly what
%   section 2b tests for separately.
% NILPOTENT_TOL: relative threshold above which a rung of the ladder
%   N^k Pi B counts as NONZERO, i.e. above which coordinate i is taken to have
%   order d_i >= k. Rung k of state i is compared against the LARGEST rung of
%   THAT SAME state (see section 2b), which is the scale-free comparison: under
%   x -> S x every rung of state i is multiplied by the same s_i, so the ratio
%   is invariant. The threshold MUST be looser than UNIT_TOL, and is tied to it
%   for that reason. An eigenvalue accepted as a unit root is only 1 to within
%   UNIT_TOL, so even a perfectly semisimple root deposits up to ~UNIT_TOL of
%   spurious mass in N and hence in rung 1 relative to rung 0; a threshold at or
%   below UNIT_TOL would read that round-off as a Jordan chain and wrongly
%   promote the state to d_i = 1 (a root at 1 + 1e-9 did exactly that at
%   1e-10), which would then report an index in the wrong units. Genuine chains
%   sit at O(1) relative - HLW's ystar has rung1/rung0 = (sigma_g/4)/sigma_ystar
%   ~ 2e-2 - so the factor of 100 leaves about two decades of margin on each
%   side.
% ZERO_L_TOL: relative threshold below which L_i^(d) counts as exactly zero,
%   i.e. below which the state carries no stochastic trend at ANY rung. The
%   comparison scale is [B*B']_ii + Omega_ii, which is in state i's own units
%   squared, exactly like L_i^(d) (N^d B transforms like B under x -> S x), so
%   the test is invariant to a rescaling of the states (no cross-state
%   comparison of magnitudes is involved).
UNIT_TOL      = 1e-6;
RESID_TOL     = 1e-8;
NILPOTENT_TOL = 100 * UNIT_TOL;
ZERO_L_TOL    = 1e-10;

% Outputs default to "undefined": every failure path below simply returns
% these, so a bad draw inside prior_function/posterior_function yields NaN
% rather than an error.
nu   = nan(numel(state_indices), 1);
d    = nan(numel(state_indices), 1);
L    = nan(n, 1);
Pi   = nan(n, n);
info = struct('n_unit_roots', 0, 'rho_stable', NaN, 'd', nan(n, 1), ...
              'rung_norms', nan(n, 0), 'max_order', NaN, ...
              'spectrum_ok', false, 'assumption1_ok', false, ...
              'assumption2_ok', false, ...
              'unit_root_defective', false, 'jordan_chain_states', [], ...
              'projector_source', 'none', 'projector_residual', NaN, ...
              'power_check_gap', NaN, 'messages', {{}});

%% ------------------------------------------------------------------
%  1. Check Assumption 1 of the note on the eigenvalues of A
% -------------------------------------------------------------------
% Assumption 1 (Rate_of_convergence_R2.pdf): every eigenvalue of A is either
% equal to 1 or lies strictly inside the unit circle. Note that it does NOT
% require the unit eigenvalue to be semisimple - Jordan blocks of any size are
% allowed, and are exactly what the order d_i of section 2b measures. Here we
% only classify the spectrum; semisimplicity is recorded separately (as
% .assumption1_ok, the old .assumption2_ok) because it is the condition under
% which every d_i is 0 and the pre-2026-09 formula applies unchanged.
[Vr, Dg, Wl] = eig(A);          % Vr: right eigenvectors (A*Vr = Vr*Dg)
lambda       = diag(Dg);        % Wl: left  eigenvectors (Wl'*A = Dg*Wl')
is_unit      = abs(lambda - 1) <= UNIT_TOL;

info.n_unit_roots = sum(is_unit);
if any(~is_unit)
    info.rho_stable = max(abs(lambda(~is_unit)));
else
    info.rho_stable = 0;
end

% Any eigenvalue that is not a unit root but sits on or outside the unit
% circle (an explosive root, or a root of unit modulus away from +1, e.g. -1
% or a complex pair) violates Assumption 1: A^j then fails to decompose into a
% projector plus a geometrically vanishing part, D_i(h) has no power-law growth,
% and h^alpha*(1 - R2_{h,i}) has no non-degenerate limit at any alpha. Return
% NaN rather than a number that would not mean anything.
if info.rho_stable >= 1 - UNIT_TOL
    info.messages{end+1} = sprintf(['Assumption 1 fails: A has a non-unit ' ...
        'eigenvalue of modulus %.6g (>= 1). nu is not defined.'], info.rho_stable);
    if verbose
        warning('tail_unreliability_index:assumption1', '%s', info.messages{end});
    end
    return
end
info.spectrum_ok = true;

% A stationary model has no unit root at all: nothing to index. This is not an
% error (the paper's unconditional R2 is the right measure there), so return
% quietly with Pi = 0 and nu = NaN.
if info.n_unit_roots == 0
    Pi   = zeros(n, n);
    L    = zeros(n, 1);
    info.assumption1_ok    = true;
    info.assumption2_ok    = true;   % deprecated alias, kept for old callers
    info.max_order         = 0;
    info.rung_norms        = zeros(n, 0);
    info.projector_source  = 'trivial';   % Pi = 0, no route needed
    info.projector_residual = 0;
    info.messages{end+1} = ['A has no unit root: every state is stationary, ' ...
        'nu is undefined (use the unconditional R2).'];
    return
end

%% ------------------------------------------------------------------
%  2. Riesz projector Pi on the unit-root GENERALIZED eigenspace
% -------------------------------------------------------------------
% When the unit eigenvalue is SEMISIMPLE the object wanted here is the ordinary
% spectral projector and A^j = Pi + A_s^j. Some estimated models make it
% defective instead; HLW is the leading example,
%
%       ystar = ystar(-1) + gg(-1)/4 + shock,      gg = gg(-1) + shock,
%
% whose trend block is [1 1/4; 0 1] - a JORDAN BLOCK at 1, because ystar
% accumulates a random walk and is therefore I(2). There the unit eigenvalue
% has algebraic multiplicity above its geometric multiplicity, lim A^j does not
% exist (the block grows linearly in j), and the eigenvector formula
% Pi = R(L'R)^{-1}L' breaks down because L'R is singular.
%
% Assumption 1 of the note admits that case on purpose, so what we compute is
% the RIESZ projector onto the GENERALIZED eigenspace of lambda = 1, which
% exists whether or not the eigenvalue is semisimple and coincides with the
% spectral projector when it is. Writing N = (A - I)*Pi, which is nilpotent and
% satisfies Pi*N = N*Pi = N,
%
%       A^j Pi = Pi + j*N + C(j,2)*N^2 + ... ,      A^j = A^j Pi + A_s^j,
%
% so for coordinate i:
%   - if e_i'N = 0 then e_i'A^j = e_i'Pi + e_i'A_s^j -> e_i'Pi geometrically,
%     D_i(h) is asymptotically affine, the order is d_i = 0 and
%     nu_i = Omega_ii / ||e_i' Pi B||^2 is the familiar I(1) index;
%   - if e_i'N ~= 0 the state accumulates a random walk and is I(d_i + 1) with
%     d_i >= 1: the binomial coefficient C(j,d) above makes the shock
%     accumulation in D_i(h) grow like h^(2d+1), and eq. (5) of the note divides
%     that growth out. This is the case the pre-2026-09 code returned NaN for.
% Section 2b builds the ladder N^k Pi B from which d_i is read off row by row.
k = info.n_unit_roots;   % algebraic multiplicity of the unit eigenvalue
projector_ok = false;

% --- Route 1 (primary): ordered real Schur + one Sylvester solve -----------
% Reorder the real Schur form so that the lambda = 1 cluster comes first,
%   U'AU = T = [T11 T12; 0 T22],   spec(T11) = {1},  spec(T22) strictly inside,
% and solve the Sylvester equation T11*Y - Y*T22 = T12. Then
%   Pi = U * [I Y; 0 0] * U'
% is the Riesz projector: [I Y; 0 0] is idempotent for ANY Y, and the Sylvester
% solution is exactly what makes it commute with T. The equation has a unique
% solution because spec(T11) and spec(T22) are disjoint, which section 1 has
% already guaranteed. This route never forms an eigenvector basis, so a
% defective eigenvalue costs it nothing.
P_schur = [];
try
    [U, Tsch] = schur(A, 'real');
    sel = abs(ordeig(Tsch) - 1) <= UNIT_TOL;
    if sum(sel) == n
        % Every eigenvalue is a unit root: the generalized eigenspace is all of
        % R^n and the projector is the identity (no Sylvester solve needed, and
        % T22 would be empty).
        P_schur = eye(n);
    else
        [U, Tsch] = ordschur(U, Tsch, sel);
        kk_       = sum(sel);
        T11 = Tsch(1:kk_, 1:kk_);
        T12 = Tsch(1:kk_, kk_+1:n);
        T22 = Tsch(kk_+1:n, kk_+1:n);
        if exist('sylvester', 'builtin') || exist('sylvester', 'file')
            Y = sylvester(T11, -T22, T12);     % solves T11*Y - Y*T22 = T12
        else
            % Portable fallback for MATLAB releases without sylvester (pre
            % R2014a): the same equation in vec form,
            %   (I (x) T11 - T22' (x) I) vec(Y) = vec(T12).
            Y = reshape((kron(eye(n-kk_), T11) - kron(T22', eye(kk_))) \ T12(:), ...
                        kk_, n-kk_);
        end
        P_schur = U * [eye(kk_), Y; zeros(n-kk_, kk_), zeros(n-kk_, n-kk_)] * U';
    end
catch
    P_schur = [];   % schur/ordschur/sylvester unavailable or failed: try route 2
end

if ~isempty(P_schur) && all(isfinite(P_schur(:)))
    res = riesz_projector_residual(P_schur, A, k);
    if res <= RESID_TOL
        Pi                      = P_schur;
        info.projector_source   = 'schur';
        info.projector_residual = res;
        projector_ok            = true;
    end
end

% --- Route 2 (fallback): the eigenvector route of the note -----------------
% Pi = R (L'R)^{-1} L', with the columns of R spanning the RIGHT unit
% eigenspace (A*R = R) and those of L the LEFT one (L'*A = L'). This is the
% route named in the note ("Computing Pi") and it is correct whenever the unit
% eigenvalue is semisimple, which covers UC3 (three independent random-walk
% trends) but not HLW. Kept as a fallback and as an independent check.
%
% Pi does not depend on how the bases R and L are chosen inside the eigenspace:
% replacing R by R*S1 and L by L*S2 leaves R(L'R)^{-1}L' unchanged, so the
% arbitrary scaling (and, for a repeated root, the arbitrary basis) that eig
% returns is harmless. Note also that Pi is in general NOT symmetric, because A
% is not normal: it is essential to keep the left and right eigenspaces apart
% and to use the LEFT action e_i'*Pi when reading off state i (note, Remark 4 -
% Pi*e_i = e_i is neither necessary nor sufficient for anything).
if ~projector_ok
    R  = Vr(:, is_unit);
    Lv = Wl(:, is_unit);
    M  = Lv' * R;
    % The conditioning of L'R is the semisimplicity guard: for a defective unit
    % eigenvalue the left and right eigenspaces are orthogonal and L'R is
    % singular. cond is used rather than rcond because M is complex-typed
    % whenever A has any complex eigenvalues (it does here: the ygap block of
    % UC3 has a complex pair) and MATLAB's rcond rejects complex input.
    if 1/cond(M) > eps^(2/3)
        Pi_eig = R * (M \ Lv');
        if norm(imag(Pi_eig), 'fro') <= RESID_TOL * max(1, norm(real(Pi_eig), 'fro'))
            Pi_eig = real(Pi_eig);   % unit eigenvalues are real, so Pi is real
            res    = riesz_projector_residual(Pi_eig, A, k);
            if res <= RESID_TOL
                Pi                      = Pi_eig;
                info.projector_source   = 'eig';
                info.projector_residual = res;
                projector_ok            = true;
            end
        end
    end
end

% --- Route 3 (cross-check, and last-resort fallback) -----------------------
% Pi = lim_{j->Inf} A^j, immediate from A^j = Pi + A_s^j in the SEMISIMPLE case
% only - with a Jordan block the limit does not exist, which is precisely why
% this route cannot rescue HLW and why route 1 is now primary. Computed by
% repeated squaring, A^(2^k), reaching J = 2^16 = 65536 in 16 products. The
% squaring count is CAPPED rather than run to machine convergence on purpose: a
% "unit" root that is numerically 1 + 1e-7 would be raised to the power J and
% drift away from 1, so we stop as soon as successive squarings agree and treat
% non-convergence as "unavailable" instead of trusting the result.
Apow = A;
gap  = Inf;
for kk_ = 1:16
    Anext = Apow * Apow;
    gap   = norm(Anext - Apow, 'fro') / max(1, norm(Anext, 'fro'));
    Apow  = Anext;
    if gap <= RESID_TOL || ~all(isfinite(Apow(:)))
        break
    end
end
power_converged = gap <= RESID_TOL && all(isfinite(Apow(:)));

if projector_ok
    % Record the agreement between the routes as a numerical health check. It
    % is only available when A^j converges, i.e. in the semisimple case.
    if power_converged
        info.power_check_gap = norm(Pi - Apow, 'fro') / max(1, norm(Pi, 'fro'));
    end
elseif power_converged
    Pi                      = Apow;
    info.projector_source   = 'power';
    info.projector_residual = riesz_projector_residual(Pi, A, k);
    info.power_check_gap    = 0;
    projector_ok            = info.projector_residual <= sqrt(RESID_TOL);
    info.messages{end+1} = ['Schur and eigenvector routes for the unit-root ' ...
        'projector failed their residual tests; used Pi = lim A^j instead.'];
    if verbose
        warning('tail_unreliability_index:projectorFallback', '%s', info.messages{end});
    end
end

if ~projector_ok
    % No route produced a genuine projector. Unlike before, a defective unit
    % eigenvalue is NOT a reason to land here (route 1 handles it), so this now
    % means a real numerical breakdown.
    Pi = nan(n, n);
    info.messages{end+1} = ['Could not compute the unit-root projector. ' ...
        'nu is not defined.'];
    if verbose
        warning('tail_unreliability_index:noProjector', '%s', info.messages{end});
    end
    return
end
%% ------------------------------------------------------------------
%  2b. The cumulation ladder N^k Pi B, and the order d_i of each state
% -------------------------------------------------------------------
% N = 0 exactly when the unit eigenvalue is semisimple (A*Pi = Pi). Otherwise N
% carries the Jordan chains: N^k x_t = N^k x_{t-1} + N^{k+1} x_{t-1} + N^k B e_t
% is a ladder of cumulations whose top rung is an exact driftless random walk
% and each of whose lower rungs cumulates the one above it. Which rung state i
% sits on is its order d_i, and that is what this section computes.
N = (A - eye(n)) * Pi;

% --- The ladder -----------------------------------------------------------
% Rung k of state i is e_i' N^k Pi B, with N^0 = Pi (note, Definition 2.1). The
% rungs are built by left-multiplying the previous one by N rather than by
% forming N^k first: N*Pi = Pi*N = N, so N*(N^{k-1} Pi B) = N^k B is exactly the
% object of eq. (4), and each step costs one n x n_shocks product.
%
% N is nilpotent on range(Pi), whose dimension is k = the algebraic multiplicity
% of the unit eigenvalue, so N^k = 0 EXACTLY and there are at most k rungs
% (k = 1, i.e. a single unit root, leaves rung 0 alone). The loop is therefore
% bounded a priori and needs no convergence test.
%
% ROW NORMS, NOT THE n x n PRODUCTS. For any M, ||e_i' M||^2 = [M M']_ii,
% because the right factor of the quadratic form, M' e_i, is the transpose of
% the left one, e_i' M, so the form collapses to a squared Euclidean length;
% equivalently, sum(M.^2, 2) is exactly diag(M*M'). (Checked against the literal
% diag(Pi*B*B'*Pi') to 8e-16 relative over random non-normal systems.) Two
% reasons to prefer it: a sum of squares cannot come out negative, and L_i^(d)
% is both a denominator and the quantity tested against zero below, so a
% slightly negative value would surface as a large negative nu instead of the
% NaN a trendless direction is supposed to get; and it costs n x n_shocks work
% instead of building an n x n matrix and then discarding every off-diagonal
% entry. The shortcut is valid for the DIAGONAL only - nu never uses the rest.
%
% Rung 0, Pi*B, is the matrix of long-run cumulative impulse responses
% lim_j dx_{t+j}/de_t'; rung d, once divided by d! in section 3, is the
% normalised long-run response lim_j j^{-d} dx_{i,t+j}/de_t' of the note's
% eq. (4). Squaring a rung BEFORE the row sum decomposes L_i^(d) across shocks
% at no extra cost: for HLW's ystar the d = 1 rung is carried entirely by the
% shock to gg, which is Remark 6's point that L is not [B*B']_ii.
rung_norms      = zeros(n, k);
M_rung          = Pi * B;                   % rung 0: N^0 = Pi
rung_norms(:,1) = sqrt(sum(M_rung.^2, 2));
for rr_ = 2:k
    M_rung             = N * M_rung;        % rung rr_-1: N^(rr_-1) * B
    rung_norms(:,rr_)  = sqrt(sum(M_rung.^2, 2));
end
info.rung_norms = rung_norms;

% --- The order d_i, read off the ladder -----------------------------------
% d_i = max{ k >= 0 : e_i' N^k B ~= 0 } (note, eq. (4)). "Nonzero" is judged
% against the LARGEST rung of the SAME state, which is the scale-free test:
% under x -> S x every rung of state i is multiplied by the same s_i, so the
% ratio is invariant, and no cross-state comparison of magnitudes is involved
% (see NILPOTENT_TOL in section 0). A row that is zero at every rung carries no
% stochastic trend in any direction and gets d_i = NaN; section 3 turns that
% into nu = NaN.
row_max = max(rung_norms, [], 2);
d_all   = nan(n, 1);
for ii = 1:n
    active = find(rung_norms(ii, :) > NILPOTENT_TOL * row_max(ii));
    if ~isempty(active)
        d_all(ii) = active(end) - 1;        % column j holds rung j-1
    end
end
info.d         = d_all;
info.max_order = max([d_all(~isnan(d_all)); 0]);

% --- Diagnostics: is the unit eigenvalue defective? -----------------------
% Row i of N tells us whether coordinate i sits on a Jordan chain at all. This
% is now a STRUCTURAL diagnostic only - it no longer decides anything, because
% eq. (5) covers the chained states too. It is not the same test as d_i > 0:
% e_i'N ~= 0 with e_i'N B = 0 is a chain that carries no shocks, and such a
% state is still handled at d_i = 0.
%
% Scale for the row test. N and Pi transform identically under a diagonal
% rescaling of the states (both become S*X*S^{-1}), so comparing row norms of
% one against the other is the natural scale-free test; the max(.,1) floor
% covers rows of Pi that are themselves zero. Note that a zero row of Pi forces
% a zero row of N - e_i'Pi = 0 gives e_i'N = e_i'A*Pi - e_i'Pi = e_i'Pi*A = 0
% using A*Pi = Pi*A - so a stationary coordinate can never be misread as
% sitting on a Jordan chain.
N_row = sqrt(sum(N.^2,  2));
P_row = sqrt(sum(Pi.^2, 2));

on_chain = N_row > NILPOTENT_TOL * max(P_row, 1);

info.unit_root_defective = any(on_chain);
info.jordan_chain_states = find(on_chain)';
info.assumption1_ok      = ~info.unit_root_defective;
info.assumption2_ok      = info.assumption1_ok;  % deprecated alias, old callers

if info.unit_root_defective
    info.messages{end+1} = sprintf(['Unit eigenvalue is defective (Jordan ' ...
        'block at 1): state(s) [%s] sit on a chain and are I(2) or higher. ' ...
        'They now receive the general-order index of eq. (5) rather than NaN; ' ...
        'their nu is in periods^(2d+1) and must not be pooled in one column ' ...
        'with the periods of the d = 0 states. Orders: [%s].'], ...
        num2str(info.jordan_chain_states), num2str(d_all(on_chain)'));
    if verbose
        warning('tail_unreliability_index:defectiveUnitRoot', '%s', info.messages{end});
    end
end

%% ------------------------------------------------------------------
%  3. Generalised long-run innovation variances L_i^(d) and the index nu_i
% -------------------------------------------------------------------
% L_i^(d_i) = || e_i' N^{d_i} B ||^2 / (d_i!)^2 (note, eq. (4)): the innovation
% variance of the TOP rung of state i's ladder. At d_i = 0 it is the old
% L_i = [Pi B B' Pi']_ii, the per-period variance of the permanent (random-walk)
% component of state i, so nothing changes for the I(1) states.
L = zeros(n, 1);
for ii = 1:n
    if ~isnan(d_all(ii))
        L(ii) = (rung_norms(ii, d_all(ii) + 1) / factorial(d_all(ii)))^2;
    end
end

% Scale against which L_i^(d) counts as zero, in state i's own units squared
% (see ZERO_L_TOL above). N^d B transforms like B under a diagonal rescaling of
% the states, so the same scale serves every order. A state with no permanent
% component has L_i equal to pure round-off, ~1e-32 against a scale of order one.
L_scale = diag(B*B') + diag(Omega);

for kk = 1:numel(state_indices)
    ii = state_indices(kk);
    if ~isnan(d_all(ii)) && L(ii) > ZERO_L_TOL * max(L_scale(ii), realmin)
        % Eq. (5) of the note: h^(2d+1) * (1 - R2_{h,i}) -> (2d+1)*Omega_ii/L_i^(d),
        % in periods^(2d+1) (quarters^(2d+1) in the applications). d = 0 gives
        % back the pre-2026-09 nu_i = Omega_ii / L_i exactly.
        d(kk)  = d_all(ii);
        nu(kk) = (2*d_all(ii) + 1) * Omega(ii, ii) / L(ii);
    else
        % No stochastic trend at any rung: a stationary state, or a
        % cointegrating direction. nu is undefined; report the paper's
        % unconditional R2 for these states instead.
        d(kk)  = NaN;
        nu(kk) = NaN;
    end
end

end


%% ========================================================================
function res = riesz_projector_residual(P, A, k)
%RIESZ_PROJECTOR_RESIDUAL  Relative violation of the defining properties of the
%   Riesz projector for the lambda = 1 cluster of A. Returns the largest of
%
%       ||P*P - P|| / ||P||        (idempotent)
%       ||A*P - P*A|| / ||P||      (commutes with A, so range(P) is A-invariant)
%       |trace(P) - k| / max(1,k)  (range(P) has the right dimension: the trace
%                                   of a projector is the dimension of its range)
%
%   A*P = P is deliberately NOT tested: it characterizes the SEMISIMPLE case and
%   fails, correctly, whenever there is a Jordan block at 1. Whether it holds is
%   determined separately in section 2b through N = (A - I)*P.
    scale = max(1, norm(P, 'fro'));
    res   = max([norm(P*P - P,   'fro') / scale, ...
                 norm(A*P - P*A, 'fro') / scale, ...
                 abs(trace(P) - k) / max(1, k)]);
end
