function R = model_assessment(D, opts)
%MODEL_ASSESSMENT  Posterior-predictive assessment with a damped Newton (Le Cam) re-estimator.
%
%   R = MODEL_ASSESSMENT(D)          D artificial datasets (default 300)
%   R = MODEL_ASSESSMENT(D, opts)    see "Options" below
%
%   REWRITTEN 2026-09-15 (was: 50,000 replications at a single parameter
%   vector). Implements the idea of an unpublished note on a one-step mode
%   update (14 Sep 2026): for each artificial dataset ytilde^(d), simulated at
%   the posterior mode theta_hat, approximate the estimate an econometrician
%   facing that dataset would have used by Newton steps from the mode with a
%   FIXED curvature matrix,
%
%       theta_{i+1} = theta_i + lambda_i * Sigma * [ grad ell(theta_i; ytilde) + grad pi(theta_i) ],
%
%   starting at theta_0 = theta_hat, where the posterior gradient is the note's
%   differenced score s(ytilde) - s(y) (because grad pi(theta_hat) = -s(y)),
%   and at later iterates the two gradients are recomputed by central
%   differences. lambda_i is chosen by backtracking on the log posterior of
%   ytilde (opts.lambdas); up to opts.newton_iter iterations (default 3).
%
%   WHY NOT THE NOTE'S PLAIN ONE-STEP (findings on HLW_2023_US_no_lags, 2026-09-15):
%   - Dynare's Sigma (oo_.posterior.optimization.Variance) is exactly the
%     inverse observed Hessian of the log posterior on the actual data at the
%     mode (verified against our own numerical Hessian, ratio 1.00), but the
%     FISHER information at theta_hat - the covariance of the likelihood scores
%     over datasets simulated from the model - disagrees with it: generalised
%     eigenvalues of observed vs Fisher information from 0.09 to 1.9 (phi, by,
%     rho, a_r much less informative in the data than the model implies). A
%     Newton step scaled by the observed Sigma therefore overshoots by up to
%     10x on model-generated data (log posterior on ytilde FELL by 100+ at
%     lambda = 1). This is the information-matrix identity failing, i.e. a
%     misspecification diagnostic in its own right; it is reported.
%   - Hence the default curvature is opts.curvature = 'fisher':
%     Sigma_F = [I_F(theta_hat) - Hessian(log prior)(theta_hat)]^-1, with I_F
%     estimated from the D simulated scores (they are needed for the first step
%     anyway, so it is free), and steps are damped. 'hessian' uses Dynare's
%     Sigma. The note's estimator is opts.newton_iter = 1, opts.lambdas = 1,
%     opts.curvature = 'hessian'.
%
%   Three sets of statistics are produced and compared:
%     (i)   at theta_hat on the ACTUAL data y                       -> R.data
%     (ii)  at theta_hat on each ytilde^(d)                         -> R.sim_mode
%     (iii) at theta_hat^(d) on ytilde^(d) (re-estimated)           -> R.sim_lecam
%   Statistics: the engine's empirical objects - TV_RMSE, DeltaR2 =
%   CorrH_Observables.^2 - TV_CorrY_H.^2 (observables x horizons), and the
%   smoother-revision regression fr (R2_emp, alpha_emp, beta_emp, R2_theory)
%   for the base states of M_.CellList_States. (i) vs (ii) is the parametric
%   bootstrap at a fixed parameter; (iii) adds the estimation step. Two-sided
%   percentile p-values are printed at opts.h_report. The validation against
%   the full mode-finder (note, last paragraph) is NOT run here.
%
%   OTHER DEVIATIONS FROM THE NOTE, AND WHY
%   - Parameterisation. Dynare's mode_compute = 4 (csminwel) optimises in
%     NATURAL units with the prior support as bounds, and Sigma is the inverse
%     Hessian in those units, so the steps are taken in natural units; the
%     line search rejects points outside the prior support (prior_bounds).
%   - Heteroskedastic scales. Dynare's driver evaluates the `scales`
%     expressions of the heteroskedastic_shocks block ONCE, at the calibrated
%     parameter values, and estimation never re-evaluates them: the kappa
%     parameters do not enter Dynare's likelihood (their posterior is the
%     prior, and their modes are prior modes), and theta_hat is the mode of the
%     objective with the FROZEN scales stored in M_.heteroskedastic_shocks.Qscale.
%     A Newton update is coherent only around that objective, so the default
%     opts.qscale = 'stored' uses the frozen matrix for every likelihood
%     evaluation and every simulation (verified: standardised posterior
%     gradient at theta_hat = 0.000 in every parameter). 'rebuild' recomputes
%     the scales from the current kappas (build_qscale). The mode condition is
%     checked and printed under either convention.
%
%   Options (struct opts, all optional):
%     .curvature       'fisher' (default) | 'hessian'
%     .newton_iter     max Newton iterations (3)
%     .lambdas         backtracking step sizes ([1 0.5 0.25 0.1])
%     .tol_step        stop when the accepted step is below this many sd (0.05)
%     .fisher_draws    minimum number of simulated score vectors used for I_F
%                      (200); if D is smaller, extra datasets are simulated for
%                      the Fisher estimate only (each costs ~0.2 s)
%     .qscale          'stored' (default) | 'rebuild'
%     .c               finite-difference constant, h_j = c*sd_j (1e-3)
%     .seed            RNG seed for the artificial data (20260915)
%     .workspace       .mat with the Dynare workspace after estimation and
%                      calib_smoother ('tmp.mat', saved by HLW_2023_US_no_lags.mod)
%     .CellList_States states to diagnose (default: the HLW .mod's list)
%     .fr_smoother_lags m grid of the f-on-r regression ([4 20])
%     .h_report        horizons printed in the tables ([1 4 8 20 40])
%     .save_file       output .mat ('model_assessment_lecam.mat'; '' = no save)
%     .model_dir       folder of the Dynare run: +<fname> package, .workspace,
%                      and where .save_file is written ('' = Code/HLW_new/hlw_level,
%                      the HLW run; UC3 passes Code/UC3) (ADDED 2026-09-17)
%     .legacy_list     state names of the cheap likelihood system ({} = the base
%                      states of CellList_States) (ADDED 2026-09-17)
%     .mask_ignore_shocks  exo names whose heteroskedastic scales are a data
%                      device (SPF: {'e_spfobs'}) and must not exclude periods
%                      from the empirical statistics; forwarded to the engine as
%                      M_.reliability_mask_ignore_shocks ({} = none) (ADDED 2026-09-19).
%                      CHANGED 2026-09-26: the artificial datasets are simulated
%                      with UNIT variance for these shocks and the observations
%                      they load on are then set to the placeholder (obs_ss)
%                      wherever the stored scale is not 1 - the treatment the
%                      actual data received. Until now the placeholder cells were
%                      simulated at the device's scale (SPF: spf_obs pre-1992 with
%                      s.d. 10000x), which made the RMSE/DeltaR2 bands of THAT
%                      observable meaningless (median simulated RMSE 199 vs 0.86
%                      in the data); nothing else was affected, since the filter
%                      gives those cells weight 1e-8 either way.
%     Models estimated OUTSIDE Dynare (user solver, see solve_model_draw; ADDED 2026-09-26):
%     .theta_hat       posterior mode (free parameters, bayestopt_.name order);
%                      default oo_.posterior.optimization.mode
%     .Sigma_H         its covariance (inverse observed Hessian or a posterior
%                      covariance); default oo_.posterior.optimization.Variance
%     .lprior          function handle theta -> log prior density (up to a constant);
%                      default Dynare's priordens on bayestopt_
%     .bounds          struct with .lb/.ub (k x 1) of the prior support used by
%                      the line search; default prior_bounds(bayestopt_, prior_trunc)
%     The workspace's M_.X_InitialDistribution (incl. .type) is kept when present
%     (else the HLW default X0 = 0, B0 = 100 I is set), and the oo_.bayestopt_*
%     smoother fields are taken from oo_ itself when bayestopt_ lacks them
%     (dr_from_abcd sets them on every solve).
%   Models without a heteroskedastic_shocks block are handled (Qscale = 1).
%
%   REQUIREMENTS: Dynare on the path (dynare_config run), ReliabilityToolbox
%   on the path (added here); the model folder is added from .model_dir.
%
%   See also: compute_state_reliability_diagnostics, apply_tv_filter,
%             build_qscale, posterior_predictive_wrapper, evaluate_at_reweighted_means.

%% ------------------------------------------------------------------ setup
if nargin < 1 || isempty(D), D = 300; end
if nargin < 2, opts = struct(); end
def = struct('curvature', 'fisher', 'newton_iter', 3, 'lambdas', [1 0.5 0.25 0.1], 'tol_step', 0.05, ...
             'fisher_draws', 200, 'qscale', 'stored', 'c', 1e-3, 'seed', 20260915, 'workspace', 'tmp.mat', ...
             'CellList_States', {{'ygap', [0 1 4 20], []; 'rrstar', [0 1 4 20], [4 20]; 'drrstar', [0 1], []}}, ...
             'fr_smoother_lags', [4 20], 'h_report', [1 4 8 20 40], 'save_file', 'model_assessment_lecam.mat', ...
             'model_dir', '', 'legacy_list', {{}}, 'mask_ignore_shocks', {{}}, ...
             'theta_hat', [], 'Sigma_H', [], 'lprior', [], 'bounds', [], 'solve_method', 'stoch_simul');
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = def.(fn{i}); end
end
if ~any(strcmp(opts.qscale, {'stored', 'rebuild'})),      error('model_assessment:opts', 'opts.qscale must be ''stored'' or ''rebuild''.'); end
if ~any(strcmp(opts.curvature, {'fisher', 'hessian'})),   error('model_assessment:opts', 'opts.curvature must be ''fisher'' or ''hessian''.'); end

% MODEL FOLDER (GENERALISED 2026-09-17): the folder holding the Dynare run -
% the +<fname> package, opts.workspace and, on exit, opts.save_file. This file
% lives in ReliabilityToolbox (MOVED there 2026-09-18 from HLW_new/hlw_level);
% the default model folder is the HLW run, Code/HLW_new/hlw_level; other models
% pass their own (UC3: Code/UC3).
here = fileparts(mfilename('fullpath'));          % ReliabilityToolbox
if isempty(opts.model_dir), model_dir = fullfile(here, '..', 'HLW_new', 'hlw_level'); else, model_dir = opts.model_dir; end
addpath(model_dir);                               % +<fname> package of the model
addpath(here);                                    % the toolbox itself (not genpath: octave_compat/ is Octave-only)
if ~exist('stoch_simul', 'file')
    error('model_assessment:dynare', 'Dynare is not on the path (addpath <dynare>/matlab; dynare_config).');
end
t_all = tic;

W = load(fullfile(model_dir, opts.workspace), 'M_', 'oo_', 'options_', 'estim_params_', 'bayestopt_', 'dataset_', 'dataset_info', 'var_list_');
M_ = W.M_; oo_ = W.oo_; options_ = W.options_; estim_params_ = W.estim_params_; bayestopt_ = W.bayestopt_;
dataset_ = W.dataset_; dataset_info = W.dataset_info; var_list_ = W.var_list_;

% the .mod's post-save block
M_.observable_names   = dataset_.name;
M_.data               = dataset_info.rawdata(~any(isnan(dataset_info.rawdata')), :);
M_.CellList_States    = opts.CellList_States;
M_.states_of_interest = expand_states_of_interest(M_.CellList_States);
M_.fr_smoother_lags   = opts.fr_smoother_lags;
M_.var_list_          = var_list_;
M_.verbose            = false;
if ~isempty(opts.mask_ignore_shocks), M_.reliability_mask_ignore_shocks = opts.mask_ignore_shocks; end   % carried into M_t via ctx.M_
if isempty(opts.solve_method), opts.solve_method = 'stoch_simul'; end                        % drivers may pass getenv(...) directly
if ~strcmpi(opts.solve_method, 'stoch_simul'), M_.solve_method = opts.solve_method; end   % 'backward': closed form for models without leads (solve_model_draw, 2026-09-28)
nstates = length(oo_.dr.state_var);
if ~isfield(M_, 'X_InitialDistribution') || ~isfield(M_.X_InitialDistribution, 'B0')   % keep a workspace's own initialisation (user-solver models, 2026-09-26)
    M_.X_InitialDistribution.X0 = zeros(nstates, 1);
    M_.X_InitialDistribution.B0 = eye(nstates) * 100;
end
if isfield(bayestopt_, 'smoother_var_list')
    oo_.bayestopt_smoother_var_list         = bayestopt_.smoother_var_list;
    oo_.bayestopt_smoother_restrict_columns = bayestopt_.smoother_restrict_columns;
    oo_.bayestopt_mf                        = bayestopt_.mf;
elseif ~isfield(oo_, 'bayestopt_smoother_var_list')
    error('model_assessment:bookkeeping', 'neither bayestopt_ nor oo_ carries the smoother bookkeeping (smoother_var_list, restrict_columns, mf).');
end
oo_.dr.qz_crit = 1 + 1e-6;

y     = M_.data;
T     = size(y, 1);
n_obs = size(y, 2);
obs_names = M_.observable_names;

% posterior mode and inverse observed Hessian (natural units, see header); a model
% estimated outside Dynare supplies them through opts (2026-09-26)
if ~isempty(opts.theta_hat), theta_hat = opts.theta_hat(:); else, theta_hat = oo_.posterior.optimization.mode(:); end
if ~isempty(opts.Sigma_H),   Sigma_H   = opts.Sigma_H;       else, Sigma_H   = oo_.posterior.optimization.Variance; end
Sigma_H   = (Sigma_H + Sigma_H') / 2;
k         = numel(theta_hat);
psd       = sqrt(diag(Sigma_H));
pnames    = bayestopt_.name;
if k ~= numel(pnames) || any(psd <= 0)
    error('model_assessment:theta', 'theta_hat must have one entry per bayestopt_.name (%d) with a positive posterior s.d. for each.', numel(pnames));
end
if ~isempty(opts.bounds)
    lb = opts.bounds.lb(:);  ub = opts.bounds.ub(:);
else
    bounds = prior_bounds(bayestopt_, options_.prior_trunc);
    lb = bounds.lb(:);  ub = bounds.ub(:);
end
hstep     = opts.c * psd;
if ~isempty(opts.lprior)
    lprior = opts.lprior;
else
    lprior = @(x) priordens(x, bayestopt_.pshape, bayestopt_.p6, bayestopt_.p7, bayestopt_.p3, bayestopt_.p4);
end

% frozen scales Dynare estimated with (variance scales, NaN = 1), T columns.
% A model without a heteroskedastic_shocks block (UC3) is homoskedastic:
% Qscale = 1 throughout, and 'rebuild' has nothing to rebuild.
if isfield(M_, 'heteroskedastic_shocks') && isfield(M_.heteroskedastic_shocks, 'Qscale') && ~isempty(M_.heteroskedastic_shocks.Qscale)
    Qstored = M_.heteroskedastic_shocks.Qscale;
    Qstored(isnan(Qstored)) = 1;
    Qstored = Qstored(:, 1:T);
else
    Qstored = ones(M_.exo_nbr, T);
    if strcmp(opts.qscale, 'rebuild'), opts.qscale = 'stored'; end
end

% cheap system for the likelihood scores: the base states of the list (was the
% HLW literal {'ygap'; 'rrstar'}); any valid state names give the same likelihood
if isempty(opts.legacy_list), legacy_list = M_.CellList_States(:, 1); else, legacy_list = opts.legacy_list(:); end

fprintf('model_assessment (damped Newton / Le Cam): D = %d, curvature = %s, newton_iter = %d, qscale = %s, k = %d, T = %d, %d observables, model %s\n', ...
        D, opts.curvature, opts.newton_iter, opts.qscale, k, T, n_obs, M_.fname);

% shared arguments for the local solvers
ctx = struct('M_', M_, 'options_', options_, 'oo_', oo_, 'estim_params_', estim_params_, 'Qstored', Qstored, ...
             'qscale', opts.qscale, 'T', T, 'lb', lb, 'ub', ub, 'legacy_list', {legacy_list});

%% ------------------------------------ likelihood caches at theta_hat +/- h_j e_j
fprintf('  building %d likelihood caches at theta_hat ... ', 2*k + 1); t0 = tic;
[cache0, cp, cm] = score_caches(theta_hat, hstep, ctx, y);
fprintf('%.1fs\n', toc(t0));
ell_y_0 = loglik(cache0, y);
s_y     = score_from_caches(cp, cm, hstep, y);
dpi_hat = prior_grad(theta_hat, hstep, lprior);
grad_post_sd = psd .* (s_y + dpi_hat);
fprintf('  log-likelihood at the mode (%s scales): %.4f; mode condition max |sd_j (s_j + dpi_j)| = %.3f\n', opts.qscale, ell_y_0, max(abs(grad_post_sd)));

%% ------------------------------------------- (i) statistics at theta_hat on the data
fprintf('  statistics at theta_hat on the actual data ... '); t0 = tic;
P0 = solve_full(theta_hat, ctx, M_.CellList_States, y);
if ~P0.ok, error('model_assessment:solve', 'the model does not solve at the posterior mode.'); end
S_data = P0.LSD;
fprintf('%.1fs\n', toc(t0));
H = size(S_data.TV_RMSE, 2);
theory_R2_hat = S_data.CorrH_Observables.^2;
data = struct('RMSE', S_data.TV_RMSE, 'DeltaR2', theory_R2_hat - S_data.TV_CorrY_H.^2, ...
              'fr', S_data.fr, 'theta', theta_hat, 'loglik', ell_y_0, 'LSD', S_data);

%% -------------------------------- pass 1: simulate the D datasets, scores at theta_hat
G = cache0;                                  % DGP at theta_hat: the cheap legacy-list system (same observables)
n_st = size(G.A, 1);  n_sh = size(G.B, 2);
try, rng(opts.seed); catch, randn('state', opts.seed); end
BurnIn = 250;
Dsim = max(D, opts.fisher_draws);            % datasets 1..D get statistics; all Dsim feed the Fisher estimate
Ys  = zeros(T, n_obs, D);
S0  = zeros(k, Dsim);                        % likelihood scores s(ytilde^(d)) at theta_hat
L0  = zeros(D, 1);                           % ell(theta_hat; ytilde^(d))
% Missing-data devices (2026-09-26, see .mask_ignore_shocks in the header): the
% exempted shocks are simulated at unit variance and the observations they load
% on are set to the placeholder obs_ss wherever the stored scale is not 1.
Qsim = G.Q;
placeholder = false(T, n_obs);
if isfield(M_, 'reliability_mask_ignore_shocks') && ~isempty(M_.reliability_mask_ignore_shocks)
    ign = find(ismember(M_.exo_names, M_.reliability_mask_ignore_shocks));
    for j = ign(:)'
        obs_j = find(abs(G.D(:, j)) > 0);                       % observables the device shock loads on
        placeholder(G.Q(j, :) ~= 1, obs_j) = true;
        Qsim(j, :) = 1;
    end
    fprintf('  missing-data device: %d placeholder cells in %d observable(s) set to obs_ss in every artificial dataset\n', sum(placeholder(:)), numel(unique(find(any(placeholder, 1)))));
end
fprintf('  simulating %d datasets (%d for the statistics) and scoring them at theta_hat ... ', Dsim, D); t0 = tic;
for d = 1:Dsim
    s = zeros(n_st, 1);  Y = zeros(T, n_obs);
    for tt = 1:BurnIn + T
        if tt > BurnIn
            e = sqrt(Qsim(:, tt - BurnIn)) .* randn(n_sh, 1);
            Y(tt - BurnIn, :) = (G.C * s + G.D * e + G.obs_ss)';
        else
            e = randn(n_sh, 1);
        end
        s = G.A * s + G.B * e;
    end
    if any(placeholder(:))
        obs_ss_mat = repmat(G.obs_ss(:)', T, 1);
        Y(placeholder) = obs_ss_mat(placeholder);
    end
    S0(:, d) = score_from_caches(cp, cm, hstep, Y);
    if d <= D
        Ys(:, :, d) = Y;
        L0(d)       = loglik(cache0, Y);
    end
end
fprintf('%.1fs\n', toc(t0));

% Fisher information at theta_hat under the model, prior curvature, curvature matrix
I_F   = cov(S0');
Hpi   = prior_hess_diag(theta_hat, 1e-2 * psd, lprior);
Sigma_F = inv(I_F - Hpi);  Sigma_F = (Sigma_F + Sigma_F') / 2;
switch opts.curvature
    case 'fisher',  Sigma = Sigma_F;
    case 'hessian', Sigma = Sigma_H;
end
ge = eig(inv(Sigma_H), I_F - Hpi);
info_check = struct('sd_hessian', psd, 'sd_fisher', sqrt(diag(Sigma_F)), 'gen_eig_obs_vs_fisher', sort(ge), ...
                    'mean_score_sd', mean(S0, 2) .* psd, 'I_F', I_F, 'Sigma_F', Sigma_F, 'Sigma_H', Sigma_H, 'n_score_draws', Dsim);
fprintf('  information check: generalised eigenvalues of observed vs Fisher information in [%.3f, %.3f] (1 = identity holds)\n', min(ge), max(ge));

%% -------------------------------- pass 2: damped Newton per dataset, then statistics
% PARFOR over datasets (2026-09-15). Every iteration is independent given the
% pass-1 objects (caches at theta_hat, S0, L0, Ys, Sigma, P0), which are
% broadcast to the workers once (~35 Sequence structs, a few hundred MB on the
% 77-state system). parfor requires that everything written inside the loop be
% a plain array indexed by d ("sliced"), not a struct field, and assigned on
% every path - hence the per-dataset temporaries below and the assembly into
% sim_mode / sim_lecam after the loop. Under Octave or without a parallel pool
% parfor runs as a plain for. The progress line prints out of order in parallel.
J = numel(S_data.fr.states);  n_m = numel(S_data.fr.m);
z2 = nan(n_obs, H, D);  z4 = nan(J, H, n_m, D);
RMSE_m = z2;  DR2_m = z2;  frR2e_m = z4;  frA_m = z4;  frB_m = z4;  frR2t_m = z4;   % (ii)  at theta_hat
RMSE_l = z2;  DR2_l = z2;  frR2e_l = z4;  frA_l = z4;  frB_l = z4;  frR2t_l = z4;   % (iii) at theta_hat^(d)
frOR_l = nan(J, n_m, D);  frRb_l = nan(J, H, D);    % (iii) [Omega_m]_ii/Omega_ii and R2_h of the base state at theta_hat^(d) (ADDED 2026-09-17)
theta_l  = nan(k, D);   iter_l  = zeros(D, 1);   lambda_l = nan(opts.newton_iter, D);
gain_l   = nan(D, 1);   gnorm_l = nan(D, 1);     step_l   = nan(D, 1);
solve_l  = false(D, 1); noasc_l = false(D, 1);
n_it = opts.newton_iter;  lambdas = opts.lambdas;  tol_step = opts.tol_step;  CLS = M_.CellList_States;
n_lag_report = 10;
fprintf('  %d datasets: damped Newton (<= %d it.), statistics at theta_hat and at theta_hat^(d)\n', D, n_it);
t0 = tic;
parfor d = 1:D
    Y = Ys(:, :, d);

    % --- damped Newton from theta_hat with fixed Sigma ------------------------
    th   = theta_hat;
    g    = S0(:, d) - s_y;                       % grad ell(theta_hat; Y) + grad pi(theta_hat)
    Lcur = L0(d) + lprior(theta_hat);
    Lstart = Lcur;
    lam_d = nan(n_it, 1);  it_d = 0;  step_d = NaN;  gnorm_d = NaN;  noasc_d = false;
    for it = 1:n_it
        dirn = Sigma * g;
        accepted = false;  lam = NaN;  th_try = th;  L_try = Lcur;
        for lam_try = lambdas
            th_c = th + lam_try * dirn;
            if any(th_c < lb | th_c > ub), continue; end
            c_try = solve_cheap(th_c, ctx, Y);
            if ~c_try.ok, continue; end
            L_c = loglik(c_try, Y) + lprior(th_c);
            if L_c > Lcur
                accepted = true;  lam = lam_try;  th_try = th_c;  L_try = L_c;
                break;
            end
        end
        if ~accepted
            if it == 1, noasc_d = true; end
            break;
        end
        step_d = sqrt((lam * dirn)' * (Sigma \ (lam * dirn)));
        th = th_try;  Lcur = L_try;
        lam_d(it) = lam;  it_d = it;
        if step_d < tol_step || it == n_it, break; end
        % gradient at the new point (2k solves): likelihood score by central differences + prior gradient.
        % FIXED 2026-09-29: a perturbation point theta +/- h_j e_j that does not solve (determinacy /
        % support boundary; ECM 2-step dataset ~260, parameter 2) used to throw inside the parfor and
        % abort the whole assessment; now the Newton iterations of that dataset stop at the current point.
        try
            [~, cp_i, cm_i] = score_caches(th, hstep, ctx, Y);
        catch
            gnorm_d = NaN;
            break;
        end
        g = score_from_caches(cp_i, cm_i, hstep, Y) + prior_grad(th, hstep, lprior);
        gnorm_d = sqrt(g' * Sigma * g);
    end
    theta_l(:, d)  = th;     iter_l(d)  = it_d;    lambda_l(:, d) = lam_d;
    gain_l(d)      = Lcur - Lstart;   gnorm_l(d) = gnorm_d;   step_l(d) = step_d;   noasc_l(d) = noasc_d;

    % --- (ii) statistics at theta_hat on ytilde^(d) ---------------------------
    Sm = stats_on_data(P0, Y);
    RMSE_m(:, :, d)     = Sm.TV_RMSE;
    DR2_m(:, :, d)      = theory_R2_hat - Sm.TV_CorrY_H.^2;
    frR2e_m(:, :, :, d) = Sm.fr.R2_emp;
    frA_m(:, :, :, d)   = Sm.fr.alpha_emp;
    frB_m(:, :, :, d)   = Sm.fr.beta_emp;
    frR2t_m(:, :, :, d) = Sm.fr.R2_theory;

    % --- (iii) statistics at theta_hat^(d) on ytilde^(d) ----------------------
    Pd = solve_full(th, ctx, CLS, Y);
    rm_d = nan(n_obs, H);  dr_d = nan(n_obs, H);  fe_d = nan(J, H, n_m);  fa_d = fe_d;  fb_d = fe_d;  ft_d = fe_d;
    or_d = nan(J, n_m);  rb_d = nan(J, H);
    if Pd.ok
        Sd = Pd.LSD;
        rm_d = Sd.TV_RMSE;
        dr_d = Sd.CorrH_Observables.^2 - Sd.TV_CorrY_H.^2;
        fe_d = Sd.fr.R2_emp;  fa_d = Sd.fr.alpha_emp;  fb_d = Sd.fr.beta_emp;  ft_d = Sd.fr.R2_theory;
        or_d = Sd.fr.omega_ratio;  rb_d = Sd.fr.R2_h_base;
    end
    solve_l(d)          = Pd.ok;
    RMSE_l(:, :, d)     = rm_d;
    DR2_l(:, :, d)      = dr_d;
    frR2e_l(:, :, :, d) = fe_d;
    frA_l(:, :, :, d)   = fa_d;
    frB_l(:, :, :, d)   = fb_d;
    frR2t_l(:, :, :, d) = ft_d;
    frOR_l(:, :, d)     = or_d;
    frRb_l(:, :, d)     = rb_d;

    if mod(d, n_lag_report) == 0 || d == D
        fprintf('    d = %3d / %d  (%.0fs elapsed)\n', d, D, toc(t0));
    end
end
sim_mode  = struct('RMSE', RMSE_m, 'DeltaR2', DR2_m, 'fr_R2_emp', frR2e_m, 'fr_alpha_emp', frA_m, 'fr_beta_emp', frB_m, 'fr_R2_theory', frR2t_m);
sim_lecam = struct('RMSE', RMSE_l, 'DeltaR2', DR2_l, 'fr_R2_emp', frR2e_l, 'fr_alpha_emp', frA_l, 'fr_beta_emp', frB_l, 'fr_R2_theory', frR2t_l, ...
                   'fr_omega_ratio', frOR_l, 'fr_R2_h_base', frRb_l, ...
                   'theta', theta_l, 'iterations', iter_l, 'lambda', lambda_l, ...
                   'logpost_gain', gain_l, 'grad_norm_final_sd', gnorm_l, 'step_norm_sd', step_l, ...
                   'solve_ok', solve_l, 'no_ascent_at_start', noasc_l);

%% ------------------------------------------------------------------ report
R = struct('D', D, 'opts', opts, 'theta_hat', theta_hat, 'Sigma', Sigma, 'param_names', {pnames}, ...
           'mode_check', struct('s_y', s_y, 'dprior', dpi_hat, 'grad_post_sd', grad_post_sd), ...
           'info_check', info_check, 'data', data, 'sim_mode', sim_mode, 'sim_lecam', sim_lecam, ...
           'obs_names', {obs_names}, 'states_of_interest', {M_.states_of_interest}, 'h', 1:H, 'fr_m', S_data.fr.m, ...
           'fr_states', {S_data.fr.states}, 'elapsed_s', toc(t_all));
print_report(R);

if ~isempty(opts.save_file)
    Rs = R;
    if isa(Rs.opts.lprior, 'function_handle'), Rs.opts.lprior = func2str(Rs.opts.lprior); end   % handles do not survive save under Octave (2026-09-26)
    save(fullfile(model_dir, opts.save_file), '-struct', 'Rs', '-v7');
    fprintf('saved %s\n', fullfile(model_dir, opts.save_file));
end
end


%% ========================================================================
%  local functions
%% ========================================================================
function c = solve_cheap(theta, ctx, ydata)
%SOLVE_CHEAP  Solve at theta, build the data-independent diffuse/TV Sequence (loglik_only engine).
c = struct('ok', false);
M_t = set_all_parameters(theta, ctx.estim_params_, ctx.M_);
M_t.verbose = false;
[info, oo_t, ~, M_t] = solve_model_draw(M_t, ctx.options_, ctx.oo_);   % stoch_simul or M_.solve_function (2026-09-26)
if info ~= 0, return; end
oo_t.dr.qz_crit = 1 + 1e-6;
M_t.heteroskedastic_shocks.Qscale = qscale_of(M_t, ctx);
M_t.loglik_only = true;
[LSD, MM] = compute_state_reliability_diagnostics(M_t.observable_names, ctx.legacy_list, oo_t, M_t, ydata);
c = struct('ok', true, 'A', MM.A, 'B', MM.B, 'C', MM.C, 'D', MM.D, 'Seq', LSD.TV_KF_Sequence, ...
           'X0', MM.X0_pruned, 'obs_ss', MM.obs_ss(:), 'Q', M_t.heteroskedastic_shocks.Qscale, 'theta', theta);
end

function P = solve_full(theta, ctx, state_list, ydata)
%SOLVE_FULL  Solve at theta and run the fast-mode engine with the full state list on ydata.
P = struct('ok', false);
M_t = set_all_parameters(theta, ctx.estim_params_, ctx.M_);
M_t.verbose = false;
[info, oo_t, ~, M_t] = solve_model_draw(M_t, ctx.options_, ctx.oo_);   % stoch_simul or M_.solve_function (2026-09-26)
if info ~= 0, return; end
oo_t.dr.qz_crit = 1 + 1e-6;
M_t.heteroskedastic_shocks.Qscale = qscale_of(M_t, ctx);
M_t.fast_mode = true;
[LSD, MM, ~, ~, ~, ~, ll] = compute_state_reliability_diagnostics(M_t.observable_names, state_list, oo_t, M_t, ydata);
P = struct('ok', true, 'LSD', LSD, 'MM', MM, 'loglik', ll, 'theta', theta, 'oo_t', oo_t, 'M_t', M_t, 'state_list', {state_list});
end

function Q = qscale_of(M_t, ctx)
switch ctx.qscale
    case 'stored',  Q = ctx.Qstored;
    case 'rebuild', Q = build_qscale(M_t, ctx.T);
end
end

function ll = loglik(c, ydata)
%LOGLIK  Diffuse/TV log-likelihood of ydata at a cached solution point.
ll = apply_tv_filter(c.A, c.C, c.Seq, ydata - c.obs_ss', c.X0, 1e-6);
end

function [c0, cp, cm] = score_caches(theta, hstep, ctx, ydata)
%SCORE_CACHES  Solution caches at theta and theta +/- h_j e_j (2k+1 solves).
k = numel(theta);
c0 = solve_cheap(theta, ctx, ydata);
if ~c0.ok, error('model_assessment:solve', 'the model does not solve at the requested point.'); end
cp = cell(k, 1);  cm = cell(k, 1);
for j = 1:k
    e = zeros(k, 1);  e(j) = hstep(j);
    cp{j} = solve_cheap(theta + e, ctx, ydata);
    cm{j} = solve_cheap(theta - e, ctx, ydata);
    if ~cp{j}.ok || ~cm{j}.ok
        error('model_assessment:solve', 'the model does not solve at a perturbation point (parameter %d).', j);
    end
end
end

function s = score_from_caches(cp, cm, hstep, ydata)
%SCORE_FROM_CACHES  Central-difference likelihood score of ydata.
k = numel(cp);  s = zeros(k, 1);
for j = 1:k
    s(j) = (loglik(cp{j}, ydata) - loglik(cm{j}, ydata)) / (2 * hstep(j));
end
end

function g = prior_grad(theta, hstep, lprior)
k = numel(theta);  g = zeros(k, 1);
for j = 1:k
    e = zeros(k, 1);  e(j) = hstep(j);
    g(j) = (lprior(theta + e) - lprior(theta - e)) / (2 * hstep(j));
end
end

function Hd = prior_hess_diag(theta, h, lprior)
%PRIOR_HESS_DIAG  Hessian of the log prior (diagonal: priors are independent).
k = numel(theta);  Hd = zeros(k);  L0 = lprior(theta);
for j = 1:k
    e = zeros(k, 1);  e(j) = h(j);
    Hd(j, j) = (lprior(theta + e) - 2 * L0 + lprior(theta - e)) / h(j)^2;
end
end

function S = stats_on_data(P, ydata)
%STATS_ON_DATA  Re-run the fast-mode engine at an already solved point on new data.
S = compute_state_reliability_diagnostics(P.M_t.observable_names, P.state_list, P.oo_t, P.M_t, ydata);
end

function print_report(R)
hs = R.opts.h_report(R.opts.h_report <= numel(R.h));
pv = @(p) 2 * min(p, 1 - p);
k  = numel(R.theta_hat);

fprintf('\n==================== damped-Newton (Le Cam) posterior-predictive assessment ====================\n');
fprintf('D = %d datasets simulated at the posterior mode; curvature = %s; qscale = %s; elapsed %.0f s\n', R.D, R.opts.curvature, R.opts.qscale, R.elapsed_s);
fprintf('re-estimation: %d datasets with no ascent direction at theta_hat, %d final-point solves failed\n', sum(R.sim_lecam.no_ascent_at_start), sum(~R.sim_lecam.solve_ok));

fprintf('\n--- information check at theta_hat: posterior sd from the observed Hessian on y vs the Fisher information under the model (%d simulated scores) ---\n', R.info_check.n_score_draws);
for j = 1:k
    fprintf('  %-10s Hessian %.4f   Fisher %.4f   ratio F/H %5.2f   mode-condition grad %6.3f sd\n', R.param_names{j}, ...
            R.info_check.sd_hessian(j), R.info_check.sd_fisher(j), R.info_check.sd_fisher(j) / R.info_check.sd_hessian(j), R.mode_check.grad_post_sd(j));
end
fprintf('  generalised eigenvalues of observed vs Fisher information: %s\n', sprintf('%.2f ', R.info_check.gen_eig_obs_vs_fisher));

fprintf('\n--- re-estimates theta_hat^(d): iterations, step sizes, and (theta^(d) - theta_hat) in Fisher-sd units ---\n');
it = R.sim_lecam.iterations;
fprintf('  iterations used: %s (of max %d)\n', sprintf('%d:%d%% ', [unique(it)'; round(100 * histc(it, unique(it))' / R.D)]), R.opts.newton_iter);
for i = 1:R.opts.newton_iter
    lam = R.sim_lecam.lambda(i, :);  lam = lam(~isnan(lam));
    if isempty(lam), continue; end
    fprintf('  iteration %d: lambda = 1 in %.0f%%, 0.5 in %.0f%%, 0.25 in %.0f%%, 0.1 in %.0f%% of %d datasets\n', i, ...
            100 * mean(lam == 1), 100 * mean(lam == 0.5), 100 * mean(lam == 0.25), 100 * mean(lam == 0.1), numel(lam));
end
fprintf('  log-posterior gain on ytilde over the iterations: median %.2f, 5%% %.2f, 95%% %.2f (k/2 = %.1f)\n', ...
        median(R.sim_lecam.logpost_gain), quantile(R.sim_lecam.logpost_gain, 0.05), quantile(R.sim_lecam.logpost_gain, 0.95), k / 2);
gn = R.sim_lecam.grad_norm_final_sd(~isnan(R.sim_lecam.grad_norm_final_sd));
if ~isempty(gn), fprintf('  Sigma-norm of the posterior gradient at the last recomputed point: median %.2f (0 = at the mode)\n', median(gn)); end
dz = (R.sim_lecam.theta - R.theta_hat) ./ R.info_check.sd_fisher;
for j = 1:k
    fprintf('  %-10s mean %7.3f   sd %6.3f\n', R.param_names{j}, mean(dz(j, :), 'omitnan'), std(dz(j, :), 'omitnan'));
end
dth = R.sim_lecam.theta - R.theta_hat;
sn = sqrt(sum(dth .* (R.Sigma \ dth), 1));
fprintf('  ||theta^(d) - theta_hat||_Sigma: median %.2f, 90%% %.2f, max %.2f  (chi_k, k = %d: median %.2f)\n', ...
        median(sn), quantile(sn, 0.9), max(sn), k, sqrt(k - 2/3));

for stat = {'RMSE', 'DeltaR2'}
    fprintf('\n--- %s: data at theta_hat | null at theta_hat: median [5%%,95%%] p | null at theta_hat^(d) (re-estimated): median [5%%,95%%] p ---\n', stat{1});
    dat = R.data.(stat{1});  A = R.sim_mode.(stat{1});  B = R.sim_lecam.(stat{1});
    for i = 1:numel(R.obs_names)
        fprintf('  %s\n', R.obs_names{i});
        for h = hs
            a = squeeze(A(i, h, :));  b = squeeze(B(i, h, :));
            fprintf('    h=%2d  %8.4f | %8.4f [%8.4f,%8.4f] p=%.2f | %8.4f [%8.4f,%8.4f] p=%.2f\n', h, dat(i, h), ...
                    median(a, 'omitnan'), quantile(a, 0.05), quantile(a, 0.95), pv(mean(a <= dat(i, h), 'omitnan')), ...
                    median(b, 'omitnan'), quantile(b, 0.05), quantile(b, 0.95), pv(mean(b <= dat(i, h), 'omitnan')));
        end
    end
end

fprintf('\n--- smoother-revision regression f = alpha + beta r (base states): data | null at theta_hat | null at theta_hat^(d) ---\n');
for j = 1:numel(R.fr_states)
    if all(isnan(R.data.fr.R2_emp(j, :, :))), fprintf('  %s: degenerate (zero filter gain), skipped\n', R.fr_states{j}); continue; end
    for im = 1:numel(R.fr_m)
        fprintf('  %s, m = %d\n', R.fr_states{j}, R.fr_m(im));
        for h = hs
            for f = {'R2_emp', 'beta_emp', 'alpha_emp'}
                dat = R.data.fr.(f{1})(j, h, im);
                a = squeeze(R.sim_mode.(['fr_' f{1}])(j, h, im, :));  b = squeeze(R.sim_lecam.(['fr_' f{1}])(j, h, im, :));
                fprintf('    h=%2d %-9s %7.3f | %7.3f [%7.3f,%7.3f] p=%.2f | %7.3f [%7.3f,%7.3f] p=%.2f\n', h, f{1}, dat, ...
                        median(a, 'omitnan'), quantile(a, 0.05), quantile(a, 0.95), pv(mean(a <= dat, 'omitnan')), ...
                        median(b, 'omitnan'), quantile(b, 0.05), quantile(b, 0.95), pv(mean(b <= dat, 'omitnan')));
            end
            bt = squeeze(R.sim_lecam.fr_R2_theory(j, h, im, :));
            fprintf('    h=%2d R2_theory %7.3f |  (fixed at theta_hat)              | %7.3f [%7.3f,%7.3f]\n', h, R.data.fr.R2_theory(j, h, im), ...
                    median(bt, 'omitnan'), quantile(bt, 0.05), quantile(bt, 0.95));
        end
    end
end

% Implied reliability (ADDED 2026-09-17): the empirical P^2 mapped into R^2_{i,h}
% by Prop. 5.1, R^2 = P (1 - rho) / (1 - rho P) with rho = [Omega_m]_ii / Omega_ii
% (fr_implied_R2). Under the null at theta_hat rho is the data's; under the
% re-estimated null it is each dataset's own. The last column is the model's
% R^2_{i,h} at theta_hat, which the implied values should scatter around if
% the model is right.
if isfield(R.data.fr, 'omega_ratio') && isfield(R.sim_lecam, 'fr_omega_ratio')
    fprintf('\n--- R^2_{i,h} implied by the empirical P^2 (Prop. 5.1): data | null at theta_hat | null at theta_hat^(d) | model R^2_{i,h} at theta_hat ---\n');
    for j = 1:numel(R.fr_states)
        if all(isnan(R.data.fr.R2_emp(j, :, :))), continue; end
        for im = 1:numel(R.fr_m)
            rho0 = R.data.fr.omega_ratio(j, im);
            fprintf('  %s, m = %d  (rho = [Omega_m]_ii/Omega_ii = %.3f)\n', R.fr_states{j}, R.fr_m(im), rho0);
            for h = hs
                dat = fr_implied_R2(R.data.fr.R2_emp(j, h, im), rho0);
                a   = fr_implied_R2(squeeze(R.sim_mode.fr_R2_emp(j, h, im, :)), rho0);
                b   = fr_implied_R2(squeeze(R.sim_lecam.fr_R2_emp(j, h, im, :)), squeeze(R.sim_lecam.fr_omega_ratio(j, im, :)));
                fprintf('    h=%2d implied R2 %7.3f | %7.3f [%7.3f,%7.3f] p=%.2f | %7.3f [%7.3f,%7.3f] p=%.2f | %7.3f\n', h, dat, ...
                        median(a, 'omitnan'), quantile(a, 0.05), quantile(a, 0.95), pv(mean(a <= dat, 'omitnan')), ...
                        median(b, 'omitnan'), quantile(b, 0.05), quantile(b, 0.95), pv(mean(b <= dat, 'omitnan')), R.data.fr.R2_h_base(j, h));
            end
        end
    end
end
fprintf('================================================================================================\n');
end
