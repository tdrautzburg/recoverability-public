% TEST_BACKWARD_LOOKING_DR  backward_looking_dr against stoch_simul (through solve_model_draw) for HLW and UC3.
%
%   ADDED 2026-09-28. Loads the Dynare workspaces tmp.mat of the HLW level model and
%   of UC3 (both purely backward-looking), compares ghx / ghu / ys of the closed-form
%   solution with stoch_simul at the stored parameters and at 20 posterior draws
%   (perturbed mode), and times both. PASS if the largest absolute difference is
%   below 1e-10. Needs Dynare's matlab/ folder on the path (setup_dynare_path).
here = fileparts(mfilename('fullpath'));  addpath(here);
setup_dynare_path;
models = {fullfile(here, '..', 'HLW_new', 'hlw_level'), 'tmp.mat'; fullfile(here, '..', 'UC3'), 'tmp.mat'; ...
          fullfile(here, '..', 'HLW_new', 'hlw_ecm'), 'tmp_ECM_flat.mat'; fullfile(here, '..', 'HLW_new', 'hlw_ecm'), 'tmp_ECM_2step.mat'; ...
          fullfile(here, '..', 'HLW_new', 'hlw_spf'), 'tmp_SPF.mat'};   % workspaces absent locally are skipped
worst = 0;
for m = 1:size(models, 1)
    mdir = models{m, 1};
    if ~exist(fullfile(mdir, models{m, 2}), 'file'), fprintf('skipping %s (no %s)\n', mdir, models{m, 2}); continue; end
    addpath(mdir);
    W = load(fullfile(mdir, models{m, 2}), 'M_', 'oo_', 'options_', 'estim_params_', 'bayestopt_', 'var_list_');
    M_ = W.M_;  oo_ = W.oo_;  options_ = W.options_;  estim_params_ = W.estim_params_;
    M_.verbose = false;  M_.var_list_ = W.var_list_;             % solve_model_draw needs it for stoch_simul
    fprintf('%s: %d endogenous, %d shocks, nfwrd %d\n', M_.fname, M_.endo_nbr, M_.exo_nbr, M_.nfwrd);
    % at the stored parameters, against the stored dr
    [dr, ok, info] = backward_looking_dr(M_, oo_);
    if ~ok, error('test_backward_looking_dr: %s (%s)', info, M_.fname); end
    d_stored = max([max(abs(dr.ghx(:) - oo_.dr.ghx(:))), max(abs(dr.ghu(:) - oo_.dr.ghu(:))), max(abs(dr.ys(:) - oo_.dr.ys(:)))]);
    fprintf('  vs stored dr: max |diff| %.2e\n', d_stored);  worst = max(worst, d_stored);
    % at perturbed draws, against stoch_simul via solve_model_draw
    theta0 = oo_.posterior.optimization.mode(:);  sd = sqrt(diag(oo_.posterior.optimization.Variance));
    try, rng(1); catch, randn('state', 1); end
    t_ss = 0;  t_bw = 0;  n_ok = 0;
    for r = 1:20
        th = theta0 + 0.5 * sd .* randn(numel(theta0), 1);
        M_r = set_all_parameters(th, estim_params_, M_);  M_r.verbose = false;
        t0 = tic;  [info_ss, oo_ss] = solve_model_draw(M_r, options_, oo_);  t_ss = t_ss + toc(t0);
        if info_ss ~= 0, continue; end
        t0 = tic;  [dr_r, ok_r] = backward_looking_dr(M_r, oo_);  t_bw = t_bw + toc(t0);
        if ~ok_r, error('test_backward_looking_dr: closed form failed at a draw (%s)', M_.fname); end
        d = max([max(abs(dr_r.ghx(:) - oo_ss.dr.ghx(:))), max(abs(dr_r.ghu(:) - oo_ss.dr.ghu(:))), max(abs(dr_r.ys(:) - oo_ss.dr.ys(:)))]);
        worst = max(worst, d);  n_ok = n_ok + 1;
    end
    fprintf('  %d draws: stoch_simul %.1f ms, closed form %.2f ms per solve; worst |diff| so far %.2e\n', n_ok, 1000 * t_ss / max(n_ok, 1), 1000 * t_bw / max(n_ok, 1), worst);
    rmpath(mdir);
end
if worst < 1e-10, fprintf('test_backward_looking_dr: PASS (max |diff| %.2e)\n', worst);
else, error('test_backward_looking_dr: FAIL (max |diff| %.2e)', worst); end
