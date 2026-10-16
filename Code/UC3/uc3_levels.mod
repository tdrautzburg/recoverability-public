@#define USE_LEGACY_DYNARE_SYNTAX = 0

@#define DO_ESTIMATE = 1

@#ifndef MH_REPLIC
@#define MH_REPLIC = 50000
@#endif
@#ifndef N_DRAWS
@#define N_DRAWS = 10000
@#endif
// MH_REPLIC: Metropolis-Hastings draws per chain (2 chains, 20% burn-in);
// N_DRAWS: prior and posterior reliability draws. The defaults are the settings
// of the full run; main_run_uc3_levels.m passes smaller values from the
// environment variables UC3_MH_REPLIC / UC3_NDRAWS, e.g.
// dynare uc3_levels.mod -DMH_REPLIC=2000 -DN_DRAWS=200

// UC3 in levels: output, inflation and unemployment with random-walk trends and
// a common AR(2) gap. UPDATED 2026-09-17 to the current toolbox pipeline (same
// structure as HLW_new/hlw_level/HLW_2023_US_no_lags.mod): estimation with
// mh_replic = 50000, calib_smoother, toolbox bookkeeping, M_.CellList_States,
// run_reliability_draws (10000 prior + 10000 posterior draws) with staged saves,
// direct engine call and the Dynare-vs-toolbox comparison figure. The lag and
// change states (l<m>NAME, d<h>NAME) are built inside the engine. Run from this
// folder: main_run_uc3_levels.m does that and then the model assessment and plots.

set_dynare_seed(1);
//////////////////////////////////////////////////////////////////////////
%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------
close all;

%----------------------------------------------------------------
% 1. Variable Definition
%----------------------------------------------------------------
var y ystar ygap pi pistar u ustar;
var yobs piobs uobs;

varexo e_pi e_y e_u;
varexo e_ystar e_ygap e_pistar e_ustar;

parameters a1 rho b1 b2 c1 c2  se_ygap se_pi se_y se_u se_ystar se_pistar se_ustar;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
a1	=	1	;
rho     =       0.5	; % this is a1+a2
b1	=	0.1	;
b2	=	-0.1;
c1	=	-0.3	;
c2	=	-0.1;
%
se_ygap	=	0.4	;
se_pi	=	0.25	;
se_y	=	0.25	;
se_u	=	0.2	;
se_ystar=	0.5	;
se_ustar=	0.1	;
se_pistar=	0.75	;
%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model(linear);
% Structural equations (transition equations for ygap and pi)
#a2=rho-a1;
ygap = a1*ygap(-1) + a2*ygap(-2) + se_ygap*e_ygap; %(IS equation)
%
%Transition equations for other states
ystar = ystar(-1) + se_ystar*e_ystar; %ystar is in 100*log(rgdp) units, hence its deltas are in quarterly
pistar = pistar(-1) + se_pistar*e_pistar;
ustar = ustar(-1) + se_ustar*e_ustar;
%
% Measurement equations, linking states to observables
y = ystar + ygap + se_y*e_y;
pi = pistar + b1*ygap + b2*ygap(-1) + se_pi*e_pi; % (Phillips curve)
u = ustar + c1*ygap + c2*ygap(-1)  + se_u*e_u;
yobs = y;
piobs = pi;
uobs = u;

end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------
% Steady state
steady_state_model;
pi=0;
pistar=0;
u=0;
ustar=0;
ygap=0;
ystar=0;
y=0;
yobs=y;
piobs=pi;
uobs=u;
end;
%
resid;
steady;
check;
model_diagnostics;
%
shocks; % all have unit se's
var e_pi; stderr 1;
var e_y; stderr 1;
var e_u; stderr 1;
var e_ustar; stderr 1;
var e_ystar; stderr 1;
var e_ygap ; stderr 1;
var e_pistar; stderr 1;
end;

%
write_latex_dynamic_model;
write_latex_static_model;
%
@#if DO_ESTIMATE==1
    %----------------------------------------------------------------
    % 5. Estimation
    %----------------------------------------------------------------
    // PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
    // PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF, UNIFORM_PDF
%
estimated_params;
%
    a1,               1.5,       0,      2,        NORMAL_PDF,      1.5,    0.5;
    rho,              0.5,       0,      0.99,       BETA_PDF,      0.5,    0.15;   % UB 0.99 (was 0.9: the mode sat on it, Hessian singular, MH could not start - run of 2026-09-17)
    b1,               0.1,       0,      2,        NORMAL_PDF,      0.1,    0.5;
    b2,               0.05,     -1,      1,        NORMAL_PDF,      0.05,    0.5;
    c1,               -0.5,     -1,      1,        NORMAL_PDF,      -0.5,    0.5;
    c2,               -0.25,     -1,      1,        NORMAL_PDF,      -0.25,    0.5;
    se_ygap      ,    0.4,     1e-4,     10,       INV_GAMMA_PDF,   0.5,    .4;
    se_ystar,         0.75,    1e-4,     10,       INV_GAMMA_PDF,   0.5,    .4;
    se_ustar   ,      0.1,     1e-4,     10,       INV_GAMMA_PDF,   0.1,    .4;
    se_pistar,        0.75,     1e-4,     10,       INV_GAMMA_PDF,   0.5,    .4;
    se_y,             0.1,     1e-4,     10,       INV_GAMMA_PDF,   0.5,    2;
    se_u,             0.1,     1e-4,     10,       INV_GAMMA_PDF,   0.2,    2;
    se_pi,            0.25,     1e-4,     10,       INV_GAMMA_PDF,   0.25,    1;

    end;

    % Observables
    varobs yobs piobs uobs;

% ReliabilityToolbox: engine, reliability_draw_wrapper, run_reliability_draws,
% expand_states_of_interest, augment_ssm. Repo layout Code/UC3 -> Code/ReliabilityToolbox;
% a copy next to this file (flat HPC layout) takes precedence.
if exist('ReliabilityToolbox', 'dir')
    addpath('ReliabilityToolbox');
else
    addpath('../ReliabilityToolbox');   % not genpath: octave_compat/ is Octave-only
end

    % Mode: global search (simulated annealing) first, then local refinement from
    % that mode, then MH: 50000 draws x 2 chains, 20% burn-in (as in the HLW .mod).
    estimation(nograph,order=1,optim=('MaxIter',100000,'Tolfun',1.0e-06), datafile='UC3_data_yobs_1960_2019',mode_compute=2,cova_compute=1,mh_replic=0,first_obs=1,lik_init=3,prefilter=0, mode_check, plot_priors=0);

    estimation(graph,order=1,optim=('MaxIter',100000,'Tolfun',1.0e-06), datafile='UC3_data_yobs_1960_2019', mode_file='uc3_levels/Output/uc3_levels_mode.mat', mode_compute=4,mh_replic=@{MH_REPLIC},mh_nblocks=2,mh_jscale=0.45,mh_drop=0.2,first_obs=1,lik_init=3,prefilter=0,filtered_vars, smoothed_state_uncertainty, mode_check, plot_priors=0);

    % Smoother at the posterior mean (no heteroskedastic_shocks block in this
    % model, so no heteroskedastic_filter option); filtered_vars gives
    % oo_.UpdatedVariables for the comparison figure below.
    calib_smoother(datafile='UC3_data_yobs_1960_2019',first_obs=1,kalman_algo =2,prefilter=0,filtered_vars);

@#endif

original_estim_params_=estim_params_;

stoch_simul(nograph);

save tmp.mat

M_.observable_names = dataset_.name ;
M_.data = dataset_info.rawdata(~any(isnan(dataset_info.rawdata')),: );
% States to diagnose: {variable, fixed-lag smoother lags m (0 = filter), change horizons h}.
% The three trends are random walks; their 4- and 20-quarter changes are the
% analogues of d4rrstar / d20rrstar in HLW.
M_.CellList_States    = {'ygap',   [0 1 4 20], [];
                         'ustar',  [0 1 4 20], [4 20];
                         'ystar',  [0 1 4 20], [4 20];
                         'pistar', [0 1 4 20], [4 20]};
M_.states_of_interest = expand_states_of_interest(M_.CellList_States);
M_.var_list_ = var_list_;
M_.verbose = false;
nstates = length(oo_.dr.state_var);
M_.X_InitialDistribution.X0 = zeros(nstates,1);
M_.X_InitialDistribution.B0 = eye(nstates)*100;
oo_.bayestopt_smoother_var_list=bayestopt_.smoother_var_list;
oo_.bayestopt_smoother_restrict_columns=bayestopt_.smoother_restrict_columns;
oo_.bayestopt_mf=bayestopt_.mf;
oo_.dr.qz_crit = 1+1e-6;

% STAGED SAVES: (1) after Dynare, (2) after the prior stage, (3) final (see the
% HLW .mod for the rationale; a restart can load_results_split('UC3_results.mat')
% and run the missing stage). save_results_split writes UC3_results.mat plus
% UC3_results_prior.mat / UC3_results_posterior.mat with the draws in chunks
% below 1 GB: a monolithic oo_ was 4.7 GB here and unreadable (MAT-5 limit).
% HLW_plots is pointed at UC3_results.mat (HLW_RESULTS_FILE) and reattaches them.
save_results_split('UC3_results.mat');

tic;
[oo_.prior_function_results,     prior_draws,     prior_info]     = run_reliability_draws('prior',     @{N_DRAWS}, M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info);
toc;
save_results_split('UC3_results.mat');    % staged save (2)
tic;
[oo_.posterior_function_results, posterior_draws, posterior_info] = run_reliability_draws('posterior', @{N_DRAWS}, M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info);
toc;

% Direct (single) run of the diagnostics engine at the posterior mean
[LatentStateDiagnostics, ModelMatrices, obs_indices, key_state_indices, TV_state_estimate, TV_surprises, TV_log_llhood]= compute_state_reliability_diagnostics(M_.observable_names, M_.CellList_States, oo_, M_, M_.data);

save_results_split('UC3_results.mat');    % staged save (3): final

% Comparison of the toolbox's diffuse filter with Dynare, one panel per state of
% the list (same block as in the HLW .mod): base states against Dynare's x_{t|t}
% and x_{t|T}; l<m>NAME against Dynare's filter shifted by m and its smoother;
% d<h>NAME against the same differences of Dynare's paths. Panels start after
% the diffuse period plus the lag/horizon of the state.
figure(100)
set(gcf, 'Position', [100, 100, 1600, 1000]);
nS = length(M_.states_of_interest);
nr = ceil(sqrt(nS)); nc = ceil(nS / nr);
T  = size(TV_state_estimate, 1);
if isfield(LatentStateDiagnostics, 'TV_KF_Sequence') && isfield(LatentStateDiagnostics.TV_KF_Sequence, 't_proper')
    t_proper = LatentStateDiagnostics.TV_KF_Sequence.t_proper;
else
    t_proper = 1;
end
shift = @(x, k) [nan(k, 1); x(1:end-k)];
for vv = 1:nS
    lab = M_.states_of_interest{vv};
    tok = regexp(lab, '^([ld])(\d+)([A-Za-z_].*)$', 'tokens', 'once');
    if isempty(tok), kind = 'base'; k = 0; base = lab; else, kind = tok{1}; k = str2double(tok{2}); base = tok{3}; end
    if ~isfield(oo_.UpdatedVariables, base), continue; end
    ss  = oo_.dr.ys(strcmp(M_.endo_names, base));
    upd = oo_.UpdatedVariables.(base)  - ss;
    smo = oo_.SmoothedVariables.(base) - ss;
    switch kind
        case 'base'
            ref = [upd, smo];
            leg = {'Dynare (t|t)', 'Dynare (t|T)', 'TV filter (t|t)'};
        case 'l'
            ref = [shift(upd, k), shift(smo, k)];
            leg = {sprintf('Dynare filter (t-%d|t-%d)', k, k), sprintf('Dynare smoother (t-%d|T)', k), sprintf('TV fixed-lag smoother (t-%d|t)', k)};
        case 'd'
            ref = [upd - shift(upd, k), smo - shift(smo, k)];
            leg = {sprintf('Dynare filter: x(t|t)-x(t-%d|t-%d)', k, k), sprintf('Dynare smoother: x(t|T)-x(t-%d|T)', k), 'TV filter (t|t)'};
    end
    t0 = min(T, t_proper + k);
    subplot(nr, nc, vv);
    plot(t0:T, [ref(t0:T, :), TV_state_estimate(t0:T, vv)], 'LineWidth', 1.2);
    xlim([t0, T]);
    legend(leg, 'location', 'northwest', 'FontSize', 5);
    title(lab, 'Interpreter', 'none');
    grid on;
end
print('-dpng', 'Dynare_Comparison_uc3_levels.png');
