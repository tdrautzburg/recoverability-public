@#define USE_LEGACY_DYNARE_SYNTAX = 0

@#define DO_ESTIMATE = 1

// UC3 in levels with a RANDOM-WALK GROWTH STATE g* in potential output (the I(2)
// variant the paper's UC3 tables were produced with: the g* column, the
// sigma_{g*} ~ IG(0.1, 0.4) row of the parameter table, and "nu ... for I(2)
// processes like y*"). ADDED 2026-09-18 alongside uc3_levels.mod, which has
// y* as a pure random walk. Differences to uc3_levels.mod:
//     ystar = ystar(-1) + gstar(-1) + se_ystar*e_ystar;   gstar = gstar(-1) + se_gstar*e_gstar;
// (gstar in quarterly log points x 100, no /4 - unlike HLW's annualised gg),
// gstar added to the diagnosed states, and separate file names (tmp_i2.mat,
// UC3_i2_results.mat) so both models' results coexist. Pipeline otherwise as in
// uc3_levels.mod; run from this folder via main_run_uc3_levels_i2.m.

set_dynare_seed(1);
//////////////////////////////////////////////////////////////////////////
%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------
close all;

%----------------------------------------------------------------
% 1. Variable Definition
%----------------------------------------------------------------
var y ystar gstar ygap pi pistar u ustar;
var yobs piobs uobs;

varexo e_pi e_y e_u;
varexo e_ystar e_gstar e_ygap e_pistar e_ustar;

parameters a1 rho b1 b2 c1 c2  se_ygap se_pi se_y se_u se_ystar se_gstar se_pistar se_ustar;

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
se_gstar=	0.05	;
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
ystar = ystar(-1) + gstar(-1) + se_ystar*e_ystar; %ystar in 100*log(rgdp) units; gstar is its quarterly trend growth (I(2) level)
gstar = gstar(-1) + se_gstar*e_gstar;
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
gstar=0;
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
var e_gstar; stderr 1;
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
%
estimated_params;
%
    a1,               1.5,       0,      2,        NORMAL_PDF,      1.5,    0.5;
    rho,              0.5,       0,      0.99,       BETA_PDF,      0.5,    0.15;
    b1,               0.1,       0,      2,        NORMAL_PDF,      0.1,    0.5;
    b2,               0.05,     -1,      1,        NORMAL_PDF,      0.05,    0.5;
    c1,               -0.5,     -1,      1,        NORMAL_PDF,      -0.5,    0.5;
    c2,               -0.25,     -1,      1,        NORMAL_PDF,      -0.25,    0.5;
    se_ygap      ,    0.4,     1e-4,     10,       INV_GAMMA_PDF,   0.5,    .4;
    se_ystar,         0.75,    1e-4,     10,       INV_GAMMA_PDF,   0.5,    .4;
    se_ustar   ,      0.1,     1e-4,     10,       INV_GAMMA_PDF,   0.1,    .4;
    se_pistar,        0.75,     1e-4,     10,       INV_GAMMA_PDF,   0.5,    .4;
    se_gstar,         0.05,     1e-4,     10,       INV_GAMMA_PDF,   0.1,    .4;   % paper's table: sigma_{g*} ~ IG(0.1, 0.4)
    se_y,             0.1,     1e-4,     10,       INV_GAMMA_PDF,   0.5,    2;
    se_u,             0.1,     1e-4,     10,       INV_GAMMA_PDF,   0.2,    2;
    se_pi,            0.25,     1e-4,     10,       INV_GAMMA_PDF,   0.25,    1;

    end;

    % Observables
    varobs yobs piobs uobs;

if exist('ReliabilityToolbox', 'dir')
    addpath(genpath('ReliabilityToolbox'));
else
    addpath(genpath('../ReliabilityToolbox'));
end

    % Mode: global search first, local refinement from that mode, then MH 50000 x 2
    estimation(nograph,order=1,optim=('MaxIter',100000,'Tolfun',1.0e-06), datafile='UC3_data_yobs_1960_2019',mode_compute=2,cova_compute=1,mh_replic=0,first_obs=1,lik_init=3,prefilter=0, mode_check, plot_priors=0);

    estimation(graph,order=1,optim=('MaxIter',100000,'Tolfun',1.0e-06), datafile='UC3_data_yobs_1960_2019', mode_file='uc3_levels_i2/Output/uc3_levels_i2_mode.mat', mode_compute=4,mh_replic=50000,mh_nblocks=2,mh_jscale=0.45,mh_drop=0.2,first_obs=1,lik_init=3,prefilter=0,filtered_vars, smoothed_state_uncertainty, mode_check, plot_priors=0);

    calib_smoother(datafile='UC3_data_yobs_1960_2019',first_obs=1,kalman_algo =2,prefilter=0,filtered_vars);

@#endif

original_estim_params_=estim_params_;

stoch_simul(nograph);

save tmp_i2.mat

M_.observable_names = dataset_.name ;
M_.data = dataset_info.rawdata(~any(isnan(dataset_info.rawdata')),: );
M_.CellList_States    = {'ygap',   [0 1 4 20], [];
                         'ustar',  [0 1 4 20], [4 20];
                         'ystar',  [0 1 4 20], [4 20];
                         'pistar', [0 1 4 20], [4 20];
                         'gstar',  [0 1 4 20], [4 20]};
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

% Staged, size-safe saves (draws chunked into UC3_i2_results_prior/_posterior.mat)
save_results_split('UC3_i2_results.mat');

tic;
[oo_.prior_function_results,     prior_draws,     prior_info]     = run_reliability_draws('prior',     10000, M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info);
toc;
save_results_split('UC3_i2_results.mat');
tic;
[oo_.posterior_function_results, posterior_draws, posterior_info] = run_reliability_draws('posterior', 10000, M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info);
toc;

[LatentStateDiagnostics, ModelMatrices, obs_indices, key_state_indices, TV_state_estimate, TV_surprises, TV_log_llhood]= compute_state_reliability_diagnostics(M_.observable_names, M_.CellList_States, oo_, M_, M_.data);

save_results_split('UC3_i2_results.mat');

% Dynare-vs-toolbox comparison, one panel per diagnosed state (see uc3_levels.mod)
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
        case 'base', ref = [upd, smo];                          leg = {'Dynare (t|t)', 'Dynare (t|T)', 'TV filter (t|t)'};
        case 'l',    ref = [shift(upd, k), shift(smo, k)];      leg = {sprintf('Dynare filter (t-%d|t-%d)', k, k), sprintf('Dynare smoother (t-%d|T)', k), sprintf('TV fixed-lag smoother (t-%d|t)', k)};
        case 'd',    ref = [upd - shift(upd, k), smo - shift(smo, k)]; leg = {sprintf('Dynare filter: x(t|t)-x(t-%d|t-%d)', k, k), sprintf('Dynare smoother: x(t|T)-x(t-%d|T)', k), 'TV filter (t|t)'};
    end
    t0 = min(T, t_proper + k);
    subplot(nr, nc, vv);
    plot(t0:T, [ref(t0:T, :), TV_state_estimate(t0:T, vv)], 'LineWidth', 1.2);
    xlim([t0, T]);
    legend(leg, 'location', 'northwest', 'FontSize', 5);
    title(lab, 'Interpreter', 'none');
    grid on;
end
print('-dpng', 'Dynare_Comparison_uc3_levels_i2.png');
