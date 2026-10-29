% RUN_SW_VARIANT  Smets-Wouters (2007) pipeline in its own run folder: estimation, reliability
% draws, model assessment, plots and tables (the UC3 pipeline of main_run_uc3_levels.m).
%
%   Called by main_run_sw_baseline.m (SW_ME = 0) and main_run_sw_me.m (SW_ME = 1). Settings,
%   read from environment variables (unset = default):
%       SW_ME           0 = baseline, 1 = calibrated measurement error          (default 0)
%       SW_MH_REPLIC    Metropolis draws per chain, 2 chains, 20% burn-in       (default 5000)
%       SW_NDRAWS       prior and posterior reliability draws                   (default 500)
%       SW_MA_D         artificial datasets of model_assessment                 (default 50)
%       SW_ME_SHARE     measurement-error variance share of each series (SW_ME = 1; default 0.10)
%       SW_ME_SD        one absolute measurement-error s.d. for all five series (overrides SW_ME_SHARE)
%       SW_RUN_DIR      run folder under Code/SmetsWouters   (default SW_baseline or SW_ME)
%   The defaults give a SHORT, illustrative chain, not a converged posterior.
%
%   Steps (all outputs in the run folder):
%     0. get_sw_data                 usmodel_data.mat, usmodel_mode.mat (downloaded once into
%                                    Code/SmetsWouters, copied into the run folder with the .mod)
%     1. dynare Smets_Wouters_2007_45.mod -DSW_ME=.. -DMH_REPLIC=.. -DN_DRAWS=..
%                                    mode + MH, calib_smoother, prior and posterior reliability
%                                    draws -> tmp.mat, SW_results{,_prior,_posterior}.mat,
%                                    Dynare's output folder Smets_Wouters_2007_45/
%     2. model_assessment(SW_MA_D)   -> model_assessment_lecam.mat
%     3. HLW_plots (main plots only, SW configuration) -> Graphs_level/, Tables_level/
%   Dynare is located by setup_dynare_path (a Dynare already on the path, else DYNARE_ROOT,
%   else the standard macOS install folders).
clear; clc; close all;
dbclear all
tstart = tic;

sw_me     = str2double(getenv('SW_ME'));         if isnan(sw_me),     sw_me     = 0;    end
mh_replic = str2double(getenv('SW_MH_REPLIC'));  if isnan(mh_replic), mh_replic = 5000; end
n_draws   = str2double(getenv('SW_NDRAWS'));     if isnan(n_draws),   n_draws   = 500;  end
run_name  = getenv('SW_RUN_DIR');
if isempty(run_name), if sw_me == 1, run_name = 'SW_ME'; else, run_name = 'SW_baseline'; end, end

code_dir = fileparts(mfilename('fullpath'));            % Code/SmetsWouters
if isempty(code_dir), code_dir = pwd; end
addpath(fullfile(code_dir, '..', 'ReliabilityToolbox'));   % not genpath: octave_compat/ is Octave-only (HLW_plots adds it there)
setup_dynare_path;   % Dynare's matlab/ folder, dynare_config, and under MATLAB its MEX folders
addpath(fullfile(code_dir, '..', 'HLW_new', 'hlw_level'));   % HLW_plots
dynare_version

% ---- 0. Data, run folder ----
run(fullfile(code_dir, 'get_sw_data.m'));
run_dir = fullfile(code_dir, run_name);
if ~exist(run_dir, 'dir'), mkdir(run_dir); end
copyfile(fullfile(code_dir, 'Smets_Wouters_2007_45.mod'), run_dir);
copyfile(fullfile(code_dir, 'usmodel_data.mat'), run_dir);
copyfile(fullfile(code_dir, 'usmodel_mode.mat'), run_dir);
if any(strcmp(strsplit(path, pathsep), code_dir)), rmpath(code_dir); end   % a +Smets_Wouters_2007_45 package of another run must not shadow this one
cd(run_dir);
fprintf('run_sw_variant: SW_ME = %d, mh_replic = %d per chain, %d + %d reliability draws, folder %s\n', sw_me, mh_replic, n_draws, n_draws, run_dir);

% ---- 1. Estimation, smoother, prior/posterior reliability draws ----
eval(sprintf('dynare Smets_Wouters_2007_45.mod -DSW_ME=%d -DMH_REPLIC=%d -DN_DRAWS=%d', sw_me, mh_replic, n_draws));
run_dir = pwd;                                          % dynare cleared the workspace
t_dynare = toc(tstart);
CLS_sw = {'ygap', [0 1 4 20], [],     'y - yf';
          'rrf',  [0 1 4 20], [4 20], ''};

% ---- 2. Posterior-predictive assessment at the posterior mode (reads tmp.mat here) ----
ma_D = str2double(getenv('SW_MA_D'));  if isnan(ma_D), ma_D = 50; end
R = model_assessment(ma_D, struct('model_dir', run_dir, 'workspace', 'tmp.mat', 'CellList_States', {CLS_sw}, ...
                                  'save_file', 'model_assessment_lecam.mat'));
t_assess = toc(tstart) - t_dynare;
% timing printed here: HLW_plots starts with clearvars and wipes these variables
fprintf('\nDynare + draw stages: %.1f h;  model_assessment: %.1f min;  now HLW_plots\n', t_dynare/3600, t_assess/60);

% ---- 3. Plots and tables (main plots only), SW configuration for HLW_plots ----
HLW_RESULTS_FILE   = fullfile(run_dir, 'SW_results.mat');
HLW_OUT_DIR        = run_dir;
DO_MAIN_PLOTS_ONLY = true;
HLW_PLOT_CONFIG = struct();
HLW_PLOT_CONFIG.PlotVar_Names     = {'ygap'; 'rrf'; 'l20rrf'; 'd20rrf'};
HLW_PLOT_CONFIG.display_to_engine = cell(0, 2);            % engine labels are the display names here
HLW_PLOT_CONFIG.time_start        = 1947.5;                % usmodel_data: 1947Q3 - 2004Q4
HLW_PLOT_CONFIG.rstar_overlay     = false;                 % no published r* for this model
HLW_PLOT_CONFIG.common_scale_base = 'rrf';                 % rrf, its 20-quarter smoother and change share the y-axis
HLW_PLOT_CONFIG.sd_scatter_var    = {'ygap', 'y-gap'};
HLW_PLOT_CONFIG.scatter_pairs     = {{'rrf', 'natural rate'}, {'ygap', 'y-gap'}; {'rrf', 'natural rate'}, {'d20rrf', '5-year change in natural rate'}};
HLW_PLOT_CONFIG.rmse_vars         = {{'rrf', 'natural rate'}; {'ygap', 'y-gap'}; {'d20rrf', '5-year change in natural rate'}};
HLW_PLOT_CONFIG.tracking_vars     = {{'rrf', 'natural rate'}; {'ygap', 'y-gap'}};
HLW_PLOT_CONFIG.param_moment_vars = {};                    % reliability-conditional parameter tables: HLW models only
HLW_PLOT_CONFIG.rstar_vs_realrate = {};                    % no observed real rate in the data set
HLW_PLOT_CONFIG.target_params     = {'crpi', 'Taylor rule inflation'; 'crr', 'interest rate smoothing'; 'cry', 'Taylor rule gap'; ...
                                     'crdy', 'Taylor rule gap change'; 'cprobp', 'Calvo prices'; 'cprobw', 'Calvo wages'; ...
                                     'chabb', 'habit'; 'csigma', 'risk aversion'; 'csadjcost', 'investment adj. cost'; ...
                                     'crhoa', 'persistence TFP'; 'crhob', 'persistence risk premium'; 'crhog', 'persistence spending'};
HLW_PLOT_CONFIG.ma_file           = fullfile(run_dir, 'model_assessment_lecam.mat');
HLW_PLOT_CONFIG.fig_prefix        = 'SW_';
HLW_plots
