% MAIN_RUN_UC3_LEVELS  UC3 (levels) pipeline: Dynare run, model assessment, plots and tables.
%
%   ADDED 2026-09-17. Run from THIS folder (Code/UC3): Dynare writes +uc3_levels/
%   and uc3_levels/ here and the chain files are resolved relative to it.
%     1. dynare uc3_levels.mod      estimation (mh_replic = 50000), calib_smoother,
%                                   10000 prior + 10000 posterior reliability draws
%                                   -> tmp.mat, UC3_results.mat (staged saves)
%     2. model_assessment(300, ...) damped-Newton posterior-predictive assessment on
%                                   tmp.mat -> model_assessment_lecam.mat (here)
%     3. HLW_plots                  with DO_MAIN_PLOTS_ONLY = true and the UC3 plot
%                                   configuration -> Graphs_level/, Tables_level/ (here);
%                                   ends with model_assessment_plots on step 2's file
%   The generic code lives in ReliabilityToolbox (model_assessment,
%   model_assessment_plots) and HLW_new/hlw_level (HLW_plots); nothing UC3-specific
%   is in them beyond the options passed below.
%
%   Run length (environment variables; unset = the full run):
%       UC3_MH_REPLIC   Metropolis draws per chain, 2 chains   (default 50000)
%       UC3_NDRAWS      prior and posterior reliability draws  (default 10000)
%       UC3_MA_D        artificial datasets of model_assessment (default 300)
%   Dynare is located by setup_dynare_path (a Dynare already on the path, else
%   DYNARE_ROOT, else the standard macOS install folders).
clear; clc; close all;
dbclear all
tstart = tic;

mh_replic = str2double(getenv('UC3_MH_REPLIC'));  if isnan(mh_replic), mh_replic = 50000; end
n_draws   = str2double(getenv('UC3_NDRAWS'));     if isnan(n_draws),   n_draws   = 10000; end

% ---- code folders: this one, ReliabilityToolbox, and hlw_level (plots); Dynare ----
uc3_dir = fileparts(mfilename('fullpath'));
if isempty(uc3_dir), uc3_dir = pwd; end
addpath(uc3_dir);
addpath(fullfile(uc3_dir, '..', 'ReliabilityToolbox'));   % not genpath: octave_compat/ is Octave-only (HLW_plots adds it there)
setup_dynare_path;   % Dynare's matlab/ folder, dynare_config, and under MATLAB its MEX folders
addpath(fullfile(uc3_dir, '..', 'HLW_new', 'hlw_level'));
dynare_version

% ---- 1. Estimation, smoother, prior/posterior reliability draws ----
fprintf('main_run_uc3_levels: mh_replic = %d per chain, %d + %d reliability draws\n', mh_replic, n_draws, n_draws);
eval(sprintf('dynare uc3_levels.mod -DMH_REPLIC=%d -DN_DRAWS=%d', mh_replic, n_draws));
t_dynare = toc(tstart);

% dynare clears the workspace: re-establish the folder and the state list
uc3_dir = fileparts(mfilename('fullpath'));
if isempty(uc3_dir), uc3_dir = pwd; end
CLS_uc3 = {'ygap',   [0 1 4 20], [];
           'ustar',  [0 1 4 20], [4 20];
           'ystar',  [0 1 4 20], [4 20];
           'pistar', [0 1 4 20], [4 20]};

% ---- 2. Posterior-predictive assessment at the posterior mode (reads tmp.mat here) ----
ma_D = str2double(getenv('UC3_MA_D'));  if isnan(ma_D), ma_D = 300; end
R = model_assessment(ma_D, struct('model_dir', uc3_dir, 'workspace', 'tmp.mat', 'CellList_States', {CLS_uc3}, 'solve_method', getenv('UC3_SOLVER')));   % UC3_SOLVER=backward: closed-form decision rule (2026-09-28)
t_assess = toc(tstart) - t_dynare;
% timing printed here: HLW_plots starts with clearvars and wipes these variables
fprintf('\nDynare + draw stages: %.1f h;  model_assessment: %.1f min;  now HLW_plots\n', t_dynare/3600, t_assess/60);

% ---- 3. Plots and tables (main plots only), UC3 configuration for HLW_plots ----
HLW_RESULTS_FILE   = fullfile(uc3_dir, 'UC3_results.mat');
HLW_OUT_DIR        = uc3_dir;
DO_MAIN_PLOTS_ONLY = true;
HLW_PLOT_CONFIG = struct();
HLW_PLOT_CONFIG.PlotVar_Names     = {'ygap'; 'ustar'; 'l20ustar'; 'ystar'; 'pistar'; 'd20ustar'};
HLW_PLOT_CONFIG.display_to_engine = cell(0, 2);            % engine labels are the display names here
HLW_PLOT_CONFIG.time_start        = 1960;                  % UC3_data_yobs_1960_2019: 1960q1 - 2019q4
HLW_PLOT_CONFIG.rstar_overlay     = false;                 % no published r* in this model
HLW_PLOT_CONFIG.common_scale_base = 'ustar';               % u* and its 20-quarter smoother share the y-axis
HLW_PLOT_CONFIG.sd_scatter_var    = {'ygap', 'y-gap'};
HLW_PLOT_CONFIG.scatter_pairs     = {{'ustar', 'u*'}, {'ygap', 'y-gap'}; {'ustar', 'u*'}, {'d20ustar', '5-year change in u*'}};
HLW_PLOT_CONFIG.rmse_vars         = {{'ustar', 'u*'}; {'ygap', 'y-gap'}; {'d20ustar', '5-year change in u*'}};
HLW_PLOT_CONFIG.tracking_vars     = {{'ustar', 'u*'}; {'ygap', 'y-gap'}};
HLW_PLOT_CONFIG.param_moment_vars = {};                    % reliability-conditional parameter tables: HLW models only
HLW_PLOT_CONFIG.target_params     = {'a1', '1st lag coeff'; 'rho', 'sum of lag coeffs'; 'b1', 'PC slope'; 'b2', 'PC lagged slope'; ...
                                     'c1', 'Okun slope'; 'c2', 'Okun lagged slope'; 'se_ygap', 'SE y-gap'; 'se_ystar', 'SE potential y'; ...
                                     'se_ustar', 'SE u*'; 'se_pistar', 'SE trend inflation'; 'se_y', 'SE y meas.'; 'se_u', 'SE u meas.'; 'se_pi', 'SE infl. meas.'};
HLW_PLOT_CONFIG.ma_file           = fullfile(uc3_dir, 'model_assessment_lecam.mat');
HLW_plots
