% MAIN_RUN_UC3_LEVELS_I2  UC3 levels with a growth state g* (I(2) potential output): Dynare, assessment, plots.
%
%   ADDED 2026-09-18. Same three stages as main_run_uc3_levels.m for
%   uc3_levels_i2.mod (the variant behind the paper's UC3 tables, with the g*
%   column). Own file names so both models' results coexist in this folder:
%     tmp_i2.mat, UC3_i2_results{,_prior,_posterior}.mat, model_assessment_lecam_i2.mat,
%     outputs under Code/UC3/i2/{Graphs_level,Tables_level}.
%   Run from Code/UC3.
clear; clc; close all;
dbclear all
tstart = tic;

dynare_dirs = {'/opt/homebrew/opt/dynare/lib/dynare/matlab', '/apps/dynare/dynare-6.0/matlab', 'C:\dynare\6.4\matlab'};
if ~exist('dynare', 'file') && ~(exist('OCTAVE_VERSION', 'builtin') > 0), od = dir('/Applications/Dynare/*/matlab/dynare_config.m'); if ~isempty(od), addpath(od(end).folder); dynare_config; end, end   % official Dynare macOS package (MATLAB MEX), newest version (2026-09-28)
if ~exist('dynare', 'file')
    for k = 1:numel(dynare_dirs)
        if exist(dynare_dirs{k}, 'dir'), addpath(dynare_dirs{k}); break; end
    end
end
if ~exist('dynare', 'file'), error('main_run_uc3_levels_i2: Dynare not found; addpath its matlab/ folder.'); end
if exist('OCTAVE_VERSION', 'builtin') && exist('dynare_config', 'file'), dynare_config; end
dynare_version

uc3_dir = fileparts(mfilename('fullpath'));
if isempty(uc3_dir), uc3_dir = pwd; end
addpath(uc3_dir);
addpath(genpath(fullfile(uc3_dir, '..', 'ReliabilityToolbox')));
if exist('setup_dynare_path', 'file'), setup_dynare_path; end   % under MATLAB adds Dynare's MEX folder if dynare_config skipped it, warns if the compiled Riccati solver is missing (2026-09-28)
addpath(fullfile(uc3_dir, '..', 'HLW_new', 'hlw_level'));

% ---- 1. Estimation, smoother, prior/posterior reliability draws ----
dynare uc3_levels_i2.mod;
t_dynare = toc(tstart);

% dynare clears the workspace: re-establish the folder and the state list
uc3_dir = fileparts(mfilename('fullpath'));
if isempty(uc3_dir), uc3_dir = pwd; end
CLS_i2 = {'ygap',   [0 1 4 20], [];
          'ustar',  [0 1 4 20], [4 20];
          'ystar',  [0 1 4 20], [4 20];
          'pistar', [0 1 4 20], [4 20];
          'gstar',  [0 1 4 20], [4 20]};

% ---- 2. Posterior-predictive assessment at the posterior mode ----
R = model_assessment(300, struct('model_dir', uc3_dir, 'workspace', 'tmp_i2.mat', 'CellList_States', {CLS_i2}, ...
                                 'save_file', 'model_assessment_lecam_i2.mat'));
t_assess = toc(tstart) - t_dynare;
fprintf('\nDynare + draw stages: %.1f h;  model_assessment: %.1f min;  now HLW_plots\n', t_dynare/3600, t_assess/60);

% ---- 3. Plots and tables (main plots only), UC3-i2 configuration ----
HLW_RESULTS_FILE   = fullfile(uc3_dir, 'UC3_i2_results.mat');
HLW_OUT_DIR        = fullfile(uc3_dir, 'i2');
if ~exist(HLW_OUT_DIR, 'dir'), mkdir(HLW_OUT_DIR); end
DO_MAIN_PLOTS_ONLY = true;
HLW_PLOT_CONFIG = struct();
HLW_PLOT_CONFIG.PlotVar_Names     = {'ygap'; 'ustar'; 'l20ustar'; 'ystar'; 'pistar'; 'gstar'; 'd20ustar'};   % paper order c, u*, y*, pi*, g* (+ l20u*, d20u*)
HLW_PLOT_CONFIG.display_to_engine = cell(0, 2);
HLW_PLOT_CONFIG.time_start        = 1960;
HLW_PLOT_CONFIG.rstar_overlay     = false;
HLW_PLOT_CONFIG.common_scale_base = 'ustar';
HLW_PLOT_CONFIG.sd_scatter_var    = {'ygap', 'y-gap'};
HLW_PLOT_CONFIG.scatter_pairs     = {{'ustar', 'u*'}, {'ygap', 'y-gap'}; {'ustar', 'u*'}, {'d20ustar', '5-year change in u*'}};
HLW_PLOT_CONFIG.rmse_vars         = {{'ustar', 'u*'}; {'ygap', 'y-gap'}; {'d20ustar', '5-year change in u*'}};
HLW_PLOT_CONFIG.tracking_vars     = {{'ustar', 'u*'}; {'ygap', 'y-gap'}};
HLW_PLOT_CONFIG.param_moment_vars = {};                    % reliability-conditional parameter tables: HLW models only
HLW_PLOT_CONFIG.target_params     = {'a1', '1st lag coeff'; 'rho', 'sum of lag coeffs'; 'b1', 'PC slope'; 'b2', 'PC lagged slope'; ...
                                     'c1', 'Okun slope'; 'c2', 'Okun lagged slope'; 'se_ygap', 'SE y-gap'; 'se_ystar', 'SE potential y'; ...
                                     'se_gstar', 'SE trend growth'; 'se_ustar', 'SE u*'; 'se_pistar', 'SE trend inflation'; ...
                                     'se_y', 'SE y meas.'; 'se_u', 'SE u meas.'; 'se_pi', 'SE infl. meas.'};
HLW_PLOT_CONFIG.ma_file           = fullfile(uc3_dir, 'model_assessment_lecam_i2.mat');
HLW_PLOT_CONFIG.param_tex_names   = {'a1', '$a_1$'; 'rho', '$\rho$'; 'b1', '$b_1$'; 'b2', '$b_2$'; 'c1', '$c_1$'; 'c2', '$c_2$'; ...
                                     'se_ygap', '$\sigma_c$'; 'se_ystar', '$\sigma_{y^*}$'; 'se_ustar', '$\sigma_{u^*}$'; ...
                                     'se_pistar', '$\sigma_{\pi^*}$'; 'se_gstar', '$\sigma_{g^*}$'; ...
                                     'se_y', '$\sigma_y$'; 'se_u', '$\sigma_u$'; 'se_pi', '$\sigma_\pi$'};   % paper notation (draft.tex, app:uc3)
HLW_PLOT_CONFIG.param_table_label = 'tab:uc3_prior_posterior';
HLW_plots
