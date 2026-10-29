%%
% Workspace variables a CALLING script may set before `HLW_plots` (they survive
% the clearvars below; everything else is cleared):
%   HLW_RESULTS_FILE   results .mat to plot (default HLW_results_no_lags.mat in the code folder)
%   HLW_OUT_DIR        where Graphs_level/ and Tables_level/ are written (default: code folder)
%   DO_MAIN_PLOTS_ONLY true: skip the FEV, scatter,
%                      RMSE-scatter and correlation blocks (default false)
%   DO_TABLES_ONLY     true: no figures at all, only the TeX tables (implies
%                      DO_MAIN_PLOTS_ONLY; ~20 s) (ADDED 2026-09-18)
%   BR_directory       folder with the Bauer-Rudebusch posterior summary (default Code/BauerRudebusch)
%   HLW_code_dir       folder of this script and the .mod (default: located via mfilename)
%   HLW_PLOT_CONFIG    struct overriding the MODEL-SPECIFIC choices (states shown, sample
%                      start, parameter list, state pairs of the scatter plots, ...); any
%                      field left out keeps its HLW default - see "cfg" below. This is
%                      what makes the script usable for other models (UC3: main_run_uc3_levels.m).
% Environment variables of the same names (HLW_RESULTS_FILE, HLW_OUT_DIR) are
% the fallback when the variable is not set, so a shell loop works too.
% NOTE the clearvars: every other variable of the caller is wiped on each call,
% including loop bookkeeping. Loop over results files with setenv, which
% survives, and derive the output folder from the loop variable alone:
%     for f = {'/a/HLW_results_x.mat', '/b/HLW_results_y.mat'}
%         setenv('HLW_RESULTS_FILE', f{1});  setenv('HLW_OUT_DIR', fileparts(f{1}));
%         DO_MAIN_PLOTS_ONLY = true;  HLW_plots
%     end
% (DO_MAIN_PLOTS_ONLY, BR_directory and HLW_code_dir survive as variables.)
clearvars -except HLW_code_dir HLW_OUT_DIR HLW_RESULTS_FILE DO_MAIN_PLOTS_ONLY DO_TABLES_ONLY BR_directory HLW_PLOT_CONFIG
close all

% ------------------------------------------------------------------------
% CONFIGURATION (CHANGED 2026-09-17: was pinned to the HPC scratch folder,
% the retired HLW_results workspace and the old assessment's
% DeltaR2_RMSE_plots.mat; now runs from this folder under MATLAB or Octave).
%   HLW_code_dir    this folder (Code/HLW_new/hlw_level): code, .mod outputs,
%                   model_assessment_lecam.mat, default results file
%   HLW_directory   OUTPUT root: Graphs_level/ and Tables_level/ go under it
%                   (every print/fopen below uses [HLW_directory '/Graphs_level'])
%   results_file    HLW_results_no_lags.mat unless HLW_RESULTS_FILE says otherwise
%   BR_directory    Bauer-Rudebusch posterior files; the BR comparison panels
%                   are skipped when the files are not there (HAVE_BR)
%   model_assessment_lecam.mat, if present, feeds model_assessment_plots at
%   the end (replaces the old DeltaR2 / RMSE fan-chart block).
% The two roles used to share the name HLW_directory, which broke repeated
% calls: the script reassigns HLW_directory to the output root, so a second
% call found the code folder pointing at the previous output folder.
% ------------------------------------------------------------------------

if ~exist('DO_MAIN_PLOTS_ONLY', 'var') || isempty(DO_MAIN_PLOTS_ONLY)
  DO_MAIN_PLOTS_ONLY = false;
end
if ~exist('DO_TABLES_ONLY', 'var') || isempty(DO_TABLES_ONLY)
  DO_TABLES_ONLY = false;
end
if DO_TABLES_ONLY, DO_MAIN_PLOTS_ONLY = true; end

% ------------------------------------------------------------------------
% MODEL-SPECIFIC CHOICES (ADDED 2026-09-17), HLW defaults; override any field
% through HLW_PLOT_CONFIG. Everything below indexes results through these, so
% the script has no hard-wired HLW state or parameter names any more.
%   PlotVar_Names      display names of the states shown / tabulated (paper convention)
%   display_to_engine  {display, engine} pairs where the engine label differs
%   time_start         decimal year of the first observation (quarterly data)
%   rstar_overlay      overlay the published HLW r* (HLW_Rstar_up_to_205q3.m) on the
%                      common-scale panel of common_scale_base
%   common_scale_base  base state whose panels share one y-axis; its l20<base> panel
%                      overlays the filter on the fixed-lag smoother
%   sd_scatter_var     {name, label}: state of the R^2-vs-s.d. and R^2-vs-half-life scatters
%   scatter_pairs      {name,label} x 2 per row: cross-state R^2 / half-life scatters
%   rmse_vars          {name,label} rows: R^2-vs-RMSE scatters and parameter-correlation tables
%   tracking_vars      {name,label} rows: filtered paths conditional on R^2 masks
%   param_moment_vars  {name,label} rows: parameter moments conditional on R^2 > 0.5
%   target_params      {param, label} rows for the correlation tables (missing ones dropped)
%   ma_file            model_assessment_lecam.mat to hand to model_assessment_plots
%                      ('' = look in the output folder, then in HLW_code_dir)
%   param_tex_names    {Dynare name, LaTeX name} rows for table_prior_posterior_params.tex
%                      (unlisted parameters keep their Dynare name) (ADDED 2026-09-18)
%   param_table_label  LaTeX label of that table ('tab:prior_posterior'); give each
%                      model its own so several tables can be \input in one paper
%   rstar_vs_realrate  {state, observable} for the filtered-r*-vs-real-rate figures
%                      (LW_rstar_vs_realrate{,_simple,_1990_bands}.png); skipped when
%                      the model lacks either name; {} disables (ADDED 2026-09-19)
%   late_start         decimal year at which the LW_KF_<state>_<year>.png panels start
%                      (1965 for HLW: after the diffuse-start band has died down); the
%                      file suffix is the year (ADDED 2026-09-27)
%   band_skip          number of initial quarters whose 95% bands are NOT drawn in the
%                      filtered-state figures (0 for HLW: the diffuse phase leaves them
%                      NaN anyway). A model with a PROPER wide initialisation (BR: s.d.
%                      100 pp at t = 0) has a finite but enormous first-quarter band
%                      that would set the y-scale; 4 hides it. Tables and the summary
%                      file are not affected (ADDED 2026-09-27)
%   common_scale_ylim  outer y-limits of the common-scale panels ([-10 10] for r*;
%                      the limits fitted to the lines are clipped to this) (ADDED 2026-09-27)
%   lag_presample_mask true: in the 20-quarter-smoother panels, do not draw the first
%                      20 smoother values (they estimate PRE-SAMPLE states, i.e. the
%                      prior); false (HLW default) keeps them (ADDED 2026-09-27)
%   fig_prefix         prefix of every draw-by-draw figure file ('LW_KF_' for the HLW
%                      line and UC3, as always; 'BR_KF_' for Bauer-Rudebusch). The
%                      MA_* assessment figures, the tables and LW_rstar_vs_realrate*
%                      keep their names (ADDED 2026-09-27)
% ------------------------------------------------------------------------
cfg = struct();
cfg.PlotVar_Names     = {'ygap';'rrstar';'l20rrstar';'ldrrstar'; 'drrstar'; 'd20rrstar'};
cfg.display_to_engine = {'ldrrstar', 'l1drrstar'};
cfg.time_start        = 1960;
cfg.rstar_overlay     = true;
cfg.common_scale_base = 'rrstar';
cfg.sd_scatter_var    = {'ygap', 'y-gap'};
cfg.scatter_pairs     = {{'rrstar', 'r*'}, {'ygap', 'y-gap'}; {'rrstar', 'r*'}, {'d20rrstar', '5-year change'}};
cfg.rmse_vars         = {{'rrstar', 'r*'}; {'ygap', 'y-gap'}; {'d20rrstar', '5-year change'}};
cfg.tracking_vars     = {{'rrstar', 'r*'}; {'ygap', 'y-gap'}};
cfg.param_moment_vars = {{'rrstar', 'r*'}; {'ygap', 'y-gap'}};
cfg.target_params     = {'rho', 'sum of lag coeffs'; 'a1', '1st lag coeffs'; 'a_r', 'r-y'; 'bpi', 'lagged inflation'; ...
                         'by', 'slope PC'; 'phi', 'COVID-loading'; 'se_ygap', 'SE y-gap'; 'se_pi', 'SE infl'; ...
                         'se_ystar', 'SE potential y'; 'se_zz', 'SE z'; 'se_gg', 'SE growth'; 'se_rr', 'SE rr'};
cfg.ma_file           = '';
cfg.param_tex_names   = cell(0, 2);
cfg.param_table_label = 'tab:prior_posterior';
cfg.rstar_vs_realrate = {'rrstar', 'rrobs'};
cfg.late_start        = 1965;
cfg.band_skip         = 0;
cfg.common_scale_ylim = [-10, 10];
cfg.lag_presample_mask = false;
cfg.fig_prefix        = 'LW_KF_';
if exist('HLW_PLOT_CONFIG', 'var') && isstruct(HLW_PLOT_CONFIG)
    fn_cfg = fieldnames(HLW_PLOT_CONFIG);
    for k_cfg = 1:numel(fn_cfg)
        if ~isfield(cfg, fn_cfg{k_cfg}), error('HLW_plots: unknown HLW_PLOT_CONFIG field "%s".', fn_cfg{k_cfg}); end
        cfg.(fn_cfg{k_cfg}) = HLW_PLOT_CONFIG.(fn_cfg{k_cfg});
    end
end

if ~exist('HLW_code_dir', 'var') || isempty(HLW_code_dir)
  HLW_code_dir = fileparts(mfilename('fullpath'));
  if isempty(HLW_code_dir), HLW_code_dir = pwd; end
end
% path setup is idempotent and cheap, so it runs on every call
addpath(HLW_code_dir);
addpath(fullfile(HLW_code_dir, '..'));                          % HLW_Rstar_up_to_205q3.m
if exist(fullfile(HLW_code_dir, 'ReliabilityToolbox'), 'dir')
    addpath(genpath(fullfile(HLW_code_dir, 'ReliabilityToolbox')));
else
    addpath(fullfile(HLW_code_dir, '..', '..', 'ReliabilityToolbox'));
end
IS_OCTAVE = exist('OCTAVE_VERSION', 'builtin') > 0;
if IS_OCTAVE
    pkg load statistics                                          % prctile, norminv, corr
    addpath(fullfile(fileparts(which('compute_state_reliability_diagnostics')), 'octave_compat'));   % contains, histogram
    warning('off', 'Octave:legend:unimplemented-location');
end

% inputs/outputs: workspace variable first, environment variable second, default third
if exist('HLW_OUT_DIR', 'var') && ~isempty(HLW_OUT_DIR), out_dir = HLW_OUT_DIR; else, out_dir = getenv('HLW_OUT_DIR'); end
if isempty(out_dir), out_dir = HLW_code_dir; end
if exist('HLW_RESULTS_FILE', 'var') && ~isempty(HLW_RESULTS_FILE), results_file = HLW_RESULTS_FILE; else, results_file = getenv('HLW_RESULTS_FILE'); end
if isempty(results_file), results_file = fullfile(HLW_code_dir, 'HLW_results_no_lags.mat'); end

% Bauer-Rudebusch comparison inputs (LW_KF_BRcomp_rrstar*.png, used by draft.tex).
% PRIMARY (2026-09-18): <repo>/data/BR_istar.mat, variable BR_istar_data (T x 5):
%   [decimal date, i* posterior mean, 2.5%, 97.5%, LW inflation expectations],
%   i* columns in DECIMALS (x100 below), expectations in PERCENT; 1971q1-2026q1.
%   BR_directory, if set, is searched for BR_istar.mat first.
% FALLBACK: the BR posterior summary (tau_filtered_draw_* / states_names, as
%   saved by BauerRudebusch/BR_computations.m; BR_diagnostics_posterior.mat or
%   BR_plots_*.mat in Code/BauerRudebusch) with PTR from Code/HLW/LONGBASE.TXT
%   or ptr_frbus.mat.
if ~exist('BR_directory', 'var') || isempty(BR_directory)
  BR_directory = fullfile(HLW_code_dir, '..', '..', 'BauerRudebusch');
end
br_istar_file = fullfile(BR_directory, 'BR_istar.mat');
if exist(br_istar_file, 'file') ~= 2, br_istar_file = fullfile(HLW_code_dir, '..', '..', '..', 'data', 'BR_istar.mat'); end
br_file = fullfile(BR_directory, 'BR_diagnostics_posterior.mat');
if exist(br_file, 'file') ~= 2
    br_cand = dir(fullfile(BR_directory, 'BR_plots_*.mat'));
    if ~isempty(br_cand), br_file = fullfile(BR_directory, br_cand(end).name); end
end
longbase_file = fullfile(HLW_code_dir, '..', '..', 'HLW', 'LONGBASE.TXT');
HAVE_BR_ISTAR = exist(br_istar_file, 'file') == 2;
HAVE_BR = HAVE_BR_ISTAR || (exist(br_file, 'file') == 2 && ...
          (exist(longbase_file, 'file') == 2 || exist(fullfile(BR_directory, 'ptr_frbus.mat'), 'file') == 2));

% States shown in the panels and the summary-table columns. PlotVar.Names are
% the DISPLAY names used in file names, titles and table headers - the paper's
% convention, so that draft.tex keeps finding LW_KF_ldrrstar.png etc.
% PlotVar.Engine are the engine's augment_ssm labels (2026-09-15: l<m>NAME =
% NAME lagged m quarters, d<h>NAME = h-quarter change) used to index
% M_.states_of_interest; the two differ only for ldrrstar = l1drrstar.
PlotVar.Names  = cfg.PlotVar_Names(:);
PlotVar.Engine = PlotVar.Names;
for k_cfg = 1:size(cfg.display_to_engine, 1)
    PlotVar.Engine = strrep(PlotVar.Engine, cfg.display_to_engine{k_cfg, 1}, cfg.display_to_engine{k_cfg, 2});
end

SummaryTableMoments = nan(7, length(PlotVar.Names ));


%%
if HAVE_BR_ISTAR
    fprintf('HLW_plots: Bauer-Rudebusch i* from %s\n', br_istar_file);
    load(br_istar_file)                                        % BR_istar_data: [date, mean, 2.5%, 97.5%, LW_infl_exp]
    BR_data.NumericDate = BR_istar_data(:, 1)';
    BR_data.istar = 100 * BR_istar_data(:, 2);                 % decimals -> percent
    BR_data.p2_5  = 100 * BR_istar_data(:, 3);
    BR_data.p97_5 = 100 * BR_istar_data(:, 4);
    PTR = BR_istar_data(:, 5);                                 % inflation expectations, percent, on the BR dates
    BR_scale = 1;
    clear BR_istar_data
elseif HAVE_BR
    fprintf('HLW_plots: Bauer-Rudebusch posterior from %s\n', br_file);
    load(br_file)                                              % tau_filtered_draw_* , states_names (was BR_diagnostics_posterior / BR_plots_06042026)
    idx1 = find(strcmp(cellstr(states_names), 'tau'));         % states_names is a string array in the MATLAB-written file
    BR_T = size(tau_filtered_draw_mean, 2);                    % BR sample starts 1971Q1
    BR_data.NumericDate = 1971+(0:BR_T-1)/4;
    if exist(longbase_file, 'file') == 2
        % PTR from FRB/US LONGBASE.TXT: header of quoted names, rows "1962Q1",v1,v2,...; 'NA' = missing
        lb_txt   = fileread(longbase_file);
        lb_lines = regexp(lb_txt, '\r?\n', 'split');
        lb_lines = lb_lines(~cellfun(@isempty, lb_lines));
        lb_names = regexp(lb_lines{1}, '"([^"]*)"', 'tokens');  lb_names = cellfun(@(c) c{1}, lb_names, 'UniformOutput', false);
        c_obs = find(strcmp(lb_names, 'OBS'), 1);  c_ptr = find(strcmp(lb_names, 'PTR'), 1);
        lb_obs = nan(numel(lb_lines) - 1, 1);  lb_ptr = lb_obs;
        for ll = 2:numel(lb_lines)
            f = strsplit(lb_lines{ll}, ',');
            q = sscanf(strrep(f{c_obs}, '"', ''), '%dQ%d');
            lb_obs(ll-1) = q(1) + (q(2) - 1) / 4;
            lb_ptr(ll-1) = str2double(f{c_ptr});               % 'NA' -> NaN
        end
        PTR = interp1(lb_obs, lb_ptr, BR_data.NumericDate(:), 'nearest');
        fprintf('HLW_plots: PTR from %s (%d of %d BR quarters matched)\n', longbase_file, sum(~isnan(PTR)), BR_T);
    else
        load(fullfile(BR_directory, 'ptr_frbus'))              % ptr, older vintage ending 2024Q4
        PTR = [ptr(:); 2.11; 2.2; 2.18; 2.162; 2.1458];        % FRB/US PTR projections appended (as before)
        PTR = PTR(1:min(end, BR_T));
        if numel(PTR) < BR_T, PTR(end+1:BR_T) = PTR(end); end
    end
    BR_scale=1; %multiple needed to get the same units
    BR_data.istar = 100*tau_filtered_draw_mean(idx1,:)';
    BR_data.p2_5 = 100*tau_filtered_draw_pc2_5(idx1,:)';
    BR_data.p97_5 = 100*tau_filtered_draw_pc97_5(idx1,:)';
    clear results* tau* var* r_star* ptr lb_*
else
    fprintf(2, ['HLW_plots: Bauer-Rudebusch i* not found (data/BR_istar.mat, or BR_diagnostics_posterior.mat / BR_plots_*.mat in %s) - ' ...
                'LW_KF_BRcomp_rrstar*.png skipped.\n'], BR_directory);
end

fprintf('HLW_plots: loading %s ...\n', results_file); t_load = tic;
% FIXED 2026-09-29: a results file written by a .mod that ran after another model variant in the
% same MATLAB session (the ECM driver) contains that session's whole workspace, including this
% script's own control variables from its previous call (out_dir, results_file, cfg with ma_file,
% PlotVar, ...). The plain load below would silently substitute them for the caller's - the 2-step
% ECM plots landed in ECM_flat/ and used the flat assessment file. Everything defined before the
% load is therefore snapshotted and restored afterwards; only variables NOT yet defined come from the file.
HLW_pre_names = who;  HLW_pre = struct();
for HLW_i = 1:numel(HLW_pre_names)
    if ~any(strcmp(HLW_pre_names{HLW_i}, {'HLW_pre', 'HLW_pre_names', 'HLW_i'})), HLW_pre.(HLW_pre_names{HLW_i}) = eval(HLW_pre_names{HLW_i}); end
end
load_results_split(results_file);      % plain `load` for monolithic files; reattaches <name>_prior/_posterior.mat draws if present
HLW_pre_names = fieldnames(HLW_pre);
for HLW_i = 1:numel(HLW_pre_names), eval([HLW_pre_names{HLW_i} ' = HLW_pre.(HLW_pre_names{HLW_i});']); end
clear HLW_pre HLW_pre_names HLW_i
fprintf('  %.0f s\n', toc(t_load));
keep_pv = ismember(PlotVar.Engine, M_.states_of_interest);                       % keep only states the run diagnosed
if ~all(keep_pv), fprintf('HLW_plots: not in this run''s states_of_interest, dropped: %s\n', strjoin(PlotVar.Names(~keep_pv)', ', ')); end
PlotVar.Names = PlotVar.Names(keep_pv);  PlotVar.Engine = PlotVar.Engine(keep_pv);
SummaryTableMoments = SummaryTableMoments(:, 1:length(PlotVar.Names));
HAVE_BR = HAVE_BR && any(strcmp(PlotVar.Names, 'rrstar'));                       % the BR comparison is an r* comparison
% base state of every column (l20rrstar -> rrstar, d20rrstar -> rrstar, l1drrstar -> drrstar):
% supplies the innovation s.d. that normalises rows 1-2 of the summary table
PlotVar.Base = regexprep(PlotVar.Engine, '^[ld]\d+', '');

if ~exist(fullfile(out_dir, 'Graphs_level'), 'dir'), mkdir(fullfile(out_dir, 'Graphs_level')); end
if ~exist(fullfile(out_dir, 'Tables_level'), 'dir'), mkdir(fullfile(out_dir, 'Tables_level')); end
HLW_directory = out_dir;      % every print/fopen below writes to [HLW_directory '/Graphs_level'] etc.

% Read the entire Excel file
%BR_data = readtable([BR_directory 'Bauer and Rudebusch tau.xlsx']);
% BR_scale=100; %this was needed to get the same units
%{
load('BR_istar.mat');
% transform i* to percent
BR_istar_data(:,2:end-1)=BR_istar_data(:,2:end-1)*100;
PTR = BR_istar_data(:,end);
BR_scale=1; %multiple needed to get the same units
BR_data.NumericDate = BR_istar_data(:,1);
BR_data.istar=BR_istar_data(:,2);
BR_data.p2_5=BR_istar_data(:,3);
BR_data.p95_5=BR_istar_data(:,4);
%}


% Published HLW r* (from 1961), overlaid on the rrstar common-scale panel only
if cfg.rstar_overlay
    HLW_Rstar_up_to_205q3;
else
    HLW = [];
end

%%

% Sample dates from the data actually estimated on (was hard-wired 1960q1-2025q3)
T_data     = size(M_.data, 1);
TimeLabels = cfg.time_start + (0:T_data-1)/4;


T = length(TimeLabels);

% Delete non-sense draws from the prior: Those with imaginary correlations
% and those for which the Blanchard-Quah condition failed
NaN_Prior = cell2mat(cellfun(@isscalar, oo_.prior_function_results(:,3), 'UniformOutput', false));
Real_Prior = cell2mat(cellfun(@isreal, oo_.prior_function_results(:,2), 'UniformOutput', false));
N_prior = sum((1-NaN_Prior).*Real_Prior);
oo_.prior_function_results(NaN_Prior==1 | Real_Prior==0, :)=[];

% max(abs(imag(squeeze(CorrH.prior_ThreeD(PlotVar.Idx{pp},:,:)))))

Error_ACF.ThreeD = cat(3, oo_.posterior_function_results{:, 3});
[nstates,nacf, ndraws] = size(Error_ACF.ThreeD);
Error_ACF.prior_ThreeD = cat(3, oo_.prior_function_results{:, 3});

% Half-life of the estimation error of state j (CHANGED 2026-09-30):
%     tau_j = min{ h >= 1 : ACF_{j,h} < 1/2 },   tau_j = Inf if no such h (within the engine's horizon)
% where ACF_{j,h} is the autocorrelation of the filtering error at lag h (cell 3 of the wrapper,
% lag h in column h). Until 2026-09-30 this was the LAST down-crossing of 0.5, and 1 when the ACF
% never fell below 0.5 - which is the opposite of the truth for a draw whose error is persistent
% over the whole horizon. Inf propagates to the medians below (a median of Inf means more than
% half the draws never cross); the TeX tables print it as $\infty$, the histograms cap it, and the
% log-half-life scatters / correlations use the finite draws only.
Error_ACF.half_life       = first_crossing_below(Error_ACF.ThreeD, 0.5);                     % nstates x ndraws
Error_ACF.prior_half_life = first_crossing_below(Error_ACF.prior_ThreeD, 0.5);               % nstates x N_prior
% UNDEFINED half-life (ADDED 2026-10-16): a draw that solved (finite CorrH, cell 2) but whose error
% ACF is NaN at every lag has no steady-state estimation error (Omega_ii = 0, e.g. the 'invertible'
% shortcut of solve_steady_state_riccati: a square, invertible model such as Smets-Wouters without
% measurement error), so its ACF is 0/0. first_crossing_below would call that "never crosses"
% (Inf); it is set to NaN here, left out of the medians below, and a state whose half-life is
% undefined in every draw is printed as '--'. Failed draws (NaN in every cell) keep Inf as before.
hl_undefined = @(acf3, corr3) squeeze(all(isnan(acf3), 2) & all(isfinite(corr3), 2));
Error_ACF.half_life(reshape(hl_undefined(Error_ACF.ThreeD, cat(3, oo_.posterior_function_results{:, 2})), size(Error_ACF.half_life)))   = NaN;
Error_ACF.prior_half_life(reshape(hl_undefined(Error_ACF.prior_ThreeD, cat(3, oo_.prior_function_results{:, 2})), size(Error_ACF.prior_half_life))) = NaN;

state_estimates.prior_KF_Std_TV = cat(3, oo_.prior_function_results{:, 13});

state_estimates.prior_KF_Std = cat(2, oo_.prior_function_results{:, 6});

state_estimates.prior_FEV = cat(3, oo_.prior_function_results{:, 8});
state_estimates.FEV = cat(3, oo_.posterior_function_results{:, 8});

CorrH.ThreeD = cat(3, oo_.posterior_function_results{:, 2});
CorrH.prior_ThreeD = cat(3, oo_.prior_function_results{:, 2});
% R^2 of eq. (3.11) = squared correlation with the truth. CHANGED 2026-09-19:
% every graph and parameter-correlation table that used to show the correlation
% now shows this R^2 (files LW_KF_R2_ptiles_*, LW_KF_R2_BayesFactor_*,
% LW_KF_scatter_R2_*, LW_KF_scatter_R2_RMSE_*); the correlation histograms
% LW_KF_CorrH<h>_* were dropped in favour of the existing LW_KF_R2H<h>_* ones.
R2.ThreeD       = CorrH.ThreeD.^2;
R2.prior_ThreeD = CorrH.prior_ThreeD.^2;

% ADDED (2026-08): eq. (3.17), the R2 that may be compared ACROSS INFORMATION
% SETS. CorrH.^2 above is eq. (3.11), which normalises every state of interest
% by ITS OWN conditional variance (3.10). For the fixed-lag smoother columns
% - l20rrstar(t) = rrstar(t-20), ldrrstar(t) = drrstar(t-1) - that "own"
% variance is the one of the EARLIER state, so the smoothed R2 was compared
% against a denominator short by m quarters of accumulated innovation
% variance, and came out BELOW the filtered R2 in Table 3 even though the
% smoothing error variance is uniformly smaller. Section 3.5 flags exactly
% this ("if each estimate is allowed to generate its own denominator, the
% ratios need not inherit that ordering") and eq. (3.17) fixes it by keeping
% the FILTER's denominator, Var[r*_{t+h} | y^t], for both:
%
%     R2^(m)_{i,h} = 1 - [Omega_m]_ii / Var[x_{i,t+h} | y^t].
%
% Cell 10 of the wrapper carries it, cell 11 the m per state and cell 12 the
% base state supplying each denominator. Rows for states that are not lagged
% copies are identical to CorrH.^2, so the filtered columns of the table are
% unchanged. CorrH itself is left in place: the h-profile, histogram, scatter
% and Bayes-factor figures below are about the FILTERED states, where the two
% measures coincide.
%
% Cells 10-12 were EMPTY placeholders between the removal of the steady-state
% data pass and the eq.-(3.17) correction (both 2026-08), and the wrapper was
% 17 wide before and after, so the width cannot tell the two vintages apart.
% Test the CONTENT instead: a draw that carries eq. (3.17) has a cell 10 the
% same size as its CorrH. That test doubles as the failed-draw filter, since a
% failed solve writes a scalar NaN there.
if size(oo_.posterior_function_results, 2) >= 12
    ok_draw = find(cellfun(@(v) isnumeric(v) && isequal(size(v), size(CorrH.ThreeD(:,:,1))), ...
                           oo_.posterior_function_results(:, 10)), 1);
else
    ok_draw = [];
end
if ~isempty(ok_draw)
    R2_h.ThreeD       = cat(3, oo_.posterior_function_results{:, 10});
    R2_h.prior_ThreeD = cat(3, oo_.prior_function_results{:,     10});
    % The map is a property of the MODEL, not of the draw, so any successful
    % draw carries the same one.
    R2_h.lag_m = oo_.posterior_function_results{ok_draw, 11};
    R2_h.base  = oo_.posterior_function_results{ok_draw, 12};
else
    % Results saved before the eq.-(3.17) correction: fall back to eq. (3.11)
    % so that the script still runs, and say so loudly.
    warning(['HLW_plots: this HLW_results predates the eq.-(3.17) correction ' ...
             '(wrapper cell 10 does not carry R2_h). Falling back to eq. (3.11) - the ' ...
             'smoothed R2 columns of the summary table are then NOT comparable to the ' ...
             'filtered ones. Re-run the .mod file to refresh them.']);
    R2_h.ThreeD       = CorrH.ThreeD.^2;
    R2_h.prior_ThreeD = CorrH.prior_ThreeD.^2;
    R2_h.lag_m        = zeros(size(CorrH.ThreeD, 1), 1);
    R2_h.base         = M_.states_of_interest;
end

param_estimates.matrix= cat(2, oo_.posterior_function_results{:, 9});
prior_param_estimates.matrix= cat(2, oo_.prior_function_results{:, 9});
% probability bands for estimates
try
    clear surprises
    surprises.ThreeD = cat(3, oo_.posterior_function_results{:, 14});
% catch
    % surprises_struct.ThreeD = cat(3, oo_.posterior_function_results{:, 11});
end
% surprises.ThreeD = cat(3, oo_.posterior_function_results{:, 11});
state_estimates.ThreeD = cat(3, oo_.posterior_function_results{:, 13});
state_estimates.mean = mean(state_estimates.ThreeD, 3);
state_estimates.median = median(state_estimates.ThreeD, 3);
state_estimates.StdOfMeans= std(state_estimates.ThreeD, 0, 3);

state_estimates.KF_Std = cat(3, oo_.posterior_function_results{:, 16});
% state_estimates.KF_Std = cat(2, oo_.posterior_function_results{:, 6});
state_estimates.Avg_KF_Var= mean(state_estimates.KF_Std.^2, 3);
state_estimates.Median_KF_Var= median(state_estimates.KF_Std.^2, 3);

% state_estimates.OverallStd = sqrt(ones(T,1)*state_estimates.Avg_KF_Var' + state_estimates.StdOfMeans.^2);
state_estimates.OverallStd = sqrt(state_estimates.Avg_KF_Var' + state_estimates.StdOfMeans.^2);
% PERCENTILE BASIS (CHANGED 2026-09-20). The Gaussian 'overall' band above is
% dominated by the few pile-up draws with enormous filter variance (HLW
% baseline: top 5% of draws = 80% of the mean variance). Every band and the
% 'Average s.d.' table row now use the 2.5/97.5 percentiles of the mixture of
% per-draw N(mean, filter variance) (mixture_band); OverallStd becomes the
% Gaussian-equivalent s.d. (95% band width / 3.92), the old definition is kept
% as OverallStd_gauss for reference.
state_estimates.OverallStd_gauss = state_estimates.OverallStd;
n_states_se = size(state_estimates.ThreeD, 2);
state_estimates.q_lo = nan(size(state_estimates.mean));  state_estimates.q_hi = state_estimates.q_lo;  state_estimates.q_med = state_estimates.q_lo;
for j_se = 1:n_states_se
    [state_estimates.q_lo(:, j_se), state_estimates.q_hi(:, j_se), state_estimates.OverallStd(:, j_se)] = ...
        mixture_band(squeeze(state_estimates.ThreeD(:, j_se, :)), squeeze(state_estimates.KF_Std(j_se, :, :)));
    state_estimates.q_med(:, j_se) = mixture_band(squeeze(state_estimates.ThreeD(:, j_se, :)), squeeze(state_estimates.KF_Std(j_se, :, :)), 0.5, 0.5);   % mixture median (ADDED 2026-09-27; saved in the summary file)
end

nu_draws = cat(2, oo_.posterior_function_results{:, 17}); % n_key_states x n_draws
nu_prior = cat(2, oo_.prior_function_results{:, 17});

% UNITS GUARD (ADDED 2026-09). Since tail_unreliability_index moved to the
% general-integration-order formula of Rate_of_convergence_R2.pdf (eq. (5)),
% nu_i is measured in quarters^(2d+1) with d = d_i the integration order of
% state i, returned alongside nu in cell 18. It is in QUARTERS - which is what
% the summary tables below print - only where d = 0. That is the case for every
% one of HLW's states_of_interest (ygap is stationary and gets NaN; rrstar and
% its lags/differences are I(1)), so the tables are unaffected; but a state of
% higher order added later would be silently mislabelled and, worse, pooled
% into a median with entries in a different unit. Fail loudly instead.
% Guarded on the width the way the eq.-(3.17) block above is: a results set
% saved by the pre-2026-09, 17-wide wrapper has no column 18, and every state
% of interest was I(1) or stationary then anyway.
if size(oo_.posterior_function_results, 2) >= 18
    nu_order_draws = cat(2, oo_.posterior_function_results{:, 18});
    % CHANGED 2026-09-18 (was an error): states with d > 0 - e.g. an I(2)
    % potential output with a random-walk growth state - have nu in
    % quarters^(2d+1) and are NOT printed in the quarters rows; their nu is set
    % to NaN here, which the tables print as '--' (the paper leaves nu blank for
    % I(2) processes). The raw values stay in cell 17 of the results.
    d_med = round(median(nu_order_draws, 2, 'omitnan'));
    if any(d_med(~isnan(d_med)) ~= 0)
        for nn = find(d_med(:)' > 0)
            fprintf(2, 'HLW_plots: %s has integration order d = %d; its nu is in quarters^%d and is left blank in the tables (median raw value %.3g).\n', ...
                    M_.states_of_interest{nn}, d_med(nn), 2*d_med(nn)+1, median(nu_draws(nn, :), 'omitnan'));
            nu_draws(nn, :) = NaN;
        end
    end
else
    nu_order_draws = zeros(size(nu_draws));
end

%% Steady-state Riccati convergence (ADDED 2026-09)
% Cell 19 is the 1x4 row [converged, iterations, rel_step, rho_cl] that
% reliability_draw_wrapper records for the steady-state Riccati solve, or
% nan(1,4) for a draw that failed before reaching it. Everything built on the
% steady-state covariance - cells {1}-{8} and {10}, hence every half-life,
% FEV, CorrH and R2 printed below - is only meaningful for draws whose first
% entry is 1, so the share is reported before the summary tables rather than
% after them.
%
% Guarded on the width the way the eq.-(3.17) and nu_order blocks above are:
% a results set saved by the pre-2026-09, 18-wide wrapper has no column 19.
%
% NOTE on the two denominators. The prior set has already had its failed draws
% deleted above (the NaN_Prior/Real_Prior mask), while the posterior set has
% not, so the two would otherwise not be comparable. Both are printed: the
% share among draws that actually REACHED the Riccati solve is the one that
% says something about the solver, and the share of all draws folds in the
% Blanchard-Kahn/indeterminacy failures that never got that far.
if size(oo_.posterior_function_results, 2) >= 19
    ric_sets = {'posterior', oo_.posterior_function_results};
    if isfield(oo_, 'prior_function_results') && ~isempty(oo_.prior_function_results) ...
            && size(oo_.prior_function_results, 2) >= 19
        ric_sets(end+1, :) = {'prior', oo_.prior_function_results};
    end

    fprintf('\n=== Steady-state Riccati convergence (results cell 19) ===\n');
    for ss = 1:size(ric_sets, 1)
        set_name = ric_sets{ss, 1};
        ric      = cat(1, ric_sets{ss, 2}{:, 19});   % n_draws x 4, NaN row = failed draw
        solved   = ~isnan(ric(:, 1));                % draws that reached the solve
        reached  = ric(:, 1) == 1;                   % NaN rows compare false, as intended
        n_all    = size(ric, 1);
        n_solved = sum(solved);
        n_conv   = sum(reached);

        if n_solved == 0
            fprintf(2, '  %-9s : no draw reached the Riccati solve (%d draws)\n', ...
                    set_name, n_all);
            continue
        end
        disp('  ')
        fprintf(['  %-9s : %d/%d converged = %.2f%% of draws reaching the solve ' ...
                 '(%.2f%% of all %d draws)\n'], ...
                set_name, n_conv, n_solved, 100*n_conv/n_solved, ...
                100*n_conv/n_all, n_all);

        if n_conv > 0
            rho   = ric(reached, 4);
            iters = ric(reached, 2);
            fprintf('%14s rho(A-KC): median %.6f, max %.6f\n', '', median(rho), max(rho));
            % Since 2026-09-16 the Riccati equation is solved directly
            % (solve_steady_state_riccati); iterations = 0 then means a direct
            % DARE solve, > 0 the fixed-point fallback, and column 3 is an
            % error estimate r/(1-rho^2) rather than the last iteration step.
            fprintf('%14s iterations: median %d, max %d (0 = direct DARE solve); error estimate: median %.1e, max %.1e\n', '', ...
                    round(median(iters)), max(iters), median(ric(reached, 3)), max(ric(reached, 3)));

            % rho at 1 is the pile-up signature: sigma_g or sigma_ystar at the
            % zero boundary leaves A-KC with a unit root, so Omega is a genuine
            % fixed point but the level direction is never contracted. These
            % draws DO converge, so the share above will not isolate them.
            n_unit = sum(rho >= 1 - 1e-8);
            if n_unit > 0
                fprintf(2, ['%14s NOTE: %d converged draw(s) have rho(A-KC) >= 1-1e-8 ' ...
                            '(pile-up: Omega is a fixed point, but the level ' ...
                            'direction is never contracted).\n'], '', n_unit);
            end
        end

        n_bad = n_solved - n_conv;
        if n_bad > 0
            bad         = solved & ~reached;
            n_nonfinite = sum(isnan(ric(bad, 4)));   % non-finite iterate => rho_cl = NaN
            fprintf(2, ['%14s WARNING: %d draw(s) did NOT converge (%d hit the ' ...
                        'iteration cap, %d went non-finite). Screen on ' ...
                        'cell 19 column 1 before taking medians across draws.\n'], ...
                    '', n_bad, n_bad - n_nonfinite, n_nonfinite);
        end
    end
    fprintf('\n');
else
    fprintf(2, ['\nSteady-state Riccati convergence: results cell 19 is absent ' ...
                '(this set is %d wide, saved by a pre-2026-09 wrapper). Re-run ' ...
                'the prior/posterior functions to record it.\n\n'], ...
            size(oo_.posterior_function_results, 2));
end

%% Start populating Summary Table - Table 1: SD

% Rows 1-3 for every column, normalised by the innovation s.d. of the column's
% BASE state (GENERALISED 2026-09-17: was two loops over names containing
% 'rrstar' and 'gap', both using rrstar's / ygap's B row; PlotVar.Base gives the
% same rows for HLW - drrstar and rrstar share their shock loading - and works
% for any model).
for pp = 1:length(PlotVar.Names)
    PlotVar.Idx{pp} = strcmp(M_.states_of_interest, PlotVar.Engine{pp});
    base_position_in_x = key_state_indices(strcmp(M_.states_of_interest, PlotVar.Base{pp}));
    if isempty(base_position_in_x), continue; end
    sd_innov = norm(ModelMatrices.B(base_position_in_x,:));
    SummaryTableMoments(1,pp) = sqrt(mean(state_estimates.Median_KF_Var(PlotVar.Idx{pp},end-80:end)))/sd_innov; % Median s.d. given parameters
    SummaryTableMoments(2,pp) = mean(state_estimates.OverallStd(end-80:end,PlotVar.Idx{pp}))/sd_innov;          % Averaged s.d. over time
    SummaryTableMoments(3,pp) = sd_innov;
end

% Summary of the filtered paths across draws (mean, median, dispersion), WITHOUT
% the per-draw arrays .ThreeD / .KF_Std / .prior_* (CHANGED 2026-09-18: with
% them the file was 3.8 GB for UC3, above MATLAB's 2 GB per-variable limit).
% Explicit .mat: Octave does not append it.
state_estimates_summary = rmfield(state_estimates, intersect(fieldnames(state_estimates), {'ThreeD', 'KF_Std', 'prior_KF_Std_TV', 'prior_KF_Std', 'prior_FEV', 'FEV'}));
state_estimates_summary.states_of_interest = M_.states_of_interest;   % column order of .mean/.OverallStd (ADDED 2026-09-19: lets Figure17_HLW_SPF.m read the file alone)
state_estimates_summary.TimeLabels         = TimeLabels;
save(fullfile(HLW_directory, 'HLW_state_estimates.mat'), 'state_estimates_summary', '-v7')
clear state_estimates_summary

%% Plot filtered series
close all;

% Hide the bands of the first cfg.band_skip quarters in every figure below (ADDED
% 2026-09-27, see the config block; 0 for HLW). The summary file has been saved
% above with the full bands; the tables use OverallStd, which is untouched.
if cfg.band_skip > 0
    n_skip = min(cfg.band_skip, size(state_estimates.q_lo, 1));
    state_estimates.q_lo(1:n_skip, :) = NaN;
    state_estimates.q_hi(1:n_skip, :) = NaN;
end

% Indices of the states that share a y-axis (names containing cfg.common_scale_base)
rrstar_indices = find(contains(PlotVar.Names, cfg.common_scale_base));

% First pass: calculate common y-axis limits for those variables
rrstar_y_min = inf;
rrstar_y_max = -inf;

for pp = rrstar_indices'
    PlotVar.Idx{pp} = strcmp(M_.states_of_interest, PlotVar.Engine{pp});

    if isequal(PlotVar.Names{pp}, ['l20' cfg.common_scale_base])
        PlotLag = 20;
        % Include regular rrstar in y-limit calculation for l20rrstar plot
        regular_pp = find(strcmp(PlotVar.Names, cfg.common_scale_base));
        PlotVar.Idx{regular_pp} = strcmp(M_.states_of_interest, PlotVar.Engine{regular_pp});
        regular_data = [state_estimates.q_lo(:,PlotVar.Idx{regular_pp}), state_estimates.mean(:,PlotVar.Idx{regular_pp}), state_estimates.q_hi(:,PlotVar.Idx{regular_pp})];
        rrstar_y_min = min(rrstar_y_min, min(regular_data(:)));
        rrstar_y_max = max(rrstar_y_max, max(regular_data(:)));
    else
        PlotLag = 0;
    end

    plot_data = [state_estimates.q_lo(:,PlotVar.Idx{pp}), state_estimates.mean(:,PlotVar.Idx{pp}), state_estimates.q_hi(:,PlotVar.Idx{pp})];
    rrstar_y_min = min(rrstar_y_min, min(plot_data(:)));
    rrstar_y_max = max(rrstar_y_max, max(plot_data(:)));
end

% Add some padding to the y-limits
y_range = rrstar_y_max - rrstar_y_min;
rrstar_y_min = rrstar_y_min - 0.05 * y_range;
rrstar_y_max = rrstar_y_max + 0.05 * y_range;
rrstar_y_min = max(rrstar_y_min, cfg.common_scale_ylim(1));   % HLW: r* panels capped at [-10, 10]; configurable since 2026-09-27
rrstar_y_max = min(rrstar_y_max, cfg.common_scale_ylim(2));

if ~DO_TABLES_ONLY   % (filtered-series and BR panels)
% Second pass: create plots
for pp = 1:length(PlotVar.Names)
    PlotVar.Idx{pp} = strcmp(M_.states_of_interest, PlotVar.Engine{pp});

    if isequal(PlotVar.Names{pp}, ['l20' cfg.common_scale_base])
        PlotLag = 20;
    else
        PlotLag = 0;
    end

    % Regular scale version
    figure("Name", PlotVar.Names{pp})
    yy_pp = [state_estimates.q_lo(:,PlotVar.Idx{pp}), state_estimates.mean(:,PlotVar.Idx{pp}), state_estimates.q_hi(:,PlotVar.Idx{pp})];
    if cfg.lag_presample_mask && PlotLag == 20
        yy_pp(1:min(PlotLag, size(yy_pp, 1)), :) = NaN;   % smoother values dated before the sample start (ADDED 2026-09-27)
    end
    my_hndl = plot(TimeLabels-PlotLag/4, yy_pp, 'Color','b', 'LineWidth',2);
    set(my_hndl, {'LineStyle'}, {'--';'-';'--'})

    if PlotLag == 20
        set(my_hndl, {'Color'}, {'r';'r';'r'});
        hold all;
        regular_pp = find(strcmp(PlotVar.Names, cfg.common_scale_base));
        % FIXED 2026-09-27: the filter x(t|t) is dated t; it used to be plotted at
        % t - 5 years like the smoother, i.e. misdated by 20 quarters in this panel
        Overall_hndl = plot(TimeLabels, [state_estimates.q_lo(:,PlotVar.Idx{regular_pp}), state_estimates.mean(:,PlotVar.Idx{regular_pp}), state_estimates.q_hi(:,PlotVar.Idx{regular_pp})], 'Color','b', 'LineWidth',2);
        set(Overall_hndl, {'LineStyle'}, {'--';'-';'--'});
        legend([Overall_hndl(2), my_hndl(2)], {'Filter';'20-qtr lag smoother'})
    end
    grid on;
    print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix PlotVar.Names{pp} '.png']);

    % Post-1965 version (ADDED 2026-09-20, full plot set only): same figure
    % restricted to 1965 onwards, y-axis fitted to the plotted lines over that
    % window so the diffuse-start band of the first years no longer sets the scale.
    if ~DO_MAIN_PLOTS_ONLY
        ax_lines = findall(gca, 'Type', 'line');
        y_from = cfg.late_start;                      % 1965 for HLW; configurable since 2026-09-27 (file suffix = the year)
        y_lo = inf; y_hi = -inf;
        for k_ln = 1:numel(ax_lines)
            xd = get(ax_lines(k_ln), 'XData');  yd = get(ax_lines(k_ln), 'YData');
            sel = xd >= y_from & isfinite(yd);
            if any(sel), y_lo = min(y_lo, min(yd(sel))); y_hi = max(y_hi, max(yd(sel))); end
        end
        if isfinite(y_lo) && isfinite(y_hi) && y_hi > y_lo
            pad = 0.05 * (y_hi - y_lo);
            xlim([y_from, max(TimeLabels)]);
            ylim([y_lo - pad, y_hi + pad]);
            print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix PlotVar.Names{pp} '_' num2str(round(cfg.late_start)) '.png']);
        end
    end

    % Common scale version for the common_scale_base variables
    if contains(PlotVar.Names{pp}, cfg.common_scale_base)
        figure("Name", [PlotVar.Names{pp} '_common_scale'])
        yy_med = state_estimates.mean(:,PlotVar.Idx{pp}) + sqrt(state_estimates.Median_KF_Var(PlotVar.Idx{pp},:))'*[norminv(0.025), -norminv(0.025)];
        if cfg.lag_presample_mask && PlotLag == 20
            yy_med(1:min(PlotLag, size(yy_med, 1)), :) = NaN;   % (yy_pp already masked above)
        end
        my_hndl_common = plot(TimeLabels-PlotLag/4, yy_pp, 'Color','b', 'LineWidth',2);
        set(my_hndl_common, {'LineStyle'}, {'--';'-';'--'})
        hold all;
        hndl_median_var = plot(TimeLabels-PlotLag/4, yy_med, 'Color','b', 'LineWidth',1);
        set(hndl_median_var , {'LineStyle'}, {'--';'--'})
        if cfg.rstar_overlay && isequal(PlotVar.Names{pp}, cfg.common_scale_base) && numel(HLW) == numel(TimeLabels) - 4
            plot(TimeLabels(5:end), HLW, '--k', 'LineWidth', 2);
        end

        if PlotLag == 20
            set(my_hndl_common, {'Color'}, {'r';'r';'r'});
            hold all;
            regular_pp = find(strcmp(PlotVar.Names, cfg.common_scale_base));
            % FIXED 2026-09-27: filter dated t, not t - 5 years (see the regular-scale panel)
            Overall_hndl_common = plot(TimeLabels, [state_estimates.q_lo(:,PlotVar.Idx{regular_pp}), state_estimates.mean(:,PlotVar.Idx{regular_pp}), state_estimates.q_hi(:,PlotVar.Idx{regular_pp})], 'Color','b', 'LineWidth',2);
            set(Overall_hndl_common, {'LineStyle'}, {'--';'-';'--'});
            legend([Overall_hndl_common(2), my_hndl_common(2)], {'Filter';'20-qtr lag smoother'})
        end

        ylim([max(rrstar_y_min, cfg.common_scale_ylim(1)), min(rrstar_y_max, cfg.common_scale_ylim(2))]);
        grid on;
        print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix PlotVar.Names{pp} '_common_scale.png']);
    end
end


% BR comparison (skipped without the Bauer-Rudebusch posterior files, see HAVE_BR)
for DO_NOM=0:double(HAVE_BR)*1
    if ~HAVE_BR, break; end
    switch DO_NOM
        case 0
            suffix = '';
            exp_mult=0; % Multiply inflation by this number
        otherwise
            suffix = '_nom';
            exp_mult=1; % Multiply inflation by this number
    end
    pp=find(strcmp(PlotVar.Names, 'rrstar'));
    figure("Name",[PlotVar.Names{pp} ' vs BR'])
    PlotVar.Idx{pp} = strcmp(M_.states_of_interest, PlotVar.Engine{pp});

    PlotLag=0;

    BR_to_LW_idx=(TimeLabels>=BR_data.NumericDate(1) & TimeLabels<=BR_data.NumericDate(end));
    my_hndl=plot(TimeLabels(BR_to_LW_idx), repmat(PTR(TimeLabels(end)>=BR_data.NumericDate),1,3)*exp_mult + [state_estimates.q_lo(BR_to_LW_idx,PlotVar.Idx{pp}), state_estimates.mean(BR_to_LW_idx,PlotVar.Idx{pp}), state_estimates.q_hi(BR_to_LW_idx,PlotVar.Idx{pp})], 'Color','b', 'LineWidth',2);
    % my_hndl=plot(TimeLabels-PlotLag/4, repmat(BR_data{:,5},1,3)*exp_mult+[state_estimates.q_lo(:,PlotVar.Idx{pp}), state_estimates.mean(:,PlotVar.Idx{pp}), state_estimates.q_hi(:,PlotVar.Idx{pp})], 'Color','b', 'LineWidth',2);
    set(my_hndl, {'LineStyle'}, {'--';'-';'--'})
    hold all;
    BR_hndl=plot(BR_data.NumericDate, BR_scale*[BR_data.istar BR_data.p2_5 BR_data.p97_5], 'Color','r', 'LineWidth',2);
    set(BR_hndl, {'LineStyle'}, {'-';'--';'--'})
    grid on;

    if DO_NOM==1
        legend([my_hndl(2), BR_hndl(1)], {'HLW+E[\pi]';'BR'})
    else
        legend([my_hndl(2), BR_hndl(1)], {'HLW';'BR'})
    end

    % BR_corr=corr(BR_scale*BR_data{:,2}-BR_data{:,5}*exp_mult, state_estimates.mean(TimeLabels>=BR_data.NumericDate(1) & TimeLabels<=BR_data.NumericDate(end),PlotVar.Idx{pp}));
    % BR1982_corr=corr(BR_scale*BR_data.MeanTau(BR_data.NumericDate>=1982)-exp_mult*BR_data{BR_data.NumericDate>=1982,5}, state_estimates.mean(TimeLabels>=1982 & TimeLabels<=BR_data.NumericDate(end),PlotVar.Idx{pp}));
    BR_corr=corr(BR_data.istar(TimeLabels(end)>=BR_data.NumericDate), PTR(TimeLabels(end)>=BR_data.NumericDate)*exp_mult+state_estimates.mean(TimeLabels>=BR_data.NumericDate(1) & TimeLabels<=BR_data.NumericDate(end),PlotVar.Idx{pp}));
    % BR1982_corr=corr(100*BR_data.MeanTau(BR_data.NumericDate>=1982)-exp_mult*BR_data{BR_data.NumericDate>=1982,5}, state_estimates.mean(TimeLabels>=1982 & TimeLabels<=BR_data.NumericDate(end),PlotVar.Idx{pp}));

    % Get the position of the axes within the figure
    ax_pos = get(gca, 'Position'); % [left, bottom, width, height] of axes

    % Position textbox inside the axes area
    % Add small margins from the axes edges
    left_margin = 0.02;   % Small margin from left edge of plot
    top_margin = 0.05;    % Small margin from top edge of plot
    box_width = 0.25;
    box_height = 0.08;


    my_string=sprintf('Overall correlation=%.2f', BR_corr);
    annotation('textbox', [ax_pos(1) + left_margin, ...
                          ax_pos(2) + ax_pos(4) - top_margin - box_height, ...
                          box_width, box_height], ...
        'String', my_string, ...
        'FitBoxToText', 'on', ...
        'BackgroundColor', 'white', ...
        'EdgeColor', 'black', ...
        'FontSize', 10);

    print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'BRcomp_' PlotVar.Names{pp} suffix '.png']);

    delete(my_hndl([1 3]))
    delete(BR_hndl([2 3]))

    print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'BRcomp_' PlotVar.Names{pp} suffix '_simple.png']);
end
end   % ~DO_TABLES_ONLY (filtered-series and BR panels)

%% Continue populating Summary Table - Table 1: now Half-life

for pp = 1:length(PlotVar.Names)          % (GENERALISED 2026-09-17: was the 'rrstar' and 'gap' name loops)
    PlotVar.Idx{pp} = strcmp(M_.states_of_interest, PlotVar.Engine{pp});
    hl_pp = Error_ACF.half_life(PlotVar.Idx{pp},:);
    SummaryTableMoments(7,pp) = median(hl_pp(~isnan(hl_pp)));   % undefined draws left out; all undefined -> NaN -> '--' (2026-10-16)
end

%% Plot histogram of half-life

if ~DO_TABLES_ONLY
for pp=1:length(PlotVar.Names)
    figure("Name",PlotVar.Names{pp})
    PlotVar.Idx{pp} = strcmp(M_.states_of_interest, PlotVar.Engine{pp});

    histogram(min(Error_ACF.half_life(PlotVar.Idx{pp},:), 50), 'Normalization', 'probability');
    grid on;
    hold on;
    prior_hist_hndl=histogram(min(Error_ACF.prior_half_life(PlotVar.Idx{pp},:), 50), 'Normalization', 'probability');

    old_ylim=ylim;
    hl_pp = Error_ACF.half_life(PlotVar.Idx{pp},:);  hl_pp = hl_pp(~isnan(hl_pp));   % undefined draws left out (2026-10-16)
    plot(min(median(hl_pp), 50)*[1,1],old_ylim,'--k', 'LineWidth', 2)   % Inf median shown at the cap; no line if undefined
    axis([0 50 old_ylim]);
    disp([PlotVar.Names{pp} ' Half-life: ' num2str(median(hl_pp))  ])

    print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'HalfLife_' PlotVar.Names{pp} '_prior.png']);
    if pp==1
        % Get the current axes
        ax = gca;

        % Get the color order and use the appropriate index (get(): a graphics
        % handle is a double in Octave, so ax.ColorOrder does not work there)
        color_order = get(ax, 'ColorOrder');
        color_index = 2; % For the first histogram, use index 1

        % Get the actual color
        actual_color = color_order(2, :);
    end
    delete(prior_hist_hndl)
    print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'HalfLife_' PlotVar.Names{pp} '.png']);
end
end   % ~DO_TABLES_ONLY (half-life histograms)
% return


%% Continue populating Summary Table - Table 1: now R^2

% CHANGED (2026-08): rows 4-5 now report eq. (3.17) (R2_h) instead of eq.
% (3.11) (CorrH.^2). The table puts the FILTER and the 5-yr fixed-lag SMOOTHER
% of r* in adjacent columns, which is a comparison across information sets, and
% Section 3.5 requires reliability to be weakly increasing in the hindsight m.
% Eq. (3.11) does not deliver that here because the smoothed estimate lives on
% an augmented LAGGED copy of the state and so brings its own, smaller
% denominator; eq. (3.17) normalises both columns by Var[r*_{t+h} | y^t]. The
% filtered columns are numerically unchanged (m = 0 makes the two identical).
%
% THIS IS ABOUT THESE TWO ROWS ONLY. Eq. (3.11) is still what every CHART in
% this script plots - the CorrH histograms and percentile fans, the Bayes-
% factor panels, and the high/low-R^2 masks that condition the r* paths. That
% is deliberate and it is safe wherever a panel compares states at the SAME
% information set: the cross-state scatters pair rrstar with ygap and with
% d20rrstar, all contemporaneous (m = 0), where the two formulas return the
% same number.
%
% The exception to keep in mind is the per-state panel loop, which runs over
% every name in PlotVar and so includes the LAGGED copies. A CorrH panel for
% l20rrstar is normalised by the conditional variance of the EARLIER state, so
% its level is not on the same footing as the rrstar panel beside it; that
% particular comparison is what rows 4-5 and the audit print below are for. The
% audit lists m state by state, so which panels this touches can be read off
% there rather than guessed from the names.
for pp = 1:length(PlotVar.Names)          % (GENERALISED 2026-09-17: was the 'rrstar' and 'gap' name loops)
    PlotVar.Idx{pp} = strcmp(M_.states_of_interest, PlotVar.Engine{pp});
    hh=4;
    SummaryTableMoments(4,pp) = median(squeeze(R2_h.ThreeD(PlotVar.Idx{pp},hh,:)));
    hh=20;
    SummaryTableMoments(5,pp) = median(squeeze(R2_h.ThreeD(PlotVar.Idx{pp},hh,:)));
    SummaryTableMoments(6,pp) = median(nu_draws(PlotVar.Idx{pp},:));
end

%% Eq. (3.17) audit: what the correction changed, and does hindsight help?
% Printed, not plotted. The first block is the like-for-like comparison of the
% two normalisations at the two horizons in the table; the second checks the
% Section 3.5 ordering draw by draw, which must hold at EVERY horizon.
fprintf('\n--- Summary-table R^2: eq. (3.17) vs eq. (3.11) (medians across draws) ---\n');
fprintf('%-12s %5s   %-24s %-24s\n', 'state', 'm', 'eq.(3.17)  h=4    h=20', 'eq.(3.11)  h=4    h=20');
for pp = 1:length(PlotVar.Names)
    PlotVar.Idx{pp} = strcmp(M_.states_of_interest, PlotVar.Engine{pp});
    fprintf('%-12s %5d   %11.4f %8.4f      %11.4f %8.4f\n', PlotVar.Names{pp}, ...
        R2_h.lag_m(PlotVar.Idx{pp}), ...
        median(squeeze(R2_h.ThreeD(PlotVar.Idx{pp},4,:))),  median(squeeze(R2_h.ThreeD(PlotVar.Idx{pp},20,:))), ...
        median(squeeze(CorrH.ThreeD(PlotVar.Idx{pp},4,:).^2)), median(squeeze(CorrH.ThreeD(PlotVar.Idx{pp},20,:).^2)));
end

% Weak monotonicity in m, draw by draw and horizon by horizon. Each smoothed
% state is checked against the FILTERED state it is a lagged copy of.
for pp = 1:length(PlotVar.Names)
    if R2_h.lag_m(PlotVar.Idx{pp}) == 0
        continue
    end
    % The base state comes from the engine's own structural map (cell 12), so
    % no name parsing is repeated here. It is skipped when the base is not
    % itself one of the reported states - there is then nothing to compare to.
    base_pp = find(strcmp(PlotVar.Names, R2_h.base{PlotVar.Idx{pp}}));
    if isempty(base_pp)
        continue
    end
    gap = squeeze(R2_h.ThreeD(PlotVar.Idx{pp},:,:) - R2_h.ThreeD(PlotVar.Idx{base_pp},:,:));
    fprintf('Hindsight gain %s vs %s: min over horizons and draws = %+.2e, median at h=20 = %+.4f\n', ...
        PlotVar.Names{pp}, PlotVar.Names{base_pp}, min(gap(:)), median(gap(20,:)));
    if min(gap(:)) < -1e-8
        warning('HLW_plots: eq. (3.17) is not weakly increasing in m for %s. Section 3.5 says it must be.', ...
            PlotVar.Names{pp});
    end
end
fprintf('\n');

%% Percentile fans and Bayes factors of R^2 (eq. 3.11) by horizon
% (the correlation histograms that used to open this block were dropped
% 2026-09-19; the R^2 histograms LW_KF_R2H<h>_* below cover them)

if ~DO_TABLES_ONLY
for pp=1:length(PlotVar.Names)
    PlotVar.Idx{pp} = strcmp(M_.states_of_interest, PlotVar.Engine{pp});

      % Percentile fans of R^2 by horizon, posterior (blue) and prior; part of
      % the main set since 2026-09-19 (the paper's ECM comparison uses them)
      figure("Name",PlotVar.Names{pp})

      my_hndl=semilogx(1:size(R2.ThreeD,2), prctile(squeeze(R2.ThreeD(PlotVar.Idx{pp},:,:)), [5 16 50 84 95], 2), 'LineWidth', 2, 'Color', 'b');
      set(gca, 'XTick', [1 2 4 10 20 50 100 250], 'XTickLabel', {'1','2','4','10','20','50','100','250'});   % plain labels (Octave prints 2x10^0 otherwise)
      set(my_hndl, {'LineStyle'}, {'-.';'--';'-';'--';'-.'});
      xlabel('quarters (log-scale)');
      ylabel('R^2');
      ylim([0 1]);
      grid on;
      hold on
      my_prior_hndl=semilogx(1:size(R2.prior_ThreeD,2), prctile(squeeze(R2.prior_ThreeD(PlotVar.Idx{pp},:,:)), [5 16 50 84 95], 2), 'LineWidth', 2, 'Color', actual_color);
      set(gca, 'XTick', [1 2 4 10 20 50 100 250], 'XTickLabel', {'1','2','4','10','20','50','100','250'});   % plain labels (Octave prints 2x10^0 otherwise)
      set(my_hndl, {'LineStyle'}, {'-.';'--';'-';'--';'-.'});
      set(my_prior_hndl, {'LineStyle'}, {'-.';'--';'-';'--';'-.'});
      grid on;

      print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'R2_ptiles_' PlotVar.Names{pp} '_prior.png']);
      % pause(0.1);
      delete(my_prior_hndl);
      print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'R2_ptiles_' PlotVar.Names{pp} '.png']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bayes Factor plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure("Name",PlotVar.Names{pp})
    my_hndl=semilogx(1:size(R2.ThreeD,2), log10(mean(squeeze(R2.ThreeD(PlotVar.Idx{pp},:,:))>0.5,2)./(mean(squeeze(R2.prior_ThreeD(PlotVar.Idx{pp},:,:))>0.5,2))), 'LineWidth', 2, 'Color', 'b');
    set(gca, 'XTick', [1 2 4 10 20 50 100 250], 'XTickLabel', {'1','2','4','10','20','50','100','250'});   % plain labels (Octave prints 2x10^0 otherwise)
    xlabel('quarters (log-scale)');
    ylabel('Bayes factor: R^2 at h > x (log10)');
    grid on;
    hold on
    threshold1 = 1/4;
    my_hndl=semilogx(1:size(R2.ThreeD,2), log10(mean(squeeze(R2.ThreeD(PlotVar.Idx{pp},:,:))>threshold1,2)./(mean(squeeze(R2.prior_ThreeD(PlotVar.Idx{pp},:,:))>threshold1,2))), '--r','LineWidth', 2);
    threshold2 = 3/4;
    my_hndl=semilogx(1:size(R2.ThreeD,2), log10(mean(squeeze(R2.ThreeD(PlotVar.Idx{pp},:,:))>threshold2,2)./(mean(squeeze(R2.prior_ThreeD(PlotVar.Idx{pp},:,:))>threshold2,2))), '--k','LineWidth', 2);
    legend({'x=0.5'; ['x=' num2str(threshold1, '%1.2f')]; ['x=' num2str(threshold2, '%1.2f')]}, 'location', 'best')
    legend('boxoff')
    axis([1 size(R2.ThreeD,2) min(ylim) max(ylim)])
    print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'R2_BayesFactor_' PlotVar.Names{pp} '.png']);
end
end   % ~DO_TABLES_ONLY (R^2 percentile fans, Bayes factors)
close all


%% Plot histogram of R2h

if ~DO_TABLES_ONLY
for pp=1:length(PlotVar.Names)
    PlotVar.Idx{pp} = strcmp(M_.states_of_interest, PlotVar.Engine{pp});

    for hh=[1 4 20 100 250]
        figure("Name",PlotVar.Names{pp})
        histogram( squeeze(R2_h.ThreeD(PlotVar.Idx{pp},hh,:)), 'Normalization', 'probability');
        grid on;
        hold on;
        old_ylim=ylim;
        plot(median(squeeze(R2_h.ThreeD(PlotVar.Idx{pp},hh,:)))*[1,1],old_ylim,'--k', 'LineWidth', 2)
        axis([0 1 old_ylim]);
        disp([PlotVar.Names{pp} ' R2: ' num2str(hh) ' :: ' num2str(median(squeeze(R2_h.ThreeD(PlotVar.Idx{pp},hh,:))))  ' :: ' num2str(median(squeeze(CorrH.ThreeD(PlotVar.Idx{pp},hh,:)))^2)  ])
        print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'R2H' num2str(hh) '_' PlotVar.Names{pp} '.png']);
    end


end
end   % ~DO_TABLES_ONLY (R2h histograms)
close all

%% Implied R^2 from the empirical proxy regression, draw by draw (ADDED 2026-09-17)
% Proposition 5.1 ("Feasible smoother regression: bound and identification"):
% the population R^2 of the FEASIBLE regression of the smoother revision
% f = xhat_{t+h|t+h+m} - xhat_{t+h|t} on the filter revision
% r = xhat_{t+h|t+h} - xhat_{t+h|t} is P^{2,m}_{i,h} = (V_h - Omega_ii)/(V_h - [Omega_m]_ii),
% while the INFEASIBLE reliability is R^2_{i,h} = 1 - Omega_ii/V_h. Eliminating
% V_h,
%     R^2_{i,h} = P (1 - rho) / (1 - rho P),    rho = [Omega_m]_ii / Omega_ii,
% (fr_implied_R2). Feeding the OLS estimate Phat^{2,m}_{i,h} (wrapper cell 20,
% fr.R2_emp) with the draw's own model-implied rho gives the DATA's estimate of
% the reliability of the filtered state, to be compared with the model's own
% R^2_{i,h} (cell 10, eq. (3.11) rows of the base states). Both are collected
% over the posterior draws, so the fans mix parameter uncertainty with the
% sampling variability of Phat; the model_assessment fans (model_assessment_plots)
% isolate the latter at the mode.
% rho comes from fr.omega_ratio where the results carry it (engine of
% 2026-09-17 or later) and otherwise from cell 6, the KF standard deviations of
% the key states: rho = (sd(l<m>NAME) / sd(NAME))^2.
if size(oo_.posterior_function_results, 2) >= 20 && any(cellfun(@isstruct, oo_.posterior_function_results(:, 20)))
    ok_fr = find(cellfun(@isstruct, oo_.posterior_function_results(:, 20)));
    fr0   = oo_.posterior_function_results{ok_fr(1), 20};
    Jfr = numel(fr0.states);  Hfr = numel(fr0.h);  n_mfr = numel(fr0.m);  Nd = numel(ok_fr);
    P_emp  = nan(Jfr, Hfr, n_mfr, Nd);  P_th = P_emp;
    rho_fr = nan(Jfr, n_mfr, Nd);
    for d = 1:Nd
        frd = oo_.posterior_function_results{ok_fr(d), 20};
        P_emp(:, :, :, d) = frd.R2_emp;
        P_th(:, :, :, d)  = frd.R2_theory;
        if isfield(frd, 'omega_ratio'), rho_fr(:, :, d) = frd.omega_ratio; end
    end
    std6  = cat(2, oo_.posterior_function_results{ok_fr, 6});      % n_key x Nd, sqrt(diag(KF_Omega))
    R2h10 = cat(3, oo_.posterior_function_results{ok_fr, 10});     % n_key x MaxH x Nd, eq. (3.17)/(3.11)
    R2_model = nan(Jfr, Hfr, Nd);
    use_cell6 = all(isnan(rho_fr(:)));                             % results predate fr.omega_ratio
    for j = 1:Jfr
        kb = strcmp(M_.states_of_interest, fr0.states{j});
        if ~any(kb), continue; end
        R2_model(j, :, :) = R2h10(kb, 1:Hfr, :);
        if use_cell6
            for im = 1:n_mfr
                kl = strcmp(M_.states_of_interest, sprintf('l%d%s', fr0.m(im), fr0.states{j}));
                if any(kl), rho_fr(j, im, :) = (std6(kl, :) ./ std6(kb, :)).^2; end
            end
        end
    end
    rho_full = repmat(reshape(rho_fr, [Jfr 1 n_mfr Nd]), [1 Hfr 1 1]);
    R2_impl  = fr_implied_R2(P_emp, rho_full);                    % J x H x n_m x Nd, from the data's Phat
    R2_check = fr_implied_R2(P_th,  rho_full);                    % identity: population P -> model R2_h exactly
    gap = R2_check - repmat(reshape(R2_model, [Jfr Hfr 1 Nd]), [1 1 n_mfr 1]);
    fprintf('\n=== Implied R^2 (Prop. 5.1) ===\nidentity check, implied(P_theory) vs model R2_h: max |gap| = %.2e over %d draws\n', max(abs(gap(:))), Nd);

    nu_med = median(nu_draws, 2);                                  % n_key x 1, quarters (d = 0 checked above)
    pcts = [5 16 50 84 95];
    fid_impl = fopen([HLW_directory '/Tables_level/table_implied_R2.tex'], 'w');
    fprintf(fid_impl, '%% R^2_{i,h} implied by the empirical proxy-regression R^2 through Prop. 5.1, R2 = P(1-rho)/(1-rho P),\n');
    fprintf(fid_impl, '%% rho = [Omega_m]_ii/Omega_ii, vs the model''s R^2_{i,h}; medians and [16,84] bands over %d posterior draws.\n', Nd);
    fprintf(fid_impl, '\\begin{tabular}{llr|cc|cc}\n');
    fprintf(fid_impl, 'state & $m$ & $h$ & implied $R^2_{i,h}$ & [16,84] & model $R^2_{i,h}$ & [16,84] \\\\\n\\hline\n');
    for j = 1:Jfr
        if all(isnan(P_emp(j, :, :, :))), continue; end        % degenerate regressor (zero filter gain)
        kb = strcmp(M_.states_of_interest, fr0.states{j});
        if ~DO_TABLES_ONLY
            figure("Name", ['implied R2 ' fr0.states{j}]);
            set(gcf, 'Position', [50 50 480*n_mfr 360]);
        end
        for im = 1:n_mfr
            qi = prctile(squeeze(R2_impl(j, :, im, :)),  pcts, 2);     % H x 5
            qm = prctile(squeeze(R2_model(j, :, :)),     pcts, 2);
            if ~DO_TABLES_ONLY
            subplot(1, n_mfr, im);
            hh_ = 1:Hfr;
            hold on;
            fill([hh_ fliplr(hh_)], [qi(:,1)' fliplr(qi(:,5)')], [0.18 0.44 0.72], 'FaceAlpha', 0.15, 'EdgeColor', 'none');
            fill([hh_ fliplr(hh_)], [qi(:,2)' fliplr(qi(:,4)')], [0.18 0.44 0.72], 'FaceAlpha', 0.30, 'EdgeColor', 'none');
            h_impl = plot(hh_, qi(:,3), '-',  'Color', [0.18 0.44 0.72], 'LineWidth', 2);
            h_mod  = plot(hh_, qm(:,3), '--', 'Color', 'k', 'LineWidth', 2);
            plot(hh_, qm(:,[2 4]), ':', 'Color', 'k', 'LineWidth', 1);
            hold off;
            xlim([0.5 Hfr+0.5]); ylim([0 1]); xticks([1 4:4:Hfr]); grid on; box on;
            xlabel('horizon h [quarters]'); ylabel('R^2_{i,h}');
            title(sprintf('%s, m = %d: median \\rho = %.2f, \\nu = %.1f qtrs', fr0.states{j}, fr0.m(im), median(rho_fr(j, im, :)), nu_med(kb)), 'Interpreter', 'tex');
            if im == 1
                legend([h_impl, h_mod], {'R^2 implied by empirical P^{2,m} (median, 68/90% bands)', 'model R^2_{i,h} (median, 68% band)'}, 'Location', 'southoutside', 'FontSize', 8);
                legend('boxoff');
            end
            end   % ~DO_TABLES_ONLY (implied-R2 panel)
            for hh = [4 8 20 40]
                if hh > Hfr, continue; end
                fprintf(fid_impl, '%s & %d & %d & %.2f & [%.2f, %.2f] & %.2f & [%.2f, %.2f] \\\\\n', strrep(fr0.states{j}, '_', '\_'), fr0.m(im), hh, ...
                        qi(hh, 3), qi(hh, 2), qi(hh, 4), qm(hh, 3), qm(hh, 2), qm(hh, 4));
            end
        end
        if ~DO_TABLES_ONLY, print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'impliedR2_' fr0.states{j} '.png']); end
        fprintf('%-8s median implied R2 at h = 4, 20, 40 (m = %d): %s ; model: %s ; nu = %.1f\n', fr0.states{j}, fr0.m(end), ...
                mat2str(prctile(squeeze(R2_impl(j, [4 20 min(40,Hfr)], end, :)), 50, 2)', 2), mat2str(prctile(squeeze(R2_model(j, [4 20 min(40,Hfr)], :)), 50, 2)', 2), nu_med(kb));
    end
    fprintf(fid_impl, '\\end{tabular}\n');
    fclose(fid_impl);
    close all
else
    fprintf(2, 'HLW_plots: results have no cell 20 (fr struct) - implied-R2 section skipped.\n');
end

%% Plot FEV
close all

if ~DO_MAIN_PLOTS_ONLY
  % the challenge here is that to preserve space, the FEC object lives on a
  % coarse grid, given by True_FEV_Horizon_GRID.
  PlotHorizonGrid = [1 2 4 8 12 24 50:50:150 250 ];
  True_FEV_Horizon_GRID = [1:100 150:50:1e3]; % THIS IS HARDWIRED IN f_extract_ABCD!!!
  for pp=1:length(PlotVar.Names)

      figure("Name",PlotVar.Names{pp})

      my_hndl=semilogx(True_FEV_Horizon_GRID(True_FEV_Horizon_GRID<=PlotHorizonGrid(end)), ...
          prctile(squeeze(state_estimates.FEV(PlotVar.Idx{pp},1:length(True_FEV_Horizon_GRID(True_FEV_Horizon_GRID<=PlotHorizonGrid(end))),:)), [ 16 50 84 ], 2)*100, 'LineWidth', 2, 'Color', 'b');
      set(my_hndl, {'LineStyle'}, {'--';'-';'--';});
      xlabel('forecast horizon [quarters, log-scale]');
      ylabel('Variance explained by initial state [%]');
      xticks(PlotHorizonGrid); % Set specific x-axis ticks
      xtickangle(0);
      grid on;
      hold on

      my_prior_hndl=semilogx(True_FEV_Horizon_GRID(True_FEV_Horizon_GRID<=PlotHorizonGrid(end)), ...
          prctile(squeeze(state_estimates.prior_FEV(PlotVar.Idx{pp},1:length(True_FEV_Horizon_GRID(True_FEV_Horizon_GRID<=PlotHorizonGrid(end))),:)), [ 16 50 84 ], 2)*100, 'LineWidth', 2, 'Color', 'r');
      set(my_prior_hndl, {'LineStyle'}, {'--';'-';'--'});
      xlabel('forecast horizon [quarters, log-scale]');
      ylabel('Variance explained by initial state [%]');
      grid on;
      hold on

      legend([my_hndl(2), my_hndl(1), my_prior_hndl(2), my_prior_hndl(1)], ...
             'Posterior (median)', 'Posterior (68% credible set)', ...
             'Prior (median)', 'Prior (68% credible set)', ...
             'Location', 'best');

      legend('boxoff')

      fig = gcf; % or your figure handle
      set(fig, 'ToolBar', 'none');     % (was fig.ToolBar = 'none'; handle dot-access is MATLAB-only)
      print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'FEV_' PlotVar.Names{pp} '_prior.png']);
      delete(my_prior_hndl);
      print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'FEV_' PlotVar.Names{pp} '.png']);
      % stop;
  end
end
%%
close all;

% RMSE of the one-step-ahead forecast errors per observable and draw. Used by
% the scatter block below AND by the parameter-correlation tables further down,
% which run regardless of DO_MAIN_PLOTS_ONLY, so it is computed outside the guard.
surprises.RMSE = sqrt(squeeze(mean(surprises.ThreeD.^2, 1,'omitnan')));

if ~DO_MAIN_PLOTS_ONLY

  var1 = cfg.sd_scatter_var;
  pp=find(strcmp(PlotVar.Names, var1{1,1}));
  if ~isempty(pp)
    hh=20;
    figure;
    scatter(squeeze(R2.ThreeD(PlotVar.Idx{pp},hh,:)),squeeze(state_estimates.KF_Std(PlotVar.Idx{pp},end,:)))
    grid on;
    xlabel('R^2 | horizon 20');
    ylabel('s.d. of state estimate');

    print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'scatter_R2_SD_' PlotVar.Names{pp} '.png']);

    figure;
    scatter(squeeze(R2.ThreeD(PlotVar.Idx{pp},hh,:)),Error_ACF.half_life(PlotVar.Idx{pp},:))
    grid on;
    xlabel('R^2 | horizon 20');
    ylabel('half-life');
    % uncomment below
    %print('-dpng', [HLW_directory ['/Graphs_level/' cfg.fig_prefix 'scatter_R2_HalfLife_' PlotVar.Names{pp} '.png']);
  end

  %%
  % hh=20;
  %
  % figure("Name",PlotVar.Names{pp})
  %
  % low_corr_idx = squeeze(CorrH.ThreeD(PlotVar.Idx{pp},hh,:))<0.25; %prctile(squeeze(CorrH.ThreeD(PlotVar.Idx{pp},hh,:)),25);
  % high_corr_idx = squeeze(CorrH.ThreeD(PlotVar.Idx{pp},hh,:))>0.5;%prctile(squeeze(CorrH.ThreeD(PlotVar.Idx{pp},hh,:)),50);
  %
  % % surprises.RMSE= mean(surprises.ThreeD.^2, 0, 3);
  % (surprises.RMSE is computed above the guard)

      % state_estimates.mean_low_corr = mean(state_estimates.ThreeD(:,:,low_corr_idx), 3);
      % state_estimates.Avg_KF_Var_low_corr = mean(state_estimates.KF_Std(:,low_corr_idx).^2, 2);
      % state_estimates.StdOfMeans_low_corr= std(state_estimates.ThreeD(:,:,low_corr_idx), 0, 3);
      % state_estimates.OverallStd_low_corr = sqrt(ones(T,1)*state_estimates.Avg_KF_Var_low_corr' + state_estimates.StdOfMeans_low_corr.^2);
      %
      % state_estimates.mean_high_corr = mean(state_estimates.ThreeD(:,:,high_corr_idx), 3);
      % state_estimates.Avg_KF_Var_high_corr = mean(state_estimates.KF_Std(:,high_corr_idx).^2, 2);
      % state_estimates.StdOfMeans_high_corr= std(state_estimates.ThreeD(:,:,high_corr_idx), 0, 3);
      % state_estimates.OverallStd_high_corr = sqrt(ones(T,1)*state_estimates.Avg_KF_Var_high_corr' + state_estimates.StdOfMeans_high_corr.^2);
      %
      %
      %
      % my_hndl_low=plot(TimeLabels-PlotLag/4, [state_estimates.q_lo_low_corr(:,PlotVar.Idx{pp}), state_estimates.mean_low_corr(:,PlotVar.Idx{pp}), state_estimates.q_hi_low_corr(:,PlotVar.Idx{pp})], 'Color','b', 'LineWidth',2);
      % set(my_hndl_low, {'LineStyle'}, {'--';'-';'--'})
      % grid on;
      % hold all
      %
      % my_hndl_high=plot(TimeLabels-PlotLag/4, [state_estimates.q_lo_high_corr(:,PlotVar.Idx{pp}), state_estimates.mean_high_corr(:,PlotVar.Idx{pp}), state_estimates.q_hi_high_corr(:,PlotVar.Idx{pp})], 'Color','r', 'LineWidth',2);
      % set(my_hndl_high, {'LineStyle'}, {'--';'-';'--'})
      % grid on;
      %
      %
      % legend([my_hndl_low(2), my_hndl_high(2)], {'corr<0.25'; 'corr>0.50'})
      %
      % print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'CorrH_Tracking_' PlotVar.Names{pp} '.png']);
      %
      % delete(my_hndl_low([1 3]))
      % delete(my_hndl_high([1 3]))
      % my_hndl=plot(TimeLabels-PlotLag/4, state_estimates.mean(:,PlotVar.Idx{pp}), '--k', 'LineWidth',2);
      % % set(my_hndl_high, {'LineStyle'}, {'--';'-';'--'})
      % grid on;
      % legend({'corr<0.25'; 'corr>0.50'; 'unconditional'},'location','best')
      % legend('boxoff');
      %
      % print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'CorrH_Tracking_' PlotVar.Names{pp} '_mean_only.png']);

  %% Correlation of reliability measures
  % close all
  %
  % var1 = {'rrstar', 'r*'};
  % var2 = {'ygap', 'y-gap'};
  %
  % pp=find(strcmp(PlotVar.Names, var1{1,1}));
  % qq=find(strcmp(PlotVar.Names, var2{1,1}));
  %
  % for hh=[1 4 10 20 50]
  %     figure;
  %     scatter(squeeze(CorrH.ThreeD(PlotVar.Idx{pp},hh,:)),squeeze(CorrH.ThreeD(PlotVar.Idx{qq},hh,:)))
  %     grid on;
  %     xlabel([var1{1,2} ' correlation with truth | horizon ' num2str(hh)]);
  %     ylabel([var2{1,2} ' correlation with truth | horizon ' num2str(hh)]);
  %
  %     print('-dpng', ['../../datastore/LaubachWilliams/Graphs_level/LW_KF_scatter_CorrH_' var1{1,1} '_' PlotVar.Names{qq} '_h' num2str(hh) '.png']);
  % end
  %
  % figure;
  % scatter((Error_ACF.half_life(PlotVar.Idx{pp},:)),(Error_ACF.half_life(PlotVar.Idx{qq},:)))
  % grid on;
  % xlabel(['half-life ' var1{1,2}]);
  % ylabel(['half-life ' var2{1,2}]);
  %
  % print('-dpng', ['../../datastore/LaubachWilliams/Graphs_level/LW_KF_scatter_HalfLife_' var1{1,1} '_' PlotVar.Names{qq} '.png']);


  %% Correlation of reliability measures
  close all

  for vv_idx = 1:size(cfg.scatter_pairs, 1)
      var1 = cfg.scatter_pairs{vv_idx, 1};
      var2 = cfg.scatter_pairs{vv_idx, 2};
      pp=find(strcmp(PlotVar.Names, var1{1,1}));
      qq=find(strcmp(PlotVar.Names, var2{1,1}));
      if isempty(pp) || isempty(qq), continue; end
      for hh=[1 4 10 20 50]
          figure;
          scatter(squeeze(R2.ThreeD(PlotVar.Idx{pp},hh,:)),squeeze(R2.ThreeD(PlotVar.Idx{qq},hh,:)))
          grid on;
          xlabel([var1{1,2} ' R^2 | horizon ' num2str(hh)]);
          ylabel([var2{1,2} ' R^2 | horizon ' num2str(hh)]);
          xlim([0 1]); ylim([0 1]);
          hold all;
          my_prior_hndl=scatter(squeeze(R2.prior_ThreeD(PlotVar.Idx{pp},hh,:)),squeeze(R2.prior_ThreeD(PlotVar.Idx{qq},hh,:)), 'MarkerEdgeColor', actual_color);

          print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'scatter_R2_' var1{1,1} '_' PlotVar.Names{qq} '_h' num2str(hh) '_prior.png']);
          delete(my_prior_hndl);
          print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'scatter_R2_' var1{1,1} '_' PlotVar.Names{qq} '_h' num2str(hh) '.png']);
      end

      figure;   % infinite half-lives (never cross 1/2) shown at the cap of 50, as in the histograms (2026-09-30)
      scatter(min(Error_ACF.half_life(PlotVar.Idx{pp},:), 50),min(Error_ACF.half_life(PlotVar.Idx{qq},:), 50))
      grid on;
      xlabel(['half-life ' var1{1,2} ' (50 = never crosses 1/2)']);
      ylabel(['half-life ' var2{1,2} ' (50 = never crosses 1/2)']);
      hold all;
      my_prior_hndl=scatter(min(Error_ACF.prior_half_life(PlotVar.Idx{pp},:), 50),min(Error_ACF.prior_half_life(PlotVar.Idx{qq},:), 50), 'MarkerEdgeColor', actual_color);
      print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'scatter_HalfLife_' var1{1,1} '_' PlotVar.Names{qq} '_prior.png']);
      delete(my_prior_hndl);
      print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'scatter_HalfLife_' var1{1,1} '_' PlotVar.Names{qq} '.png']);
  end



  %% RMSE plots
  for vv_idx = 1:numel(cfg.rmse_vars)
      var1 = cfg.rmse_vars{vv_idx};
      pp=find(strcmp(PlotVar.Names, var1{1,1}));
      if isempty(pp), continue; end

      % R^2 vs 1-step-ahead RMSE (the title's corr is the correlation across draws of the two)
      hh=1;
      figure;
      for ii=1:2
          subplot(1,2,ii)
          scatter(squeeze(R2.ThreeD(PlotVar.Idx{pp},hh,:)),surprises.RMSE(ii,:))
          title([M_.observable_names{ii} ': corr=' num2str(corr(squeeze(R2.ThreeD(PlotVar.Idx{pp},hh,:)),surprises.RMSE(ii,:)'), '%1.2f')])
          xlabel('R^2 of the state estimate | horizon 1')
          ylabel('RMSE')
          grid on;
          print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'scatter_R2_RMSE_' var1{1,1}  '.png']);
      end
      figure;
      for ii=1:2
          subplot(1,2,ii)
          % finite half-lives only: Inf (never crosses 1/2) has no log (2026-09-30)
          hl_pp = Error_ACF.half_life(PlotVar.Idx{pp},:);  fin_pp = isfinite(hl_pp);  pos_pp = fin_pp & hl_pp>1;
          scatter(log(hl_pp(fin_pp)),surprises.RMSE(ii,fin_pp))
          full_corr = NaN;  pos_HL_corr = NaN;                                   % NaN when fewer than 3 draws qualify
          if sum(fin_pp) >= 3, full_corr=corr(log(hl_pp(fin_pp))',surprises.RMSE(ii,fin_pp)'); end
          if sum(pos_pp) >= 3, pos_HL_corr=corr(log(hl_pp(pos_pp))',surprises.RMSE(ii,pos_pp)'); end

          title([M_.observable_names{ii} ': corr=' num2str(full_corr, '%1.2f') '. corr|HL>1=' num2str(pos_HL_corr, '%1.2f.')])
          xlabel('half-life [log]')
          ylabel('RMSE')
          grid on;
          print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'scatter_HalfLife_RMSE_' var1{1,1}  '.png']);
      end
  end
end

%% Tracking measures(=filter estimates) given restricted reliability
if ~DO_TABLES_ONLY
for vv = 1:numel(cfg.tracking_vars)
    var1 = cfg.tracking_vars{vv};
    if ~any(strcmp(PlotVar.Names, var1{1,1})), continue; end
    for hh=4:4:20


        pp=find(strcmp(PlotVar.Names, var1{1,1}));
        close all
        figure("Name",PlotVar.Names{pp})

        low_corr_idx = squeeze(CorrH.ThreeD(PlotVar.Idx{pp},hh,:)).^2<0.25; %prctile(squeeze(CorrH.ThreeD(PlotVar.Idx{pp},hh,:)),25);
        high_corr_idx = squeeze(CorrH.ThreeD(PlotVar.Idx{pp},hh,:)).^2>0.5;%prctile(squeeze(CorrH.ThreeD(PlotVar.Idx{pp},hh,:)),50);

        % surprises.RMSE= mean(surprises.ThreeD.^2, 0, 3);
        % surprises.RMSE = sqrt(squeeze(mean(surprises.ThreeD.^2, 1)));

        state_estimates.mean_low_corr = mean(state_estimates.ThreeD(:,:,low_corr_idx), 3);
        state_estimates.Avg_KF_Var_low_corr = mean(state_estimates.KF_Std(:,:,low_corr_idx).^2, 3);
        state_estimates.StdOfMeans_low_corr= std(state_estimates.ThreeD(:,:,low_corr_idx), 0, 3);
        % state_estimates.OverallStd_low_corr = sqrt(ones(T,1)*state_estimates.Avg_KF_Var_low_corr' + state_estimates.StdOfMeans_low_corr.^2);
        state_estimates.OverallStd_low_corr = sqrt(state_estimates.Avg_KF_Var_low_corr' + state_estimates.StdOfMeans_low_corr.^2);
        [state_estimates.q_lo_low_corr, state_estimates.q_hi_low_corr] = deal(nan(size(state_estimates.mean)));
        [state_estimates.q_lo_low_corr(:,PlotVar.Idx{pp}), state_estimates.q_hi_low_corr(:,PlotVar.Idx{pp})] = mixture_band(squeeze(state_estimates.ThreeD(:,PlotVar.Idx{pp},low_corr_idx)), squeeze(state_estimates.KF_Std(PlotVar.Idx{pp},:,low_corr_idx)));

        state_estimates.mean_high_corr = mean(state_estimates.ThreeD(:,:,high_corr_idx), 3);
        state_estimates.Avg_KF_Var_high_corr = mean(state_estimates.KF_Std(:,:,high_corr_idx).^2, 3);
        state_estimates.StdOfMeans_high_corr= std(state_estimates.ThreeD(:,:,high_corr_idx), 0, 3);
        state_estimates.OverallStd_high_corr = sqrt(state_estimates.Avg_KF_Var_high_corr' + state_estimates.StdOfMeans_high_corr.^2);
        [state_estimates.q_lo_high_corr, state_estimates.q_hi_high_corr] = deal(nan(size(state_estimates.mean)));
        [state_estimates.q_lo_high_corr(:,PlotVar.Idx{pp}), state_estimates.q_hi_high_corr(:,PlotVar.Idx{pp})] = mixture_band(squeeze(state_estimates.ThreeD(:,PlotVar.Idx{pp},high_corr_idx)), squeeze(state_estimates.KF_Std(PlotVar.Idx{pp},:,high_corr_idx)));


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot with three sets of lines
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        my_hndl_low=plot(TimeLabels-PlotLag/4, [state_estimates.q_lo_low_corr(:,PlotVar.Idx{pp}), state_estimates.mean_low_corr(:,PlotVar.Idx{pp}), state_estimates.q_hi_low_corr(:,PlotVar.Idx{pp})], 'Color','b', 'LineWidth',2);
        set(my_hndl_low, {'LineStyle'}, {'--';'-';'--'})
        grid on;
        hold all

        my_hndl_high=plot(TimeLabels-PlotLag/4, [state_estimates.q_lo_high_corr(:,PlotVar.Idx{pp}), state_estimates.mean_high_corr(:,PlotVar.Idx{pp}), state_estimates.q_hi_high_corr(:,PlotVar.Idx{pp})], 'Color','r', 'LineWidth',2);
        set(my_hndl_high, {'LineStyle'}, {'--';'-';'--'})
        grid on;

        % my_hndl=plot(TimeLabels-PlotLag/4, state_estimates.mean(:,PlotVar.Idx{pp}), '--k', 'LineWidth',2);
        my_hndl=plot(TimeLabels-PlotLag/4, [state_estimates.q_lo(:,PlotVar.Idx{pp}), state_estimates.mean(:,PlotVar.Idx{pp}), state_estimates.q_hi(:,PlotVar.Idx{pp})], 'Color','k', 'LineWidth',2);
        set(my_hndl_high, {'LineStyle'}, {'--';'-';'--'})
        grid on;

        % AutoUpdate off: lines are deleted and re-plotted below with the legend in place (Octave's auto-update callback errors on that)
        legend([my_hndl_low(2), my_hndl_high(2), my_hndl(2)], {'R-sq<0.25'; 'R-sq>0.50'; 'unconditional'}, 'AutoUpdate', 'off')

        print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'CorrHH' num2str(hh) '_Tracking_' PlotVar.Names{pp} '.png']);

        delete(my_hndl_low([1 3]))
        delete(my_hndl_high([1 3]))
        delete(my_hndl([1 3]))
        my_hndl=plot(TimeLabels-PlotLag/4, state_estimates.mean(:,PlotVar.Idx{pp}), '--k', 'LineWidth',2);
        % set(my_hndl_high, {'LineStyle'}, {'--';'-';'--'})
        grid on;
        legend({'R-sq<0.25'; 'R-sq>0.50';  'unconditional'},'location','best')
        legend('boxoff');

        print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'CorrH' num2str(hh) '_Tracking_' PlotVar.Names{pp} '_mean_only.png']);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot with two sets of lines
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        close all
        figure("Name",PlotVar.Names{pp})
        my_hndl_high=plot(TimeLabels-PlotLag/4, [state_estimates.q_lo_high_corr(:,PlotVar.Idx{pp}), state_estimates.mean_high_corr(:,PlotVar.Idx{pp}), state_estimates.q_hi_high_corr(:,PlotVar.Idx{pp})], 'Color','r', 'LineWidth',2);
        set(my_hndl_high, {'LineStyle'}, {'--';'-';'--'})
        grid on; hold('on');

        % my_hndl=plot(TimeLabels-PlotLag/4, state_estimates.mean(:,PlotVar.Idx{pp}), '--k', 'LineWidth',2);
        my_hndl=plot(TimeLabels-PlotLag/4, [state_estimates.q_lo(:,PlotVar.Idx{pp}), state_estimates.mean(:,PlotVar.Idx{pp}), state_estimates.q_hi(:,PlotVar.Idx{pp})], 'Color','k', 'LineWidth',2);
        set(my_hndl_high, {'LineStyle'}, {'--';'-';'--'})
        grid on;

        legend([my_hndl_high(2), my_hndl(2)], {'R-sq>0.50'; 'unconditional'}, 'AutoUpdate', 'off')

        print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'CorrHH' num2str(hh) '_SimpleTrack_' PlotVar.Names{pp} '.png']);

        % (my_hndl_low belongs to the previous figure and is already deleted; deleting it again errors in Octave)
        delete(my_hndl_high([1 3]))
        delete(my_hndl([1 3]))
        my_hndl=plot(TimeLabels-PlotLag/4, state_estimates.mean(:,PlotVar.Idx{pp}), '--k', 'LineWidth',2);
        % set(my_hndl_high, {'LineStyle'}, {'--';'-';'--'})
        grid on;
        % two lines only (FIXED 2026-09-18: the three-entry legend mislabelled the red line as R-sq<0.25)
        legend([my_hndl_high(2), my_hndl], {'R-sq>0.50'; 'unconditional'},'location','best')
        legend('boxoff');

        print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'CorrH' num2str(hh) '_SimpleTrack_' PlotVar.Names{pp} '_mean_only.png']);

        % post-1990 versions: with bands (as before) and mean only (ADDED 2026-09-18, paper Fig. tracking:reliability panel b)
        figure("Name",PlotVar.Names{pp})
        IDX=(TimeLabels-PlotLag/4>=1990);
        my_hndl_high=plot(TimeLabels(IDX)-PlotLag/4, [state_estimates.q_lo_high_corr(IDX,PlotVar.Idx{pp}), state_estimates.mean_high_corr(IDX,PlotVar.Idx{pp}), state_estimates.q_hi_high_corr(IDX,PlotVar.Idx{pp})], 'Color','r', 'LineWidth',2);
        set(my_hndl_high, {'LineStyle'}, {'--';'-';'--'})
        grid on; hold all;

        % my_hndl=plot(TimeLabels-PlotLag/4, state_estimates.mean(:,PlotVar.Idx{pp}), '--k', 'LineWidth',2);
        my_hndl=plot(TimeLabels(IDX)-PlotLag/4, [state_estimates.q_lo(IDX,PlotVar.Idx{pp}), state_estimates.mean(IDX,PlotVar.Idx{pp}), state_estimates.q_hi(IDX,PlotVar.Idx{pp})], 'Color','k', 'LineWidth',2);
        set(my_hndl_high, {'LineStyle'}, {'--';'-';'--'})
        grid on;

        legend([my_hndl_high(2), my_hndl(2)], {'R-sq>0.50'; 'unconditional'}, 'AutoUpdate', 'off')

        print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'CorrHH' num2str(hh) '_SimpleTrack_1990_' PlotVar.Names{pp} '.png']);

        delete(my_hndl_high([1 3]))
        delete(my_hndl([1 3]))
        my_hndl=plot(TimeLabels(IDX)-PlotLag/4, state_estimates.mean(IDX,PlotVar.Idx{pp}), '--k', 'LineWidth',2);
        grid on;
        legend([my_hndl_high(2), my_hndl], {'R-sq>0.50'; 'unconditional'},'location','best')
        legend('boxoff');

        print('-dpng', [HLW_directory '/Graphs_level/' cfg.fig_prefix 'CorrH' num2str(hh) '_SimpleTrack_1990_' PlotVar.Names{pp} '_mean_only.png']);
    end
end
end   % ~DO_TABLES_ONLY (tracking panels)

%% Filtered r* vs the real rate in the data (ADDED 2026-09-19)
% Used to be the stand-alone Miscellaneous_Archive/JWright_plot.m, run by hand on
% the plot workspace, so no pipeline mode produced the paper's
% LW_rstar_vs_realrate_1990_bands.png. Three files: full sample with the 95%
% band, full sample mean only, 1990- with the band and fixed axes. The band is
% OverallStd (filter uncertainty plus dispersion of the means across draws).
if ~DO_TABLES_ONLY && ~isempty(cfg.rstar_vs_realrate)
    rv_state = cfg.rstar_vs_realrate{1};  rv_obs = cfg.rstar_vs_realrate{2};
    pp_rv = find(strcmp(PlotVar.Names, rv_state), 1);
    nn_rv = find(strcmp(M_.observable_names, rv_obs), 1);
    if isempty(pp_rv) || isempty(nn_rv)
        fprintf('HLW_plots: r*-vs-real-rate figures skipped (%s / %s not in this model)\n', rv_state, rv_obs);
    else
        if strcmp(rv_state, 'rrstar'), rv_label = 'r^*'; else, rv_label = rv_state; end
        rv_mean = state_estimates.mean(:, PlotVar.Idx{pp_rv});
        rv_lo   = state_estimates.q_lo(:, PlotVar.Idx{pp_rv});
        rv_hi   = state_estimates.q_hi(:, PlotVar.Idx{pp_rv});
        rv_real = M_.data(:, nn_rv);
        T_rv    = min([numel(rv_mean), numel(rv_real), numel(TimeLabels)]);
        t_rv    = reshape(TimeLabels(1:T_rv), 1, []);
        rv_mean = reshape(rv_mean(1:T_rv), 1, []);  rv_lo = reshape(rv_lo(1:T_rv), 1, []);  rv_hi = reshape(rv_hi(1:T_rv), 1, []);  rv_real = reshape(rv_real(1:T_rv), 1, []);
        % 20-quarter-lag smoother of the same state (ADDED 2026-09-30): posterior mean of x(t|t+20),
        % which the engine stores as l20<state> dated t+20, hence plotted at TimeLabels - 5 years like
        % the l20 panels above (same presample mask); black dashed line. Skipped if not diagnosed.
        pp_rv20 = find(strcmp(PlotVar.Names, ['l20' rv_state]), 1);
        if ~isempty(pp_rv20)
            PlotVar.Idx{pp_rv20} = strcmp(M_.states_of_interest, PlotVar.Engine{pp_rv20});
            rv_sm  = state_estimates.mean(:, PlotVar.Idx{pp_rv20});
            if cfg.lag_presample_mask, rv_sm(1:min(20, numel(rv_sm))) = NaN; end
            T_rv20 = min(numel(rv_sm), numel(TimeLabels));
            t_rv20 = reshape(TimeLabels(1:T_rv20), 1, []) - 5;  rv_sm = reshape(rv_sm(1:T_rv20), 1, []);
            rv_leg3 = ['Smoothed ' rv_label ' (20-qtr lag, posterior mean)'];
        end
        rv_files = {'LW_rstar_vs_realrate.png', 'LW_rstar_vs_realrate_simple.png', 'LW_rstar_vs_realrate_1990_bands.png'};
        for rv_k = 1:3
            close all
            figure('Units', 'inches', 'Position', [1 1 10 5]);
            if rv_k ~= 2
                fill([t_rv fliplr(t_rv)], [rv_hi, fliplr(rv_lo)], [0.7 0.7 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                hold on;
                rv_leg1 = ['Filtered ' rv_label ' (posterior mean, 95% band)'];
            else
                rv_leg1 = ['Filtered ' rv_label ' (posterior mean)'];
            end
            h_rv1 = plot(t_rv, rv_mean, 'b-', 'LineWidth', 2); hold on;
            h_rv2 = plot(t_rv, rv_real, 'r-', 'LineWidth', 1.5);
            rv_h = [h_rv1, h_rv2];  rv_leg = {rv_leg1, 'Real rate (data)'};
            if ~isempty(pp_rv20)
                h_rv3 = plot(t_rv20, rv_sm, 'k--', 'LineWidth', 1.5);
                rv_h = [rv_h, h_rv3];  rv_leg = [rv_leg, {rv_leg3}];
            end
            grid on; xlabel('Year', 'FontSize', 16); ylabel('Percent', 'FontSize', 16); set(gca, 'FontSize', 16);
            if rv_k == 3, xlim([1990, max(t_rv)]); ylim([-7, 9]); else, xlim([min(t_rv), max(t_rv)]); end
            if rv_k == 3, rv_loc = 'southwest'; else, rv_loc = 'best'; end   % 1990 panel: the lower-left corner is empty in every HLW variant; 'best' put the 3-entry legend on the lines
            legend(rv_h, rv_leg, 'Location', rv_loc, 'FontSize', 16, 'Interpreter', 'tex');
            legend('boxoff');
            print('-dpng', [HLW_directory '/Graphs_level/' rv_files{rv_k}]);
        end
        fprintf('HLW_plots: r*-vs-real-rate figures written (%s vs %s)\n', rv_state, rv_obs);
    end
end

%% Correlations of reliability measures with parameters
target_params= {'rho', 'sum of lag coeffs'; ...
                'a1', '1st lag coeffs'; ...
                'a_r', 'r-y'; ...
                'bpi', 'lagged inflation'; ...
                'by', 'slope PC'; ...
                'phi', 'COVID-loading'; ...
                'se_ygap', 'SE y-gap'
                'se_pi', 'SE infl'   ; ...
                'se_ystar', 'SE potential y'; ...
                'se_zz', 'SE z'   ; ...
                'se_gg', 'SE growth'   ; ...
                'se_rr', 'SE rr'   ; ...
                };

% The list above is the HLW default kept for reference; the one actually used is
% cfg.target_params restricted to the parameters this model has (GENERALISED 2026-09-17)
target_params   = cfg.target_params(ismember(cfg.target_params(:,1), M_.param_names), :);
n_target_params = size(target_params,1);

Horizons_Grid = [1 4 10 20 50];
n_horizons = length(Horizons_Grid);

% Create column labels
horizons_labels = arrayfun(@(x) sprintf('H=%d', x), Horizons_Grid, 'UniformOutput', false);
halflife_labels = {'unconditional', 'half-life$>$1'};
rmse_labels = {['1-step RMSE ' M_.observable_names{1}], ['1-step RMSE ' M_.observable_names{2}]};   % first two observables (was 'output growth', 'inflation')

% Function to format numbers to 2 significant digits
format_2sig = @(x) sprintf('%1.2f', x);

% Calculate RMSE correlations (independent of vv_idx)
corr_matrix_RMSE = nan(n_target_params, 2);
for tt=1:size(target_params, 1)
    pp_param = find(strcmp(M_.param_names, target_params{tt,1}));
    corr_matrix_RMSE(tt,1) = corr(surprises.RMSE(1,:)',param_estimates.matrix(pp_param,:)');
    corr_matrix_RMSE(tt,2) = corr(surprises.RMSE(2,:)',param_estimates.matrix(pp_param,:)');
end

% Display and write RMSE table (since it's independent of variables)
fprintf('\n=== RMSE Correlations ===\n');
row_labels = target_params(:,2);

% Format the matrix for display (array2table is MATLAB-only; unused below, so skipped where absent)
try
    corr_matrix_RMSE_formatted = array2table(corr_matrix_RMSE, 'RowNames', row_labels, 'VariableNames', rmse_labels);
catch
end
% Apply formatting to displayed values
for i = 1:size(corr_matrix_RMSE, 1)
    for j = 1:size(corr_matrix_RMSE, 2)
        fprintf('%s: %s = %s\n', row_labels{i}, rmse_labels{j}, format_2sig(corr_matrix_RMSE(i,j)));
    end
end

% Write RMSE table to TeX file
tex_file_rmse = fopen([HLW_directory '/Tables_level/table_rmse_correlations.tex'], 'w');
%fprintf(tex_file_rmse, '\\documentclass{article}\n\\usepackage{booktabs}\n\\begin{document}\n\n');
% fprintf(tex_file_rmse, '\\begin{table}[htbp]\n\\centering\n');
% fprintf(tex_file_rmse, '\\caption{RMSE Correlations}\n');
fprintf(tex_file_rmse, '\\begin{tabular}{lcc}\n');
% fprintf(tex_file_rmse, '\\toprule\n');
fprintf(tex_file_rmse, 'Parameter & %s & %s \\\\\n\\midrule\n', strrep(rmse_labels{1}, '_', '\_'), strrep(rmse_labels{2}, '_', '\_'));

for i = 1:n_target_params
    fprintf(tex_file_rmse, '%s & %s & %s \\\\\n', target_params{i,2}, ...
        format_2sig(corr_matrix_RMSE(i,1)), format_2sig(corr_matrix_RMSE(i,2)));
end
% fprintf(tex_file_rmse, '\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');
fprintf(tex_file_rmse, '\\end{tabular}\n');
% fprintf(tex_file_rmse, '\\end{document}\n');
fclose(tex_file_rmse);

for vv_idx = 1:numel(cfg.rmse_vars)
    var1 = cfg.rmse_vars{vv_idx};
    pp=find(strcmp(PlotVar.Names, var1{1,1}));
    if isempty(pp), continue; end

    corr_matrix_CorrH = nan(n_target_params, n_horizons);
    corr_matrix_HalfLife = nan(n_target_params, 2);

    for tt=1:size(target_params, 1)
        pp_param = find(strcmp(M_.param_names, target_params{tt,1}));
        for hh=Horizons_Grid
            corr_matrix_CorrH(tt,Horizons_Grid==hh) = corr(squeeze(R2.ThreeD(PlotVar.Idx{pp},hh,:)),param_estimates.matrix(pp_param,:)');   % R^2 (not the correlation) vs parameters since 2026-09-19
        end
        hl_pp = Error_ACF.half_life(PlotVar.Idx{pp},:);  fin_pp = isfinite(hl_pp);  pos_pp = fin_pp & hl_pp>1;   % finite half-lives only (2026-09-30)
        full_corr = NaN;  pos_HL_corr = NaN;                                                                 % NaN when fewer than 3 draws qualify (SPF r*: every draw has half-life 1)
        if sum(fin_pp) >= 3, full_corr=corr(log(hl_pp(fin_pp))',param_estimates.matrix(pp_param,fin_pp)'); end
        if sum(pos_pp) >= 3, pos_HL_corr=corr(log(hl_pp(pos_pp))',param_estimates.matrix(pp_param,pos_pp)'); end

        corr_matrix_HalfLife(tt,:) = [full_corr, pos_HL_corr];
    end

    % Create tables with labels
    row_labels = target_params(:,2);

    % Create and display tables
    fprintf('\n=== %s ===\n', var1{1,2});

    % Table 1: Correlation by Horizons
    fprintf('\nCorrelation by Horizons:\n');
    fprintf('%-20s', 'Parameter');
    for i = 1:n_horizons
        fprintf('%10s', horizons_labels{i});
    end
    fprintf('\n');
    for i = 1:n_target_params
        fprintf('%-20s', row_labels{i});
        for j = 1:n_horizons
            fprintf('%10s', format_2sig(corr_matrix_CorrH(i,j)));
        end
        fprintf('\n');
    end

    % Table 2: Half-Life Correlations
    fprintf('\nHalf-Life Correlations:\n');
    fprintf('%-20s%15s%15s\n', 'Parameter', 'unconditional', 'half-life>1');
    for i = 1:n_target_params
        fprintf('%-20s%15s%15s\n', row_labels{i}, ...
            format_2sig(corr_matrix_HalfLife(i,1)), format_2sig(corr_matrix_HalfLife(i,2)));
    end

    % Create variable name for filenames (remove special characters)
    var_name_clean = strrep(var1{1,1}, '*', 'star');

    % Write Table 1 to separate TeX file
    filename_corrh = sprintf([HLW_directory '/Tables_level/table_corrh_%s.tex'], var_name_clean);
    tex_file_corrh = fopen(filename_corrh, 'w');
    % fprintf(tex_file_corrh, '\\documentclass{article}\n\\usepackage{booktabs}\n\\begin{document}\n\n');
    % fprintf(tex_file_corrh, '\\begin{table}[htbp]\n\\centering\n');
    % fprintf(tex_file_corrh, '\\caption{Correlation by Horizons: %s}\n', var1{1,2});
    fprintf(tex_file_corrh, '\\begin{tabular}{l%s}\n', repmat('c', 1, n_horizons));
    % fprintf(tex_file_corrh, '\\toprule\n');
    fprintf(tex_file_corrh, 'Parameter');
    for i = 1:n_horizons
        fprintf(tex_file_corrh, ' & H=%d', Horizons_Grid(i));
    end
    fprintf(tex_file_corrh, ' \\\\\n\\midrule\n');

    for i = 1:n_target_params
        fprintf(tex_file_corrh, '%s', target_params{i,2});
        for j = 1:n_horizons
            fprintf(tex_file_corrh, ' & %s', format_2sig(corr_matrix_CorrH(i,j)));
        end
        fprintf(tex_file_corrh, ' \\\\\n');
    end
    % fprintf(tex_file_corrh, '\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');
    fprintf(tex_file_rmse, '\\end{tabular}\n');
    % fprintf(tex_file_corrh, '\\end{document}\n');
    fclose(tex_file_corrh);

    % Write Table 2 to separate TeX file
    filename_halflife = sprintf([HLW_directory '/Tables_level/table_halflife_%s.tex'], var_name_clean);
    tex_file_halflife = fopen(filename_halflife, 'w');
    % fprintf(tex_file_halflife, '\\documentclass{article}\n\\usepackage{booktabs}\n\\begin{document}\n\n');
    % fprintf(tex_file_halflife, '\\begin{table}[htbp]\n\\centering\n');
    % fprintf(tex_file_halflife, '\\caption{Half-Life Correlations: %s}\n', var1{1,2});
    fprintf(tex_file_halflife, '\\begin{tabular}{lcc}\n');
    % fprintf(tex_file_halflife, '\\toprule\n');
    fprintf(tex_file_halflife, 'Parameter & unconditional & half-life$>$1 \\\\\n\\midrule\n');

    for i = 1:n_target_params
        fprintf(tex_file_halflife, '%s & %s & %s \\\\\n', target_params{i,2}, ...
            format_2sig(corr_matrix_HalfLife(i,1)), format_2sig(corr_matrix_HalfLife(i,2)));
    end
    % fprintf(tex_file_halflife, '\\bottomrule\n\\end{tabular}\n\\end{table}\n\n');
    fprintf(tex_file_rmse, '\\end{tabular}\n');
    % fprintf(tex_file_halflife, '\\end{document}\n');
    fclose(tex_file_halflife);

    fprintf('TeX tables written to: %s and %s\n', filename_corrh, filename_halflife);
end

fprintf('RMSE TeX table written to: table_rmse_correlations.tex\n');

%%
if isfield(oo_, 'posterior_mean')          % Dynare's posterior summary; absent for a model estimated outside Dynare (BR, 2026-09-26)
    prior_posterior_table(M_, oo_, bayestopt_, [HLW_directory '/Tables_level/table_prior_posterior_params.tex'], ...
                          struct('tex_names', {cfg.param_tex_names}, 'label', cfg.param_table_label))
else
    fprintf(2, 'HLW_plots: table_prior_posterior_params.tex skipped - oo_.posterior_mean not present (model estimated outside Dynare).\n');
end

%%
%% Write Summary Table 1

% Open file for writing
% fid = fopen([HLW_directory '/Tables_level/table_' summary '.tex'], 'w');
fid = fopen([HLW_directory '/Tables_level/table_summary.tex'], 'w');

% Provenance of the R^2 rows (LaTeX comment, does not render).
fprintf(fid, ['%% R^2 rows: eq. (3.17), R2^(m)_{i,h} = 1 - [Omega_m]_ii / Var[x_{i,t+h}|y^t].\n' ...
              '%% Filter (m=0) and fixed-lag smoother (m>0) share the denominator, the conditional\n' ...
              '%% variance (3.10) of the CONTEMPORANEOUS state, so reliability is weakly increasing\n' ...
              '%% in the hindsight m (Section 3.5). For m=0 this is eq. (3.11) = CorrH^2, unchanged.\n']);

% Write table header
fprintf(fid, '\\begin{tabular}{l|');
for pp = 1:length(PlotVar.Names)
    fprintf(fid, 'c');
end
fprintf(fid, '}\n');

% Write column headers
fprintf(fid, ' & ');
for pp = 1:length(PlotVar.Names)
    fprintf(fid, '%s', PlotVar.Names{pp});
    if pp < length(PlotVar.Names)
        fprintf(fid, ' & ');
    else
        fprintf(fid, ' \\\\\n');
    end
end
fprintf(fid, '\\hline\n');

% Row labels
row_labels = {
    'Median s.d. of state estimate \\ given parameters/innovation s.d.', ...
    'Average s.d. of state estimate/innovation s.d.', ...
    'Innovation s.d. (post mean)', ...
    'Median 4-quarter ahead $R^2$', ...
    'Median 20-quarter ahead $R^2$', ...
    'Median $\nu$ [quarters]', ...
    'Median half-life of estimation errors [quarters]'
};

% Write data rows
for rr = 1:7
    fprintf(fid, '%s & ', row_labels{rr});
    for pp = 1:length(PlotVar.Names)
        val = SummaryTableMoments(rr, pp);

        % Format based on row type
        if isnan(val)
            % not defined for this state (e.g. nu for a white-noise change state)
            fprintf(fid, '--');
        elseif rr <= 2
            % Rows 1-2: the state estimate's s.d. AS A MULTIPLE of the state's
            % own innovation s.d. - a ratio, and therefore dimensionless. These
            % used to be printed with a "pp" suffix, which belongs to row 3
            % alone; the two were sharing one branch.
            fprintf(fid, '%.1f', val);
        elseif rr == 3
            % Row 3: the innovation s.d. itself, genuinely in percentage points
            fprintf(fid, '%.1fpp', val);
        elseif rr <= 5
            % Rows 4-5: R^2, format as 0.XX
            fprintf(fid, '%.2f', val);
        elseif rr == 6
            % Row 6: nu in quarters, one decimal (CHANGED 2026-09-17, was %.2f)
            fprintf(fid, '%.1f', val);
        elseif isinf(val)
            % Row 7: more than half the draws never cross 1/2 within the horizon (2026-09-30)
            fprintf(fid, '$\\infty$');
        else
            % Row 7: format as integer
            fprintf(fid, '%d', round(val));
        end

        if pp < length(PlotVar.Names)
            fprintf(fid, ' & ');
        else
            fprintf(fid, ' \\\\\n');
        end
    end
end
% Close table
fprintf(fid, '\\end{tabular}\n');

% Close file
fclose(fid);

% Undo rescaling
OldSummaryTableMoments=SummaryTableMoments([1:2, 4:end],:);
OldSummaryTableMoments(1,:)=SummaryTableMoments(1,:).*SummaryTableMoments(3,:);
OldSummaryTableMoments(2,:)=SummaryTableMoments(2,:).*SummaryTableMoments(3,:);

fid = fopen([HLW_directory '/Tables_level/table_summary_old_format.tex'], 'w');

% Provenance of the R^2 rows: as in table_summary.tex, eq. (3.17).
fprintf(fid, ['%% R^2 rows: eq. (3.17), R2^(m)_{i,h} = 1 - [Omega_m]_ii / Var[x_{i,t+h}|y^t];\n' ...
              '%% common denominator across information sets, see Section 3.5.\n']);

% Write table header
fprintf(fid, '\\begin{tabular}{l|');
for pp = 1:length(PlotVar.Names)
    fprintf(fid, 'c');
end
fprintf(fid, '}\n');

% Write column headers
fprintf(fid, ' & ');
for pp = 1:length(PlotVar.Names)
    fprintf(fid, '%s', PlotVar.Names{pp});
    if pp < length(PlotVar.Names)
        fprintf(fid, ' & ');
    else
        fprintf(fid, ' \\\\\n');
    end
end
fprintf(fid, '\\hline\n');

% Row labels
row_labels = {
    'Median s.d. of state estimate given parameters', ...
    'Average s.d. of state estimate', ...
    'Median 4-quarter ahead $R^2$', ...
    'Median 20-quarter ahead $R^2$', ...
    'Median $\nu$ [quarters]', ...
    'Median half-life of estimation errors [quarters]'
};

% Write data rows
for rr = 1:6
    fprintf(fid, '%s & ', row_labels{rr});
    for pp = 1:length(PlotVar.Names)
        val = OldSummaryTableMoments(rr, pp);

        % Format based on row type
        if isnan(val)
            fprintf(fid, '--');
        elseif rr <= 2
            % Rows 1-2: format as X.Xpp
            fprintf(fid, '%.1fpp', val);
        elseif rr <= 4
            % Rows 3-4: R^2, format as 0.XX
            fprintf(fid, '%.2f', val);
        elseif rr == 5
            % Row 5: nu in quarters, one decimal (CHANGED 2026-09-17, was %.2f)
            fprintf(fid, '%.1f', val);
        elseif isinf(val)
            fprintf(fid, '$\\infty$');                 % half-life: median draw never crosses 1/2 (2026-09-30)
        else
            % Row 6: format as integer
            fprintf(fid, '%d', round(val));
        end

        if pp < length(PlotVar.Names)
            fprintf(fid, ' & ');
        else
            fprintf(fid, ' \\\\\n');
        end
    end
end

% Close table
fprintf(fid, '\\end{tabular}\n');

% Close file
fclose(fid);


%%
%% Parameter mean and s.d.: prior vs posterior, unconditional and conditional on R^2 > 0.5

n_params = size(param_estimates.matrix, 1);

% Prefer TeX-formatted names if available
if isfield(M_, 'param_names_tex') && ~isempty(M_.param_names_tex)
    param_names_display = M_.param_names_tex;
else
    param_names_display = M_.param_names;
end

% Unconditional moments (these do not depend on variable/horizon)
post_uncond_mean  = mean(param_estimates.matrix, 2);
post_uncond_std   = std( param_estimates.matrix, 0, 2);
prior_uncond_mean = mean(prior_param_estimates.matrix, 2);
prior_uncond_std  = std( prior_param_estimates.matrix, 0, 2);

fmt_cell = @(m,s) sprintf('%.3f (%.3f)', m, s);

for vv = 1:numel(cfg.param_moment_vars)
    var1 = cfg.param_moment_vars{vv};
    pp = find(strcmp(PlotVar.Names, var1{1,1}));
    if isempty(pp), continue; end
    PlotVar.Idx{pp} = strcmp(M_.states_of_interest, PlotVar.Engine{pp});

    for hh = 4:4:20

        % Indices of draws with high R^2 at horizon hh
        post_high_idx  = squeeze(CorrH.ThreeD(PlotVar.Idx{pp}, hh, :)).^2  > 0.5;
        prior_high_idx = squeeze(CorrH.prior_ThreeD(PlotVar.Idx{pp}, hh, :)).^2 > 0.5;

        n_post_tot  = size(param_estimates.matrix, 2);
        n_prior_tot = size(prior_param_estimates.matrix, 2);
        n_post_hi   = sum(post_high_idx);
        n_prior_hi  = sum(prior_high_idx);

        % Conditional moments
        if n_post_hi > 1
            post_cond_mean = mean(param_estimates.matrix(:, post_high_idx), 2);
            post_cond_std  = std( param_estimates.matrix(:, post_high_idx), 0, 2);
        else
            post_cond_mean = nan(n_params,1);
            post_cond_std  = nan(n_params,1);
        end

        if n_prior_hi > 1
            prior_cond_mean = mean(prior_param_estimates.matrix(:, prior_high_idx), 2);
            prior_cond_std  = std( prior_param_estimates.matrix(:, prior_high_idx), 0, 2);
        else
            prior_cond_mean = nan(n_params,1);
            prior_cond_std  = nan(n_params,1);
        end

        % ---------- Console output ----------
        fprintf('\n=== Parameter moments: %s, horizon h=%d ===\n', var1{1,2}, hh);
        fprintf('Posterior draws with R^2>0.5: %d / %d\n', n_post_hi, n_post_tot);
        fprintf('Prior    draws with R^2>0.5: %d / %d\n', n_prior_hi, n_prior_tot);
        fprintf('%-20s | %-18s %-18s | %-18s %-18s\n', ...
            'Parameter', 'Prior (uncond.)', 'Prior (R^2>0.5)', 'Post (uncond.)', 'Post (R^2>0.5)');
        fprintf('%s\n', repmat('-', 1, 100));
        for ii = 1:n_params
            fprintf('%-20s | %-18s %-18s | %-18s %-18s\n', ...
                M_.param_names{ii}, ...
                fmt_cell(prior_uncond_mean(ii), prior_uncond_std(ii)), ...
                fmt_cell(prior_cond_mean(ii),   prior_cond_std(ii)),   ...
                fmt_cell(post_uncond_mean(ii),  post_uncond_std(ii)),  ...
                fmt_cell(post_cond_mean(ii),    post_cond_std(ii)));
        end

        % ---------- TeX output ----------
        var_name_clean = strrep(var1{1,1}, '*', 'star');
        fname = sprintf('%s/Tables_level/table_param_moments_%s_h%d.tex', ...
                        HLW_directory, var_name_clean, hh);
        fid = fopen(fname, 'w');
        fprintf(fid, '%% Auto-generated by the post-estimation script.\n');
        fprintf(fid, '%% Conditioning variable: %s, horizon h=%d, threshold R^2>0.5.\n', var1{1,2}, hh);
        fprintf(fid, '%% N_post_high=%d/%d, N_prior_high=%d/%d.\n', ...
                n_post_hi, n_post_tot, n_prior_hi, n_prior_tot);
        fprintf(fid, '\\begin{tabular}{lcccc}\n');
        fprintf(fid, ' & \\multicolumn{2}{c}{Prior} & \\multicolumn{2}{c}{Posterior} \\\\\n');
        fprintf(fid, '\\cmidrule(lr){2-3} \\cmidrule(lr){4-5}\n');
        fprintf(fid, 'Parameter & Unconditional & $R^{2}_{h=%d}>0.5$ & Unconditional & $R^{2}_{h=%d}>0.5$ \\\\\n', hh, hh);
        fprintf(fid, '\\hline\n');
        for ii = 1:n_params
            fprintf(fid, '%s & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) & %.3f (%.3f) \\\\\n', ...
                param_names_display{ii}, ...
                prior_uncond_mean(ii), prior_uncond_std(ii), ...
                prior_cond_mean(ii),   prior_cond_std(ii),   ...
                post_uncond_mean(ii),  post_uncond_std(ii),  ...
                post_cond_mean(ii),    post_cond_std(ii));
        end
        fprintf(fid, '\\end{tabular}\n');
        fclose(fid);

        fprintf('TeX table written to: %s\n', fname);
    end
end




close all
%% Posterior-predictive fan charts and tables (CHANGED 2026-09-17)
% Until 2026-09-17 this block drew the DeltaR2(h) and RMSE(h) fans from
% DeltaR2_RMSE_plots.mat, the output of the OLD model assessment (50,000
% replications at one parameter vector, no re-estimation). The current
% assessment is model_assessment.m (damped-Newton re-estimated mode); its
% result struct R (model_assessment_lecam.mat) is plotted and tabulated by
% model_assessment_plots: RMSE and DeltaR2 fans per observable, the f-on-r
% proxy-regression panels, the implied-R2 chart at the mode, and the p-value
% tables. Output goes to the same Graphs_level / Tables_level folders with
% prefix MA_.
ma_file = cfg.ma_file;
if isempty(ma_file)
    ma_file = fullfile(HLW_directory, 'model_assessment_lecam.mat');          % output folder first (per-model runs) ...
    if exist(ma_file, 'file') ~= 2, ma_file = fullfile(HLW_code_dir, 'model_assessment_lecam.mat'); end   % ... then the code folder
end
if exist(ma_file, 'file') == 2
    model_assessment_plots(ma_file, struct('graph_dir', fullfile(HLW_directory, 'Graphs_level'), ...
                                           'table_dir', fullfile(HLW_directory, 'Tables_level'), ...
                                           'figures', ~DO_TABLES_ONLY));
else
    fprintf(2, 'HLW_plots: %s not found - run R = model_assessment(300) first; posterior-predictive fans skipped.\n', ma_file);
end

