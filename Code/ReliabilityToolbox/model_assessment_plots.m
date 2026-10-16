function out = model_assessment_plots(R, opts)
%MODEL_ASSESSMENT_PLOTS  Figures and TeX tables for the posterior-predictive assessment.
%
%   model_assessment_plots()            loads model_assessment_lecam.mat from this folder
%   model_assessment_plots(R)           R = output struct of model_assessment
%   model_assessment_plots(file)        path of a saved R
%   model_assessment_plots(R, opts)     see options
%
%   ADDED 2026-09-17. model_assessment.m prints a report and saves R but draws
%   nothing; the fan charts at the end of HLW_plots.m consumed the OLD
%   assessment's DeltaR2_RMSE_plots.mat (50,000 replications at one parameter
%   vector). This function replaces them and reads R. Everything is
%   "data vs the model's sampling distribution": the statistic on the actual
%   data at the posterior mode theta_hat, against its distribution over the D
%   artificial datasets simulated at theta_hat, evaluated
%     (ii)  at theta_hat            -> R.sim_mode   (parametric bootstrap, fixed theta)
%     (iii) at theta_hat^(d)        -> R.sim_lecam  (re-estimated, damped Newton)
%   With an R of model_assessment_is (posterior-mean estimator, 2026-09-29) the DGP is the
%   posterior mean theta_bar: R.data is the posterior-mean estimator on the data, R.sim_mode the
%   statistic at theta_bar on every dataset, R.sim_lecam the posterior-mean estimator on every
%   dataset, and R.data_dgp the population objects at theta_bar (blue reference line, "model"
%   column of the tables). The drivers pass the matching labels (label_ii, label_ii_tex,
%   label_iii, label_iii_tex, label_data, label_pop).
%   Two-sided percentile p-values p = 2 min(q, 1-q), q = share of the null <= data.
%
%   FIGURES (PNG, opts.graph_dir, prefix opts.prefix):
%     <prefix>RMSE_<obs>.png, <prefix>DeltaR2_<obs>.png
%           one per observable: data (black), band at theta_hat (blue),
%           band at theta_hat^(d) (red), medians dashed; DeltaR2 = model R2 -
%           empirical R2 of the h-step forecast of the observable, zero line.
%     <prefix>legend.png
%           the ONE legend shared by all MA_ figures (no legends inside panels)
%     <prefix>P2_<state>.png, <prefix>DeltaP2_<state>.png, <prefix>beta_<state>.png, <prefix>alpha_<state>.png
%           smoother-revision (proxy) regression f = alpha + beta r, one row of
%           panels (one per smoother lag m) with a common y-scale across m:
%           Phat^{2,m} with the population P^{2,m} at theta_hat; Delta_{i,h,m} =
%           Phat^{2,m} - P^{2,m} against the same difference under both nulls;
%           betahat (reference 1); alphahat (reference 0). (Until 2026-09-18
%           these were one 3 x n_m figure <prefix>fr_<state>.png.)
%     <prefix>impliedR2_<state>.png
%           R^2_{i,h} implied by the empirical P^2 through Prop. 5.1,
%           R^2 = P (1-rho)/(1-rho P), rho = [Omega_m]_ii/Omega_ii (fr_implied_R2):
%           data (black) with both bands, against the model's own R^2_{i,h} at
%           theta_hat (blue dashed). Needs fr.omega_ratio / fr.R2_h_base in R
%           (model_assessment run on or after 2026-09-17); skipped otherwise.
%   TABLES (TeX tabulars, opts.table_dir):
%     <prefix>table_RMSE.tex, <prefix>table_DeltaR2.tex
%           observable x horizon: data | median, p at theta_hat | median, p at theta_hat^(d)
%     <prefix>table_fr.tex
%           state, m, horizon: Phat^2, betahat, alphahat and implied R^2 with
%           their p-values under both nulls, plus the model R^2_{i,h}
%     <prefix>table_DeltaP2.tex
%           state, m, horizon: Delta_{i,h,m} on the data, median and band of the
%           two nulls, p-values (ADDED 2026-09-18)
%
%   Options: .graph_dir ('Graphs_level'), .table_dir ('Tables_level') - relative
%   paths are taken from this file's folder; .prefix ('MA_'); .band (90);
%   .h_report (R.opts.h_report); .states (R.fr_states, degenerate ones skipped);
%   .close (true: close each figure after printing).
%
%   Runs under MATLAB and Octave (fill/print only, no exportgraphics/tiledlayout).
%
%   See also: model_assessment, fr_implied_R2, HLW_plots.

% This file lives in ReliabilityToolbox (MOVED 2026-09-18). Relative output
% folders are taken from the folder of the R file (the model folder); with no
% argument the HLW run's file is used.
here = fileparts(mfilename('fullpath'));          % ReliabilityToolbox
if nargin < 1 || isempty(R), R = fullfile(here, '..', 'HLW_new', 'hlw_level', 'model_assessment_lecam.mat'); end
if ischar(R), rdir = fileparts(R); if isempty(rdir), rdir = pwd; end; R = load(R); else, rdir = pwd; end
if nargin < 2, opts = struct(); end
def = struct('graph_dir', 'Graphs_level', 'table_dir', 'Tables_level', 'prefix', 'MA_', 'band', 90, ...
             'h_report', R.opts.h_report, 'states', {R.fr_states}, 'close', true, 'figures', true, ...   % .figures false: tables only (ADDED 2026-09-18)
             'label_iii', 're-estimated null', 'label_iii_tex', '\hat\theta^{(d)}', ...                 % naming of null (iii) (ADDED 2026-09-28: model_assessment_is
             'label_data', 'data (statistic at the posterior mode)', ...                               %  passes 'posterior-mean null', '\bar\theta^{(d)}',
             'label_pop', 'population value at the mode', ...                                          %  'population value at the posterior mean' (2026-09-29)
             'label_ii', 'null at the posterior mode', 'label_ii_tex', '\hat\theta');                  %  null (ii): statistics at the DGP point ('null at the posterior mean', '\bar\theta')
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = def.(fn{i}); end
end
gdir = opts.graph_dir;  if ~is_absolute(gdir), gdir = fullfile(rdir, gdir); end
tdir = opts.table_dir;  if ~is_absolute(tdir), tdir = fullfile(rdir, tdir); end
if ~exist(gdir, 'dir'), mkdir(gdir); end
if ~exist(tdir, 'dir'), mkdir(tdir); end

h  = R.h(:)';
H  = numel(h);
hs = opts.h_report(opts.h_report <= H);
qlo = (100 - opts.band) / 200;  qhi = 1 - qlo;
c_mode  = [0.18 0.44 0.72];      % blue : null (ii), statistic at the DGP point (theta_hat, or theta_bar for model_assessment_is)
c_lecam = [0.80 0.30 0.10];      % red  : null (iii), re-estimated on every dataset
pv = @(null, dat) 2 * min(mean(null <= dat, 'omitnan'), 1 - mean(null <= dat, 'omitnan'));
out = struct('files', {{}});
% the DGP point, for the table comments: model_assessment_is (2026-09-29) records R.dgp ('mean' | 'mode'); model_assessment's is the mode
if isfield(R, 'dgp') && strcmpi(R.dgp, 'mean'), dgp_str = 'the posterior mean'; else, dgp_str = 'the posterior mode'; end

%% ---------------------------------------------------- common legend (one file per model)
% CHANGED 2026-09-18: no legends inside the panels; <prefix>legend.png carries
% the one legend shared by every MA_ figure (bands, null medians, data, and the
% population object at theta_hat - P^{2,m} in the P2 panels, R^2_{i,h} in the
% implied-R2 panels).
if opts.figures
    fig = figure('Name', 'legend', 'Position', [50 50 640 260]);  ax = axes('Position', [0 0 1 1]);  hold on;
    fill([0 1 1 0], [0 0 1 1], c_mode,  'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', sprintf('%d%% band of the %s', opts.band, opts.label_ii));
    plot([0 1], [0 0], '--', 'Color', c_mode,  'LineWidth', 1.0, 'DisplayName', ['median of the ' opts.label_ii]);
    fill([0 1 1 0], [0 0 1 1], c_lecam, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', sprintf('%d%% band of the %s', opts.band, opts.label_iii));
    plot([0 1], [0 0], '--', 'Color', c_lecam, 'LineWidth', 1.0, 'DisplayName', ['median of the ' opts.label_iii]);
    plot([0 1], [0 0], 'ko-', 'LineWidth', 1.4, 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'DisplayName', opts.label_data);
    plot([0 1], [0 0], 'b--', 'LineWidth', 1.6, 'DisplayName', [opts.label_pop ' (P^{2,m}_{i,h}; model R^2_{i,h})']);
    plot([0 1], [0 0], 'k:', 'LineWidth', 0.8, 'DisplayName', 'reference (0 for \Delta and \alpha, 1 for \beta)');
    lg = legend('Location', 'north', 'FontSize', 11);  legend('boxoff');
    % The drawn objects are moved out of view by the axes limits instead of Visible 'off': MATLAB
    % draws legend entries of invisible objects at reduced opacity (the "greyed-out" legends of
    % 2026-09-28), Octave does not. (CHANGED 2026-09-29)
    set(ax, 'Visible', 'off', 'XLim', [5 6], 'YLim', [5 6]);  set(lg, 'Visible', 'on');
    out = save_fig(out, fig, fullfile(gdir, sprintf('%slegend.png', opts.prefix)), opts.close);

    % ADDED 2026-09-19: wide two-row variants for placing under a row of panels,
    % without the reference-line entry and with the population entry shortened:
    %   <prefix>legend_2rows.png        six entries (with the blue population line)
    %   <prefix>legend_2rows_nopop.png  five entries (RMSE / DeltaR2 panels have no population line)
    for with_pop = [true false]
        fig = figure('Name', 'legend 2 rows', 'Position', [50 50 1450 110]);  ax = axes('Position', [0 0 1 1]);  hold on;
        fill([0 1 1 0], [0 0 1 1], c_mode,  'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', sprintf('%d%% band of the %s', opts.band, opts.label_ii));
        plot([0 1], [0 0], '--', 'Color', c_mode,  'LineWidth', 1.0, 'DisplayName', ['median of the ' opts.label_ii]);
        fill([0 1 1 0], [0 0 1 1], c_lecam, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', sprintf('%d%% band of the %s', opts.band, opts.label_iii));
        plot([0 1], [0 0], '--', 'Color', c_lecam, 'LineWidth', 1.0, 'DisplayName', ['median of the ' opts.label_iii]);
        plot([0 1], [0 0], 'ko-', 'LineWidth', 1.4, 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'DisplayName', opts.label_data);
        if with_pop
            plot([0 1], [0 0], 'b--', 'LineWidth', 1.6, 'DisplayName', opts.label_pop);
            fname = sprintf('%slegend_2rows.png', opts.prefix);
        else
            fname = sprintf('%slegend_2rows_nopop.png', opts.prefix);
        end
        lg = legend('Location', 'north', 'FontSize', 11, 'NumColumns', 3);  legend('boxoff');
        set(ax, 'Visible', 'off', 'XLim', [5 6], 'YLim', [5 6]);  set(lg, 'Visible', 'on');   % see the note above (MATLAB dims entries of invisible objects)
        out = save_fig(out, fig, fullfile(gdir, fname), opts.close);
    end
end

%% ---------------------------------------------------- RMSE and DeltaR2 per observable
for stat = {'RMSE', 'DeltaR2'}
    st = stat{1};
    dat = R.data.(st);  A = R.sim_mode.(st);  B = R.sim_lecam.(st);
    for i = 1:numel(R.obs_names)
        if ~opts.figures, break; end
        fig = figure('Name', sprintf('%s %s', st, R.obs_names{i}), 'Position', [50 50 560 380]);  hold on;
        band(h, squeeze(A(i, :, :)), c_mode,  qlo, qhi, sprintf('%d%% band, %s', opts.band, opts.label_ii));
        band(h, squeeze(B(i, :, :)), c_lecam, qlo, qhi, sprintf('%d%% band, %s', opts.band, opts.label_iii));
        plot(h, dat(i, :), 'ko-', 'LineWidth', 1.6, 'MarkerFaceColor', 'k', 'MarkerSize', 3, 'DisplayName', opts.label_data);
        if strcmp(st, 'DeltaR2'), plot(h, zeros(size(h)), 'k-', 'LineWidth', 0.6, 'HandleVisibility', 'off'); end
        xlim([0.5, H + 0.5]);  set(gca, 'XTick', [1 4:4:H]);
        xlabel('horizon h [quarters]');
        if strcmp(st, 'RMSE'), ylabel('RMSE of the h-step forecast'); else, ylabel('\DeltaR^2(h) = model R^2 - empirical R^2'); end
        title(sprintf('%s: %s', st, R.obs_names{i}), 'Interpreter', 'none');
        set(gca, 'FontSize', 11);  grid on;  box on;  hold off;
        out = save_fig(out, fig, fullfile(gdir, sprintf('%s%s_%s.png', opts.prefix, st, R.obs_names{i})), opts.close);
    end
    % TeX table: observable x horizon
    f = fullfile(tdir, sprintf('%stable_%s.tex', opts.prefix, st));
    fid = fopen(f, 'w');
    fprintf(fid, '%% %s of the h-step forecasts of the observables: %s vs the model''s sampling\n', st, opts.label_data);
    fprintf(fid, '%% distribution over D = %d datasets simulated at %s; two-sided percentile p-values. model_assessment_plots, %s.\n', R.D, dgp_str, datestr(now, 'yyyy-mm-dd'));
    fprintf(fid, '\\begin{tabular}{ll%s}\n', repmat('c', 1, 5));
    fprintf(fid, ' & & & \\multicolumn{2}{c}{null at $%s$} & \\multicolumn{2}{c}{null at $%s$} \\\\\n', opts.label_ii_tex, opts.label_iii_tex);
    fprintf(fid, 'observable & $h$ & data & median & $p$ & median & $p$ \\\\\n\\hline\n');
    for i = 1:numel(R.obs_names)
        for hh = hs
            a = squeeze(A(i, hh, :));  b = squeeze(B(i, hh, :));
            fprintf(fid, '%s & %d & %.3f & %.3f & %.2f & %.3f & %.2f \\\\\n', tex_escape(R.obs_names{i}), hh, dat(i, hh), ...
                    median(a, 'omitnan'), pv(a, dat(i, hh)), median(b, 'omitnan'), pv(b, dat(i, hh)));
        end
    end
    fprintf(fid, '\\end{tabular}\n');
    fclose(fid);
    out.files{end+1} = f;
end

%% ------------------------------------------------ proxy regression f = alpha + beta r
% Four figures per base state (CHANGED 2026-09-18: was one 3 x n_m panel figure),
% each with the SAME y-scale across the m columns, tightened to the plotted
% objects (+5%) and, where the object is bounded, clipped to [0,1]:
%   <prefix>P2_<state>        Phat^{2,m}_{i,h} (data) with both null bands and the
%                             population P^{2,m}_{i,h} at theta_hat (blue dashed)
%   <prefix>DeltaP2_<state>   Delta_{i,h,m} = Phat^{2,m}_{i,h} - P^{2,m}_{i,h}: data
%                             against the same difference under the nulls (P at
%                             theta_hat for null (ii), at theta_hat^(d) for null (iii))
%   <prefix>fr_coef_<state>   betahat (row 1, reference 1) and alphahat (row 2, reference 0)
%   <prefix>impliedR2_<state> R^2_{i,h} implied by Phat (Prop. 5.1) with the model's R^2_{i,h}
% Legend: shaded = 90% band of the null, dashed in the band's colour = median of
% that null, black = data, blue dashed = population object at theta_hat.
have_implied = isfield(R.data.fr, 'omega_ratio') && isfield(R.sim_lecam, 'fr_omega_ratio') && isfield(R.data.fr, 'R2_h_base');
n_m = numel(R.fr_m);
ftab = fullfile(tdir, sprintf('%stable_fr.tex', opts.prefix));
fid = fopen(ftab, 'w');
fprintf(fid, '%% Smoother-revision (proxy) regression f = alpha + beta r, f = xhat(t+h|t+h+m) - xhat(t+h|t), r = xhat(t+h|t+h) - xhat(t+h|t):\n');
fprintf(fid, '%% OLS on the data [%s] vs the model''s sampling distribution (D = %d datasets simulated at %s); p = two-sided percentile p-value,\n', opts.label_data, R.D, dgp_str);
fprintf(fid, '%% first under the %s, second under the %s. Implied R2 = P(1-rho)/(1-rho P), Prop. 5.1; "model" = population R2 at %s.\n', opts.label_ii, opts.label_iii, dgp_str);
fprintf(fid, '\\begin{tabular}{lrr|ccc|ccc|ccc|cccc}\n');
fprintf(fid, ' & & & \\multicolumn{3}{c|}{$\\hat P^{2,m}_{i,h}$} & \\multicolumn{3}{c|}{$\\hat\\beta$} & \\multicolumn{3}{c|}{$\\hat\\alpha$} & \\multicolumn{4}{c}{implied $R^2_{i,h}$} \\\\\n');
fprintf(fid, 'state & $m$ & $h$ & data & $p$ & $p^{(d)}$ & data & $p$ & $p^{(d)}$ & data & $p$ & $p^{(d)}$ & data & $p$ & $p^{(d)}$ & model \\\\\n\\hline\n');
fdel = fullfile(tdir, sprintf('%stable_DeltaP2.tex', opts.prefix));
fid_d = fopen(fdel, 'w');
fprintf(fid_d, '%% Delta_{i,h,m} = Phat^{2,m}_{i,h} (OLS on the data [%s]) - P^{2,m}_{i,h} (population at %s, the DGP), vs the model''s\n', opts.label_data, dgp_str);
fprintf(fid_d, '%% sampling distribution of the same difference: %s (P at the DGP), %s (P at the re-estimated point). D = %d.\n', opts.label_ii, opts.label_iii, R.D);
fprintf(fid_d, '\\begin{tabular}{lrr|c|cc|cc}\n');
fprintf(fid_d, ' & & & & \\multicolumn{2}{c|}{null at $%s$} & \\multicolumn{2}{c}{null at $%s$} \\\\\n', opts.label_ii_tex, opts.label_iii_tex);
fprintf(fid_d, 'state & $m$ & $h$ & $\\Delta$ data & median [%d\\%%] & $p$ & median [%d\\%%] & $p$ \\\\\n\\hline\n', opts.band, opts.band);
for j = 1:numel(R.fr_states)
    if ~any(strcmp(opts.states, R.fr_states{j})), continue; end
    if all(isnan(R.data.fr.R2_emp(j, :, :))), continue; end     % degenerate regressor (zero filter gain)
    sj = R.fr_states{j};

    % ---- objects per m: data, null (ii) at theta_hat, null (iii) re-estimated -----
    P_dat = cell(n_m, 1);  P_A = P_dat;  P_B = P_dat;  P_th = P_dat;
    D_dat = P_dat;  D_A = P_dat;  D_B = P_dat;
    I_dat = P_dat;  I_A = P_dat;  I_B = P_dat;  rho0 = nan(n_m, 1);
    % Population objects of the DGP (theta_hat) for null (ii). With model_assessment they coincide with
    % the data estimator's (R.data is the statistic at theta_hat); with model_assessment_is R.data holds
    % the POSTERIOR-MEAN estimator, whose population P and omega ratio differ from those at theta_hat,
    % and R.data_mode carries the theta_hat objects (ADDED 2026-09-28).
    % (RENAMED 2026-09-29: model_assessment_is stores them as R.data_dgp; the datasets are simulated at the
    % posterior MEAN there, so the DGP objects are the population values at theta-bar. The blue
    % reference line and the tables' "model" column show the DGP's population values.)
    if isfield(R, 'data_dgp'), fr_dgp = R.data_dgp.fr; elseif isfield(R, 'data_mode'), fr_dgp = R.data_mode.fr; else, fr_dgp = R.data.fr; end
    for im = 1:n_m
        P_dat{im} = squeeze(R.data.fr.R2_emp(j, :, im));                    % 1 x H
        P_th{im}  = squeeze(R.data.fr.R2_theory(j, :, im));                 % population P of the data estimator
        P_dgp     = squeeze(fr_dgp.R2_theory(j, :, im));                    % population P at theta_hat (the DGP)
        P_A{im}   = reshape(R.sim_mode.fr_R2_emp(j, :, im, :), H, []);      % H x D
        P_B{im}   = reshape(R.sim_lecam.fr_R2_emp(j, :, im, :), H, []);
        D_dat{im} = P_dat{im} - P_th{im};
        D_A{im}   = P_A{im} - P_dgp(:);                                      % P at theta_hat under null (ii)
        D_B{im}   = P_B{im} - reshape(R.sim_lecam.fr_R2_theory(j, :, im, :), H, []);   % P at theta_hat^(d) under null (iii)
        if have_implied
            rho0(im)  = R.data.fr.omega_ratio(j, im);
            I_dat{im} = fr_implied_R2(P_dat{im}, rho0(im));
            I_A{im}   = fr_implied_R2(P_A{im}, fr_dgp.omega_ratio(j, im));
            rho_d     = reshape(R.sim_lecam.fr_omega_ratio(j, im, :), 1, []);
            I_B{im}   = fr_implied_R2(P_B{im}, repmat(rho_d, H, 1));
        end
    end

    if opts.figures
        % Each figure: one row of n_m panels sharing a y-scale, no in-panel legend
        % (see <prefix>legend.png), panel 460 x 380 px, 11pt axes.
        pw = 460;  ph = 380;
        % ---- Phat^2 ----
        yl = tight_limits([P_dat; P_th; band_edges(P_A, qlo, qhi); band_edges(P_B, qlo, qhi)], 0.05, [0 1]);
        fig = figure('Name', sprintf('proxy R2, %s', sj), 'Position', [50 50 pw * n_m, ph]);
        for im = 1:n_m
            subplot(1, n_m, im);  hold on;
            band(h, P_A{im}, c_mode,  qlo, qhi, '');
            band(h, P_B{im}, c_lecam, qlo, qhi, '');
            plot(h, squeeze(fr_dgp.R2_theory(j, :, im)), 'b--', 'LineWidth', 1.6);   % population P at the DGP point
            plot(h, P_dat{im}, 'ko-', 'LineWidth', 1.4, 'MarkerFaceColor', 'k', 'MarkerSize', 3);
            panel_format(H, yl, 'horizon h', 'P^{2,m}_{i,h}', im == 1);
            title(sprintf('%s: P^{2,m}_{i,h}, m = %d', sj, R.fr_m(im)), 'Interpreter', 'tex', 'FontSize', 12);
            hold off;
        end
        out = save_fig(out, fig, fullfile(gdir, sprintf('%sP2_%s.png', opts.prefix, sj)), opts.close);

        % ---- Delta = Phat^2 - P^2 ----
        yl = tight_limits([D_dat; band_edges(D_A, qlo, qhi); band_edges(D_B, qlo, qhi); {0}], 0.05, []);
        fig = figure('Name', sprintf('Delta P2, %s', sj), 'Position', [50 50 pw * n_m, ph]);
        for im = 1:n_m
            subplot(1, n_m, im);  hold on;
            band(h, D_A{im}, c_mode,  qlo, qhi, '');
            band(h, D_B{im}, c_lecam, qlo, qhi, '');
            plot(h, zeros(size(h)), 'k:', 'LineWidth', 0.8);
            plot(h, D_dat{im}, 'ko-', 'LineWidth', 1.4, 'MarkerFaceColor', 'k', 'MarkerSize', 3);
            panel_format(H, yl, 'horizon h', '\Delta_{i,h,m}', im == 1);
            title(sprintf('%s: \\Delta_{i,h,m}, m = %d', sj, R.fr_m(im)), 'Interpreter', 'tex', 'FontSize', 12);
            hold off;
        end
        out = save_fig(out, fig, fullfile(gdir, sprintf('%sDeltaP2_%s.png', opts.prefix, sj)), opts.close);

        % ---- beta and alpha, one figure each ----
        coef = {'beta_emp', '\beta', 1, 'beta'; 'alpha_emp', '\alpha', 0, 'alpha'};
        for r = 1:2
            fld = coef{r, 1};
            C_dat = cell(n_m, 1);  C_A = C_dat;  C_B = C_dat;
            for im = 1:n_m
                C_dat{im} = squeeze(R.data.fr.(fld)(j, :, im));
                C_A{im}   = reshape(R.sim_mode.(['fr_' fld])(j, :, im, :), H, []);
                C_B{im}   = reshape(R.sim_lecam.(['fr_' fld])(j, :, im, :), H, []);
            end
            yl = tight_limits([C_dat; band_edges(C_A, qlo, qhi); band_edges(C_B, qlo, qhi); {coef{r, 3}}], 0.05, []);
            fig = figure('Name', sprintf('f-on-r %s, %s', coef{r, 4}, sj), 'Position', [50 50 pw * n_m, ph]);
            for im = 1:n_m
                subplot(1, n_m, im);  hold on;
                band(h, C_A{im}, c_mode,  qlo, qhi, '');
                band(h, C_B{im}, c_lecam, qlo, qhi, '');
                plot(h, coef{r, 3} * ones(size(h)), 'k:', 'LineWidth', 0.8);
                plot(h, C_dat{im}, 'ko-', 'LineWidth', 1.4, 'MarkerFaceColor', 'k', 'MarkerSize', 3);
                panel_format(H, yl, 'horizon h', ['OLS ' coef{r, 2}], im == 1);
                title(sprintf('%s: %s, m = %d', sj, coef{r, 2}, R.fr_m(im)), 'Interpreter', 'tex', 'FontSize', 12);
                hold off;
            end
            out = save_fig(out, fig, fullfile(gdir, sprintf('%s%s_%s.png', opts.prefix, coef{r, 4}, sj)), opts.close);
        end

        % ---- implied R2 vs the model's R2_{i,h} ----
        if have_implied
            yl = tight_limits([I_dat; {fr_dgp.R2_h_base(j, :)}; band_edges(I_A, qlo, qhi); band_edges(I_B, qlo, qhi)], 0.05, [0 1]);
            fig = figure('Name', sprintf('implied R2, %s', sj), 'Position', [50 50 pw * n_m, ph]);
            for im = 1:n_m
                subplot(1, n_m, im);  hold on;
                band(h, I_A{im}, c_mode,  qlo, qhi, '');
                band(h, I_B{im}, c_lecam, qlo, qhi, '');
                plot(h, fr_dgp.R2_h_base(j, :), 'b--', 'LineWidth', 1.6);          % model R^2 at the DGP point
                plot(h, I_dat{im}, 'ko-', 'LineWidth', 1.4, 'MarkerFaceColor', 'k', 'MarkerSize', 3);
                panel_format(H, yl, 'horizon h', 'R^2_{i,h}', im == 1);
                title(sprintf('%s: implied R^2, m = %d (\\rho = %.2f)', sj, R.fr_m(im), rho0(im)), 'Interpreter', 'tex', 'FontSize', 12);
                hold off;
            end
            out = save_fig(out, fig, fullfile(gdir, sprintf('%simpliedR2_%s.png', opts.prefix, sj)), opts.close);
        end
    end   % opts.figures

    % ---- table rows ----
    coef_all = {'R2_emp', 'beta_emp', 'alpha_emp'};
    for im = 1:n_m
        for hh = hs
            cells = cell(1, 4);
            for r = 1:3
                fld = coef_all{r};
                d0 = R.data.fr.(fld)(j, hh, im);
                a = squeeze(R.sim_mode.(['fr_' fld])(j, hh, im, :));  b = squeeze(R.sim_lecam.(['fr_' fld])(j, hh, im, :));
                cells{r} = sprintf('%.3f & %.2f & %.2f', d0, pv(a, d0), pv(b, d0));
            end
            if have_implied
                d0 = I_dat{im}(hh);  a = I_A{im}(hh, :)';  b = I_B{im}(hh, :)';
                cells{4} = sprintf('%.3f & %.2f & %.2f & %.3f', d0, pv(a, d0), pv(b, d0), fr_dgp.R2_h_base(j, hh));   % model R^2 at the DGP point
            else
                cells{4} = '-- & -- & -- & --';
            end
            fprintf(fid, '%s & %d & %d & %s & %s & %s & %s \\\\\n', tex_escape(sj), R.fr_m(im), hh, cells{1}, cells{2}, cells{3}, cells{4});
            a = D_A{im}(hh, :)';  b = D_B{im}(hh, :)';  d0 = D_dat{im}(hh);
            fprintf(fid_d, '%s & %d & %d & %.3f & %.3f [%.3f, %.3f] & %.2f & %.3f [%.3f, %.3f] & %.2f \\\\\n', tex_escape(sj), R.fr_m(im), hh, d0, ...
                    median(a, 'omitnan'), quantile(a, qlo), quantile(a, qhi), pv(a, d0), median(b, 'omitnan'), quantile(b, qlo), quantile(b, qhi), pv(b, d0));
        end
    end
end
fprintf(fid, '\\end{tabular}\n');    fclose(fid);    out.files{end+1} = ftab;
fprintf(fid_d, '\\end{tabular}\n');  fclose(fid_d);  out.files{end+1} = fdel;
fprintf('model_assessment_plots: %d files written to %s and %s\n', numel(out.files), gdir, tdir);
if ~have_implied
    fprintf('  (implied-R2 figure and columns skipped: R predates fr.omega_ratio / fr.R2_h_base - re-run model_assessment)\n');
end
end


%% ========================================================================
function band(h, sims, col, qlo, qhi, label, label_median)
% shaded quantile band over the 2nd dimension of sims (H x D) with its median
% (the median is a legend entry when label_median is given)
sims = reshape(sims, numel(h), []);
lo  = quantile(sims, qlo, 2)';  hi = quantile(sims, qhi, 2)';  md = median(sims, 2, 'omitnan')';
ok  = ~isnan(lo) & ~isnan(hi);
if any(ok)
    fill([h(ok), fliplr(h(ok))], [lo(ok), fliplr(hi(ok))], col, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', label);
end
if nargin >= 7 && ~isempty(label_median)
    plot(h, md, '--', 'Color', col, 'LineWidth', 1.0, 'DisplayName', label_median);
else
    plot(h, md, '--', 'Color', col, 'LineWidth', 1.0, 'HandleVisibility', 'off');
end
end

function panel_format(H, yl, xlab, ylab, first)
% common axes formatting of the one-row MA_ figures
xlim([0.5, H + 0.5]);  ylim(yl);  set(gca, 'XTick', [1 4:4:H], 'FontSize', 11);  grid on;  box on;
xlabel(xlab, 'FontSize', 11);
if first, ylabel(ylab, 'FontSize', 12); end
end

function E = band_edges(S, qlo, qhi)
% lower and upper quantile curves of each H x D array in the cell S, as a cell of vectors
E = cell(2 * numel(S), 1);
for k = 1:numel(S)
    E{2*k-1} = quantile(S{k}, qlo, 2)';
    E{2*k}   = quantile(S{k}, qhi, 2)';
end
end

function yl = tight_limits(series, pad, clip)
% common y-limits for a set of curves (cell of numeric arrays): range + pad*range, clipped to clip = [lo hi] if given
v = [];
for k = 1:numel(series), s = series{k}(:); v = [v; s(isfinite(s))]; end
if isempty(v), yl = [0 1]; return; end
lo = min(v);  hi = max(v);  r = max(hi - lo, 1e-3);
yl = [lo - pad * r, hi + pad * r];
if ~isempty(clip), yl = [max(yl(1), clip(1)), min(yl(2), clip(2))]; end
if yl(2) <= yl(1), yl = [yl(1) - 1e-3, yl(1) + 1e-3]; end
end

function out = save_fig(out, fig, file, do_close)
print(fig, '-dpng', '-r110', file);
out.files{end+1} = file;
if do_close, close(fig); end
end

function s = tex_escape(s)
s = strrep(s, '_', '\_');
end

function tf = is_absolute(p)
tf = ~isempty(p) && (p(1) == '/' || p(1) == '\' || (numel(p) > 1 && p(2) == ':'));
end
