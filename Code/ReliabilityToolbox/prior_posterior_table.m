function prior_posterior_table(M_, oo_, bayestopt_, filename, opts)
% PRIOR_POSTERIOR_TABLE Creates a LaTeX table comparing prior and posterior
%
% Inputs:
%   M_         - Dynare model structure
%   oo_        - Dynare results structure
%   bayestopt_ - Dynare Bayesian estimation options structure
%   filename   - (optional) output filename (default: 'prior_posterior.tex')
%   opts       - (optional, ADDED 2026-09-18) struct with any of
%                .tex_names  n x 2 cell {Dynare name, LaTeX name}: parameters listed
%                            here are printed with the LaTeX name (e.g. 'se_ygap',
%                            '$\sigma_c$'); others keep the escaped Dynare name
%                .label      LaTeX label of the table (default 'tab:prior_posterior')
%                .caption    caption (default 'Prior and Posterior Distributions')
%
% Example:
%   prior_posterior_table(M_, oo_, bayestopt_)
%   prior_posterior_table(M_, oo_, bayestopt_, 'results/estimates.tex', ...
%       struct('tex_names', {{'a1', '$a_1$'; 'se_ygap', '$\sigma_c$'}}, 'label', 'tab:uc3_prior_posterior'))

if nargin < 4 || isempty(filename)
    filename = 'prior_posterior.tex';
end
if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'tex_names') || isempty(opts.tex_names), opts.tex_names = cell(0, 2); end
if ~isfield(opts, 'label')   || isempty(opts.label),   opts.label   = 'tab:prior_posterior'; end
if ~isfield(opts, 'caption') || isempty(opts.caption), opts.caption = 'Prior and Posterior Distributions'; end

% Get parameter information
param_names = M_.param_names;
% Rows = the ESTIMATED parameters in bayestopt_ order (pshape/p1/p2 are indexed that way). FIXED
% 2026-09-29: the loop ran over fieldnames(oo_.posterior_mean.parameters), which in the 2-step ECM run
% (calibrated ECM_coeff, 24 estimated parameters) held 25 fields and indexed bayestopt_ past its end.
% Names without posterior moments (not estimated in this run) are skipped.
post_names = fieldnames(oo_.posterior_mean.parameters);
if isfield(bayestopt_, 'name') && ~isempty(bayestopt_.name)
    est_param_names = bayestopt_.name(:);
    keep = ismember(est_param_names, post_names);
    est_param_names = est_param_names(keep);  bay_index = find(keep);
else
    est_param_names = post_names;  bay_index = 1:numel(post_names);
end
n_estimated = length(est_param_names);

% Open file
fid = fopen(filename, 'w');

% Write LaTeX preamble
fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{%s}\n', opts.caption);
fprintf(fid, '\\label{%s}\n', opts.label);
fprintf(fid, '\\begin{tabular}{lcccccc}\n');
fprintf(fid, '\\hline\\hline\n');
fprintf(fid, 'Parameter & Prior & Prior & Prior & Posterior & Posterior & Posterior \\\\\n');
fprintf(fid, '          & Distribution & Mean & Std & Mean & Std & 90\\%% HPD \\\\\n');
fprintf(fid, '\\hline\n');

% Loop through estimated parameters
for i = 1:n_estimated
    % Get parameter name
    pname = est_param_names{i};
    k_tex = find(strcmp(opts.tex_names(:, 1), pname), 1);
    if isempty(k_tex)
        pname_latex = strrep(pname, '_', '\_'); % Escape underscores for LaTeX
    else
        pname_latex = opts.tex_names{k_tex, 2};
    end

    % Prior information (bayestopt_ index of this parameter)
    ib = bay_index(i);
    prior_dist = get_prior_name(bayestopt_.pshape(ib));
    prior_mean = bayestopt_.p1(ib);  % Prior mean
    prior_std = bayestopt_.p2(ib);   % Prior std

    % Posterior information
    post_mean = oo_.posterior_mean.parameters.(pname);
    post_std = oo_.posterior_std.parameters.(pname);
    hpd_lower = oo_.posterior_hpdinf.parameters.(pname);
    hpd_upper = oo_.posterior_hpdsup.parameters.(pname);

    % Write row
    fprintf(fid, '%s & %s & %.4f & %.4f & %.4f & %.4f & [%.4f, %.4f] \\\\\n', ...
        pname_latex, prior_dist, prior_mean, prior_std, post_mean, post_std, ...
        hpd_lower, hpd_upper);
end

% Close table
fprintf(fid, '\\hline\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');

fclose(fid);

fprintf('LaTeX table saved to: %s\n', filename);

end
