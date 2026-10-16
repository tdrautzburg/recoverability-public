function [estim_params_, bayestopt_, M_] = estim_params_from_names(param_names, M_)
%ESTIM_PARAMS_FROM_NAMES  Minimal estim_params_/bayestopt_ for a model whose draws are its parameter vector.
%
%   [estim_params_, bayestopt_]     = ESTIM_PARAMS_FROM_NAMES(param_names)
%   [estim_params_, bayestopt_, M_] = ESTIM_PARAMS_FROM_NAMES(param_names, M_)
%
%   The third output is M_ with the fields set_all_parameters reads
%   unconditionally added when missing: M_.Sigma_e, M_.Correlation_matrix
%   (identity, q shocks), M_.Skew_e (zeros(0,4): no skewed shocks), M_.H (zeros)
%   and M_.Correlation_matrix_ME (identity, p observables) - Dynare 7 reads H,
%   Correlation_matrix_ME and Skew_e before testing whether any measurement
%   error or skewness is estimated. dr_from_abcd sets the same
%   fields on every solve, so only the first call of a per-draw routine on a
%   never-solved M_ needs this.
%
%   ADDED 2026-09-26. The per-draw routines write a draw xparam1 into the model
%   with Dynare's set_all_parameters(xparam1, estim_params_, M_), which needs the
%   estim_params_ bookkeeping of an estimated_params block. A model estimated
%   OUTSIDE Dynare (user solver, see solve_model_draw) has none, so this builds the
%   smallest pair of structures for which
%
%       set_all_parameters(x, estim_params_, M_)  does  M_.params(idx) = x
%
%   with idx the positions of param_names in M_.param_names (1:k when M_ is not
%   given or has no param_names, i.e. the draw IS M_.params), no shock standard
%   deviations, no measurement errors, no correlations. bayestopt_.name is
%   param_names, which is what run_reliability_draws checks the draws against and
%   what HLW_plots uses to label cell 9 of the wrapper output; the prior fields
%   (pshape = 0, bounds -Inf/Inf) are placeholders - the prior of such a model
%   lives outside Dynare, so prior draws must be handed to run_reliability_draws
%   through opts.draws.
%
%   See also: solve_model_draw, run_reliability_draws, dr_from_abcd.

param_names = cellstr(param_names(:));
k = numel(param_names);
if nargin >= 2 && isfield(M_, 'param_names') && ~isempty(M_.param_names)
    [tf, idx] = ismember(param_names, cellstr(M_.param_names));
    if ~all(tf)
        error('estim_params_from_names:unknown', 'parameter(s) not in M_.param_names: %s', strjoin(param_names(~tf), ', '));
    end
else
    idx = (1:k)';
end
idx = idx(:);

estim_params_ = struct();
estim_params_.nvx       = 0;
estim_params_.nvn       = 0;
estim_params_.ncx       = 0;
estim_params_.ncn       = 0;
estim_params_.nsx       = 0;          % Dynare >= 6 (skewness) - ignored by older versions
estim_params_.nendoinit = 0;          % Dynare >= 7 - ignored by older versions
estim_params_.np        = k;
estim_params_.var_exo   = zeros(0, 10);
estim_params_.var_endo  = zeros(0, 10);
estim_params_.corrx     = zeros(0, 11);
estim_params_.corrn     = zeros(0, 11);
estim_params_.param_vals = [idx, nan(k, 5)];   % column 1 is all set_all_parameters reads

bayestopt_ = struct();
bayestopt_.name   = param_names;
bayestopt_.pshape = zeros(k, 1);
bayestopt_.p1     = nan(k, 1);
bayestopt_.p2     = nan(k, 1);
bayestopt_.p3     = -inf(k, 1);
bayestopt_.p4     =  inf(k, 1);
bayestopt_.p5     = nan(k, 1);
bayestopt_.p6     = nan(k, 1);
bayestopt_.p7     = nan(k, 1);
bayestopt_.lb     = -inf(k, 1);
bayestopt_.ub     =  inf(k, 1);

if nargin >= 2
    if isfield(M_, 'exo_names')
        q = numel(M_.exo_names);
        if ~isfield(M_, 'Sigma_e'),               M_.Sigma_e = eye(q);            end
        if ~isfield(M_, 'Correlation_matrix'),    M_.Correlation_matrix = eye(q); end
        if ~isfield(M_, 'Skew_e'),                M_.Skew_e = zeros(0, 4);        end   % no skewed shocks (Dynare 7 csn_update_specification)
    end
    if isfield(M_, 'observable_names')
        p = numel(M_.observable_names);
        if ~isfield(M_, 'H'),                     M_.H = zeros(p);                   end
        if ~isfield(M_, 'Correlation_matrix_ME'), M_.Correlation_matrix_ME = eye(p); end
    end
else
    M_ = [];
end
end
