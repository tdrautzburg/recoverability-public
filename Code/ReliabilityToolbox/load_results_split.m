function S = load_results_split(fname)
%LOAD_RESULTS_SPLIT  Load a results file written by save_results_split (or a plain one).
%
%   S = LOAD_RESULTS_SPLIT('X.mat') returns a struct of the variables in X.mat;
%   if X_prior.mat / X_posterior.mat exist next to it, their chunked draws are
%   concatenated and reattached as S.oo_.prior_function_results /
%   S.oo_.posterior_function_results. A monolithic results file (pre-2026-09-17,
%   e.g. HLW_results_no_lags.mat) loads unchanged.
%
%   LOAD_RESULTS_SPLIT('X.mat') with no output assigns the variables into the
%   CALLER's workspace, like `load X.mat` (used by HLW_plots).
%
%   See also: save_results_split.

S = load(fname);
[fdir, fbase] = fileparts(fname);
for set_name = {'prior', 'posterior'}
    side = fullfile(fdir, [fbase '_' set_name{1} '.mat']);
    if exist(side, 'file') == 2
        C = load(side);
        parts = cell(C.n_parts, 1);
        for p = 1:C.n_parts
            parts{p} = C.(sprintf('results_part_%d', p));
        end
        R = vertcat(parts{:});
        if size(R, 1) ~= C.n_draws || size(R, 2) ~= C.n_cells
            error('load_results_split:size', '%s: reassembled %d x %d, header says %d x %d.', side, size(R, 1), size(R, 2), C.n_draws, C.n_cells);
        end
        if ~isfield(S, 'oo_'), S.oo_ = struct(); end
        S.oo_.([set_name{1} '_function_results']) = R;
        fprintf('load_results_split: %s -> oo_.%s_function_results (%d draws)\n', side, set_name{1}, C.n_draws);
    end
end

if nargout == 0
    fn = fieldnames(S);
    for k = 1:numel(fn)
        assignin('caller', fn{k}, S.(fn{k}));
    end
    clear S
end
end
