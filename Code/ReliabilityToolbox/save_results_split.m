function save_results_split(fname, max_part_bytes)
%SAVE_RESULTS_SPLIT  Save the caller's workspace with the reliability draws split into size-safe files.
%
%   SAVE_RESULTS_SPLIT('X.mat')  saves every variable of the calling workspace
%   (like `save X.mat`) EXCEPT that oo_.prior_function_results and
%   oo_.posterior_function_results are removed from oo_ and written to
%   X_prior.mat and X_posterior.mat as chunked cell arrays
%       results_part_1, ..., results_part_<n_parts>   (vertical slices of draws)
%       n_parts, n_draws, n_cells, set_name
%   each part below max_part_bytes (default 1 GB). All files are MATLAB v7.
%   LOAD_RESULTS_SPLIT reverses this and hands back oo_ with the two fields
%   reattached, so downstream code (HLW_plots) is unchanged.
%
%   WHY (2026-09-17). Octave's `save` writes MAT-5 (v7) files, whose per-element
%   size field is 32 bits: a variable above 4 GB is written CORRUPT without any
%   warning (UC3 run of 2026-09-17: oo_ = 4.7 GB with 22 diagnosed states x
%   20 000 draws -> `load: invalid element type`, six hours of draws lost), and
%   MATLAB refuses to read variables above 2 GB from v7 files (HLW's 3.3 GB
%   oo_). The draws are the only large objects - the per-draw cells carry the
%   filtered state paths and KF standard deviations (cells 13, 14, 16) - so
%   splitting them by draws keeps every variable small, and every other
%   variable is checked against the 2 GB limit before writing.
%
%   See also: load_results_split, run_reliability_draws, HLW_plots.

if nargin < 2 || isempty(max_part_bytes), max_part_bytes = 1e9; end
limit_bytes = 2^31 - 1;

% ---- collect the caller's variables (as `save` would) ----
W = evalin('caller', 'whos');
names = {W.name};
names = names(~strcmp(names, 'ans'));
S = struct();
for k = 1:numel(names)
    S.(names{k}) = evalin('caller', names{k});
end

[fdir, fbase] = fileparts(fname);
if isempty(fdir), fdir = pwd; end

% ---- split the two result sets out of oo_ ----
for set_name = {'prior', 'posterior'}
    fld = [set_name{1} '_function_results'];
    if isfield(S, 'oo_') && isfield(S.oo_, fld) && ~isempty(S.oo_.(fld))
        R = S.oo_.(fld);
        S.oo_ = rmfield(S.oo_, fld);
        write_chunked(fullfile(fdir, [fbase '_' set_name{1} '.mat']), R, set_name{1}, max_part_bytes);
    end
end

% ---- guard: nothing left above the v7 per-variable limit ----
fn = fieldnames(S);
for k = 1:numel(fn)
    tmp = S.(fn{k});  w = whos('tmp');
    if w.bytes > limit_bytes
        error('save_results_split:toolarge', ...
              'variable %s is %.2f GB, above the 2 GB MATLAB v7 per-variable limit; it would be written corrupt by Octave above 4 GB.', fn{k}, w.bytes / 1e9);
    end
end
save(fullfile(fdir, [fbase '.mat']), '-struct', 'S', '-v7');
fprintf('save_results_split: wrote %s (%d variables)\n', fullfile(fdir, [fbase '.mat']), numel(fn));
end


function write_chunked(file, R, set_name, max_part_bytes)
n_draws = size(R, 1);  n_cells = size(R, 2);
probe = R(1:min(n_draws, 50), :);  w = whos('probe');
bytes_per_draw = max(w.bytes / size(probe, 1), 1);
per_part = max(1, floor(max_part_bytes / bytes_per_draw));
n_parts  = ceil(n_draws / per_part);
C = struct('n_parts', n_parts, 'n_draws', n_draws, 'n_cells', n_cells, 'set_name', set_name);
for p = 1:n_parts
    rows = (p - 1) * per_part + 1 : min(p * per_part, n_draws);
    C.(sprintf('results_part_%d', p)) = R(rows, :);
end
save(file, '-struct', 'C', '-v7');
fprintf('save_results_split: wrote %s (%d draws in %d part(s), ~%.0f MB each)\n', file, n_draws, n_parts, per_part * bytes_per_draw / 1e6);
end
