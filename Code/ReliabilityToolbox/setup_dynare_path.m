function info = setup_dynare_path(opts)
%SETUP_DYNARE_PATH  Put Dynare's matlab/ folder (and, under MATLAB, its compiled MEX files) on the path.
%
%   info = SETUP_DYNARE_PATH()          info = SETUP_DYNARE_PATH(struct('root', '/Applications/Dynare/7.2-arm64'))
%
%   ADDED 2026-09-28. The Homebrew Dynare (/opt/homebrew/opt/dynare) ships the m-files
%   and OCTAVE-only MEX files; under MATLAB its kalman_steady_state / mjdgges MEX are
%   absent, so the toolbox falls back to its own QZ Riccati solver (measured: the
%   population engine ran at 103 ms under MATLAB vs 66 ms under Octave with the MEX).
%   The official Dynare macOS package (dynare.org/download, the arm64 .pkg) installs
%   /Applications/Dynare/<version>-arm64/ with matlab/ and mex/matlab/<MATLAB version range>/.
%   Search order for the Dynare root:
%     1. opts.root, then the environment variable DYNARE_ROOT (a folder holding matlab/)
%     2. under MATLAB: the newest /Applications/Dynare/*/ (official package)
%     3. /opt/homebrew/opt/dynare/lib/dynare (Homebrew; complete under Octave, no MEX under MATLAB)
%   Then addpath(<root>/matlab), dynare_config, and under MATLAB, if kalman_steady_state
%   is still not a MEX file (dynare_config adds only the mex/matlab folder matching its
%   own MATLAB version range; a newer MATLAB is not matched although the MEX load), every
%   <root>/mex/matlab/* folder is added. A warning names the fallback when no MEX is found.
%   Returns info.root, info.matlab_dir, info.is_octave, info.mex_kalman_steady_state,
%   info.mex_mjdgges, info.added (path folders added here).
%   Nothing is done twice: an existing Dynare path is kept (but the MEX check still runs).

if nargin < 1, opts = struct(); end
is_octave = exist('OCTAVE_VERSION', 'builtin') > 0;
info = struct('root', '', 'matlab_dir', '', 'is_octave', is_octave, 'mex_kalman_steady_state', false, 'mex_mjdgges', false, 'added', {{}});
root = '';
if isfield(opts, 'root') && ~isempty(opts.root), root = opts.root; end
if isempty(root), root = getenv('DYNARE_ROOT'); end
if isempty(root) && ~is_octave
    cands = dir('/Applications/Dynare');
    cands = cands([cands.isdir] & ~ismember({cands.name}, {'.', '..'}));
    if ~isempty(cands)
        names = sort({cands.name});  root = fullfile('/Applications/Dynare', names{end});   % newest version string last
    end
end
if isempty(root) && exist('/opt/homebrew/opt/dynare/lib/dynare/matlab', 'dir'), root = '/opt/homebrew/opt/dynare/lib/dynare'; end
if isempty(root)
    if exist('dynare_config', 'file'), root = fileparts(fileparts(which('dynare_config')));
    else, error('setup_dynare_path:notfound', 'no Dynare installation found (set DYNARE_ROOT or install Dynare).'); end
end
info.root = root;  info.matlab_dir = fullfile(root, 'matlab');
if ~exist('dynare_config', 'file') || ~strcmp(fileparts(which('dynare_config')), info.matlab_dir)
    if ~exist(info.matlab_dir, 'dir'), error('setup_dynare_path:matlab', '%s has no matlab/ folder.', root); end
    addpath(info.matlab_dir);  info.added{end+1} = info.matlab_dir;
    dynare_config;
end
% MEX files under MATLAB: add the mex/matlab/* folders that dynare_config skipped
if ~is_octave && ~(exist('kalman_steady_state', 'file') == 3)
    mexdirs = dir(fullfile(root, 'mex', 'matlab'));
    mexdirs = mexdirs([mexdirs.isdir] & ~ismember({mexdirs.name}, {'.', '..'}));
    for i = 1:numel(mexdirs)
        d = fullfile(root, 'mex', 'matlab', mexdirs(i).name);
        if ~isempty(dir(fullfile(d, ['*.' mexext])))
            addpath(d);  info.added{end+1} = d;
        end
    end
end
info.mex_kalman_steady_state = exist('kalman_steady_state', 'file') == 3;
info.mex_mjdgges             = exist('mjdgges', 'file') == 3;
if ~info.mex_kalman_steady_state
    warning('setup_dynare_path:mex', ['Dynare''s kalman_steady_state MEX is not available in %s (%s): the steady-state Riccati ' ...
            'equations use the toolbox''s own QZ solver (slower). Under MATLAB install the official Dynare macOS package ' ...
            '(dynare.org/download) or set DYNARE_ROOT to it.'], root, iff(is_octave, 'Octave', 'MATLAB'));
end
end

function s = iff(c, a, b)
if c, s = a; else, s = b; end
end
