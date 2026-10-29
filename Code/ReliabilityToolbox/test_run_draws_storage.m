function test_run_draws_storage()
%TEST_RUN_DRAWS_STORAGE  Numeric storage of run_reliability_draws (ADDED 2026-10-27).
%
%   Fake wrappers (function handles) produce rows with the cell layout of reliability_draw_wrapper:
%   numeric matrices, a label cell, a struct with numeric/char/logical/cell fields, empty cells,
%   and failed draws written as the real wrapper writes them (scalar NaNs, {}, the draw itself,
%   nan(1,4), []). Checks:
%     1. storage 'numeric' and 'cell' return the same results (class, size and value of every
%        element, recursively), serial and with parfor (a plain loop without a pool / under
%        Octave), with the FIRST draw failing and with a row whose shape departs from the
%        template (stored as an exception);
%     2. a storage that does not fit into the available memory stops with
%        run_reliability_draws:memory before allocating;
%     3. when the wrapper throws an error half-way, the preallocated storage is released.
%        MATLAB keeps freed memory for reuse instead of returning it to the operating system,
%        so process memory does not fall back after the first failed call; the test therefore
%        repeats the failed call three times: with a leak each call would add its storage
%        again (about 470 MB), with the storage released the third call adds (almost) nothing.
%   Needs Dynare on the path (set_all_parameters); prints a pass message or errors.

here = fileparts(mfilename('fullpath'));
addpath(here);
if ~exist('set_all_parameters', 'file'), setup_dynare_path; end
warning('off', 'run_reliability_draws:memoryunknown');

% minimal Dynare structures: k estimated structural parameters, nothing else
k = 3;
M_ = struct('params', zeros(k, 1), 'param_names', {{'a'; 'b'; 'c'}}, 'H', 0, 'Correlation_matrix_ME', 1, ...
            'Sigma_e', zeros(0), 'Correlation_matrix', zeros(0), 'exo_nbr', 0, 'Skew_e', zeros(0, 4));   % Skew_e: Dynare 7
estim_params_ = struct('nvx', 0, 'nvn', 0, 'ncx', 0, 'ncn', 0, 'nsx', 0, 'nendoinit', 0, 'np', k, ...
                       'param_vals', [(1:k)', zeros(k, 9)]);
bayestopt_ = struct('name', {{'a'; 'b'; 'c'}});
options_ = struct('prior_trunc', 1e-10, 'mh_drop', 0.2);
oo_ = struct();

% ---- 1. numeric vs cell storage ----
if exist('rng', 'file'), rng(20261027); else, randn('state', 20261027); end %#ok<RAND>
n = 60;
D = randn(k, n);
D(1, 1) = -1;                      % the first draw fails: the template must come from a later one
D(1, 7) = abs(D(1, 7));  D(2, 7) = 3;   % draw 7 solves, but cell 2 gets an unexpected shape -> exception
for par = [false, true]
    base = struct('draws', D, 'function', @fake_wrapper, 'parallel', par, 'verbose', false);
    Rc = run_reliability_draws('prior', n, M_, options_, oo_, estim_params_, bayestopt_, [], [], setfield(base, 'storage', 'cell')); %#ok<SFLD>
    for cs = [7, 500]                 % 7: nine chunks, the last one partial; 500: one chunk
        o = base;  o.storage = 'numeric';  o.chunk_size = cs;
        Rn = run_reliability_draws('prior', n, M_, options_, oo_, estim_params_, bayestopt_, [], [], o);
        assert(isequal(size(Rc), size(Rn)), 'case 1 (parallel %d, chunk %d): sizes differ', par, cs);
        for i = 1:n
            for j = 1:size(Rc, 2)
                assert(equal_strict(Rc{i, j}, Rn{i, j}), 'case 1 (parallel %d, chunk %d): row %d, cell %d differs', par, cs, i, j);
            end
        end
    end
end
n_failed = sum(D(1, :) < 0);
assert(n_failed >= 5, 'case 1: only %d failed draws', n_failed);

% ---- 2. insufficient memory: error before allocating ----
avail = available_memory_bytes();
if isfinite(avail)
    n_big = ceil(4 * avail / 8 / 1e6) + 1;       % rows of 1e6 doubles: about 4x the available memory
    Dbig = abs(randn(k, n_big));
    t2 = tic;
    try
        run_reliability_draws('prior', n_big, M_, options_, oo_, estim_params_, bayestopt_, [], [], ...
                              struct('draws', Dbig, 'function', @fake_big_wrapper, 'parallel', false, 'verbose', false));
        error('test_run_draws_storage:nomemerror', 'case 2: expected run_reliability_draws:memory');
    catch err
        assert(strcmp(err.identifier, 'run_reliability_draws:memory'), 'case 2: wrong error %s (%s)', err.identifier, err.message);
    end
    assert(toc(t2) < 30, 'case 2: the memory error came late (%.1f s)', toc(t2));
    clear Dbig
else
    fprintf('test_run_draws_storage: available memory unknown on this platform; case 2 skipped.\n');
end

% ---- 3. storage released after an error inside the loop ----
m0 = process_memory_mb();
if isfinite(m0)
    n3 = 150;                               % 150 rows of 1e6 doubles (8 MB each), chunks of 50
    D3 = abs(randn(k, n3));  D3(3, 60) = -99; % draw 60 throws: chunk 1 complete (400 MB), chunk 2 partly filled
    m_err = NaN(1, 3);  m_after = NaN(1, 3);
    for rep = 1:3
        try
            run_reliability_draws('prior', n3, M_, options_, oo_, estim_params_, bayestopt_, [], [], ...
                                  struct('draws', D3, 'function', @fake_big_wrapper, 'parallel', false, 'verbose', false, ...
                                         'check_memory', false, 'chunk_size', 50));
            error('test_run_draws_storage:noerror', 'case 3: the wrapper error was not raised');
        catch err
            assert(strcmp(err.identifier, 'fake:fail'), 'case 3: wrong error %s (%s)', err.identifier, err.message);
            m_err(rep) = str2double(err.message);
        end
        m_after(rep) = process_memory_mb();
    end
    fprintf('test_run_draws_storage: process memory %.0f MB before; at the error / after it: %s MB\n', m0, ...
            strjoin(arrayfun(@(a, b) sprintf('%.0f / %.0f', a, b), m_err, m_after, 'UniformOutput', false), ', '));
    assert(m_err(1) - m0 > 300, 'case 3: storage did not grow as expected (%.0f MB)', m_err(1) - m0);
    assert(m_after(3) - m_after(2) < 150, 'case 3: the third failed call added %.0f MB: storage not released', m_after(3) - m_after(2));
else
    fprintf('test_run_draws_storage: process memory unknown on this platform; case 3 skipped.\n');
end

fprintf('test_run_draws_storage: all checks passed (%d draws, %d failed, serial and parfor).\n', n, n_failed);
end


function row = fake_wrapper(x, M_, varargin)
% rows with the layout of reliability_draw_wrapper; x(1) < 0 = failed draw
row = cell(1, 8);
if x(1) < 0
    row{1} = NaN;  row{2} = NaN;  row{3} = x;  row{4} = {};  row{5} = NaN;  row{6} = nan(1, 4);  % row{7}, row{8} stay []
    return
end
row{1} = x(1) * ones(1, 5);
row{2} = reshape(1:12, 3, 4) * x(2);
if x(2) > 2, row{2} = [1 2; 3 4]; end              % departs from the template -> exception
row{3} = M_.params;
row{4} = {'aa'; 'bbb'};
row{5} = struct('R2', x(:)', 'lbl', {{'x', 'yy'}}, 'flag', logical([1 0]), 's', 'abc', 'e', [], 'nan', [NaN 1]);
row{6} = [1, x(3), -0, 4];
row{7} = [];
row{8} = x(3);
end

function row = fake_big_wrapper(x, varargin)
% 8 MB rows; x(3) = -99 throws, with the process memory (MB) in the message
if x(3) == -99, error('fake:fail', '%.1f', process_memory_mb()); end
row = {x(1) * ones(1e6, 1), [x(2), 1]};            % cell 2 non-scalar: a draw that solved
end

function tf = equal_strict(a, b)
% same class, size and value (NaN == NaN), recursively through cells and structs
tf = strcmp(class(a), class(b)) && isequal(size(a), size(b));
if ~tf, return; end
if iscell(a)
    for e = 1:numel(a), if ~equal_strict(a{e}, b{e}), tf = false; return; end, end
elseif isstruct(a)
    tf = isequal(fieldnames(a), fieldnames(b));
    if ~tf, return; end
    fl = fieldnames(a);
    for e = 1:numel(a)
        for q = 1:numel(fl)
            if ~equal_strict(a(e).(fl{q}), b(e).(fl{q})), tf = false; return; end
        end
    end
else
    tf = isequaln(a, b);
end
end

function mb = process_memory_mb()
% memory of this MATLAB/Octave process in MB (Windows: memory(); macOS/Linux: resident set size)
mb = NaN;
try
    if ispc && exist('memory') > 0                                          %#ok<EXIST>
        u = memory;  mb = u.MemUsedMATLAB / 2^20;
    else
        if exist('OCTAVE_VERSION', 'builtin'), pid = getpid(); else, pid = feature('getpid'); end
        [status, out] = system(sprintf('ps -o rss= -p %d', pid));
        if status == 0, mb = str2double(strtrim(out)) / 1024; end
    end
catch
    mb = NaN;
end
end
