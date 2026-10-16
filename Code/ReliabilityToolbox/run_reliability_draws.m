function [results, draws, info] = run_reliability_draws(type, n_draws, M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info, opts)
%RUN_RELIABILITY_DRAWS  prior_function / posterior_function done by hand, with parfor over draws.
%
%   [results, draws, info] = RUN_RELIABILITY_DRAWS('prior',     n, M_, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info)
%   [results, draws, info] = RUN_RELIABILITY_DRAWS('posterior', n, ...)
%   [...]                   = RUN_RELIABILITY_DRAWS(..., opts)
%
%   Replaces Dynare's
%       prior_function(function='reliability_draw_wrapper', sampling_draws=n);
%       posterior_function(function='reliability_draw_wrapper', sampling_draws=n);
%   which run the wrapper in a SERIAL loop inside execute_prior_posterior_function
%   (and, for the posterior, need the MATLAB-written *_mh_history_0.mat that
%   Octave cannot read). This function draws the parameter vectors itself and
%   evaluates the wrapper in a PARFOR over draws - each draw is independent
%   (set_all_parameters -> stoch_simul -> engine), so the speed-up is close to
%   the number of workers. Under Octave, or with no parallel pool, parfor runs
%   as a plain for.
%
%   DRAWS
%     'prior'      Dynare's own prior sampler: dprior(bayestopt_, options_.prior_trunc).draws(n)
%                  (prior_draw as a fallback on older Dynare).
%     'posterior'  read directly from the MH chain files
%                  <M_.dname>/metropolis/<M_.fname>_mh<j>_blck<b>.mat: the x2
%                  blocks of every file are concatenated per chain in file
%                  order, the first mh_drop share of EACH chain is discarded
%                  (as Dynare's CutSample does), the chains are pooled and n
%                  rows are drawn at random without replacement (with
%                  replacement, with a warning, if n exceeds the pool). Set
%                  opts.mh_dir when the metropolis folder is not at
%                  <M_.dname>/metropolis relative to the current directory.
%     opts.draws   ADDED 2026-09-26: a k x n_avail matrix of parameter vectors
%                  supplied by the caller (rows in bayestopt_.name order)
%                  REPLACES the sampler for either type: a model estimated
%                  outside Dynare (user solver, see solve_model_draw) hands its
%                  own prior/posterior draws in this way. n_draws of the
%                  columns are used - all of them (n_draws := n_avail, with a
%                  warning if more were asked for) or a random subset without
%                  replacement when n_draws < n_avail.
%   Draws are reproducible through opts.seed.
%
%   OUTPUTS
%     results   n x n_cells cell, one wrapper output row per draw - the same
%               layout as oo_.<type>_function_results, so downstream scripts
%               can do  oo_.posterior_function_results = results;
%     draws     k x n matrix of the parameter vectors (bayestopt_.name order)
%     info      struct: .type .n_draws .n_failed (rows whose cell 2 is a scalar
%               NaN, i.e. the solver failed) .seed .mh_drop .mh_files
%               .elapsed_s .n_cells .function
%
%   OPTIONS (struct, all optional)
%     .function   wrapper name ('reliability_draw_wrapper')
%     .draws      k x n_avail caller-supplied parameter draws (see DRAWS); [] = sample
%     .seed       RNG seed (20260915); the wrapper itself is deterministic
%     .mh_drop    share of each chain discarded (options_.mh_drop)
%     .mh_dir     metropolis folder (fullfile(M_.dname, 'metropolis'))
%     .parallel   true (default): parfor; false: for (useful for debugging)
%     .verbose    true: progress line per 50 draws in the serial path
%
%   The .mod's post-save block (M_.observable_names, M_.data, M_.CellList_States,
%   M_.states_of_interest, M_.var_list_, X_InitialDistribution, the
%   oo_.bayestopt_* fields, oo_.dr.qz_crit) must have been run first, exactly
%   as for Dynare's prior_function / posterior_function.
%
%   See also: reliability_draw_wrapper, execute_prior_posterior_function,
%             manual_posterior_functions.

if nargin < 10, opts = struct(); end
def = struct('function', 'reliability_draw_wrapper', 'seed', 20260915, 'mh_drop', [], 'mh_dir', '', ...
             'parallel', true, 'verbose', true, 'draws', []);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}), opts.(fn{i}) = def.(fn{i}); end
end
if ~any(strcmpi(type, {'prior', 'posterior'}))
    error('run_reliability_draws:type', 'type must be ''prior'' or ''posterior''.');
end
type = lower(type);
if ~exist(opts.function, 'file')
    error('run_reliability_draws:function', 'wrapper %s.m is not on the path.', opts.function);
end
f = str2func(opts.function);
t0 = tic;

try, rng(opts.seed); catch, rand('seed', opts.seed); randn('seed', opts.seed); end

%% ------------------------------------------------------------ parameter draws
info = struct('type', type, 'n_draws', n_draws, 'seed', opts.seed, 'function', opts.function, 'mh_drop', NaN, 'mh_files', {{}}, 'source', 'sampler');
if ~isempty(opts.draws)
    % caller-supplied draws (ADDED 2026-09-26): no sampling, no MH files
    info.source = 'user';
    n_avail = size(opts.draws, 2);
    if n_draws >= n_avail
        if n_draws > n_avail
            warning('run_reliability_draws:fewer', 'n_draws = %d exceeds the %d supplied draws; using all of them.', n_draws, n_avail);
        end
        n_draws = n_avail;
        draws   = opts.draws;
    else
        draws = opts.draws(:, randperm(n_avail, n_draws));
    end
    info.n_draws = n_draws;
    type_src = 'user';
else
    type_src = type;
end
switch type_src
    case 'user'
        % draws already set
    case 'prior'
        if exist('dprior', 'class') || exist('dprior', 'file')
            Prior = dprior(bayestopt_, options_.prior_trunc);
            draws = Prior.draws(n_draws);
        else
            prior_draw(bayestopt_, options_.prior_trunc);   % initialisation call (older Dynare)
            draws = zeros(numel(bayestopt_.name), n_draws);
            for i = 1:n_draws
                draws(:, i) = prior_draw();
            end
        end
    case 'posterior'
        if isempty(opts.mh_drop), opts.mh_drop = options_.mh_drop; end
        if isempty(opts.mh_dir),  opts.mh_dir  = fullfile(M_.dname, 'metropolis'); end
        info.mh_drop = opts.mh_drop;
        files = dir(fullfile(opts.mh_dir, [M_.fname '_mh*_blck*.mat']));
        if isempty(files)
            error('run_reliability_draws:mh', 'no MH chain files %s_mh*_blck*.mat in %s.', M_.fname, opts.mh_dir);
        end
        tok = regexp({files.name}, '_mh(\d+)_blck(\d+)\.mat$', 'tokens', 'once');
        mhno = cellfun(@(t) str2double(t{1}), tok);
        blck = cellfun(@(t) str2double(t{2}), tok);
        pool = [];
        for b = unique(blck)
            sel = find(blck == b);
            [~, ord] = sort(mhno(sel));
            chain = [];
            for s = sel(ord)
                S = load(fullfile(opts.mh_dir, files(s).name), 'x2');
                chain = [chain; S.x2];                                  %#ok<AGROW>
                info.mh_files{end+1} = files(s).name;                   %#ok<AGROW>
            end
            chain = chain(floor(opts.mh_drop * size(chain, 1)) + 1 : end, :);   % CutSample: drop the burn-in of each chain
            pool = [pool; chain];                                       %#ok<AGROW>
        end
        if n_draws <= size(pool, 1)
            idx = randperm(size(pool, 1), n_draws);
        else
            warning('run_reliability_draws:replacement', ...
                'n_draws = %d exceeds the %d retained MH draws; sampling with replacement.', n_draws, size(pool, 1));
            idx = randi(size(pool, 1), n_draws, 1);
        end
        draws = pool(idx, :)';
end
k = size(draws, 1);
if k ~= numel(bayestopt_.name)
    error('run_reliability_draws:dim', 'draws have %d parameters, bayestopt_ names %d.', k, numel(bayestopt_.name));
end

%% ---------------------------------------------------- output width from draw 1
M_1  = set_all_parameters(draws(:, 1), estim_params_, M_);
row1 = f(draws(:, 1), M_1, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info);
n_cells = numel(row1);
results = cell(n_draws, n_cells);
results(1, :) = row1;

%% ---------------------------------------------------------- the draw loop
if opts.parallel
    parfor i = 2:n_draws
        M_i = set_all_parameters(draws(:, i), estim_params_, M_);
        results(i, :) = f(draws(:, i), M_i, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info);
    end
else
    for i = 2:n_draws
        M_i = set_all_parameters(draws(:, i), estim_params_, M_);
        results(i, :) = f(draws(:, i), M_i, options_, oo_, estim_params_, bayestopt_, dataset_, dataset_info);
        if opts.verbose && mod(i, 50) == 0
            fprintf('\n  run_reliability_draws: %d / %d draws (%.0fs)\n', i, n_draws, toc(t0));
        end
    end
end

info.n_cells   = n_cells;
info.n_failed  = sum(cellfun(@(c) isscalar(c) && isnumeric(c), results(:, 2)));
info.elapsed_s = toc(t0);
if opts.verbose
    fprintf('\nrun_reliability_draws (%s): %d draws, %d failed solves, %d cells, %.0f s\n', ...
            type, n_draws, info.n_failed, n_cells, info.elapsed_s);
end
end
