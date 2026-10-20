% PROFILE_POSTERIOR_DRAW  Time and profile single posterior draws of UC3-i2 (diagnostic for a slow machine).
%
%   TEMPORARY (added 2026-10-20): a debugging aid for slow posterior draws on one test machine.
%   Not part of the replication workflow. DELETE THIS FILE once the debugging is complete.
%
%   Run from Code/UC3 of recoverability-public AFTER main_run_uc3_levels_i2 has passed the prior
%   stage (it needs UC3_i2_results.mat and uc3_levels_i2/metropolis/*_blck*.mat), e.g.
%       matlab -batch "profile_posterior_draw"            (default threads)
%       matlab -singleCompThread -batch "profile_posterior_draw"
%   Prints: versions and threads, whether the 2026-10-20 solver fix is on the path, the time of
%   five posterior draws with the Riccati method each used, and the 25 most expensive functions
%   of one draw (MATLAB profiler, self time). Writes nothing.

uc3_dir = pwd;
addpath(fullfile(uc3_dir, '..', 'ReliabilityToolbox'));      % not genpath (octave_compat stays off)
setup_dynare_path;

fprintf('\n==== environment ====\n');
fprintf('MATLAB %s on %s; maxNumCompThreads = %d\n', version, computer, maxNumCompThreads);
dynare_version
ssr = which('solve_steady_state_riccati');
txt = fileread(ssr);
fprintf('solve_steady_state_riccati: %s\n  2026-10-20 ordeig fix present: %d\n', ssr, contains(txt, 'FIXED 2026-10-20'));
fprintf('kalman_steady_state MEX: %d (3 = MEX file)\n', exist('kalman_steady_state', 'file'));
pinfo = ver('parallel');
fprintf('Parallel Computing Toolbox installed: %d\n', ~isempty(pinfo));

fprintf('\n==== workspace and one MCMC chain ====\n');
W = load(fullfile(uc3_dir, 'UC3_i2_results.mat'), 'M_', 'options_', 'oo_', 'estim_params_', 'bayestopt_', 'dataset_', 'dataset_info');
mh = dir(fullfile(uc3_dir, 'uc3_levels_i2', 'metropolis', 'uc3_levels_i2_mh*_blck1.mat'));
C = load(fullfile(mh(end).folder, mh(end).name), 'x2');
fprintf('chain file %s, %d draws\n', mh(end).name, size(C.x2, 1));
rows = size(C.x2, 1) - (0:4);

fprintf('\n==== five posterior draws ====\n');
for k = 1:numel(rows)
    x = C.x2(rows(k), :)';
    M_k = set_all_parameters(x, W.estim_params_, W.M_);
    t0 = tic;
    out = reliability_draw_wrapper(x, M_k, W.options_, W.oo_, W.estim_params_, W.bayestopt_, W.dataset_, W.dataset_info);
    t = toc(t0);
    rf = out{19};                                            % [converged, iterations, rel_step, rho_cl]
    fprintf('draw %d: %.2f s; Riccati converged %d, iterations %d (0 = direct solve), rho(A-KC) %.4f\n', k, t, rf(1), rf(2), rf(4));
end

fprintf('\n==== profile of one posterior draw (top 25 by self time) ====\n');
x = C.x2(rows(1), :)';
M_k = set_all_parameters(x, W.estim_params_, W.M_);
profile clear
profile on
reliability_draw_wrapper(x, M_k, W.options_, W.oo_, W.estim_params_, W.bayestopt_, W.dataset_, W.dataset_info);
profile off
P = profile('info');
FT = P.FunctionTable;
self_t = zeros(1, numel(FT));
for j = 1:numel(FT)
    ch = FT(j).Children;
    if isempty(ch), self_t(j) = FT(j).TotalTime; else, self_t(j) = FT(j).TotalTime - sum([ch.TotalTime]); end
end
[~, ord] = sort(self_t, 'descend');
fprintf('%-60s %10s %10s %8s\n', 'function', 'self [s]', 'total [s]', 'calls');
for j = ord(1:min(25, numel(ord)))
    fprintf('%-60s %10.3f %10.3f %8d\n', FT(j).FunctionName(1:min(end, 60)), self_t(j), FT(j).TotalTime, FT(j).NumCalls);
end
