M_slow = M_;  M_slow.fast_mode = false;  M_slow.verbose = false;
M_fast = M_;  M_fast.fast_mode = true;   M_fast.verbose = false;

tic; S_slow = compute_state_reliability_diagnostics( ...
    M_.observable_names, M_.states_of_interest, oo_, M_slow, M_.data); t_slow = toc;
tic; S_fast = compute_state_reliability_diagnostics( ...
    M_.observable_names, M_.states_of_interest, oo_, M_fast, M_.data); t_fast = toc;

% steady_CorrY_H / steady_RMSE were removed with the steady-state data pass
% (2026-08); the TV objects are the only empirical statistics now.
check = {'CorrH_Observables','TV_CorrY_H','TV_RMSE'};
for k = 1:numel(check)
    a = S_slow.(check{k});  b = S_fast.(check{k});
    same = isequaln(a, b);              % isequaln so NaN == NaN
    fprintf('%-20s identical: %d   max|diff|: %.3g\n', check{k}, same, ...
        max(abs(a(:) - b(:)), [], 'omitnan'));
end
fprintf('slow %.1fs  fast %.1fs  speedup %.1fx\n', t_slow, t_fast, t_slow/t_fast);

%% see what is demanding (this revealed the 4a fix)
M_fast = M_;  M_fast.fast_mode = true;  M_fast.verbose = false;

profile clear
profile on -historysize 5000000
S = compute_state_reliability_diagnostics( ...
    M_.observable_names, M_.states_of_interest, oo_, M_fast, M_.data);
profile off

p = profile('info');
[~, ix] = sort([p.FunctionTable.TotalTime], 'descend');
fprintf('%10s %8s  %s\n','TotalTime','Calls','Function');
for k = ix(1:min(20, numel(ix)))
    fprintf('%10.3f %8d  %s\n', p.FunctionTable(k).TotalTime, ...
        p.FunctionTable(k).NumCalls, p.FunctionTable(k).FunctionName);
end