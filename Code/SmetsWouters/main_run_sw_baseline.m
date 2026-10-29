% MAIN_RUN_SW_BASELINE  Smets-Wouters (2007) baseline: estimation (short MH chain), reliability
% draws, model assessment, plots and tables, without measurement error.
%
%   Outputs in Code/SmetsWouters/SW_baseline/. Chain length and draw counts: see
%   run_sw_variant.m (environment variables SW_MH_REPLIC, SW_NDRAWS, SW_MA_D).
setenv('SW_ME', '0');
run(fullfile(fileparts(mfilename('fullpath')), 'run_sw_variant.m'));
