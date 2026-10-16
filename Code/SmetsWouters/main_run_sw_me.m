% MAIN_RUN_SW_ME  Smets-Wouters (2007) re-estimated with CALIBRATED measurement error on dy, dc,
% dinve, dw and labobs: estimation (short MH chain), reliability draws, model assessment, plots
% and tables.
%
%   The measurement-error standard deviations are fixed, not estimated: by default
%   sqrt(SW_ME_SHARE) times the sample s.d. of each series, SW_ME_SHARE = 0.10; SW_ME_SD sets one
%   absolute s.d. for all five instead (see Smets_Wouters_2007_45.mod). Outputs in
%   Code/SmetsWouters/SW_ME/. Chain length and draw counts: see run_sw_variant.m.
setenv('SW_ME', '1');
run(fullfile(fileparts(mfilename('fullpath')), 'run_sw_variant.m'));
