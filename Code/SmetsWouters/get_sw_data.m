% GET_SW_DATA  Download the Smets-Wouters (2007) data and published mode into Code/SmetsWouters.
%
%   Files (from Johannes Pfeifer's DSGE_mod replication of Smets and Wouters 2007,
%   github.com/JohannesPfeifer/DSGE_mod/tree/master/Smets_Wouters_2007):
%       usmodel_data.mat   observables dy dc dinve labobs pinfobs dw robs, 1947Q3-2004Q4
%       usmodel_mode.mat   published posterior mode, the starting point of the mode search
%   Files already present are not downloaded again. run_sw_variant.m calls this script.
%   Without internet access, download the two files from the URLs below by hand
%   and place them in this folder.

sw_dir = fileparts(mfilename('fullpath'));
if isempty(sw_dir), sw_dir = pwd; end
sw_url = 'https://raw.githubusercontent.com/JohannesPfeifer/DSGE_mod/master/Smets_Wouters_2007/';
sw_files = {'usmodel_data.mat', 'usmodel_mode.mat'};
for sw_k = 1:numel(sw_files)
    sw_target = fullfile(sw_dir, sw_files{sw_k});
    if exist(sw_target, 'file')
        fprintf('get_sw_data: %s present\n', sw_target);
        continue
    end
    fprintf('get_sw_data: downloading %s%s\n', sw_url, sw_files{sw_k});
    if exist('OCTAVE_VERSION', 'builtin')
        urlwrite([sw_url sw_files{sw_k}], sw_target);
    else
        websave(sw_target, [sw_url sw_files{sw_k}]);
    end
end
sw_check = load(fullfile(sw_dir, 'usmodel_data.mat'));
fprintf('get_sw_data: usmodel_data.mat has %d quarters of %s\n', numel(sw_check.dy), strjoin(intersect(fieldnames(sw_check), {'dy', 'dc', 'dinve', 'labobs', 'pinfobs', 'dw', 'robs'})', ' '));
clear sw_dir sw_url sw_files sw_k sw_target sw_check
