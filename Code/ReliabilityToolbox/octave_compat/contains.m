function tf = contains(str, pattern, varargin)
%CONTAINS  Octave stand-in for MATLAB's contains (ADDED 2026-09-17).
%   tf = CONTAINS(str, pattern) is true where the char/cellstr str contains
%   any of the char/cellstr patterns. 'IgnoreCase', true is honoured.
%   Lives in ReliabilityToolbox/octave_compat, which the plotting scripts add
%   to the path ONLY under Octave (MATLAB has the real function).
ignore_case = false;
for k = 1:2:numel(varargin)
    if strcmpi(varargin{k}, 'IgnoreCase'), ignore_case = logical(varargin{k+1}); end
end
if ischar(str), str = {str}; was_char = true; else, was_char = false; end
if ischar(pattern), pattern = {pattern}; end
if ignore_case, str = lower(str); pattern = lower(pattern); end
tf = false(size(str));
for i = 1:numel(str)
    for p = 1:numel(pattern)
        if ~isempty(strfind(str{i}, pattern{p})), tf(i) = true; break; end
    end
end
if was_char, tf = tf(1); end
end
