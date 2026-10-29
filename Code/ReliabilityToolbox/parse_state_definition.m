function [coefs, vars] = parse_state_definition(def, endo_names)
%PARSE_STATE_DEFINITION  Read the algebraic definition of a state of interest as a linear combination of model variables.
%
%   [coefs, vars] = PARSE_STATE_DEFINITION(def)              syntax check only
%   [coefs, vars] = PARSE_STATE_DEFINITION(def, endo_names)  also checks that every variable exists
%
%   ADDED 2026-09-30. def is the optional 4th column of M_.CellList_States, e.g.
%       { 'ygap', [0 1 4 20], [], 'y - yf' }
%   a sum of terms  [+|-] [number*] variable  such as 'y - yf', '-r + 0.25*pinf', '2*a'.
%   Only constant numeric coefficients are allowed (no parameter names, no parentheses):
%   the engine turns the definition into a row of the decision rule, coef' * [ghx; ghu]
%   of the named variables, so the combination is exact for any parameter draw.
%   Returns the coefficients (column) and the variable names (cellstr); repeated
%   variables are summed.
%
%   See also: expand_states_of_interest, compute_state_reliability_diagnostics.

if ~ischar(def) || isempty(strtrim(def))
    error('parse_state_definition:empty', 'a state definition must be a non-empty character vector.');
end
s = regexprep(def, '\s+', '');
if any(s == '(') || any(s == ')')
    error('parse_state_definition:syntax', 'definition "%s": parentheses are not supported (write the sum out).', def);
end
% split into signed terms by scanning: a '+' or '-' starts a new term unless it is the sign of an
% exponent, i.e. it follows an 'e'/'E' that itself follows a digit or '.' (as in 1e-3)
cut = false(1, numel(s));
for p = 2:numel(s)
    if any(s(p) == '+-')
        is_exp = p >= 3 && any(s(p-1) == 'eE') && (isstrprop(s(p-2), 'digit') || s(p-2) == '.');
        cut(p) = ~is_exp;
    end
end
starts = [1, find(cut)];  ends = [find(cut) - 1, numel(s)];
vars = {};  coefs = [];
for k = 1:numel(starts)
    t = s(starts(k):ends(k));  sg = 1;
    if t(1) == '+', t = t(2:end); elseif t(1) == '-', sg = -1; t = t(2:end); end
    if isempty(t), error('parse_state_definition:syntax', 'definition "%s": dangling sign.', def); end
    parts = strsplit(t, '*');
    isname = cellfun(@(x) ~isempty(regexp(x, '^[A-Za-z_]\w*$', 'once')), parts);
    if numel(parts) == 1 && isname
        c = 1;  v = parts{1};
    elseif numel(parts) == 2 && sum(isname) == 1
        v = parts{isname};  c = str2double(parts{~isname});
        if isnan(c), error('parse_state_definition:syntax', 'definition "%s": term "%s" needs a numeric coefficient.', def, t); end
    else
        error('parse_state_definition:syntax', 'definition "%s": term "%s" is not [number*]variable.', def, t);
    end
    j = find(strcmp(vars, v), 1);
    if isempty(j), vars{end+1} = v;  coefs(end+1) = sg * c; else, coefs(j) = coefs(j) + sg * c; end   %#ok<AGROW>
end
if isempty(vars), error('parse_state_definition:syntax', 'definition "%s" names no variable.', def); end
coefs = coefs(:);  vars = vars(:);
if nargin >= 2 && ~isempty(endo_names)
    miss = vars(~ismember(vars, endo_names));
    if ~isempty(miss)
        error('parse_state_definition:unknownVariable', 'definition "%s": not endogenous variable(s) of the model: %s.', def, strjoin(miss', ', '));
    end
end
end
