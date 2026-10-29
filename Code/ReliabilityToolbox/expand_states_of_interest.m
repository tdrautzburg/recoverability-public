function [labels, sp] = expand_states_of_interest(s)
%EXPAND_STATES_OF_INTEREST  Diagnosed-state names, in engine order, from M_.CellList_States.
%
%   labels       = EXPAND_STATES_OF_INTEREST(M_.CellList_States)
%   [labels, sp] = EXPAND_STATES_OF_INTEREST(...)
%
%   M_.CellList_States is the J x 3 cell (a J x 2 cell is read as lags only)
%
%       { 'ygap',    [0 1 4 20], []
%         'rrstar',  [0 1 4 20], [4 20]
%         'drrstar', [0 1],      []    }
%
%   column 1: model variable; column 2: fixed-lag smoother lags m (0 = the
%   filtered state itself); column 3: change horizons h; OPTIONAL column 4
%   (ADDED 2026-09-30): the algebraic definition of a state that is not a
%   model variable, as a linear combination of endogenous variables with
%   numeric coefficients, e.g.
%
%       { 'ygap', [0 1 4 20], [], 'y - yf';
%         'rrf',  [0],        [], ''       }
%
%   ('' or [] for states that are model variables; parse_state_definition
%   gives the syntax). The engine builds the combination's row of the
%   decision rule from the named variables' rows, so lags and changes of a
%   defined state work as for any other. `labels` is the
%   n_key x 1 cellstr of the states compute_state_reliability_diagnostics
%   diagnoses for that list, in the ROW ORDER of its outputs (CorrH, R2_h,
%   nu, TV_state_estimate columns, wrapper cells): per row of the list, the
%   base name if 0 is among its lags, then 'l<m><NAME>' for each m > 0
%   ascending, then 'd<h><NAME>' for each h ascending - the labels the .mod
%   files used to declare by hand (l20rrstar, d20rrstar).
%
%   Intended use in the .mod file, so that every script that indexes results
%   with strcmp(M_.states_of_interest, name) keeps working:
%
%       M_.CellList_States   = { 'ygap', [0 1 4 20], []; 'rrstar', [0 1 4 20], [4 20] };
%       M_.states_of_interest = expand_states_of_interest(M_.CellList_States);
%
%   reliability_draw_wrapper passes M_.CellList_States to the engine when the
%   field exists and M_.states_of_interest otherwise. The engine itself
%   accepts either form as its second argument and returns `labels` as
%   LatentStateDiagnostics.state_labels.
%
%   A legacy cellstr of names (any shape) is accepted and returned as a
%   column, unchanged (each name with lags {0}, no changes).
%
%   sp  struct with the parsed list: .base (J x 1 names), .lags, .changes
%       (J x 1 cells of row vectors), .defs (J x 1 cellstr, the definition of
%       column 4 or '' - ADDED 2026-09-30), .is_spec (false for a legacy cellstr),
%       and per label .lag_m (n_key x 1, the m of each label; 0 for a base
%       or change state) and .base_of (n_key x 1, the base name of each lag
%       label, the label itself otherwise - the eq.-(3.17) denominator
%       state).
%
%   Validation is strict: names must be non-empty char rows, lags finite
%   integers >= 0, changes finite integers >= 1; repeated names are merged.
%
%   See also: compute_state_reliability_diagnostics, augment_ssm,
%             reliability_draw_wrapper.

    [base, lags, changes, is_spec, defs] = parse_state_spec(s);
    J = numel(base);

    labels = {};  lag_m = [];  base_of = {};
    for j = 1:J
        if any(lags{j} == 0)
            labels{end+1} = base{j};   lag_m(end+1) = 0;  base_of{end+1} = base{j};                     %#ok<AGROW>
        end
        for k = lags{j}(lags{j} > 0)
            labels{end+1} = sprintf('l%d%s', k, base{j});  lag_m(end+1) = k;  base_of{end+1} = base{j}; %#ok<AGROW>
        end
        for h = changes{j}
            lab = sprintf('d%d%s', h, base{j});
            labels{end+1} = lab;  lag_m(end+1) = 0;  base_of{end+1} = lab;                                %#ok<AGROW>
        end
    end
    labels = labels(:);
    if isempty(labels)
        error('expand_states_of_interest:empty', 'the list names no state to diagnose.');
    end

    sp = struct('base', {base}, 'lags', {lags}, 'changes', {changes}, 'defs', {defs}, 'is_spec', is_spec, ...
                'lag_m', lag_m(:), 'base_of', {base_of(:)});
end


%% ========================================================================
function [base, lags, changes, is_spec, defs] = parse_state_spec(s)
%PARSE_STATE_SPEC  Read the list in either of its two forms.
%   Legacy: a cellstr (row or column) of variable names -> each with lags
%   {0} (the filtered state) and no changes; is_spec = false.
%   Spec:   a J x 2, J x 3 or J x 4 cell { name, lags[, changes[, definition]] }; is_spec = true.
%   defs: J x 1 cellstr, the definition of column 4 ('' when absent).
    if ~iscell(s) || isempty(s)
        error('expand_states_of_interest:badSpec', ...
            'the list must be a non-empty cell array (names, or a J x 3 cell {name, lags, changes}).');
    end
    % A cell made of character vectors only is the legacy list of names,
    % whatever its shape; anything else must be a J x 2, J x 3 or J x 4 spec.
    allchar = all(cellfun(@ischar, s(:)));
    if allchar
        is_spec = false;
    elseif size(s, 2) >= 2 && size(s, 2) <= 4
        is_spec = true;
    else
        error('expand_states_of_interest:badSpec', ...
            ['the list must be a cellstr of names or a J x 3 (J x 4) cell ' ...
             '{name, smoother lags, change horizons[, definition]}; got a %d x %d cell.'], size(s, 1), size(s, 2));
    end

    if ~is_spec
        names = s(:);
        if ~all(cellfun(@(x) ischar(x) && isrow(x) && ~isempty(x), names))
            error('expand_states_of_interest:badSpec', ...
                'the list (legacy form) must contain non-empty character row vectors.');
        end
        base    = names;
        lags    = repmat({0}, numel(names), 1);
        changes = repmat({zeros(1, 0)}, numel(names), 1);
        defs    = repmat({''}, numel(names), 1);
        if numel(unique(base)) ~= numel(base)
            error('expand_states_of_interest:badSpec', 'the list repeats a name.');
        end
        return
    end

    base = {};  lags = {};  changes = {};  defs = {};
    for r = 1:size(s, 1)
        nm = s{r, 1};
        if ~ischar(nm) || ~isrow(nm) || isempty(nm)
            error('expand_states_of_interest:badSpec', ...
                'list entry {%d,1} must be a non-empty character row vector.', r);
        end
        v = check_list(s{r, 2}, r, 2, 0);
        if size(s, 2) >= 3
            w = check_list(s{r, 3}, r, 3, 1);
        else
            w = zeros(1, 0);
        end
        d = '';
        if size(s, 2) == 4 && ~isempty(s{r, 4})
            d = s{r, 4};
            if ~ischar(d) || ~isrow(d)
                error('expand_states_of_interest:badSpec', 'list entry {%d,4} (definition) must be a character row vector or empty.', r);
            end
            parse_state_definition(d);                             % syntax check; the engine checks the variable names
        end
        % reshape: Octave's unique returns a 0xN empty for an empty row
        k = find(strcmp(base, nm), 1);
        if isempty(k)
            base{end+1, 1}    = nm;                                %#ok<AGROW>
            lags{end+1, 1}    = reshape(unique(v), 1, []);         %#ok<AGROW>
            changes{end+1, 1} = reshape(unique(w), 1, []);         %#ok<AGROW>
            defs{end+1, 1}    = d;                                 %#ok<AGROW>
        else
            lags{k}    = reshape(unique([lags{k}, v]), 1, []);
            changes{k} = reshape(unique([changes{k}, w]), 1, []);
            if ~isempty(d)
                if ~isempty(defs{k}) && ~strcmp(regexprep(defs{k}, '\s+', ''), regexprep(d, '\s+', ''))
                    error('expand_states_of_interest:badSpec', 'the list defines "%s" twice, differently.', nm);
                end
                defs{k} = d;
            end
        end
    end
end

function v = check_list(v, r, col, vmin)
    if isempty(v)
        v = zeros(1, 0);
        return
    end
    if ~isnumeric(v) || ~isvector(v) || ~isreal(v)
        error('expand_states_of_interest:badSpec', ...
            'list entry {%d,%d} must be a real numeric vector or empty.', r, col);
    end
    v = double(reshape(v, 1, []));
    if any(~isfinite(v)) || any(mod(v, 1) ~= 0) || any(v < vmin)
        error('expand_states_of_interest:badSpec', ...
            'list entry {%d,%d} must contain finite integers >= %d.', r, col, vmin);
    end
end
