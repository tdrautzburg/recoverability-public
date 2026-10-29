function [Aa, Ba, Ca, Da, stateNamesAug] = augment_ssm(A, B, C, D, stateNames, spec)
%AUGMENT_SSM  Augment a linear state space system with lags and h-period changes.
%
%   [Aa, Ba, Ca, Da, stateNamesAug] = ...
%       AUGMENT_SSM(A, B, C, D, stateNames, spec)
%
%   Input system.  Note the timing convention: the observation equation
%   loads on x_{t-1}, not on x_t.
%
%       x_t = A x_{t-1} + B e_t
%       y_t = C x_{t-1} + D e_t
%
%   with A (n x n), B (n x ne), C (ny x n), D (ny x ne), and stateNames a
%   cell of n unique labels, one per row of x_t.
%
%   spec is an m-by-3 cell array
%
%       { 'StateName1', [1 4 20], [20]
%         'StateName2', [],       [1 4]
%         'StateName3', [2],      []    }
%
%   naming a state of x_t in column 1, the requested lags in column 2, and
%   the requested h-period change horizons in column 3.  The row above asks
%   for lags 1, 4 and 20 of StateName1 together with its 20-period change;
%   for the 1- and 4-period changes of StateName2 and no lags of it in
%   particular; and for lag 2 of StateName3 with no change.
%
%   An m-by-2 spec is also accepted and read as lags only, which is the
%   original lag-augmentation interface.  Either list may be empty, a row
%   with both lists empty contributes nothing, and an empty spec returns
%   the inputs unchanged.
%
%   Appended states and their labels (CHANGED 2026-09-15 from 'Lagk_NAME' /
%   'dh_NAME' to the convention of the project's .mod files, l20rrstar,
%   d20rrstar):
%
%       'lkNAME'  =  NAME_{t-k}
%       'dhNAME'  =  NAME_t - NAME_{t-h}
%
%   so on quarterly data d4NAME is the four-quarter (year-on-year) change
%   and d1NAME the one-quarter change. Reading the base back out of a label
%   is unambiguous because a variable name cannot start with a digit.
%
%   Output ordering.  The augmented state vector is
%
%       xa_t = [ x_t ; z_t ; d_t ]
%
%   the original states, then every lag state, then every change state.
%   Within the lag and change blocks, base states appear in order of first
%   appearance in spec; within a base state, lags and horizons ascend.
%   Repeated rows for the same state are merged and duplicate requests
%   collapsed.
%
%   Which lags actually appear.  A lag is realised as a unit shift of the
%   lag before it, so the lag block holds every lag from 1 up to the
%   deepest one needed by that state -- whether or not each was asked for:
%
%       maxLag_i = max( requested lags_i , requested horizons_i - 1 )
%
%   Asking for lag 20 therefore brings lags 1..19 along, and asking for the
%   4-period change brings lags 1..3 along even when column 2 is empty.
%
%   Construction.  Write z^{(i,k)}_t = x_{i,t-k} and
%   d^{(i,h)}_t = x_{i,t} - x_{i,t-h}.  The lag states follow
%
%       z^{(i,1)}_t = x_{i,t-1}           (select row i of x_{t-1})
%       z^{(i,k)}_t = z^{(i,k-1)}_{t-1}   (unit shift, k >= 2)
%
%   A change state must be a linear function of x_{t-1} and e_t to qualify
%   as a state, so substitute the transition equation and re-date the
%   lagged level, x_{i,t-h} = x_{i,(t-1)-(h-1)} = z^{(i,h-1)}_{t-1}:
%
%       d^{(i,h)}_t = (e_i'A) x_{t-1} + (e_i'B) e_t - z^{(i,h-1)}_{t-1}
%
%   with the subtracted term reading x_{i,t-1} out of the original block
%   when h = 1.  Hence
%
%       Aa = [ A   0   0 ]        Ba = [ B  ]      Ca = [ C  0  0 ]
%            [ S   L   0 ]             [ 0  ]      Da = D
%            [ Sd  Ld  0 ]             [ Bd ]
%
%   where S holds the selection rows, L the shifts, Sd the rows e_i'A less
%   e_i' when h = 1, Ld the -1 entries picking Lag(h-1)_NAME, and Bd the
%   rows e_i'B.
%
%   An h-period change needs lag h-1, NOT lag h.  d4NAME leans on
%   l1NAME..l3NAME; d20NAME on l1NAME..l19NAME.
%
%   Properties worth knowing before passing the result to a Kalman filter:
%
%     - Aa is block lower triangular with diagonal blocks A, L and 0, and L
%       is nilpotent, so eig(Aa) = eig(A) together with one zero per
%       appended state.  The augmentation adds no dynamics and preserves
%       stability.
%     - The lag rows of Ba are zero, but the change rows are not: they
%       equal e_i'B, because x_{i,t} responds to the date-t shock.
%     - Every appended state is an exact linear identity, so Ba*Ba' is
%       singular by construction.  That is expected and harmless for a
%       filter that does not invert it, but the initial state covariance
%       must be set consistently, the appended block being a deterministic
%       function of the history of x.
%
%   Reading the columns of Ca.  Ca places zeros on the appended block,
%   leaving y_t = C x_{t-1} + D e_t unchanged in content, but mind the
%   off-by-one when filling those columns in.  Because the measurement
%   equation loads on xa_{t-1}, the labels are exact for the state vector
%   and shifted one period for the columns of Ca:
%
%       column lkNAME loads y_t on NAME_{t-1-k};
%       column dhNAME loads y_t on NAME_{t-1} - NAME_{t-1-h}.
%
%   To load y_t on NAME_{t-j} for j >= 2, use column l(j-1)NAME; j = 1
%   is the corresponding column of the original C block.
%
%   A change state is a linear combination of states already present, so it
%   adds no information to the filter.  It is worth carrying explicitly
%   when you want the change filtered, smoothed, reported, or available as
%   an observable loading.
%
%   Sparsity of A and B is preserved.  The orientation (row or column) of
%   stateNames is preserved, as is its element type (char or string).

narginchk(6, 6);

% ---------------------------------------------------------------- shapes
n  = size(A, 1);
ne = size(B, 2);
ny = size(C, 1);

if ~isnumeric(A) || ndims(A) ~= 2 || size(A, 2) ~= n %#ok<ISMAT>
    error('augment_ssm:ANotSquare', 'A must be a square numeric matrix.');
end
if ~isnumeric(B) || size(B, 1) ~= n
    error('augment_ssm:BRows', ...
        'B must be numeric with %d rows to conform with A.', n);
end
if ~isnumeric(C) || size(C, 2) ~= n
    error('augment_ssm:CCols', ...
        'C must be numeric with %d columns to conform with A.', n);
end
if ~isnumeric(D) || size(D, 1) ~= ny
    error('augment_ssm:DRows', ...
        'D must be numeric with %d rows to conform with C.', ny);
end
if size(D, 2) ~= ne
    error('augment_ssm:DCols', ...
        'D must have %d columns to conform with B.', ne);
end

% ----------------------------------------------------------- state names
if ~iscell(stateNames)
    error('augment_ssm:NamesNotCell', 'stateNames must be a cell array.');
end
if numel(stateNames) ~= n
    error('augment_ssm:NamesLength', ...
        'stateNames has %d entries but A implies %d states.', ...
        numel(stateNames), n);
end

names   = cell(n, 1);   % canonical char copies, used for all matching
nString = 0;
for i = 1:n
    s = stateNames{i};
    if isa(s, 'string') && isscalar(s)
        nString  = nString + 1;
        names{i} = char(s);
    elseif ischar(s) && isrow(s)
        names{i} = s;
    else
        error('augment_ssm:BadName', ...
            ['stateNames{%d} must be a character row vector or a string ' ...
             'scalar.'], i);
    end
    if isempty(names{i})
        error('augment_ssm:EmptyName', 'stateNames{%d} is empty.', i);
    end
end
if nString ~= 0 && nString ~= n
    error('augment_ssm:MixedNameTypes', ...
        'stateNames must be homogeneous: either all char or all string.');
end
useString = (nString == n);

if numel(unique(names)) ~= n
    error('augment_ssm:DuplicateNames', ...
        ['stateNames must be unique; matching a lag or change request ' ...
         'would be ambiguous.']);
end

% stateNames is documented as one label per row of x_t, i.e. a column
% cell, but a row cell is accepted and echoed back in the same shape.
isRowCell = (size(stateNames, 1) == 1) && (n > 1);

% --------------------------------------------------- trivial passthrough
if isempty(spec)
    Aa = A;  Ba = B;  Ca = C;  Da = D;  stateNamesAug = stateNames;
    return
end

% ------------------------------------------------------- parse the spec
if ~iscell(spec) || ndims(spec) ~= 2 %#ok<ISMAT>
    error('augment_ssm:BadSpec', ...
        'spec must be a cell array: { name, lags, horizons; ... }.');
end
nCol = size(spec, 2);
if nCol ~= 2 && nCol ~= 3
    error('augment_ssm:BadSpecWidth', ...
        ['spec must have 2 columns { name, lags } or 3 columns ' ...
         '{ name, lags, horizons }; it has %d.'], nCol);
end

m       = size(spec, 1);
baseIdx = zeros(0, 1);   % index into x_t of each base state, spec order
lagReq  = cell(0, 1);    % requested lags for that base state
difReq  = cell(0, 1);    % requested change horizons for that base state

for r = 1:m
    s = spec{r, 1};
    if isa(s, 'string') && isscalar(s)
        s = char(s);
    end
    if ~ischar(s) || ~isrow(s) || isempty(s)
        error('augment_ssm:BadSpecName', ...
            ['spec{%d,1} must be a non-empty character row vector or ' ...
             'string scalar.'], r);
    end

    hit = find(strcmp(names, s));
    if isempty(hit)
        error('augment_ssm:UnknownState', ...
            'spec{%d,1} = ''%s'' is not the name of a state in x_t.', r, s);
    end
    % names are unique, so hit is a scalar.

    v = localCheckList(spec{r, 2}, r, 2, 'lags');
    if nCol == 3
        w = localCheckList(spec{r, 3}, r, 3, 'change horizons');
    else
        w = zeros(1, 0);
    end

    k = find(baseIdx == hit);
    if isempty(k)
        baseIdx(end+1, 1) = hit;   %#ok<AGROW>
        lagReq{end+1, 1}  = v;     %#ok<AGROW>
        difReq{end+1, 1}  = w;     %#ok<AGROW>
    else
        lagReq{k} = [lagReq{k}, v];
        difReq{k} = [difReq{k}, w];
    end
end

J = numel(baseIdx);

% Deepest lag each base state needs: the deepest requested outright, and
% one less than the deepest change horizon, since dh leans on Lag(h-1).
maxLag = zeros(J, 1);
for j = 1:J
    % unique sorts ascending and collapses repeats; the reshape keeps an
    % empty list a 1x0 row (Octave returns a 0xN empty for an empty input,
    % which breaks the concatenation below).
    lagReq{j} = reshape(unique(lagReq{j}), 1, []);
    difReq{j} = reshape(unique(difReq{j}), 1, []);
    maxLag(j) = max([lagReq{j}, difReq{j} - 1, 0]);
end

off = [0; cumsum(maxLag(1:end-1))];          % lag block offsets
qL  = sum(maxLag);

nDif = zeros(J, 1);
for j = 1:J
    nDif(j) = numel(difReq{j});
end
offD = [0; cumsum(nDif(1:end-1))];           % change block offsets
qD   = sum(nDif);

if qL == 0 && qD == 0
    Aa = A;  Ba = B;  Ca = C;  Da = D;  stateNamesAug = stateNames;
    return
end

% ---------------------------------------------------------------- labels
lagNames  = cell(qL, 1);
diffNames = cell(qD, 1);
for j = 1:J
    base = names{baseIdx(j)};
    for k = 1:maxLag(j)
        lagNames{off(j) + k} = sprintf('l%d%s', k, base);
    end
    hs = difReq{j};
    for p = 1:numel(hs)
        diffNames{offD(j) + p} = sprintf('d%d%s', hs(p), base);
    end
end

allNames = [names; lagNames; diffNames];
if numel(unique(allNames)) ~= numel(allNames)
    [u, ~, ic] = unique(allNames);
    cnt = accumarray(ic(:), 1);
    dup = u(cnt > 1);
    error('augment_ssm:NameCollision', ...
        ['Generated labels collide with existing state names: %s. Rename ' ...
         'the existing states or drop the request.'], ...
        strjoin(reshape(dup, 1, []), ', '));
end

% -------------------------------------------------------------- matrices
% S  (qL x n)  selects the lag-1 states out of x_{t-1}
% L  (qL x qL) shifts lag k-1 into lag k within each base-state block
% Sd (qD x n)  holds e_i'A, less e_i' when h = 1
% Ld (qD x qL) holds the -1 picking l(h-1)NAME when h >= 2
% Bd (qD x ne) holds e_i'B
if issparse(A)
    S  = sparse(qL, n);   L  = sparse(qL, qL);
    Sd = sparse(qD, n);   Ld = sparse(qD, qL);
else
    S  = zeros(qL, n);    L  = zeros(qL, qL);
    Sd = zeros(qD, n);    Ld = zeros(qD, qL);
end
if issparse(B)
    Bd = sparse(qD, ne);
else
    Bd = zeros(qD, ne);
end

for j = 1:J
    i = baseIdx(j);
    o = off(j);

    if maxLag(j) >= 1
        S(o + 1, i) = 1;
        for k = 2:maxLag(j)
            L(o + k, o + k - 1) = 1;
        end
    end

    hs = difReq{j};
    for p = 1:numel(hs)
        h = hs(p);
        rw = offD(j) + p;

        % d_t = x_{i,t} - x_{i,t-h}
        %     = (e_i'A) x_{t-1} + (e_i'B) e_t - x_{i,(t-1)-(h-1)}
        Sd(rw, :) = A(i, :);
        Bd(rw, :) = B(i, :);

        if h == 1
            Sd(rw, i) = Sd(rw, i) - 1;            % x_{i,t-1}, original block
        else
            Ld(rw, o + h - 1) = Ld(rw, o + h - 1) - 1;   % l(h-1)NAME
        end
    end
end

Aa = [A,  zeros(n,  qL), zeros(n,  qD)
      S,  L,             zeros(qL, qD)
      Sd, Ld,            zeros(qD, qD)];
Ba = [B; zeros(qL, ne); Bd];
Ca = [C, zeros(ny, qL + qD)];
Da = D;

% ---------------------------------------------------------------- names
appended = [lagNames; diffNames];
if useString
    appended = cellfun(@string, appended, 'UniformOutput', false);
end
stateNamesAug = [reshape(stateNames, [], 1); appended];
if isRowCell
    stateNamesAug = reshape(stateNamesAug, 1, []);
end

end

% ------------------------------------------------------------------------
function v = localCheckList(v, r, col, what)
%LOCALCHECKLIST  Validate one list of lags or horizons, returned as a row.

if isempty(v)
    v = zeros(1, 0);
    return
end
if ~isnumeric(v) || ~isvector(v)
    error('augment_ssm:BadList', ...
        'spec{%d,%d} must be a numeric vector of %s, or empty.', r, col, what);
end
if ~isreal(v)
    error('augment_ssm:ComplexList', ...
        'spec{%d,%d} must be real.', r, col);
end
v = double(reshape(v, 1, []));
if any(~isfinite(v)) || any(mod(v, 1) ~= 0) || any(v < 1)
    error('augment_ssm:BadListValues', ...
        ['spec{%d,%d} must contain finite integers >= 1; lag 0 is the ' ...
         'contemporaneous state and a 0-period change is identically ' ...
         'zero.'], r, col);
end

end
