%TEST_AUGMENT_SSM  Verification script for augment_ssm (see MERGE_NOTES.md, section 9).
%   Prints 'All checks passed.' or errors out at the first failed assertion.
%   Runs under Octave and MATLAB; needs nothing but augment_ssm.m on the path.

function test_augment_ssm()

    % ------------------------------------------------------------ system
    % 5 states, 3 shocks, 2 observables; stable A with a mild unit-root-free
    % spectrum so the simulated identities are checked on O(1) numbers.
    try
        rng(20260915);
    catch
        rand('state', 20260915); randn('state', 20260915);
    end
    n = 5; ne = 3; ny = 2;
    A = randn(n); A = 0.8 * A / max(abs(eig(A)));
    B = randn(n, ne);
    C = randn(ny, n);
    D = randn(ny, ne);
    names = {'y'; 'pi'; 'u'; 'g'; 'r'};

    spec = {'pi', [1 4 20], [20]
            'u',  [],       [1 4]
            'g',  [2],      []   };

    [Aa, Ba, Ca, Da, namesA] = augment_ssm(A, B, C, D, names, spec);

    % --------------------------------------------------- labels / sizes
    qL = 20 + 3 + 2;          % maxLag: pi 20, u max(4-1)=3, g 2
    qD = 1 + 2 + 0;
    na = n + qL + qD;
    assert(isequal(size(Aa), [na na]), 'Aa size');
    assert(isequal(size(Ba), [na ne]), 'Ba size');
    assert(isequal(size(Ca), [ny na]), 'Ca size');
    assert(isequal(Da, D), 'Da must equal D');
    assert(iscolumn(namesA) && numel(namesA) == na, 'names shape');
    assert(isequal(namesA(1:n), names), 'original names must come first');

    expLag = [arrayfun(@(k) sprintf('l%dpi', k), 1:20, 'UniformOutput', false), ...
              arrayfun(@(k) sprintf('l%du',  k), 1:3,  'UniformOutput', false), ...
              arrayfun(@(k) sprintf('l%dg',  k), 1:2,  'UniformOutput', false)]';
    expDif = {'d20pi'; 'd1u'; 'd4u'};
    assert(isequal(namesA(n+1:n+qL), expLag), 'lag labels / order');
    assert(isequal(namesA(n+qL+1:end), expDif), 'change labels / order');

    idx = @(lab) find(strcmp(namesA, lab));

    % ------------------------------------------------- block structure
    assert(isequal(Aa(1:n, 1:n), A), 'A block preserved');
    assert(all(all(Aa(1:n, n+1:end) == 0)), 'appended states must not feed back into x_t');
    assert(all(all(Aa(n+qL+1:end, n+qL+1:end) == 0)), 'change rows must not load on change states');
    assert(all(all(Ca(:, n+1:end) == 0)), 'Ca zero on the appended block');
    assert(isequal(Ca(:, 1:n), C), 'C block preserved');
    assert(all(all(Ba(n+1:n+qL, :) == 0)), 'lag rows of Ba must be zero');
    assert(isequal(Ba(idx('d20pi'), :), B(2, :)), 'd20pi row of Ba = B(pi,:)');
    assert(isequal(Ba(idx('d1u'),   :), B(3, :)), 'd1u row of Ba = B(u,:)');
    assert(isequal(Ba(idx('d4u'),   :), B(3, :)), 'd4u row of Ba = B(u,:)');

    % the hand-traced wiring: d20_pi carries -1 on Lag19_pi, 0 on Lag20_pi
    assert(Aa(idx('d20pi'), idx('l19pi')) == -1, 'd20pi must load -1 on l19pi');
    assert(Aa(idx('d20pi'), idx('l20pi')) ==  0, 'd20pi must not touch l20pi');
    assert(idx('d20pi') == 31 && idx('l19pi') == 24, 'offsets differ from the hand trace');
    assert(isequal(Aa(idx('d20pi'), 1:n), A(2, :)), 'd20pi loads e_i''A on x_{t-1}');
    % h = 1 change reads the level out of the original block
    assert(isequal(Aa(idx('d1u'), 1:n), A(3, :) - (1:n == 3)), 'd1u = e_u''A - e_u'' on x_{t-1}');
    assert(all(Aa(idx('d1u'), n+1:end) == 0), 'd1u must not use the lag block');
    assert(Aa(idx('d4u'), idx('l3u')) == -1, 'd4u must load -1 on l3u');
    % lag chain wiring
    assert(Aa(idx('l1pi'), 2) == 1 && Aa(idx('l2pi'), idx('l1pi')) == 1, 'pi lag chain');
    assert(Aa(idx('l1g'),  4) == 1 && Aa(idx('l2g'),  idx('l1g'))  == 1, 'g lag chain');

    % ------------------------------------------- eigenvalues preserved
    ea = sort(abs(eig(Aa))); e0 = sort(abs(eig(A)));
    assert(max(abs(ea(end-n+1:end) - e0)) < 1e-8, 'eig(Aa) must contain eig(A)');
    assert(max(ea(1:na-n)) < 1e-8, 'the remaining eigenvalues of Aa must be zero');

    % ------------------------------------ the two identities, by simulation
    T = 60;
    e = randn(ne, T);
    x  = zeros(n, T);  xa = zeros(na, T);  y = zeros(ny, T);  ya = zeros(ny, T);
    xprev = zeros(n, 1); xaprev = zeros(na, 1);
    for t = 1:T
        x(:, t)  = A  * xprev  + B  * e(:, t);
        y(:, t)  = C  * xprev  + D  * e(:, t);
        xa(:, t) = Aa * xaprev + Ba * e(:, t);
        ya(:, t) = Ca * xaprev + Da * e(:, t);
        xprev = x(:, t); xaprev = xa(:, t);
    end
    assert(max(abs(xa(1:n, :) - x), [], 'all') < 1e-12, 'x block reproduced');
    assert(max(abs(ya - y), [], 'all') < 1e-12, 'y reproduced');
    lagged = @(row, k) [zeros(1, k), x(row, 1:T-k)];       % NAME_{t-k}, zero for t-k <= 0
    for k = 1:20, assert(max(abs(xa(idx(sprintf('l%dpi', k)), :) - lagged(2, k))) < 1e-12, 'Lag_pi identity'); end
    for k = 1:3,  assert(max(abs(xa(idx(sprintf('l%du', k)), :) - lagged(3, k))) < 1e-12, 'Lag_u identity');  end
    for k = 1:2,  assert(max(abs(xa(idx(sprintf('l%dg', k)), :) - lagged(4, k))) < 1e-12, 'Lag_g identity');  end
    assert(max(abs(xa(idx('d20pi'), :) - (x(2, :) - lagged(2, 20)))) < 1e-12, 'd20pi identity');
    assert(max(abs(xa(idx('d1u'),   :) - (x(3, :) - lagged(3, 1))))  < 1e-12, 'd1u identity');
    assert(max(abs(xa(idx('d4u'),   :) - (x(3, :) - lagged(3, 4))))  < 1e-12, 'd4u identity');

    % ------------------------------------------------ 2-column spec
    [A2, B2, C2, D2, N2] = augment_ssm(A, B, C, D, names, spec(:, 1:2));
    [A3, B3, C3, D3, N3] = augment_ssm(A, B, C, D, names, [spec(:, 1:2), {[]; []; []}]);
    assert(isequal(A2, A3) && isequal(B2, B3) && isequal(C2, C3) && isequal(D2, D3) && isequal(N2, N3), ...
        '2-column spec must equal 3-column spec with empty horizons');
    assert(size(A2, 1) == n + 20 + 0 + 2, '2-column spec: only requested lags (pi 20, g 2)');

    % ------------------------------------------------- edge cases
    [A0, B0, C0, D0, N0] = augment_ssm(A, B, C, D, names, {});
    assert(isequal(A0, A) && isequal(B0, B) && isequal(C0, C) && isequal(D0, D) && isequal(N0, names), 'empty spec passthrough');
    [A0, ~, ~, ~, N0] = augment_ssm(A, B, C, D, names, {'y', [], []});
    assert(isequal(A0, A) && isequal(N0, names), 'no-op row passthrough');
    [Am, ~, ~, ~, Nm] = augment_ssm(A, B, C, D, names, {'pi', [1 4 20], [20]; 'u', [], [1 4]; 'g', [2], []; 'pi', [4 2], [20]});
    assert(isequal(Am, Aa) && isequal(Nm, namesA), 'repeated rows must merge to the same system');
    [~, ~, ~, ~, Nr] = augment_ssm(A, B, C, D, names', spec);
    assert(isrow(Nr) && isequal(Nr', namesA), 'row-cell stateNames echoed as a row');

    % ------------------------------------------------- rejections
    mustFail(@() augment_ssm(A, B, C, D, names, {'zzz', [1], []}),      'augment_ssm:UnknownState');
    mustFail(@() augment_ssm(A, B, C, D, names, {'pi', [0], []}),       'augment_ssm:BadListValues');
    mustFail(@() augment_ssm(A, B, C, D, names, {'pi', [], [0]}),       'augment_ssm:BadListValues');
    mustFail(@() augment_ssm(A, B, C, D, names, {'pi', [1.5], []}),     'augment_ssm:BadListValues');
    mustFail(@() augment_ssm(A, B, C, D, names, {'pi', [Inf], []}),     'augment_ssm:BadListValues');
    mustFail(@() augment_ssm(A, B, C, D, names, {'pi', [1+2i], []}),    'augment_ssm:ComplexList');
    mustFail(@() augment_ssm(A, B, C, D, names(1:4), spec),             'augment_ssm:NamesLength');
    mustFail(@() augment_ssm(A, B, C, D, {'y';'pi';'u';'g';'pi'}, spec), 'augment_ssm:DuplicateNames');
    mustFail(@() augment_ssm(A, B, C, D, {'y';'pi';'u';'g';'l1pi'}, spec), 'augment_ssm:NameCollision');
    mustFail(@() augment_ssm(A, B, C, D, names, 'pi'),                  'augment_ssm:BadSpec');
    mustFail(@() augment_ssm(A, B, C, D, names, {'pi'}),                'augment_ssm:BadSpecWidth');
    mustFail(@() augment_ssm(A, B, C, D, names, {'pi', [1], [], []}),   'augment_ssm:BadSpecWidth');
    mustFail(@() augment_ssm(A(1:4, :), B, C, D, names, spec),          'augment_ssm:ANotSquare');

    disp('All checks passed.');
end

function mustFail(f, id)
    ok = false;
    try
        f();
    catch err
        ok = strcmp(err.identifier, id);
        if ~ok
            error('test_augment_ssm:wrongError', 'expected %s, got %s (%s)', id, err.identifier, err.message);
        end
    end
    if ~ok
        error('test_augment_ssm:noError', 'expected error %s was not raised', id);
    end
end
