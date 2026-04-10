function [lambda, X, info] = nepN(A, target, opts)
%NEPN Solve the nonlinear eigenvalue problem A(lambda)*x = 0.
%
%   [lambda, X, info] = nepN(A, target)
%   [lambda, X, info] = nepN(A, target, opts)
%
% INPUT
%   A      : function handle, A(lambda) must return an N-by-N matrix
%
%   target :
%       - complex scalar lambda0
%           -> refine one eigenpair near lambda0
%
%       - 1x4 vector [reMin reMax imMin imMax]
%           -> search for eigenvalues in that rectangle of the complex plane
%
%   opts   : optional struct
%
%       % For one eigenpair
%       opts.x0          initial vector (optional)
%       opts.c0          normalization vector c in c'*x = 1 (optional)
%
%       % Derivative info (optional but recommended if available)
%       opts.dA          function handle dA(lambda), returns N-by-N matrix
%       opts.dA_times_x  function handle dA_times_x(lambda,x), returns A'(lambda)*x
%
%       % Newton parameters
%       opts.maxIt       = 30
%       opts.tol         = 1e-10   % target relative residual for single solve
%       opts.acceptTol   = 1e-8    % acceptance threshold in box search
%       opts.fdRelStep   = 1e-6    % finite-difference step for derivative
%       opts.lineSearch  = true
%
%       % Box search parameters
%       opts.nx          = 31
%       opts.ny          = 31
%       opts.maxSeeds    = 40
%       opts.dupTol      = 1e-8
%
% OUTPUT
%   lambda : scalar or column vector of eigenvalues found
%   X      : eigenvector(s), one per column
%   info   : diagnostics
%
% NOTES
%   1) For a single eigenpair, this solves the bordered Newton system
%
%         [ A(lambda)   A'(lambda)*x ] [dx     ] = [ -A(lambda)*x ]
%         [   c'             0       ] [dlambda]   [      0       ]
%
%      with the normalization c'*x = 1.
%
%   2) For box search, this is a heuristic:
%      - sample the box
%      - find points where sigma_min(A(lambda)) is small
%      - use those as seeds for Newton
%
%      It often works well for isolated eigenvalues, but it is NOT a
%      mathematically guaranteed "find all eigenvalues in the box" method.
%
%   3) Repeated / defective eigenvalues may require a more specialized method.
%
% EXAMPLE 1: one eigenpair near a guess
%   A = @(w) [w^2+1, 1,      0;
%             2,     exp(w)-3, 1;
%             0,     1,      w-2];
%   [lam, x, info] = nepN(A, 0.5+0.2i);
%
% EXAMPLE 2: search a box
%   opts.nx = 51; opts.ny = 51; opts.maxSeeds = 80;
%   [lam, X, info] = nepN(A, [-4 4 -4 4], opts);
%
% EXAMPLE 3: provide derivative
%   dA = @(w) [2*w, 0,      0;
%              0,   exp(w), 0;
%              0,   0,      1];
%   opts.dA = dA;
%   [lam, x, info] = nepN(A, 0.5+0.2i, opts);

    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    opts = default_opts(opts);

    if ~isa(A, 'function_handle')
        error('A must be a function handle.');
    end

    % Distinguish single solve vs box search
    if isscalar(target)
        % ---- single eigenpair near a guess ----
        lambda0 = target;
        M0 = A(lambda0);
        check_square_matrix(M0);

        n = size(M0,1);

        if isfield(opts, 'x0') && ~isempty(opts.x0)
            x0 = opts.x0(:);
            if numel(x0) ~= n
                error('opts.x0 has wrong length.');
            end
        else
            x0 = initial_vector_from_svd(M0);
        end

        if isfield(opts, 'c0') && ~isempty(opts.c0)
            c0 = opts.c0(:);
            if numel(c0) ~= n
                error('opts.c0 has wrong length.');
            end
        else
            c0 = x0;
        end

        [x0, c0] = normalize_for_newton(x0, c0);

        [lambda, x, out] = refine_one(A, lambda0, x0, c0, opts);

        X = x;
        info = out;
        return
    end

    % ---- box search in C ----
    if numel(target) ~= 4
        error('target must be either a scalar guess or [reMin reMax imMin imMax].');
    end

    box = target(:).';
    if ~(box(1) < box(2) && box(3) < box(4))
        error('Invalid search box.');
    end

    [seeds, scoreGrid, Xgrid, Ygrid] = make_seeds(A, box, opts);

    lambda = zeros(0,1);
    X = zeros(0,0);
    relres = zeros(0,1);
    rawres = zeros(0,1);
    histories = {};
    usedSeeds = zeros(0,1);

    for k = 1:numel(seeds)
        z0 = seeds(k);

        % initial vector from smallest singular vector at the seed
        M0 = A(z0);
        x0 = initial_vector_from_svd(M0);
        c0 = x0;
        [x0, c0] = normalize_for_newton(x0, c0);

        [z, x, out] = refine_one(A, z0, x0, c0, opts);

        if ~out.ok
            continue
        end

        if ~in_box(z, box, 100*eps)
            continue
        end

        if out.relres > opts.acceptTol
            continue
        end

        j = find(abs(lambda - z) <= opts.dupTol*(1 + abs(z)), 1);

        if isempty(j)
            lambda(end+1,1) = z; %#ok<AGROW>
            if isempty(X)
                X = x;
            else
                X(:,end+1) = x; %#ok<AGROW>
            end
            relres(end+1,1) = out.relres; %#ok<AGROW>
            rawres(end+1,1) = out.rawres; %#ok<AGROW>
            histories{end+1,1} = out.history; %#ok<AGROW>
            usedSeeds(end+1,1) = z0; %#ok<AGROW>
        else
            if out.relres < relres(j)
                lambda(j) = z;
                X(:,j) = x;
                relres(j) = out.relres;
                rawres(j) = out.rawres;
                histories{j} = out.history;
                usedSeeds(j) = z0;
            end
        end
    end

    if ~isempty(lambda)
        [~, ord] = sortrows([real(lambda), imag(lambda)]);
        lambda = lambda(ord);
        X = X(:,ord);
        relres = relres(ord);
        rawres = rawres(ord);
        histories = histories(ord);
        usedSeeds = usedSeeds(ord);
    else
        X = zeros(size(A(mean(box(1:2)) + 1i*mean(box(3:4))),1), 0);
    end

    info.seeds = seeds;
    info.usedSeeds = usedSeeds;
    info.relres = relres;
    info.rawres = rawres;
    info.histories = histories;
    info.scoreGrid = scoreGrid;
    info.xgrid = Xgrid;
    info.ygrid = Ygrid;
    info.ok = ~isempty(lambda);
end


% ========================================================================
% Helpers
% ========================================================================

function opts = default_opts(opts)
    if ~isfield(opts, 'x0'),         opts.x0 = []; end
    if ~isfield(opts, 'c0'),         opts.c0 = []; end
    if ~isfield(opts, 'dA'),         opts.dA = []; end
    if ~isfield(opts, 'dA_times_x'), opts.dA_times_x = []; end

    if ~isfield(opts, 'maxIt'),      opts.maxIt = 30; end
    if ~isfield(opts, 'tol'),        opts.tol = 1e-10; end
    if ~isfield(opts, 'acceptTol'),  opts.acceptTol = 1e-8; end
    if ~isfield(opts, 'fdRelStep'),  opts.fdRelStep = 1e-6; end
    if ~isfield(opts, 'lineSearch'), opts.lineSearch = true; end

    if ~isfield(opts, 'nx'),         opts.nx = 31; end
    if ~isfield(opts, 'ny'),         opts.ny = 31; end
    if ~isfield(opts, 'maxSeeds'),   opts.maxSeeds = 40; end
    if ~isfield(opts, 'dupTol'),     opts.dupTol = 1e-8; end
end


function check_square_matrix(M)
    if ndims(M) ~= 2 || size(M,1) ~= size(M,2)
        error('A(lambda) must return a square matrix.');
    end
end


function x = initial_vector_from_svd(M)
    check_square_matrix(M);
    [~, ~, V] = svd(M, 'econ');
    x = V(:,end);
    if norm(x) == 0
        x = zeros(size(M,1),1);
        x(1) = 1;
    else
        x = x / norm(x);
    end
    x = fix_phase(x);
end


function [x, c] = normalize_for_newton(x, c)
    x = x(:);
    c = c(:);

    if norm(x) == 0
        error('Initial vector x0 must be nonzero.');
    end
    if norm(c) == 0
        error('Normalization vector c0 must be nonzero.');
    end

    x = x / norm(x);
    c = c / norm(c);

    alpha = c' * x;
    if abs(alpha) < 1e-14
        c = x;
        alpha = c' * x;
    end
    x = x / alpha;
end


function dAx = derivative_times_x(A, lambda, x, opts)
    % Prefer dA_times_x, then dA, then finite differences
    if ~isempty(opts.dA_times_x)
        dAx = opts.dA_times_x(lambda, x);
        return
    end

    if ~isempty(opts.dA)
        dAx = opts.dA(lambda) * x;
        return
    end

    h = opts.fdRelStep * (1 + abs(lambda));
    if h == 0
        h = opts.fdRelStep;
    end

    dAx = (A(lambda + h)*x - A(lambda - h)*x) / (2*h);
end


function rr = relative_residual(M, x, r)
    den = norm(M, 'fro') * norm(x) + eps;
    rr = norm(r) / den;
end


function x = final_normalize(x)
    if norm(x) == 0
        return
    end
    x = x / norm(x);
    x = fix_phase(x);
end


function x = fix_phase(x)
    [~, j] = max(abs(x));
    if abs(x(j)) > 0
        x = x * exp(-1i * angle(x(j)));
    end
end


function [lambda, x, out] = refine_one(A, lambda0, x0, c0, opts)
    lambda = lambda0;
    x = x0;
    c = c0;

    hist.lambda = zeros(opts.maxIt+1,1);
    hist.relres = zeros(opts.maxIt+1,1);
    hist.rawres = zeros(opts.maxIt+1,1);

    ok = false;
    itUsed = 0;

    for it = 1:opts.maxIt
        M = A(lambda);
        check_square_matrix(M);

        r = M*x;
        rr = relative_residual(M, x, r);

        hist.lambda(it) = lambda;
        hist.relres(it) = rr;
        hist.rawres(it) = norm(r);
        itUsed = it;

        if rr < opts.tol
            ok = true;
            break
        end

        dAx = derivative_times_x(A, lambda, x, opts);

        % bordered Newton system
        B = [M, dAx; c', 0];
        rhs = [-r; 0];
         
        delta = B \ rhs;
        if any(~isfinite(delta))
            break
        end

        dx = delta(1:end-1);
        dl = delta(end);

        if opts.lineSearch
            [lambdaNew, xNew, rrNew, accepted] = ...
                take_step(A, lambda, x, c, dx, dl, rr);

            if ~accepted
                break
            end

            lambda = lambdaNew;
            x = xNew;

            % optional early stop on tiny update
            if abs(dl) < opts.tol*(1 + abs(lambda)) && rrNew < sqrt(opts.tol)
                % allow one more iteration to clean up residual
            end
        else
            lambda = lambda + dl;
            x = x + dx;

            alpha = c' * x;
            if abs(alpha) < 1e-14
                break
            end
            x = x / alpha;
        end
    end

    M = A(lambda);
    r = M*x;
    rr = relative_residual(M, x, r);

    itUsed = itUsed + 1;
    if itUsed > numel(hist.lambda)
        itUsed = numel(hist.lambda);
    end
    hist.lambda(itUsed) = lambda;
    hist.relres(itUsed) = rr;
    hist.rawres(itUsed) = norm(r);

    x = final_normalize(x);

    out.ok = (rr < opts.acceptTol);
    out.relres = rr;
    out.rawres = norm(A(lambda)*x);
    out.history.lambda = hist.lambda(1:itUsed);
    out.history.relres = hist.relres(1:itUsed);
    out.history.rawres = hist.rawres(1:itUsed);
    out.nIter = itUsed - 1;
end


function [lambdaBest, xBest, rrBest, accepted] = take_step(A, lambda, x, c, dx, dl, rr0)
    alphas = [1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64];

    accepted = false;
    rrBest = inf;
    lambdaBest = lambda;
    xBest = x;

    for a = alphas
        lambdaTry = lambda + a*dl;
        xTry = x + a*dx;

        beta = c' * xTry;
        if abs(beta) < 1e-14
            continue
        end
        xTry = xTry / beta;

        Mtry = A(lambdaTry);
        rTry = Mtry * xTry;
        rrTry = relative_residual(Mtry, xTry, rTry);

        if rrTry < rrBest
            rrBest = rrTry;
            lambdaBest = lambdaTry;
            xBest = xTry;
        end

        if rrTry < rr0
            accepted = true;
            return
        end
    end

    % fallback: accept best trial if finite and not catastrophic
    if isfinite(rrBest) && rrBest < 10*max(rr0, eps)
        accepted = true;
    end
end


function [seeds, scoreGrid, Xgrid, Ygrid] = make_seeds(A, box, opts)
    xr = linspace(box(1), box(2), opts.nx);
    yi = linspace(box(3), box(4), opts.ny);
    [Xgrid, Ygrid] = meshgrid(xr, yi);

    scoreGrid = inf(size(Xgrid));

    for k = 1:numel(Xgrid)
        z = Xgrid(k) + 1i*Ygrid(k);
        try
            M = A(z);
            check_square_matrix(M);
            s = svd(M);
            scoreGrid(k) = s(end) / max(1, s(1));
        catch
            scoreGrid(k) = inf;
        end
    end

    idxLocal = [];
    [ny, nx] = size(scoreGrid);

    for iy = 2:ny-1
        for ix = 2:nx-1
            c = scoreGrid(iy,ix);
            if ~isfinite(c)
                continue
            end
            block = scoreGrid(iy-1:iy+1, ix-1:ix+1);
            block(2,2) = inf;
            if c <= min(block(:))
                idxLocal(end+1,1) = sub2ind(size(scoreGrid), iy, ix); %#ok<AGROW>
            end
        end
    end

    [~, idxGlobal] = sort(scoreGrid(:), 'ascend');
    idxGlobal = idxGlobal(1:min(numel(idxGlobal), opts.maxSeeds));

    idx = unique([idxLocal; idxGlobal], 'stable');
    [~, ord] = sort(scoreGrid(idx), 'ascend');
    idx = idx(ord);

    zcand = Xgrid(idx) + 1i*Ygrid(idx);

    dx = (box(2) - box(1)) / max(opts.nx - 1, 1);
    dy = (box(4) - box(3)) / max(opts.ny - 1, 1);
    sep = 0.75 * max(min(abs(dx), abs(dy)), eps);

    seeds = zeros(0,1);

    for k = 1:numel(zcand)
        z = zcand(k);
        if isempty(seeds) || all(abs(seeds - z) > sep)
            seeds(end+1,1) = z; %#ok<AGROW>
            if numel(seeds) >= opts.maxSeeds
                break
            end
        end
    end
end


function tf = in_box(z, box, tol)
    tf = real(z) >= box(1) - tol && real(z) <= box(2) + tol && ...
         imag(z) >= box(3) - tol && imag(z) <= box(4) + tol;
end
