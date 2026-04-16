function [Lambda, X] = beyn_integral(A, n, c, r)
%BEYN_INTEGRAL

tol = 1e-6;
kmax = n;

not_converged = true;

A0 = zeros(n, 0);
A1 = zeros(n, 0);
k = 0;

max_norm = 0;

while not_converged
    V = randn(n, 5);
    k = k + size(V, 2);

    % Approximate A0: integral( A(z) \ V )
    [A0_new, mn] = trap_circle(@(z) A(z) \ V, c, r);
    max_norm = max(max_norm, mn);
    A0 = [A0, A0_new];

    % Approximate A1: integral( z * A(z) \ V )
    A1_new = trap_circle(@(z) z * (A(z) \ V), c, r);
    A1 = [A1, A1_new];

    [U0, S0, V0] = svd(A0, 'econ');
    
    % Rank detection with truncation
    rk = sum(diag(S0) > max_norm * tol);
    
    not_converged = (rk == k) && (k < kmax);
end

U0 = U0(:, 1:rk); S0 = S0(1:rk, 1:rk); V0 = V0(:, 1:rk);
M = U0' * A1 * V0 / S0;

[X, Lambda] = eig(M);
X = U0 * X;

Lambda = diag(Lambda);

end

