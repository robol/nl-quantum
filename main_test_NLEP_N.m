clear
close all
clc
%%
A = @(w) [w^2+1, 1, 0;
          2, exp(w)-3, 1;
          0, 1, w-2];

[lam, x, info] = nepN(A, 0.5 + 0.2i);
%%
opts.x0 = randn(3,1) + 1i*randn(3,1);
[lam, x, info] = nepN(A, 0.5 + 0.2i, opts);
%%
opts.nx = 61;
opts.ny = 61;
opts.maxSeeds = 80;
opts.acceptTol = 1e-8;

[lam, X, info] = nepN(A, [-4 4 -4 4], opts);
lam
%%
% for ii = 1:length(lam)
% norm(A(lam(ii))*X(:,ii));
% end
%%
imagesc(info.xgrid(1,:), info.ygrid(:,1), log10(info.scoreGrid));
axis xy
colorbar
hold on
plot(real(lam), imag(lam), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('Re(\lambda)');
ylabel('Im(\lambda)');
title('log10(relative smallest singular value of A(\lambda))');
%%
clear
A = @(w) [w^2+1, 1, 0;
          2, exp(w)-3, 1;
          0, 1, w-2];

dA = @(w) [2*w, 0, 0;
           0, exp(w), 0;
           0, 0, 1];

opts.dA = dA;
[lam_1, x_1, info] = nepN(A, 0.5 + 0.2i, opts);
lam_1
% norm(A(lam_1)*x_1)
dA = @(w) [2*w, 0, 0;
           0, exp(w), 0;
           0, 0, 1];

opts.dA = dA;
[lam_1, x_1, info] = nepN(A, 1.0583 - 1j*0.0000, opts);
lam_1
%%
clear
A = @(w) [w^2+1, 1, 0;
          2, exp(w)-3, 1;
          0, 1, w-2];

opts = struct();
opts.dA_times_x = @(w,x) [2*w, 0, 0;
                          0, exp(w), 0;
                          0, 0, 1] * x;
[lam_2, x_2, info] = nepN(A, 0.5 + 0.2i, opts);
lam_2
% norm(A(lam_2)*x_2)
%%
