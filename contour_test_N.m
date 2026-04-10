%% Initialization
clear;
close all;
clc;

addpath(genpath(pwd));

%% Physical constants
mu0  = 4*pi*1e-7;           % Vacuum permeability [H/m]
eps0 = 8.8541878188e-12;    % Vacuum permittivity [F/m]
c    = 1 / sqrt(mu0*eps0);  % Speed of light in vacuum [m/s]

%% Geometry definition
l       = 0.1;      % Stick length [m]
dmax    = 0.12;     % Maximum x-offset [m]
nSticks = 25;       % Number of sticks

d = linspace(0, dmax, nSticks);

% Base geometry for one vertical stick:
% node 1 = [0; 0; 0]
% node 2 = [0; 0; l]
NN0 = [0 0 0; ...
       0 0 l].';

% Replicate the base stick geometry
NN = repmat(NN0, 1, nSticks);

% Assign x-offsets to each stick
NN(1,1:2:end) = d;
NN(1,2:2:end) = d;

nNodes = size(NN, 2);

% Connectivity matrix G:
% each column identifies the two nodes of one stick
G = [1:2:2*nSticks-1; ...
     2:2:2*nSticks];

% Stick barycenters
bar_stick = 0.5 * (NN(:, G(1,:)) + NN(:, G(2,:)));

% Geometric / numerical parameters
radius    = 1e-3 * ones(1, nSticks);  % Stick radius [m]
np_gauss  = 3;                        % Number of Gauss points
N.thread  = 10;                       % Number of threads
norm_eps  = 0;                        % Unused in this script, kept for compatibility

%% Geometry plot
figure;
plot3( ...
    [NN(1,G(1,:)); NN(1,G(2,:))], ...
    [NN(2,G(1,:)); NN(2,G(2,:))], ...
    [NN(3,G(1,:)); NN(3,G(2,:))], ...
    '.-', 'Color', 'k' ...
);
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
title('Geometry');
view(0, 0);

%% Incidence matrix assembly
val_neg_g = G(1,:);
val_pos_g = G(2,:);

Matrix_G = sparse(1:nSticks, val_neg_g, -ones(nSticks,1), nSticks, nNodes);
Matrix_G = Matrix_G + ...
           sparse(1:nSticks, val_pos_g,  ones(nSticks,1), nSticks, nNodes);

%% Capacitive element construction
% For each node, store the list of connected edges.
ne_cap_max = max(full(sum(abs(Matrix_G), 1)));

Cap_Elem = zeros(1 + ne_cap_max, nNodes);
bar_cap  = zeros(3, nNodes);

for k = 1:nNodes
    idE = find(Matrix_G(:, k));   % Edges connected to node k
    bar_cap(:, k) = mean(bar_stick(:, idE), 2);

    nE = length(idE);
    Cap_Elem(1, k) = nE;
    Cap_Elem(2:nE+1, k) = idE;
end

ne_cap_tot = sum(Cap_Elem(1,:)); %#ok<NASGU> % Total number of capacitive connections
A_ = Matrix_G.';

disp('...done');
disp('------------------------------------');

%% Assembly of static P and L matrices
P0 = P_stick_assemble_matlab_vec( ...
    nSticks, nNodes, NN, G, radius, np_gauss + 20, ...
    ne_cap_max, Cap_Elem, N.thread);

L0 = L_stick_assemble_matlab_vec( ...
    nSticks, nNodes, NN, G, radius, np_gauss + 1, 1, N.thread);

%% Distance matrices for retardation effects
% dist_P = distance_points(NN.');      % Alternative option
dist_P = distance_points(bar_cap.');
dist_L = distance_points(bar_stick.');

%% Frequency-dependent operators
P = @(w) P0 .* exp(-1j * dist_P * w / c);
L = @(w) L0 .* exp(-1j * dist_L * w / c);

%% Reference angular frequency
w0 = sqrt(A_(:,1)' * P0 * A_(:,1)) / sqrt(L0(1,1));
wavel = c / (2*pi*w0); %#ok<NASGU> % Reference wavelength [m]

%% Generalized eigenvalue problem (static case)
scal_ = w0^2;

% Solve:
%   (A_' * P0 * A_) u = lambda * (scal_ * L0) u
%
% Eigenvalues are then rescaled to recover normalized frequencies.
[U, Lambda] = eigs(A_' * P0 * A_, scal_ * L0, 20, 'smallestabs');
Lambda = diag(Lambda);

% Normalized eigenfrequencies
ww_eig = sqrt(Lambda * scal_) ./ w0

%% Residual check for the generalized eigenvalue problem
resEIG = zeros(length(Lambda), 1);

for ii = 1:length(Lambda)
    ui = U(:,ii) / norm(U(:,ii));
    resEIG(ii) = norm((scal_ * L0 * Lambda(ii) - A_' * P0 * A_) * ui);
end

resEIG

%% Plot static eigenvectors
for kk = 1:length(Lambda)
    figure(100 + kk);
    hold on;
    plot(abs(U(:,kk)),  '.-');
    plot(real(U(:,kk)), '.-');
    plot(imag(U(:,kk)), '.-');
    legend('Abs', 'Re', 'Im');
    title(sprintf('Static eigenvector %d', kk));
    xlabel('DoF index');
    ylabel('Amplitude');
    grid on;
    drawnow;
end

%% Nonlinear eigenvalue problem (NEP) setup
opts = [];
opts.nx        = nSticks;
opts.ny        = 25;
opts.maxSeeds  = opts.nx * opts.ny;
opts.acceptTol = 1e-6;
opts.tol       = 1e-12;

scal = w0;
ret  = 1;   % Retardation / propagation enabled (1 = on, 0 = off)

% Frequency-domain nonlinear operator:
% A(w) = -A_'*P(w*scal*ret)*A_ + w^2*scal^2*L(w*scal*ret)
A = @(w) -A_' * P(w * scal * ret) * A_ + ...
          w.^2 * scal^2 * L(w * scal * ret);

% Derivative of A(w) with respect to w
opts.dA = @(w) ...
    -ret * A_' * ((-1j * scal * dist_P / c) .* P(w * scal)) * A_ + ...
     2 * w * scal^2 * L(w * scal * ret) + ...
     ret * scal^2 * w^2 * ((-1j * scal * dist_L / c) .* L(w * scal));

%% Solve nonlinear eigenvalue problem
warning off;
[lambw, X, info] = nepN(A, [0 1.2 -0.2 0.2], opts); %#ok<ASGLU>

% Alternative search regions:
% [lambw, X, info] = nepN(A, [0 0.7 -0.2 0.2], opts);
% [lambw, X, info] = nepN(A, [0 0.5 -0.2 0.2], opts);
% [lambw, X, info] = nepN(A, [0 0.3 -0.2 0.2], opts);

% ww_nlep = lambw * scal (true eigenfrequencies)

% Normalized nonlinear eigenfrequencies
lambw 

%% Residual check for the nonlinear eigenvalue problem
resNEP = zeros(length(lambw), 1);

for ii = 1:length(lambw)
    resNEP(ii) = norm(A(lambw(ii)) * X(:,ii));
end

resNEP

%% Plot nonlinear eigenvectors
for kk = 1:length(lambw)
    figure(1000 + kk);
    hold on;
    plot(abs(X(:,kk)),  '.-');
    plot(real(X(:,kk)), '.-');
    plot(imag(X(:,kk)), '.-');
    legend('Abs', 'Re', 'Im');
    title(sprintf('Nonlinear eigenvector %d', kk));
    xlabel('DoF index');
    ylabel('Amplitude');
    grid on;
    drawnow;
end