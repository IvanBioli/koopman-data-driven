%% PARAMETERS
clear
addpath(genpath(pwd))
saving = true;

delta_t = 0.5;          % Time step for discretization

N = 20;                 % Hyperbolic cross approximation order
if N == 20
    M = 100 -1;
    L = 10;
else
    M = 3 * 100 - 1;                 % Data points per dimension
    L = 20;
end

%% TRAPEZOIDAL QUADRATURE NODES DEFINITION
% Plotting the Hermite functions to understand where to truncate the
% trapezoidal quadrature rule in x2 (which is applied in [-L, L])
n_points = 1000;
x = linspace(-L,L,n_points);
deg = 0:N;
y = hermite_fun(x,deg);
figure()
for i = deg
    plot(x, y(:,i+1));
    hold on
end

% Creating the quadrature grid and the corresponding weights
% For x1 we consider the interval [-pi,pi] since the Fourier basis is
% periodic with period 2*pi
x1_grid = linspace(-pi, pi, M+1);
w1 = 2 * pi * ones(M+1, 1) / M; w1(1) = w1(1) / 2; w1(end) = w1(end) / 2;
% For x2 we truncate the quadrature to the interval [-L,L]
x2_grid = linspace(-L, L, M+1);
w2 = 2 * L * ones(M+1, 1) / M; w2(1) = w2(1) / 2; w2(end) = w2(end) / 2;
% Meshing the quadrature nodes and weight
[x0_1,x0_2] = meshgrid(x1_grid,x2_grid); 
x0_1 = x0_1(:); x0_2 = x0_2(:); 
x0 = [x0_1,x0_2];                   % One row per datapoint
w = w1 * w2'; w = w(:);             % Weight vector

%% EIGENVALUES COMPUTATION
x1 = pendulum_step(x0, delta_t);            % One step of the pendulum
psi = hyperbolic_approximant(N);            % Hyperbolic cross approximation dictionary
% Computing the matrices outside the functions to compute them only once
psi_0 = psi_matrix(psi, x0);
psi_1 = psi_matrix(psi, x1);

%% EDMD for eigenpairs computation
[lambdas_EDMD, KFun_EDMD] = EDMD(x0, x1, w, psi, psi_0, psi_1);

%% Comparison of DMD, EDMD and ResDMD
epsilon = 0.25;
[lambdas_ResDMD, KFun_ResDMD] = ResDMD(x0, x1, w, psi, epsilon, psi_0, psi_1);
[Q,lambdas_DMD] = Snapshot_DMD(psi_0, psi_1);

%% Plots
fig = figure();
% Drawing the unit circle
theta = linspace(0, 2*pi, 1000);
plot(exp(1i*theta), '-r', 'DisplayName', 'Unit circle')
hold on

% Plotting the eigenvalues computed using ResDMD
plot_eigenvalues(setdiff(lambdas_EDMD, lambdas_ResDMD), 'm.', 'MarkerSize', 10, 'DisplayName', 'EDMD')
plot_eigenvalues(lambdas_ResDMD, 'bx', 'Markersize', 10, 'DisplayName', "ResDMD")

% Plotting the eigenvalues computed using DMD
plot_eigenvalues(lambdas_DMD, '+', 'color', [0.4660 0.6740 0.1880], 'MarkerSize', 10, 'DisplayName', 'DMD')
axis square
axis equal
legend('Interpreter','latex','Location','bestoutside')
title("$K = " + num2str(length(lambdas_EDMD)) + "$", 'Interpreter','latex', 'FontSize', 20)

if saving
    saveas(fig, "figures/pendulum/pendulum_DMD__N"+num2str(length(lambdas_EDMD)), 'epsc')
    saveas(fig, "figures/pendulum/pendulum_DMD_N"+num2str(length(lambdas_EDMD)), 'png')
end
%% Phase portrait plots
args = [0.4932, 0.9765, 1.4452, 1.8951];
eigen_phase_portraits(args, psi, psi_0, psi_1, w, saving)

%% ResDMD to estimate pseudospectrum
epsilon = 0.25;
N1 = 500; N2 = N1; a = 1.5;
grid = complexgrid(-a, a, N1, -a, a, N2);
[tau, ~] = ResDMD_pseudospectrum(x0, x1, w, psi, epsilon, grid, psi_0, psi_1);
%% Pseudospectrum plot
fig = figure();
% Drawing the unit circle
theta = linspace(0, 2*pi, 1000);
plot(exp(1i*theta), '-r')
hold on
% Plotting the eigenvalues
%plot_eigenvalues(setdiff(lambdas_EDMD, lambdas_ResDMD), 'm.')
%plot_eigenvalues(lambdas_ResDMD, 'bx')
plot_eigenvalues(lambdas_EDMD, 'm.')
% Drawing the estimated epsilon pseudospectrum
epsilon_vals = [0.25];
for e = epsilon_vals
    mask = tau < e;
    contour(real(grid), imag(grid), 2 * mask * e, 1, 'k', 'ShowText','on','LineWidth',1.5);
end
axis square
title("$K = " + num2str(length(lambdas_EDMD)) + "$", 'Interpreter','latex', 'FontSize', 20)
if saving
    saveas(fig, "figures/pendulum/pendulum_N"+num2str(length(lambdas_EDMD)), 'epsc')
    saveas(fig, "figures/pendulum/pendulum_N"+num2str(length(lambdas_EDMD)), 'png')
end

save workspaces\pendulum_N1256.mat tau  