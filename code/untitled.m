%% PARAMETERS
clear
addpath(genpath(pwd))
save = true;

delta_t = 0.5;          % Time step for discretization

N = 100;                 % Hyperbolic cross approximation order
M = 3 * 100 - 1;                 % Data points per dimension

%% TRAPEZOIDAL QUADRATURE NODES DEFINITION
% Plotting the Hermite functions to understand where to truncate the
% trapezoidal quadrature rule in x2 (which is applied in [-L, L])
L = 20;
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
%% ResDMD to control spectral pollution
epsilon = 0.5;
%[lambdas_ResDMD, KFun_ResDMD] = ResDMD(x0, x1, w, psi, epsilon, psi_0, psi_1);
%% ResDMD to estimate pseudospectrum
N1 = 1000; N2 = N1; a = 1.5;
grid = complexgrid(-a, a, N1, -a, a, N2);
[tau, ~] = ResDMD_pseudospectrum(x0, x1, w, psi, epsilon, grid, psi_0, psi_1);
%% PLOTS
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
title("$N_K = " + num2str(length(lambdas_EDMD)) + "$", 'Interpreter','latex', 'FontSize', 14)
if save
    saveas(fig, "figures/pendulum/pendulum_N"+num2str(length(lambdas_EDMD)), 'eps')
    saveas(fig, "figures/pendulum/pendulum_N"+num2str(length(lambdas_EDMD)), 'png')
end
%%
thetas = angle(lambdas_EDMD);
%%
n = 1000; z = zeros(n); f = zeros(n);
[X1,X2] = meshgrid(linspace(-4,4,n));
for i = 1:n
    if mod(i, 100) == 0
        fprintf('i = %d\n',i);
    end
    for j = 1:n
        z(i,j) = X1(i,j) + 1i*X2(i,j);
        tmp = KFun_EDMD([X1(i,j), X2(i,j)]);
        f(i,j) = tmp(5);
    end
end
PhasePlot(z,f);