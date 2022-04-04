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

%% EDMD and KEDMD for the same Kernel
%kernel = @(x, y) psi(x) * psi(y)';
kernel = @(x, y) exp(-norm(x-y)^2 / 2);

% Matrices to be used by EDMD
A = kernel_gramian(kernel, x1, x0, w);
G = kernel_gramian(kernel, x0, x0, w); G = (G + G')/2;

% EDMD
%[lambdas, KFun] = EDMD(x0, x1, w, psi, psi_0, psi_1);
% ResDMD to remove spectral pollution
%[lambdas_res, KFun_res] = ResDMD(x0, x1, w, psi, epsilon, psi_0, psi_1);

fig = figure();
plot_eigenvalues(setdiff(lambdas_EDMD, lambdas_ResDMD), 'g.', 'MarkerSize', 10, 'DisplayName', 'EDMD')
hold on
plot_eigenvalues(lambdas_ResDMD, 'bx', 'MarkerSize', 10, 'DisplayName', 'ResDMD')

% KEDMD
best_err = 1;
best_eta = 0;
%{
for eta = [logspace(-1, -20, 20), 0]
    [lambdas, KFun] = KEDMD(x0, x1, w, kernel, G, A, eta);
    err = norm(sort(lambdas) - sort(lambdas_res)) / norm(lambdas);
    if err > best_err
        err = best_err;
        best_eta = eta;
    end
    fprintf('eta = %e:\t\t %e\n', eta, err)
end
%}
best_eta = 1e-16;
[lambdas, KFun] = KEDMD(x0, x1, w, kernel, G, A, best_eta);
plot_eigenvalues(lambdas, 'r+', 'MarkerSize', 10, 'DisplayName', 'KEDMD')
legend()

if saving
    saveas(fig, "figures/pendulum/kernelized/KEDMD_", 'epsc')
    saveas(fig, "figures/pendulum/kernelized/KEDMD_", 'png')
end
%theta = linspace(0, 2*pi, 1000);
%hold on
%plot(cos(theta), sin(theta), 'r')
axis square
axis equal
%% EDMD and KEDMD RBF kernel
sigma = sqrt(2);
kernel = @(x, y) exp(-norm(x-y)^2 / sigma^2);
keySet = keys(quadratures);
results = containers.Map;

kernel = @(x, y) fun_dict(x) * fun_dict(y)';

for k = keySet
    key = string(k);
    disp(key)
    x0 = quadratures(key).x0;
    w = quadratures(key).w;
    x1 = F(x0);

    % Matrices to be used by EDMD
    A = kernel_gramian(kernel, x1, x0, w);
    G = kernel_gramian(kernel, x0, x0, w); G = (G + G')/2;
    
    % EDMD
    [lambdas, KFun] = EDMD(x0, x1, w, fun_dict);
    % ResDMD to remove spectral pollution
    [lambdas_res, KFun_res] = ResDMD(x0, x1, w, fun_dict, epsilon);

    fig = figure();
    plot_eigenvalues(setdiff(lambdas, lambdas_res), 'g.', 'MarkerSize', 10, 'DisplayName', 'EDMD')
    hold on
    plot_eigenvalues(lambdas_res, 'bx', 'MarkerSize', 10, 'DisplayName', 'ResDMD')

    % KEDMD
    best_err = 1;
    best_eta = 0;
    for eta = [logspace(-4, -20, 17), 0]
        [lambdas, KFun] = KEDMD(x0, x1, w, kernel, G, A, eta);
        err = norm(sort(lambdas) - sort(lambdas_res)) / norm(lambdas);
        if err > best_err
            err = best_err;
            best_eta = eta;
        end
        fprintf('eta = %e:\t\t %e\n', eta, err)
    end

    [lambdas, KFun] = KEDMD(x0, x1, w, kernel, G, A, best_eta);
    plot_eigenvalues(lambdas, 'r+', 'MarkerSize', 10, 'DisplayName', 'KEDMD')
    sgtitle(key);
    legend()
    
    if saving
        saveas(fig, "figures/gauss_map/kernelized/KEDMD_"+key, 'epsc')
        saveas(fig, "figures/gauss_map/kernelized/KEDMD_"+key, 'png')
    end
    %theta = linspace(0, 2*pi, 1000);
    %hold on
    %plot(cos(theta), sin(theta), 'r')
    axis square
    axis equal
end