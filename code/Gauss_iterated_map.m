%% PARAMETERS DEFINITION
clear
addpath(genpath(pwd))
save = false;
% Definition of the iteration function
gauss_map = @(x, alpha, beta) exp(-alpha * x.^2) + beta;
alpha = 2;
beta = -1-exp(-alpha);
F = @(x) gauss_map(x, alpha, beta);

% Dictionary creation
N = 40;             % Size of the dictionary
fun_dict = @(x) legendreP(0:N-1,2*x+1).*sqrt(2*(0:N-1) + 1);

M = 100;            % Number of quadrature points
epsilon = 0.01;     % Tolerance for ResDMD

%% QUADRATURE RULES
flag = true;        % Flag for computing also for quadrature rules other than Gauss-Legendre
quadratures = quadrature_nodes_weights(M, flag);

%% EDMD and ResDMD
keySet = keys(quadratures);
results = containers.Map;

for k = keySet
    key = string(k);
    x0 = quadratures(key).x0;
    w = quadratures(key).w;
    x1 = F(x0);
    
    % EDMD
    [lambdas, KFun] = EDMD(x0, x1, w, fun_dict);
    % ResDMD to remove spectral pollution
    [lambdas_res, KFun_res] = ResDMD(x0, x1, w, fun_dict, epsilon);

    fig = figure();
    plot_eigenvalues(setdiff(lambdas, lambdas_res), 'r.')
    hold on
    plot_eigenvalues(lambdas_res, 'bx')
    sgtitle(key);
    axis square
    if save
        saveas(fig, "figures/gauss_map/ResDMD_"+key, 'eps')
        saveas(fig, "figures/gauss_map/ResDMD_"+key, 'png')
    end
    
    % Saving the results
    results_struct = struct( ...
        'lambdas', lambdas, ...
        'KFun', KFun, ...
        'lambdas_res', lambdas_res, ...
        'KFun_res', KFun_res);
    results(key) = results_struct;
    clear results_struct
end

%% PSEUDOSPECTRUM APPROXIMATION
N1 = 1000; N2 = 1000;
zoom_shift = 0.6;

lambdas = results('Gauss-Legendre').lambdas;
a1 = min(real(lambdas)) - zoom_shift; b1 = max(real(lambdas)) + zoom_shift;
a2 = min(imag(lambdas)) - zoom_shift; b2 = max(imag(lambdas)) + zoom_shift;

grid = complexgrid(a1, b1, N1, a2, b2, N2);

x0 = quadratures('Gauss-Legendre').x0;
w = quadratures('Gauss-Legendre').w;
x1 = F(x0);

[tau, ~] = ResDMD_pseudospectrum(x0, x1, w, fun_dict, epsilon, grid);

% Plotting the eigenvalues and the pseudospectra
fig = figure();
plot_eigenvalues(setdiff(results('Gauss-Legendre').lambdas, results('Gauss-Legendre').lambdas_res), 'r.')
hold on
plot_eigenvalues(results('Gauss-Legendre').lambdas_res, 'bx')
epsilon_vals = [0.3, 0.1, 0.01, 0.001];
for e = epsilon_vals
    mask = tau < e;
    contour(real(grid), imag(grid), 2 * mask * e, 1, 'k', 'ShowText','on','LineWidth',1.5);
end
axis square

if save
    saveas(fig, 'figures/gauss_map/pseudospectra_contour', 'eps')
    saveas(fig, 'figures/gauss_map/pseudospectra_contour', 'png')
end
%% CONVERGENCE OF THE GALERKIN MATRIX

Kfun_dot_fun = @(x) fun_dict(x)' * fun_dict(F(x));
Galerkin_exact = integral(Kfun_dot_fun, -1, 0,'ArrayValued', true);
%%
keySet = {'Gauss-Legendre', 'Trapezoidal', 'Riemann sum', 'Montecarlo'};
valueSet = {[], [], [], []};
errors = containers.Map(keySet, valueSet);
Ms = round(logspace(1,4,100));

for index = 1:length(Ms)
    M = Ms(index);
    fprintf('Iteration %d, M = %d\n', index, M);
    quadratures = quadrature_nodes_weights(M, true);
    
    for k = keySet
        key = string(k);
        x0 = quadratures(key).x0;
        x1 = F(x0);
        w = quadratures(key).w;
        % Computing the matrices psi_0 and psi_1
        % Computing the matrices psi_0 and psi_1
        psi_0 = psi_matrix(fun_dict, x0);
        psi_1 = psi_matrix(fun_dict, x1);

        A = psi_0' * (w.* psi_1);
        errors(key) = [errors(key), max(abs(A - Galerkin_exact), [],'all')];
    end
end

%%
Ms = round(logspace(1,4,100));
fig = figure();
for k = keySet
    key = string(k);
    loglog(Ms, errors(key), 'LineWidth', 1.5, 'Displayname', key);
    hold on
end
legend('Location','east', 'AutoUpdate','off')

Ms = Ms(Ms > 1000);
legend_index = round(length(Ms)/2);
% O(M^(-1/2))
start = errors('Montecarlo'); start = start(Ms > 1000); start = start(1);
y = 10 * start * 1./ Ms.^(1/2);
loglog(Ms, y, 'k.', 'MarkerSize', 4)
text(Ms(legend_index), y(legend_index)*2 , '$\mathcal{O}(M^{-1/2})$', 'Interpreter','latex')
% O(M^(-1))
start = errors('Riemann sum'); start = start(Ms > 1000); start = start(1);
y = 10 * start * 1./ Ms;
loglog(Ms, y, 'k.', 'MarkerSize', 4)
text(Ms(legend_index), y(legend_index) , '$\mathcal{O}(M^{-1})$', 'Interpreter','latex')
% O(M^(-2))
start = errors('Trapezoidal'); start = start(Ms > 1000); start = start(1);
y = 6000 * start * 1./ Ms.^2;
loglog(Ms, y, 'k.', 'MarkerSize', 4)
text(Ms(legend_index), y(legend_index) , '$\mathcal{O}(M^{-2})$', 'Interpreter','latex')

if save
    saveas(fig, 'figures/gauss_map/Galerkin_convergence', 'eps')
    saveas(fig, 'figures/gauss_map/Galerkin_convergence', 'png')
end