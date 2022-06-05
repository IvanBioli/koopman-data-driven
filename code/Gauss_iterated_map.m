%% RESULTS FOR THE GAUSS ITERATED MAP EXAMPLE - NON KERNELIZED ALGORITHMS
%% PARAMETERS DEFINITION
clearvars -except loading
rng(0)
addpath(genpath(pwd))
saving = true;
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
flag = true;        % Flag for computing the quadrature nodes also for quadrature rules other than Gauss-Legendre
quadratures = quadrature_nodes_weights(M, flag);

%% DMD
keySet = keys(quadratures);
results_DMD = containers.Map;

for k = keySet
    % Getting the quadrature points and computing the snapshot pairs
    key = string(k);
    x0 = quadratures(key).x0;
    x1 = F(x0);
    % Evaluating the dictionary at the quadrature points
    psi_0 = psi_matrix(fun_dict, x0); psi_0 = psi_0';
    psi_1 = psi_matrix(fun_dict, x1); psi_1 = psi_1';
    
    % Applying DMD
    [Q,lambdas] = Snapshot_DMD(psi_0, psi_1, true);
    % Plotting the results for DMD
    fig = figure('visible', 'off');
    plot_eigenvalues(lambdas, 'go')
    title(key, 'FontSize', 20);
    axis square
    if saving
        saveas(fig, "figures/gauss_map/DMD_"+key, 'epsc')
        saveas(fig, "figures/gauss_map/DMD_"+key, 'png')
    end
    
    % Saving the results
    results_struct = struct( ...
        'lambdas', lambdas);
    results_DMD(key) = results_struct;
    clear results_struct
end

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

    % Plotting the eigenvalues computed with EDMD and ResDMD
    fig = figure();
    plot_eigenvalues(setdiff(lambdas, lambdas_res), 'm.', 'MarkerSize', 10)
    hold on
    plot_eigenvalues(lambdas_res, 'bx', 'MarkerSize', 10)
    
    % Plotting the eigenvalues computed using DMD
    lambdas_DMD = results_DMD(key).lambdas;
    plot_eigenvalues(lambdas_DMD, 'go', 'MarkerSize', 10)
    
    title(key, 'FontSize', 20);
    axis square
    if saving
        saveas(fig, "figures/gauss_map/ResDMD_"+key, 'epsc')
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
figure()
epsilon_vals = [0.3, 0.1, 0.01, 0.001];
N = 1000;
zoom_shift = 0.6;
% Computing the area of interest for the pseudospectrum approximation
lambdas = results('Gauss-Legendre').lambdas;
a1 = min(real(lambdas)) - zoom_shift; b1 = max(real(lambdas)) + zoom_shift;
a2 = min(imag(lambdas)) - zoom_shift; b2 = max(imag(lambdas)) + zoom_shift;

% Getting the quadrature points and computing the snapshot pairs
x0 = quadratures('Gauss-Legendre').x0;
w = quadratures('Gauss-Legendre').w;
x1 = F(x0);

% Pseudospectrum approximation based on the proposition proved in the
% project (Proposition 3.2 of the report)
opts.npts = N;
opts.ax = [a1, b1, a2, b2];
opts.levels = log10(epsilon_vals);
opts.no_graphics = 1;
[grid,sigs] = ResDMD_pseudospectrum_v2(x0, x1, w, fun_dict, opts); 

%% Plotting the eigenvalues and the pseudospectra contours
fig = figure();
plot_eigenvalues(setdiff(results('Gauss-Legendre').lambdas, results('Gauss-Legendre').lambdas_res), 'm.')
hold on
plot_eigenvalues(results('Gauss-Legendre').lambdas_res, 'bx')

for e = epsilon_vals
    mask = sigs < e;
    contour(real(grid), imag(grid), 2 * mask * e, 1, 'k', 'ShowText','on','LineWidth',1.5);
end
axis square

if saving
    saveas(fig, 'figures/gauss_map/pseudospectra_contour', 'epsc')
    saveas(fig, 'figures/gauss_map/pseudospectra_contour', 'png')
end
%% CONVERGENCE OF THE GALERKIN MATRIX
Kfun_dot_fun = @(x) fun_dict(x)' * fun_dict(F(x));

keySet = {'Gauss-Legendre', 'Trapezoidal', 'Riemann sum', 'Montecarlo'};
valueSet = {[], [], [], []};
errors = containers.Map(keySet, valueSet);
Ms = round(logspace(1,4,100));
if loading
    load workspaces\galerkin_convergence.mat
else
    % Computing the exact matrix
    Galerkin_exact = integral(Kfun_dot_fun, -1, 0,'ArrayValued', true);
    % Computing the error for different quadrature rules and increasing
    % number of datapoints
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
            psi_0 = psi_matrix(fun_dict, x0);
            psi_1 = psi_matrix(fun_dict, x1);
            % Computing the matrix A and the error
            A = psi_0' * (w.* psi_1);
            errors(key) = [errors(key), max(abs(A - Galerkin_exact), [],'all')];
        end
    end
end

%% PLOT OF THE CONVERGENCE
Ms = round(logspace(1,4,100));
fig = figure();
for k = keySet
    key = string(k);
    loglog(Ms, errors(key), 'LineWidth', 1.5, 'Displayname', key);
    hold on
end
legend('Location','east', 'AutoUpdate','off', 'Interpreter','latex', FontSize=12)

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
xlabel('M', 'Interpreter','latex');
ylabel('$\max_{1\leq j,k \leq K} |(\Psi_0^* W \Psi_1)_{jk} - \langle\mathcal{K}\psi_k, \psi_j\rangle|$', 'Interpreter','latex');
if saving
    saveas(fig, 'figures/gauss_map/Galerkin_convergence', 'epsc')
    saveas(fig, 'figures/gauss_map/Galerkin_convergence', 'png')
end