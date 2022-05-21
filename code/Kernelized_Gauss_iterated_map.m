%% PARAMETERS DEFINITION
clearvars -except grid sigs
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

epsilon = 0.1;     % Tolerance for ResDMD

%% EDMD and KEDMD for the same Kernel
M = N;
quadratures = quadrature_nodes_weights(M, false);
keySet = {'Gauss-Legendre'};
kernel = @(x, y) fun_dict(x) * fun_dict(y)';

for k = keySet
    key = string(k);
    disp(key)
    x0 = quadratures(key).x0;
    w = quadratures(key).w;
    x1 = F(x0);

    % Matrices to be used by EDMD
    %A = kernel_gramian(kernel, x1, x0, w);
    %G = kernel_gramian(kernel, x0, x0, w); G = (G + G')/2;
    
    % EDMD
    [lambdas, ~] = EDMD(x0, x1, w, fun_dict);
    % ResDMD to remove spectral pollution
    [lambdas_res, ~] = ResDMD(x0, x1, w, fun_dict, epsilon);

    fig = figure();
    plot_eigenvalues(setdiff(lambdas, lambdas_res), 'g.', 'MarkerSize', 10, 'DisplayName', 'EDMD')
    hold on
    plot_eigenvalues(lambdas_res, 'bx', 'MarkerSize', 10, 'DisplayName', 'ResDMD')

    [lambdas, KFun] = KEDMD(x0, x1, w, kernel, 0, G, A, 0);
    plot_eigenvalues(lambdas, 'ro', 'MarkerSize', 10, 'DisplayName', 'KEDMD')
    %title(key, 'FontSize', 20);
    legend('Location','northeast')
    axis square
    axis equal
    if saving
        saveas(fig, "figures/gauss_map/kernelized/KEDMD_"+key, 'epsc')
        saveas(fig, "figures/gauss_map/kernelized/KEDMD_"+key, 'png')
    end
end

%% EDMD and KEDMD: RBF kernel and polynomial kernel
epsilon = 0.01;     % Tolerance for ResDMD
M = 1000;
quadratures = quadrature_nodes_weights(M, true);
keySet = {'Montecarlo'};
K_dominants = 40;

% Kernels
gamma = 1./mean((quadratures('Montecarlo').x0-mean(quadratures('Montecarlo').x0)).^2)^2;
kernels = {@(x, y) exp(-gamma*(x-y)^2)};

for i = 1:length(kernels)
    kernel = kernels{i};
    for k = keySet
        key = string(k);
        x0 = quadratures(key).x0;
        w = quadratures(key).w;
        x1 = F(x0);
        % Matrices to be used by EDMD
        A = kernel_gramian(kernel, x1, x0, w);
        G = kernel_gramian(kernel, x0, x0, w); G = (G + G')/2;
        fprintf('cond(A) = %e,\t rank(A) = %d\n', cond(A), rank(A))
        fprintf('cond(G) = %e,\t rank(G) = %d\n', cond(G), rank(G))
        % Adding regularization
        eta = 1e-16; % Regularization parameter
        fig = figure();
        epsilon_vals = [0.3, 0.1, 0.01, 0.001];
        for e = epsilon_vals
            mask = sigs < e;
            contour(real(grid), imag(grid), 2 * mask * e, 1, 'k', 'ShowText','on','LineWidth',1.5);
            hold on
        end
        axis square
        [lambdas_1step, KFun] = KEDMD(x0, x1, w, kernel, K_dominants,G,A,eta);
        plot_eigenvalues(lambdas_1step, 'g.', 'MarkerSize', 10, 'DisplayName', 'KEDMD-Step 1')
        hold on
        %KResDMD
        quadratures = quadrature_nodes_weights(M, true);
        x0 = quadratures('Gauss-Legendre').x0;
        w = quadratures('Gauss-Legendre').w;
        x1 = F(x0);
        [lambdas_2step, ~] = EDMD(x0, x1, w, KFun);
        [lambdas_res_2step, ~] = ResDMD(x0, x1, w, KFun, epsilon);
        plot_eigenvalues(setdiff(lambdas_2step, lambdas_res_2step), 'r.', 'MarkerSize', 10, 'DisplayName', 'Step 2')
        plot_eigenvalues(lambdas_res_2step, 'bx', 'MarkerSize', 10, 'DisplayName', 'KResDMD')
        axis square
        axis equal
        legend('','', '', '', 'KEDMD-Step 1','Step 2','KResDMD', 'Location','bestoutside')
        if saving
            saveas(fig, "figures/gauss_map/kernelized/KEDMD_2step_rbf_"+key, 'epsc')
            saveas(fig, "figures/gauss_map/kernelized/KEDMD_2step_rbf_"+key, 'png')
        end
    end
end