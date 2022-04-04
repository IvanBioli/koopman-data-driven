%% PARAMETERS DEFINITION
clear
addpath(genpath(pwd))
saving = false;
% Definition of the iteration function
gauss_map = @(x, alpha, beta) exp(-alpha * x.^2) + beta;
alpha = 2;
beta = -1-exp(-alpha);
F = @(x) gauss_map(x, alpha, beta);

% Dictionary creation
N = 40;             % Size of the dictionary
fun_dict = @(x) legendreP(0:N-1,2*x+1).*sqrt(2*(0:N-1) + 1);

M = 40;            % Number of quadrature points
epsilon = 0.01;     % Tolerance for ResDMD

%% QUADRATURE RULES
flag = false;        % Flag for computing also for quadrature rules other than Gauss-Legendre
quadratures = quadrature_nodes_weights(M, flag);

%% EDMD and KEDMD for the same Kernel
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
    for eta = [logspace(-1, -20, 20), 0]
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
    title(key, 'FontSize', 20);
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

%% EDMD and KEDMD: RBF kernel and polynomial kernel

keySet = keys(quadratures);
results = containers.Map;
n = 40;

kernels = {@(x, y) exp(-norm(x-y)^2 / 2), ... 
           @(x, y) (1 + x' * y)^n; 
           };

for i = 1:length(kernels)
    kernel = kernels{i};
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
        title(key, 'FontSize', 20);
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
end