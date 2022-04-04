%% PARAMETERS DEFINITION
clear
addpath(genpath(pwd))
saving = false;
matrix = 'loss';
include_no_reorth = false;

n = 400;            % Size of the matrix
it = 50;            % Number of iterations
k = 6;              % Eigenpairs to monitor

% Definition of the matrix
if strcmp(matrix, 'random')
    A = randn(n,n);
    A = (A + A') / 2;
    d = eig(A);
    d = sort(d, 'descend');
elseif strcmp(matrix, 'loss')
    rng(0);
    e = sort([5; randn(n-1,1)]); 
    A = spdiags(e, 0, n, n); 
    A(1,end) = 1000;
    d = eig(full(A));
    d = sort(d, 'descend');
else
    beta = 1;
    alpha = 0.95;
    d = [beta; alpha.^(1:n-1)'];
    A = spdiags(d, 0, n, n);
end
Afun = @(x) A*x;
x = rand(n,1);      % Initial vector

% Definition of the DMD input
V = zeros(n, it+1);
V(:,1) = x;
for i = 2:it+1
    V(:, i) = Afun(V(:, i-1));
end

%% ARNOLDI
% ARNOLDI WITHOUT REORTHOGONALIZATION
if include_no_reorth
    [U, H] = arnoldi(x, Afun, it, 0);
    
    arnoldi_error = zeros(it, k);
    arnoldi_error((1:it)' < (1:k)) = nan;
    for i = 1:it
        eigv = eig(H(1:i, 1:i));
        eigv = sort(eigv,'descend');
        eigv = eigv(1:min(i,k));
        arnoldi_error(i, 1:min(i,k)) = abs(eigv - d(1:min(i,k)))';
    end
end
% ARNOLDI WITH REORTHOGONALIZATION
[U, H] = arnoldi(x, Afun, it, 0.7);

arnoldi_reorth_error = zeros(it, k);
arnoldi_reorth_error((1:it)' < (1:k)) = nan;
for i = 1:it
    eigv = eig(H(1:i, 1:i));
    eigv = sort(eigv,'descend');
    eigv = eigv(1:min(i,k));
    arnoldi_reorth_error(i, 1:min(i,k)) = abs(eigv - d(1:min(i,k)))';
end

%% ARNOLDI - DMD
DMD_error = zeros(it, k);
DMD_error((1:it)' < (1:k)) = nan;
for i = 1:it
    [Q,lambdas] = DMD(V(:,1:i+1));
    lambdas = sort(lambdas, 'descend');
    DMD_error(i, 1:min(i,k)) = abs(lambdas(1:min(i,k)) - d(1:min(i,k)))';
end

%% SVD - DMD
DMD_SVD_error = zeros(it, k);
DMD_SVD_error((1:it)' < (1:k)) = nan;
for i = 1:it
    [Q,lambdas] = DMD_SVD(V(:,1:i+1));
    lambdas = sort(lambdas, 'descend');
    DMD_SVD_error(i, 1:min(i,k)) = abs(lambdas(1:min(i,k)) - d(1:min(i,k)))';
end

%% PLOTS
its = 1:it;
k_plot = 3;

fig = figure();
%colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560]};
colors = {'m', 'b', 'r', 'k'};
for i = 1:k_plot
    if include_no_reorth
        semilogy(its, arnoldi_error(:,i), 'color', colors{1})
        hold on
    end
    semilogy(its, arnoldi_reorth_error(:,i), 'color', colors{2})
    hold on
    semilogy(its, DMD_error(:,i), 'color', colors{3})
    semilogy(its, DMD_SVD_error(:,i), 'color', colors{4})
end
legend('Arnoldi NO reorth.', 'Arnoldi WITH reorth.', 'Arnoldi-based DMD', 'SVD-based DMD', 'Location', 'best')

if saving
    saveas(fig, "figures/Arnoldi_vs_DMD", 'epsc')
    saveas(fig, "figures/Arnoldi_vs_DMD", 'png')
end