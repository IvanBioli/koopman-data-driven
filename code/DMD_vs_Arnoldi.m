%% COMPARISON OF DMD (ARNOLDI BASED AND SVD-BASED) WITH ARNOLDI
%% PARAMETERS DEFINITION
clearvars -except loading
rng(0);
addpath(genpath(pwd))
saving = true;

n = 400;            % Size of the matrix
it = 60;            % Number of iterations
k = 10;              % Eigenpairs to monitor

% Definition of the matrix
beta = 1;
gamma = 0.99;
alpha = 0.8;
d = [beta; gamma; alpha.^(1:n-2)'];
A = spdiags(d, 0, n, n);
[Q,~] = qr(randn(n)); 
A = Q'*A*Q;
Afun = @(x) A*x;
x = rand(n,1); x = x/norm(x);      % Initial vector

% Definition of the DMD input
V = zeros(n, it+1);
V(:,1) = x;
for i = 2:it+1
    V(:, i) = Afun(V(:, i-1));
end

%% ARNOLDI
[U, H] = arnoldi(x, Afun, it, 0.7);

arnoldi_reorth_error = zeros(it, k);
arnoldi_reorth_error((1:it)' < (1:k)) = nan;
for i = 1:it
    eigv = eig(H(1:i, 1:i));
    eigv = sort(eigv,'descend');
    eigv = eigv(1:min(i,k));
    arnoldi_reorth_error(i, 1:min(i,k)) = abs(eigv - d(1:min(i,k)))';
end

%% ARNOLDI-BASED DMD
DMD_error = zeros(it, k);
DMD_error((1:it)' < (1:k)) = nan;
for i = 1:it
    [Q,lambdas] = DMD(V(:,1:i+1));
    lambdas = sort(lambdas, 'descend');
    DMD_error(i, 1:min(i,k)) = abs(lambdas(1:min(i,k)) - d(1:min(i,k)))';
end

%% SVD-BASED DMD
DMD_SVD_error = zeros(it, k);
DMD_SVD_error((1:it)' < (1:k)) = nan;
for i = 1:it
    [Q,lambdas] = DMD_SVD(V(:,1:i+1));
    lambdas = sort(lambdas, 'descend');
    DMD_SVD_error(i, 1:min(i,k)) = abs(lambdas(1:min(i,k)) - d(1:min(i,k)))';
end

%% PLOTS
its = 1:it;
k_plot = 5;

fig = figure();
colors = {'b', 'r', 'k'};
for i = 1:k_plot
    semilogy(its, arnoldi_reorth_error(:,i), 'color', colors{1})
    hold on
    semilogy(its, DMD_error(:,i), 'color', colors{2})
    semilogy(its, DMD_SVD_error(:,i), 'color', colors{3})
end
legend('Arnoldi', 'Arnoldi-based DMD', 'SVD-based DMD', 'Location', 'northeast', 'Interpreter', 'latex', 'Fontsize', 11)
if saving
    saveas(fig, "figures/Arnoldi_vs_DMD", 'epsc')
    saveas(fig, "figures/Arnoldi_vs_DMD", 'png')
end