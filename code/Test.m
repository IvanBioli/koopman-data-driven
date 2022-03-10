%% Data generation
n = 1000;
N = 100;

A = rand(n); A = A / norm(A);
evalues = eigs(A,N-1);

V = zeros(n, N); 
V(:,1) = rand(n,1);
for i = 2:N
    V(:,i) = A * V(:, i-1);
end

%% Testing DMD
[Q,lambdas] = DMD(V);
lambdas

abs(lambdas(1:9) - evalues(1:9))

%% Testing DMD
[Q,lambdas] = DMD_SVD(V);

abs(lambdas(1:9) - evalues(1:9))

%% Testing EDMD