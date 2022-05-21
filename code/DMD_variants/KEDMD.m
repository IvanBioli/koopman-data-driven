function [lambdas, varargout] = KEDMD(x0, x1, w, k, K_dominants, varargin)
%KEDMD Summary of this function goes here
%   Detailed explanation goes here
M = size(x0,1); % Number of datapoints
tol = 1e-3;

if length(varargin) > 2     % If matrices and eta in input
    G = varargin{1};
    A = varargin{2};
    eta = varargin{3};
else                        % Computing the matrices with the kernel trick
    eta = 1e-16;
    A = kernel_gramian(k, x1, x0, w);
    G = kernel_gramian(k, x0, x0, w); 
end

% Regularization
A = A + eta * norm(A, 'fro') * eye(size(A,1));
G = G + eta * norm(G, 'fro') * eye(size(G,1));
G = (G + G')/2;

% Computing the matrix K 
[Q,S] = eig(G);
S = sqrt(S);
K = pinv(S, tol) * Q' * A * Q * pinv(S, tol);

% Solving for approximate eigenvectors only
if sum(strcmp('valuesOnly', varargin)) > 0
    if K_dominants == 0
        lambdas = eig(K);
    else
        error('Use K_dominant ~= 0 only for K-ResDMD');
    end
    return
end

% Solving for approximate eigenvectors and eigenfunctions
if sum(strcmp('modes', varargin)) > 0
    if K_dominants == 0
        [Xi, lambdas, U] = eig(K, 'vector');
    else
        error('Use K_dominant ~= 0 only for K-ResDMD');
    end
else
    if K_dominants == 0
        [Xi, lambdas] = eig(K, 'vector');
    else
        [Xi, lambdas] = eigs(K,K_dominants,'largestabs');
        lambdas = diag(lambdas);
    end
end
KFun = cell(M,1);
if K_dominants ~= 0
    Xi = orth(Xi);
end
for i = 1:size(Xi,2)
    Mat = Q * pinv(S) * Xi;
    KFun = @(x) arrayfun(@(t)k(x, t), x0') * Mat;
end
varargout{1} = KFun;

% Solving for approximate eigenmodes
if sum(strcmp('modes', varargin)) > 0
    B = pinv(S) * Q' * spdiags(w, 0, M, M)* x0;
    V = (U' * B).';
    varargout{2} = V;
end

end
