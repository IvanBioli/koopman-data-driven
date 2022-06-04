function [lambdas, varargout] = KEDMD(x0, x1, w, k, K_dominants, varargin)
% [LAMBDAS, VARARGOUT] = KEDMD(X0, X1, W, K, K_DOMINANTS, VARARGIN)
%   Kernelized Extended Dynamic Mode Decomposition
% INPUT:
%   - x0: matrix of the initial snapshots. Contains one snapshot per column.
%   - x1: matrix of the final snapshots. Contains one snapshot per column.
%   - w: weights vector
%   - k: kernel
%   - K_dominants: number of dominants eigenpairs to retain. If equal to 
%     zero, retains all the approximated eigenpairs. 
%   - varargin (optional):
%       - A, G: gramian matrices (must be varargin{1} and varargin{2})
%       - eta: regularization parameter (must be varargin{3})
%       - 'valuesOnly': compute approximate eigenvalues only
%       - 'modes': compute also approximate eigenmodes
% OUTPUT:
%   - lambdas: eigenvalues approximations
%   - varargout (optional, according to what in input in varargin):
%       - varargout{1}: approximate eigenfunctions
%       - varargout{2}: approximate eigenmodes

M = size(x0,1); % Number of datapoints
tol = 1e-3;

if length(varargin) > 2     % If matrices G, A and eta in input
    G = varargin{1};
    A = varargin{2};
    eta = varargin{3};
% If not in input, compute the matrices with the kernel trick and set the
% regularization parameter 
else                        
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

% Solving for approximate eigenvalues only
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
