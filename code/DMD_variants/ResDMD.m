function [lambdas, KFun] = ResDMD(x0, x1, w, fun_dict, epsilon, psi_0, psi_1)
% [LAMBDAS, KFUN] = RESDMD(X0, X1, W, FUN_DICT, EPSILON, PSI_0, PSI_1)
%   Residual Dynamic Mode Decomposition
% INPUT:
%   - x0: matrix of the initial snapshots. Contains one snapshot per column.
%   - x1: matrix of the final snapshots. Contains one snapshot per column.
%   - w: weights vector
%   - fun_dict: dictionary of observables 
%   - epsilon: tolerance up to which eigenpairs are retained according to
%     the residual
%   - psi_0 (optional): matrix of the evaluations of the dictionay at the
%     initial snapshots
%   - psi_1 (optional): matrix of the evaluations of the dictionay at the
%     final snapshots 
% OUTPUT:
%   - lambdas: eigenvalues approximations
%   - KFun: eigenfunctions approximations

% Computing the matrices psi_0 and psi_1 if not in input
if nargin < 6
    psi_0 = psi_matrix(fun_dict, x0);
    psi_1 = psi_matrix(fun_dict, x1);
end

% Solving the generalized eigenvalue problem
A = psi_0' * (w .* psi_0); A = (A+A')/2;
B = psi_0' * (w .* psi_1);
[V, lambdas] = eig(B,A, 'vector');

% Discarting eigenpairs with residual higher than epsilon
C = psi_1' * (w .*psi_1); C = (C+C')/2;
res_vec = zeros(length(lambdas),1);
for j = 1:length(lambdas)
    lambda = lambdas(j); 
    g = V(:,j);
    res_vec(j) = (g' * (C - lambda * B' - conj(lambda) * B + abs(lambda)^2 * A)* g) / (g' * A * g);
end
lambdas(res_vec > epsilon^2) = [];
V(:, res_vec > epsilon^2) = [];

% Computing the eigenfunctions from the eigenvectors of K
KFun = @(x) fun_dict(x) * V;

end