function [lambdas, KFun] = EDMD(x0, x1, w, fun_dict, psi_0, psi_1)
%EDMD Summary of this function goes here
%   Detailed explanation goes here

if nargin < 5
    % Computing the matrices psi_0 and psi_1
    psi_0 = psi_matrix(fun_dict, x0);
    psi_1 = psi_matrix(fun_dict, x1);
end

% Solving the generalized eigenvalue problem
A = psi_0' * (w .* psi_0); A = (A+A')/2;
B = psi_0' * (w .* psi_1);
[V, lambdas] = eig(B,A, 'vector');

% Computing the eigenfunctions from the eigenvectors of K
KFun = @(x) fun_dict(x) * V; 

end