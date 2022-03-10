function [lambdas, KFun] = ResDMD(x0, x1, w, fun_dict, epsilon, psi_0, psi_1)
%ResDMD Summary of this function goes here
%   Detailed explanation goes here

if nargin < 6
    % Computing the matrices psi_0 and psi_1
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