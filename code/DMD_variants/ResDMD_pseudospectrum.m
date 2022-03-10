function [tau, mask] = ResDMD_pseudospectrum(x0, x1, w, fun_dict, epsilon, grid, psi_0, psi_1)
%ResDMD_pseudospectrum Summary of this function goes here
%   Detailed explanation goes here

if nargin < 7
    % Computing the matrices psi_0 and psi_1
    psi_0 = psi_matrix(fun_dict, x0);
    psi_1 = psi_matrix(fun_dict, x1);
end
tic
A = psi_0' * (w .* psi_0); A = (A+A')/2;
B = psi_0' * (w .* psi_1);
C = psi_1' * (w .*psi_1); C = (C+C')/2;

% Solving the generalized eigenvalue problem for each point in the grid
[N1, N2] = size(grid);
tau = zeros(N1, N2);
for i = 1:N1
    toc
    tic
    disp(i)
    for j = 1:N2
        lambda = grid(i,j);
        D = C - lambda * B' - conj(lambda) * B + abs(lambda)^2 * A;
        D = (D+D')/2;
        [~,t] = eigs(D,A,1,'smallestabs','Tolerance',epsilon * 1e-2);
        tau(i,j) = t;
    end
end

tau = sqrt(tau);
mask = tau < epsilon;

end