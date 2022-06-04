function [tau, mask] = ResDMD_pseudospectrum(x0, x1, w, fun_dict, epsilon, grid, psi_0, psi_1)
% [TAU, MASK] = RESDMD_PSEUDOSPECTRUM(X0, X1, W, FUN_DICT, EPSILON, GRID, PSI_0, PSI_1)
%   Residual Dynamic Mode Decomposition for pseudospectrum approximation
% INPUT:
%   - x0: matrix of the initial snapshots. Contains one snapshot per column.
%   - x1: matrix of the final snapshots. Contains one snapshot per column.
%   - w: weights vector
%   - fun_dict: dictionary of observables 
%   - epsilon: value of epsilon for which the epsilon-pseudospectrum is
%     computed
%   - grid: grid of points at which the residual is evaluated
%   - psi_0 (optional): matrix of the evaluations of the dictionay at the
%     initial snapshots
%   - psi_1 (optional): matrix of the evaluations of the dictionay at the
%     final snapshots 
% OUTPUT:
%   - tau: matrix of residuals the grid points
%   - mask: mask matrix indicating which grid points lie inside the
%     approximation of the epsilon-pseudospectrum

% Waitbar
f = waitbar(0,'1','Name','Approximating pseudospectrum...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% Computing the matrices psi_0 and psi_1 if not in input
if nargin < 7
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
    % Update waitbar at each line
    waitbar(i/N1,f,sprintf('Line %d of %d',i, N1))
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