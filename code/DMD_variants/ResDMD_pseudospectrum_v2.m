function [grid,sigs] = ResDMD_pseudospectrum_v2(x0, x1, w, fun_dict, opts, psi_0, psi_1)
%ResDMD_pseudospectrum Summary of this function goes here
%   Detailed explanation goes here

if nargin < 6
    % Computing the matrices psi_0 and psi_1
    psi_0 = psi_matrix(fun_dict, x0);
    psi_1 = psi_matrix(fun_dict, x1);
end
[U,S,V] = svd(sqrt(w) .* psi_0, "econ");
Q = [U, null(U')];
M = Q' * (sqrt(w) .* psi_1 * V * pinv(S));
% Computing the pseudospectrum
[x,y,sigs] = eigtool(M, opts);
% Returning the grid
[X,Y] = meshgrid(x,y);
grid = X + 1i * Y;




end