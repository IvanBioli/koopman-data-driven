function [grid,sigs] = ResDMD_pseudospectrum_v2(x0, x1, w, fun_dict, opts, psi_0, psi_1)
% [TAU, MASK] = RESDMD_PSEUDOSPECTRUM_V2(X0, X1, W, FUN_DICT, OPTS, PSI_0, PSI_1)
%   Residual Dynamic Mode Decomposition for pseudospectrum approximation
%   Version based on Proposition 3.2 in the report. Relies on eigtool.
% INPUT:
%   - x0: matrix of the initial snapshots. Contains one snapshot per column.
%   - x1: matrix of the final snapshots. Contains one snapshot per column.
%   - w: weights vector
%   - fun_dict: dictionary of observables 
%   - opts: eigtool options
%   - psi_0 (optional): matrix of the evaluations of the dictionay at the
%     initial snapshots
%   - psi_1 (optional): matrix of the evaluations of the dictionay at the
%     final snapshots 
% OUTPUT:
%   - grid: matrix of the grid points
%   - sigs: matrix of residuals the grid points

% Computing the matrices psi_0 and psi_1 if not in input
if nargin < 6
    psi_0 = psi_matrix(fun_dict, x0);
    psi_1 = psi_matrix(fun_dict, x1);
end

% Computing the matrices in Proposition 3.2 of the report
[U,S,V] = svd(sqrt(w) .* psi_0, "econ");
Q = [U, null(U')];
M = Q' * (sqrt(w) .* psi_1 * V * pinv(S));

% Computing the pseudospectrum using eigtool
[x,y,sigs] = eigtool(M, opts);

% Computing the grid
[X,Y] = meshgrid(x,y);
grid = X + 1i * Y;
end