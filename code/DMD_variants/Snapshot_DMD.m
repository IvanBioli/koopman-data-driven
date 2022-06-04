function [Q,lambdas] = Snapshot_DMD(X, Y, tsvd_flag)
% [Q, LAMBDAS] = SNAPSHOT_DMD(X, Y, TSVD_FLAG) 
% SVD-based Dynamic Mode Decomposition for snapshot pairs
% INPUT:
%   - X: matrix of the initial snapshots. Contains one snapshot per column.
%   - Y: matrix of the final snapshots. Contains one snapshot per column.
%   - tsvd_flag: true to use TSVD
% OUTPUT:
%   - Q: Ritz vectors
%   - lambdas: Ritz values

if nargin < 3
    tsvd_flag = false;
end

[U,S,W] = svd(X,'econ');
if tsvd_flag
   k = rank(X(:, 1:end-1));
   U = U(:, 1:k);
   S = S(1:k, 1:k);
   W = W(:, 1:k);
end
A_tilde = U' * Y * W * pinv(S); 

% Computation of the Ritz values and vectors
[Q, lambdas] = eig(A_tilde, 'vector');
Q = U * Q;
% Sorting Ritz values and vectors 
[~, ind] = sort(abs(lambdas),'descend');      
lambdas = lambdas(ind);
Q = Q(:, ind);

end

