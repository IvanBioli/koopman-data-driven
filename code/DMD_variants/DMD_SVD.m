function [Q,lambdas] = DMD_SVD(V, tsvd_flag)
% [Q, LAMBDAS] = DMD_SVD(V, TSVD_FLAG) 
% SVD-based Dynamic Mode Decomposition for snapshot sequences
% INPUT:
%   - V: snapshot sequence matrix. Contains one snapshot per column.
%   - tsvd_flag: true to use TSVD
% OUTPUT:
%   - Q: Ritz vectors
%   - lambdas: Ritz values

if nargin < 2
    tsvd_flag = false;
end

% Compute the SVD of V_1^(N-1)
[U,S,W] = svd(V(:, 1:end-1),'econ');
% Truncate the SVD if requested
if tsvd_flag
   k = rank(V(:, 1:end-1));
   U = U(:, 1:k);
   S = S(1:k, 1:k);
   W = W(:, 1:k);
end
S_tilde = U' * V(:, 2:end) * W * pinv(S); 

% Computation of the Ritz values and vectors
[Q, lambdas] = eig(S_tilde, 'vector');
Q = U * Q;
% Sorting Ritz values and vectors 
[~, ind]=sort(abs(lambdas),'descend');      
lambdas = lambdas(ind);
Q = Q(:, ind);

end

