function [Q,lambdas] = DMD_SVD(V, tsvd_flag)
%DMD_SVD Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    tsvd_flag = false;
end

[U,S,W] = svd(V(:, 1:end-1),'econ');
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

