function [Q,lambdas] = DMD_SVD(V)
%DMD_SVD Summary of this function goes here
%   Detailed explanation goes here
    
[U,S,W] = svd(V(:, 1:end-1),'econ'); % ADD TSVD???
S_tilde = U' * V(:, 2:end) * W * diag(1./diag(S)); 

% Computation of the Ritz values and vectors
[Q, lambdas] = eig(S_tilde, 'vector');
Q = U * Q;
% Sorting Ritz values and vectors 
[~, ind]=sort(abs(lambdas),'descend');      
lambdas = lambdas(ind);
Q = Q(:, ind);

end

