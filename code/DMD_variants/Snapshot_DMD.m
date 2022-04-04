function [Q,lambdas] = Snapshot_DMD(X, Y)
%Snapshot_DMD Summary of this function goes here
%   Detailed explanation goes here

[U,S,W] = svd(X,'econ'); % ADD TSVD???
A_tilde = U' * Y * W * pinv(S); 

% Computation of the Ritz values and vectors
[Q, lambdas] = eig(A_tilde, 'vector');
Q = U * Q;
% Sorting Ritz values and vectors 
[~, ind] = sort(abs(lambdas),'descend');      
lambdas = lambdas(ind);
Q = Q(:, ind);

end

