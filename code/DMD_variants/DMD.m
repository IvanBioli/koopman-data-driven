function [Q,lambdas] = DMD(V)
%DMD Summary of this function goes here
%   Detailed explanation goes here

N = size(V,2);          % Number of snapshots
v_N = V(:, end);        % Last snapshot
V = V(:, 1:end-1);    % First N-1 snapshots

a = pinv(V) * v_N;      % Last column of the companion matrix
S = spdiags(ones(N-1,1), -1, N-1, N-1);
S(:, end) = a;

% Computation of the Ritz values and vectors
[Q, lambdas] = eig(full(S), 'vector');
Q = V * Q;
% Sorting Ritz values and vectors 
[~, ind]=sort(abs(lambdas),'descend');      
lambdas = lambdas(ind);
Q = Q(:, ind);

end

