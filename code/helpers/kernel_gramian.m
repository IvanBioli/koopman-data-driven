function G = kernel_gramian(k, X, Y, w, eta)
% G = KERNEL_GRAMIAN(K, X, Y, W, ETA)
%   Computes the gramian matrix G using the kernel trick with kernel k. 
% INPUT:
%   - k: kernel
%   - X: initial datapoints
%   - Y: final datapoints
%   - w (optional): weights associated to the datapoints 
%   - eta (optional): regularization parameter
% OUTPUT:
%   - G: Gramian matrix

G = zeros(size(X,1), size(Y,1));   % One row per data point
for i = 1:size(X,1)
    for j = 1:size(Y,1)
        G(i,j) = k(X(i,:)', Y(j,:)');
    end
end

if nargin > 3
    G = sqrt(w) .* G .* sqrt(w');
end

% Adding regularization
if nargin > 4
    G = G + eta * norm(G,'fro') * eye(size(G,1));
end

