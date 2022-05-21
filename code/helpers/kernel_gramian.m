function G = kernel_gramian(k, X, Y, w, eta)
%KERNEL_GRAMIAN Summary of this function goes here
%   Detailed explanation goes here

G = zeros(size(X,1), size(Y,1));   % One row per data point
for i = 1:size(X,1)
    for j = 1:size(Y,1)
        G(i,j) = k(X(i,:)', Y(j,:)');
    end
end

if nargin > 3
    G = sqrt(w) .* G .* sqrt(w');
end
if nargin > 4
    G = G + eta * norm(G,'fro') * eye(size(G,1));
end

