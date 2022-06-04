function [U, H] = arnoldi(x, Afun, it, reorth_tol)
% [U, H] = ARNOLDI(X, AFUN, IT, REORTH_TOL)
%   Implementation of the Arnoldi method

% Default reorth_tol = 0.7
if nargin < 4
    reorth_tol = 0.7;
end

% Allocating memory for U and H
n = size(x,1);
it = min(it, n-1); % Number of iterations must be lower than n-1
U = zeros(n, it+1);
H = zeros(it+1, it);

% Iterations
U(:,1) = x / norm(x);
for k = 1:it
    w = Afun(U(:, k));
    H(1:k, k) = U(:, 1:k)' * w;
    u_tilde = w - U(:, 1:k) * H(1:k, k);
    H(k+1, k) = norm(u_tilde);
    
    % Reorthogonalization to solve loss of orthogonality
    if H(k+1, k) < reorth_tol * norm(w)
        h_hat = U(:, 1:k)' * u_tilde;
        u_tilde = u_tilde - U(:, 1:k) * h_hat;
        H(1:k, k) = H(1:k, k) + h_hat;
        H(k+1, k) = norm(u_tilde);
    end

    U(:, k+1) = u_tilde / H(k+1,k);
end

end

