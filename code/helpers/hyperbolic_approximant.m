function psi = hyperbolic_approximant(N, show)
%HYPERBOLIC_APPROXIMANT Summary of this function goes here
%   Detailed explanation goes here
if nargin < 2
    show = false;
end

% Finding the hyperbolic cross approximation indices
indices = [];
for i = 0:N
    for j = 0:N
        if max(i,1) * max(j,1) <= N
        %if (1+i) * (1+j) <= N + 1
            indices = [indices; i, j];
        end
    end
end
indices = [indices; indices(:,1), -indices(:,2)];
indices = unique(indices,'rows');
if show
    plot(indices(:,1), indices(:,2), 'k.')
end

psi = @(x) hermite_fun(x(2), indices(:,1)) .* (1/sqrt(2*pi) * exp(1i*x(1)*indices(:,2).'));
end

