function psi = hyperbolic_approximant(N, show)
% PSI = HYPERBOLIC_APPROXIMANT(N, SHOW)
%   Returns the dictionary obtained using hyperbolic cross approximation of
%   order N, for the nonlinear pendulum problem. The Fourier basis is used
%   in x_1 and Hermite functions are used in x_2.

if nargin < 2
    show = false;
end

% Finding the hyperbolic cross approximation indices
indices = [];
for i = 0:N
    for j = 0:N
        if max(i,1) * max(j,1) <= N
            indices = [indices; i, j];
        end
    end
end
indices = [indices; indices(:,1), -indices(:,2)];
indices = unique(indices,'rows');
if show
    plot(indices(:,1), indices(:,2), 'k.')
end
% Evaluation 
psi = @(x) hermite_fun(x(2), indices(:,1)) .* (1/sqrt(2*pi) * exp(1i*x(1)*indices(:,2).'));
end

