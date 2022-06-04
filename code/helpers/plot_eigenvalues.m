function plot_eigenvalues(lambdas, varargin)
% PLOT_EIGENVALUES(LAMBDAS, VARARGIN) 
%   Plots the eigenvalues lambdas on the complex plane

    plot(complex(lambdas), varargin{:})
    xlabel('$\Re(\lambda)$', 'Interpreter','latex');
    ylabel('$\Im(\lambda)$', 'Interpreter','latex');
end