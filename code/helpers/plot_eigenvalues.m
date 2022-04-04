function plot_eigenvalues(lambdas, varargin)
    plot(complex(lambdas), varargin{:})
    xlabel('$\Re(\lambda)$', 'Interpreter','latex');
    ylabel('$\Im(\lambda)$', 'Interpreter','latex');
end