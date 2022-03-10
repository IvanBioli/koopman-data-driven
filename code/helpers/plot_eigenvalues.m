function plot_eigenvalues(lambdas, marker)
    plot(complex(lambdas), marker)
    xlabel('Re(\lambda)');
    ylabel('Im(\lambda)');
end