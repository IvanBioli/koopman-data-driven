function  eigen_phase_portraits(args, psi, psi_0, psi_1, w, save)
% EIGEN_PHASE_PORTRAITS(ARGS, PSI, PSI_0, PSI_1, W, SAVE)
%  Plots the phase portraits of the approximate eigenfunctions corresponding 
%  to the eigenvalues in e^(1i*args)
% INPUT:
%   - args: argument of the complex eigenvalues
%   - psi: dictionary of observables
%   - psi_0: matrix of the evaluations of the dictionary at the initial
%       datapoints
%   - psi_1: matrix of the evaluations of the dictionary at the final
%       datapoints
%   - w: weights associated to the datapoints
%   - save: true to save the figures

% Mesh for the plots
n = 1000; 
[X1,X2] = meshgrid(linspace(-3, 3, n), linspace(-4,4,n));
z = X1 + 1i * X2;

% Computing the matrices for the generalized eigenvalue problem
A = psi_0' * (w .* psi_0); A = (A+A')/2;
B = psi_0' * (w .* psi_1);
C = psi_1' * (w .*psi_1); C = (C+C')/2;

% Solving the generalized eigenvalue problem for each lambda and obtaining
% the phase portraits
for theta = args
        lambda = exp(1i*theta);
        D = C - lambda * B' - conj(lambda) * B + abs(lambda)^2 * A;
        D = (D+D')/2;
        [g,t] = eigs(D,A,1,'smallestabs');
        assert(t < 0.05)            % Assertion error if the residual is above 0.05
        f = @(x) psi(x) * g;        % Function of which the phase portrait is plotted
        fig = figure();

        vals = zeros(n,n);
        for i = 1:n
            for j = 1:n
            vals(i,j) = f([X1(i,j), X2(i,j)]);
            end
        end
        PhasePlot(z,-vals,'m');
        axis on
        set(gca,'TickLabelInterpreter','latex')
        xlabel('$x_1$', 'Interpreter','latex','FontSize', 20)
        ylabel('$x_2$', 'Interpreter','latex','FontSize', 20)
        title("$\lambda = \exp{(" + num2str(theta) + "i)}$", 'Interpreter','latex','FontSize', 20)
        if save
            set(gca,'LooseInset',get(gca,'TightInset'));
            exportgraphics(gca,"figures/pendulum/phase_portrait_"+num2str(theta* 1e4)+".png")
            exportgraphics(gca,"figures/pendulum/phase_portrait_"+num2str(theta* 1e4)+".png")
        end
end
end

