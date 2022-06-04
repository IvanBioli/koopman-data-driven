function quadratures = quadrature_nodes_weights(M, flag)
% QUADRATURES = QUADRATURE_NODES_WEIGHTS(M, FLAG) 
%   Computes M quadrature nodes and weights in the interval [-1,0]. 
%   If flag is false, only Gauss-Legendre nodes and weights are computed, 
%   otherwise also according to the trapezoidal quadrature rule, Riemann 
%   sum quadrature rule and Montecarlo quadrature rule.

quadratures = containers.Map;

% Gauss-Legendre
[data.x0,data.w] = lgwt(M,-1,0);
quadratures('Gauss-Legendre') = data;

if flag
    % Trapezoidal
    data.x0 = linspace(-1,0,M)';
    data.w = ones(M,1) / (M - 1); data.w(1) = data.w(1) / 2; data.w(end) = data.w(1);
    quadratures('Trapezoidal') = data;

    % Riemann sum
    data.x0 = linspace(-1,0,M)';
    data.w = ones(M,1) / (M - 1);
    quadratures('Riemann sum') = data;
    
    % Montecarlo
    data.x0 = rand(M,1) - 1;
    data.w = ones(M,1) / M;
    quadratures('Montecarlo') = data;
end

end