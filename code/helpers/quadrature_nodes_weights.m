function quadratures = quadrature_nodes_weights(M, flag)

quadratures = containers.Map;

% Gauss-Legendre
[data.x0,data.w] = lgwt(M,-1,0);
quadratures('Gauss-Legendre') = data;

if flag
    % Trapezoidal
    data.x0 = linspace(-1,0,M);
    data.w = ones(M,1) / (M - 1); data.w(1) = data.w(1) / 2; data.w(end) = data.w(1);
    quadratures('Trapezoidal') = data;

    % Riemann sum
    data.x0 = linspace(-1,0,M);
    data.w = ones(M,1) / (M - 1);
    quadratures('Riemann sum') = data;
    
    % Montecarlo
    data.x0 = rand(M,1) - 1;
    data.w = ones(M,1) / M;
    quadratures('Montecarlo') = data;
end

end