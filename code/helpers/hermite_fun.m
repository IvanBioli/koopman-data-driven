function y = hermite_fun(x,k_vec)
% Y = HERMITE_FUN(X, K_VEC) 
%   Computes the hermite functions of orders in k_vec, at the points x

kmax = max(k_vec);
epsi=zeros(kmax+1, length(x));
epsi(1,:)=(pi)^(-1/4)*exp(-0.5*(x.^2));                  % Zero Order (k=0) hermite function
epsi(2,:)=(((pi)^(-1/4)*exp(-0.5*(x.^2))).*x)*sqrt(2);   % First Order (k=1) hermite function
if kmax > 2
    for n=3:kmax+1
        m=n-1;
        epsi(n,:)=(x.*(sqrt(2/m)).*(epsi(n-1,:)))-((sqrt((m-1)/m)).*(epsi(n-2,:)));
    end
end
y=epsi(k_vec+1,:)';
end