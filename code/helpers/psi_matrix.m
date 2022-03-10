function psi = psi_matrix(fun_dict, x)
%PSI_MATRIX Summary of this function goes here
%   Detailed explanation goes here

% CHANGED FOR LEGENDRE POLYNOMIALS
%psi = fun_dict(x);
psi = zeros(size(x,1), length(fun_dict(x(1,:))));
for i = 1:length(x)
    psi(i,:) = fun_dict(x(i,:));
end

end

