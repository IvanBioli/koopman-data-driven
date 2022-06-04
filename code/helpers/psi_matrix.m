function psi = psi_matrix(fun_dict, x)
% PSI = PSI_MATRIX(FUN_DICT, X)
%   Computes the evaluation matrix psi of the dictionary fun_dict at the
%   datapoints x

psi = zeros(size(x,1), length(fun_dict(x(1,:))));
for i = 1:length(x)
    psi(i,:) = fun_dict(x(i,:));
end

end

