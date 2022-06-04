function M = complexgrid(a1, b1, N1, a2, b2, N2)
% M = COMPLEXGRID(A1, B1, N1, A2, B2, N2)
%   Computes a grid in the complex plane, with 
%       - N1 points between a1 and b1 along the real axis
%       - N2 points between a2 and b2 along the imaginary axis

x = linspace(a1, b1, N1);
y = linspace(a2, b2, N2);
[X,Y] = meshgrid(x,y);
M = X + 1i * Y;

end

