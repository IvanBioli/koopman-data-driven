function M = complexgrid(a1, b1, N1, a2, b2, N2)
%COMPLEXGRID Summary of this function goes here
%   Detailed explanation goes here

x = linspace(a1, b1, N1);
y = linspace(a2, b2, N2);
[X,Y] = meshgrid(x,y);
M = X + 1i * Y;

end

