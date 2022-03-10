function boundary = get_boundary(mask)
%GET_BOUNDARY Summary of this function goes here
%   Detailed explanation goes here

[n,m] = size(mask);
boundary = mask;
boundary(2:n-1, 2:m-1) = boundary(2:n-1, 2:m-1) - ... 
                               boundary(2:n-1, 2:m-1)...
                            .* boundary(3:n, 2:m-1)...
                            .* boundary(1:n-2, 2:m-1)...
                            .* boundary(2:n-1, 3:m)...
                            .* boundary(2:n-1, 1:n-2);
