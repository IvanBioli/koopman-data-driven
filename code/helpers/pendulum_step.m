function x1 = pendulum_step(x0, delta_t)
%PENDULUM_STEP Summary of this function goes here
%   Detailed explanation goes here
xdot = @(t,x) [x(2); -sin(x(1))];
M = size(x0,1); % x0 must have one data point per row
x1 = zeros(M,2);
for i = 1:M
    [~,res] = ode45(xdot, [0 delta_t], x0(i,:)); 
    x1(i,:) = res(end,:);
end
end

