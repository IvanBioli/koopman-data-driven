function x1 = pendulum_step(x0, delta_t)
% X1 = PENDULUM_STEP(X0, DELTA_T)
%   Computes one step of the nonlinear pendulum with time step size delta_t
%   and starting datapoints x0

xdot = @(t,x) [x(2); -sin(x(1))];
M = size(x0,1); % x0 must have one data point per row
x1 = zeros(M,2);
for i = 1:M
    [~,res] = ode45(xdot, [0 delta_t], x0(i,:)); 
    x1(i,:) = res(end,:);
end
end

