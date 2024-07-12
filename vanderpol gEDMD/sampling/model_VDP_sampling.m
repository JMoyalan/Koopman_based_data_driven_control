function F = model_VDP_sampling(t,x,u)
% x_dot = F + G*u 
%       = [                   x2 ] + [0]*u
%         [ (1 - x1^2) * x2 - x1 ]   [1]

F = zeros(2,1);
F(1) = x(2);
F(2) = (1-x(1)^2)*x(2)- x(1) + u;

% F = zeros(2,1);
% F(1) = x(1);
% F(2) = x(2) - x(1)^2;