function F = model_VDP_control(t,y,x,a,c)
% Control
u = double( subs(c, x, y) / subs(a, x, y) );

% dynamical system
F = zeros(2,1);
F(1) = y(2);
F(2) = (1-y(1)^2)*y(2)- y(1) + u;