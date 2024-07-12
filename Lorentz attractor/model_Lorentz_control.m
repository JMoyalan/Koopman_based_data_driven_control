function F = model_Lorenz_control(t,y,x,a,c,param)
rho = param.rho;
sigma = param.sigma;
beta = param.beta;

u = subs(c, x, y) / subs(a, x, y);

F = zeros(3,1);
F(1) = sigma*(y(2)-y(1));
F(2) = y(1)*(rho-y(3)) - y(2) + u;
F(3) = y(1)*y(2) - beta*y(3);
