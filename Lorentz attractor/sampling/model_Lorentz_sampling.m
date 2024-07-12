function F = model_Lorenz(t,x,u,param)
rho = param.rho;
sigma = param.sigma;
beta = param.beta;

F = zeros(3,1);
F(1) = sigma*(x(2)-x(1));
F(2) = x(1)*(rho-x(3)) - x(2) + u;
F(3) = x(1)*x(2) - beta*x(3);