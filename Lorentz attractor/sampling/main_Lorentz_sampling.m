% ========================================================================
% Sampling for Inverted Pendulum
% x_dot = F(x) + G(x)*u
% ========================================================================

clear; close all; clc;

% Inputs
dist = 'uniform'; % {'uniform','normal'}
xs = [-5; -5; -5];
xe = [5; 5; 5];
mean = [0; 0; 0; 0];
stdev = [sqrt(2); sqrt(2); sqrt(2); sqrt(2)];
ns = 1e4;
dt = 1e-3;
tspan = [0 : dt : 2*dt];
reltol = 1e-6;
abstol = 1e-6;
absxlim = [5+0.01; 5+0.01; 5+0.01];
sim2pts = 'yes'; % {'yes','no'}

% parameters
rho = 28;
sigma = 10;
beta = 8/3;
param.rho = rho;
param.sigma = sigma;
param.beta = beta;

%% random initial points
if strcmp(dist,'uniform')
    for i1 = 1 : length(xs)
        y(:,i1) = xs(i1) + (xe(i1) - xs(i1))*rand(ns,1);
    end
 elseif strcmp(dist,'normal')
    for i1 = 1 : length(xs)
        y(:,i1) = mean(i1) + stdev(i1).*randn(ns,1);
    end
end

%% sampling
% -----------------------------------------------------------------------
% Time-series data are collected for the following four cases:
% Case1. simluated from (1) with u = 0
% Case2. simulated from (1) with u = 1
% -----------------------------------------------------------------------
% ODE options
odeopts = odeset('RelTol',reltol,'AbsTol',abstol,'Events',@(t,x) eventfun(t,x,absxlim));

% Simulation for Case 1, generating data matrices, X1 and Y1
u = 0; 
[X1, Y1] = datMat(ns,y,u,tspan,odeopts,sim2pts,param);
fprintf(sprintf('Case 1 simulated. Number of samples = %d\n',size(X1,2)));

% Simluation for Case 2, generating data matrices, X2 and Y2
u = 1; 
[X2,Y2] = datMat(ns,y,u,tspan,odeopts,sim2pts,param);
fprintf(sprintf('Case 2 simulated. Number of samples = %d\n',size(X2,2)));

% Plot figures
figure; scatter(X1(1,:),X1(2,:),'b');
hold on; scatter(Y1(1,:),Y1(2,:),'r'); title('Sample Case 1');
figure; scatter(X2(1,:),X2(2,:),'b');
hold on; scatter(Y2(1,:),Y2(2,:),'r'); title('Sample Case 2');

% Save simluated data
save('sample_Lorenz.mat', ...
    'X1', 'X2', 'Y1', 'Y2', 'dt', 'tspan', 'dist', 'xs', 'xe', 'mean', ...
    'stdev', 'ns', 'dt', 'tspan', 'reltol', 'abstol', 'absxlim', ...
    'sim2pts', 'param');

%% simulation function part
function [X,Y] = datMat(ns,y,u,tspan,odeopts,sim2pts,param)
X = []; Y = [];
nx = size(y,2);
hbar = parfor_progressbar(ns, 'sampling...');
parfor i1 = 1 : ns
    init = zeros(nx,1);
    for i2 = 1 : nx
        init(i2,1) = y(i1,i2);
    end
    [t,x,te,ye,ie] = ode45(@(t,x) model_Lorentz_sampling(t,x,u,param), tspan, init, odeopts);
    if strcmp(sim2pts,'yes')
        if isempty(ie)
            t = t(1:2);
            x = x(1:2,:);
        else
            if length(t) >= 3
                t = t(1:2);
                x = x(1:2,:);
            else
                t = t(1);
                x = x(1,:);
            end
        end
    elseif strcmp(sim2pts,'no')
        if ~isempty(ie)
            t = t(1:end-1);
            x = x(1:end-1,:);
        end
    end
    if length(t) >= 2
        X = [X, x(1:end-1,:)'];
        Y = [Y, x(2:end,:)'];
    end
    hbar.iterate(1);
end
close(hbar);
end