% ========================================================================
% Sampling for Van der Pol oscillator
% ========================================================================

clear all; close all; clc;

% Inputs
dist = 'uniform'; % {'uniform','normal'}
xs = [-5; -5];
xe = [5; 5];
mean = [0; 0];
stdev = [sqrt(2); sqrt(2)];
ns = 1e4;
dt = 1e-2;
tspan = [0 : dt : 2];
reltol = 1e-6;
abstol = 1e-6;
absxlim = [5.01; 5.01];
sim2pts = 'yes'; % {'yes','no'}
actevent = 'yes';

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
[Xe1,Ye1,X1,dX1] = datMat(ns,y,u,tspan,odeopts,sim2pts,actevent);
fprintf(sprintf('Case 1 simulated. Number of samples = %d\n',size(X1,2)));

% Simluation for Case 2, generating data matrices, X2 and Y2
u = 1; 
[Xe2,Ye2,X2,dX2] = datMat(ns,y,u,tspan,odeopts,sim2pts,actevent);
fprintf(sprintf('Case 2 simulated. Number of samples = %d\n',size(X2,2)));

% Plot figures
figure; scatter(Xe1(1,:),Xe1(2,:),'b');
hold on; scatter(Ye1(1,:),Ye1(2,:),'r'); title('Sample Case 1');
figure; scatter(Xe2(1,:),Xe2(2,:),'b');
hold on; scatter(Ye2(1,:),Ye2(2,:),'r'); title('Sample Case 2');
figure; scatter(dX1(1,:),dX1(2,:),'r'); title('Sample Case 1');
figure; scatter(dX2(1,:),dX2(2,:),'r'); title('Sample Case 2');
% Save simluated data
save('samples.mat', ...
    'Xe1', 'Xe2', 'Ye1', 'Ye2', 'X1','X2', 'dX1', 'dX2', 'dt', 'tspan', 'dist', 'xs', 'xe', 'mean', ...
    'stdev', 'ns', 'dt', 'tspan', 'reltol', 'abstol', 'absxlim', 'sim2pts');

%% simulation function part
function [Xe,Ye,X,dX] = datMat(ns,y,u,tspan,odeopts,sim2pts,actevent)
Xe = []; Ye = []; X = []; dX = []; 
nx = size(y,2);
hbar = parfor_progressbar(ns, 'sampling...');
parfor i1 = 1 : ns
    init = zeros(nx,1);
    for i2 = 1 : nx
        init(i2,1) = y(i1,i2);
    end
    if strcmp(actevent,'yes')
        [t,x,te,ye,ie] = ode45(@(t,x) model_VDP_sampling(t,x,u), tspan, init, odeopts);
    else
        [t,x] = ode45(@(t,x) model_VDP_sampling(t,x,u), tspan, init);
        ie = [];
    end
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
    
    % compute derivatives
    dx = [];
    for i3 = 1 : size(x,1)
        dx(:,i3) = feval(@(x) model_VDP_sampling(0,x,u), x(i3,:)');
    end
    
    if length(t) >= 2
        Xe = [Xe, x(1:end-1,:)'];
        Ye = [Ye, x(2:end,:)'];
        X = [X, x'];
        dX = [dX, dx];
    end
    hbar.iterate(1);
end
close(hbar);
end