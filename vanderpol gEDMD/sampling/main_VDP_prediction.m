clear all; close all; clc;

% -----------------------------------
% Load result
% -----------------------------------
load('../gEDMD/result_gEDMD.mat','L1','L2','U','C_x_Psi','cp','Psi','norms');

% -----------------------------------
% Model ODE simulation
% -----------------------------------
dt = 1e-2;
ts = 0;
tf = 2;
tspan = [ts:dt:tf];
nt = length(tspan);
init = [4; 1];

u = 0;
[t1,x1] = ode45(@(t,x) model_VDP_sampling(t,x,u), tspan, init);

u = 1;
[t2,x2] = ode45(@(t,x) model_VDP_sampling(t,x,u), tspan, init);

% -----------------------------------
% Prediction (u=0)
% -----------------------------------
norms = [];
cp = [];
[t3, x3] = ode45(@(t,x) model_kmc_cont(t,x,Psi,L1,C_x_Psi,cp), tspan, init);  % approximate

% plot
for i1 = 1 : size(x1,2)
    figure;
    plot(t1,x1(:,i1),'k','LineWidth',1.5); hold on;
    plot(t2,x3(:,i1),'--r','LineWidth',1.2);
end

% -----------------------------------
% Prediction (u=1)
% -----------------------------------
norms = [];
cp = [];
[t4, x4] = ode45(@(t,x) model_kmc_cont(t,x,Psi,L2,C_x_Psi,cp), tspan, init);  % approximate

% plot
for i1 = 1 : size(x1,2)
    figure;
    plot(t1,x1(:,i1),'k','LineWidth',1.5); hold on;
    plot(t2,x2(:,i1),'--r','LineWidth',1.2);
end