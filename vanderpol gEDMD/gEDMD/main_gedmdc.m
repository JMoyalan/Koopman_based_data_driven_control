% ========================================================================
% KOOPMAN GENERATOR ESTIMATION USING gEDMD
% ========================================================================

clear; close all; format compact; clc;

addpath('../../SOS Toolbox');
addpath('../../SOS Toolbox/multipoly/multipoly');
%% 1. Load sample points
fprintf(sprintf('[%s] Start loading data samples\n', datetime('now')));
load('../sampling/samples.mat');
fprintf(sprintf('[%s] Finish loading %d data samples\n', datetime('now'), size(X1,2)));

%% Define parameters for OCP
% variable definition
x = sym('x',[2,1],'real');
u = sym('u',[0,1],'real');
dx = sym('dx',[2,1],'real');

% model, dxdt = fx + gx*u
fx = [x(2); (1-x(1)^2)*x(2)- x(1)];
gx = [0; 1];

% local controller design using LQR
A = double(subs(jacobian(fx,x),x,[0;0]));
B = gx;
Q = eye(2);
R = 1;
N = 0;
[K,S,e] = lqr(A,B,Q,R,N);

% constant parameters, alpha and b(x)
alpha = 4;
bx = x'*S*x;

%% 2. Estimate Koopman generator
% 2-1. create dictionary functions
opts.polybasis.enable = 1;
opts.polybasis.interval = [-1 1]; % in case of returnlegendre polynomial, set interval
opts.polybasis.round = 1e-6; % options for numerical precisions of coefficients
opts.polybasis.errnrm = 1e-6; % options for numerical precisions of coefficients
opts.polybasis.type = 'Legendre'; % type of basis {Monomial,Hermite,Legendre}
opts.polybasis.order = 9; % order of polynomial basis
opts.polybasis.noconst = 0;
opts.sinubasis.enable = 0; % use sinusoidal basis functions
opts.sinubasis.n = 10; % number of different frequencies of sinusoidal functions
opts.tpsrbf.enable = 0; % use Thin Plate Spline Radial Basis Functions (RBF)
opts.tpsrbf.nc = 300; % number of center points of RBFs
opts.grbf.enable =0;
opts.grbf.nc = 50;
opts.grbf.epsilon = 0.1;
opts.OCP.bx = bx;
opts.OCP.deg_a = 0 : 1;
opts.OCP.deg_c = 0 : 6;
opts.OCP.deg_s = 0 : 7;
deg_a = opts.OCP.deg_a(end);
deg_c = opts.OCP.deg_c(end);
% -----------------------------------------------------------------------
% Create orthogonal dictionary functions, and coefficient matries of
% a(x), c(x), a(x)*b(x), and b(x)*c(x) in terms of Psi(x)
% -----------------------------------------------------------------------
basis_name = sprintf('%d_%d_%d_%d_%s_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d.mat', ...
    length(x),length(u),length(dx), ...
    opts.polybasis.enable, opts.polybasis.type,opts.polybasis.order, opts.polybasis.noconst, ...
    opts.sinubasis.enable, opts.sinubasis.n, ...
    opts.tpsrbf.enable, opts.tpsrbf.nc, ...
    opts.grbf.enable, opts.grbf.nc, ...
    opts.OCP.deg_a(end), opts.OCP.deg_c(end), opts.OCP.deg_s(end));

if isfile(basis_name)
    load(basis_name,'C_s_Psi','C_a_Psi','C_c_Psi','C_ab_Psi','C_bc_Psi','C_x_Psi','Psi','Psi_poly','dPsi','c');
else
    fprintf(sprintf('[%s] Started generating basis functions\n', datetime('now')));
%     [C_x_Psi, Psi, dPsi, c] = FindBasis(x,dx,opts);
    [C_s_Psi,C_a_Psi,C_c_Psi,C_ab_Psi,C_bc_Psi,C_x_Psi,Psi,Psi_poly,dPsi,c] = FindBasis(x,dx,opts);
    fprintf(sprintf('[%s] Finished generating basis functions\n', datetime('now')));
    save(basis_name,'C_s_Psi','C_a_Psi','C_c_Psi','C_ab_Psi','C_bc_Psi','C_x_Psi','Psi','Psi_poly','dPsi','c');
end

% 2-2. estimate Koopman generator
gedmdopt.type = 'gedmd'; % {gedmd, edmd}
gedmdopt.method = 1; % {1}: mini-batch progress; normal for-loop o.w.
gedmdopt.batch = 100;
gedmdopt.sinubasis_enable = opts.sinubasis.enable;
gedmdopt.sinubasis_n = opts.sinubasis.n;
gedmdopt.tpsrbf_enable = opts.tpsrbf.enable;
gedmdopt.tpsrbf_nc = opts.tpsrbf.nc;
gedmdopt.tpsrbf_ctps = []; % custom RBF center point inputs; put [] to calculate new center points
gedmdopt.grbf_enable = opts.grbf.enable;
gedmdopt.grbf_nc = opts.grbf.nc;
gedmdopt.grbf_cgrbf = [];
gedmdopt.sparse = 0; % enable iterative sparsification
gedmdopt.normal.X = 1; % enable normalization of X
gedmdopt.normal.dX = 1; % enable normalization of dX

U = [];
fprintf(sprintf('[%s] Started gEDMD\n', datetime('now')));
[L1, normPsi, normdPsi, cp] = gEDMDc(X1,dX1,U,Psi,dPsi,gedmdopt);
fprintf(sprintf('[%s] Finished gEDMD\n', datetime('now')));

fprintf(sprintf('[%s] Started gEDMD\n', datetime('now')));
[L2, normPsi, normdPsi, cp] = gEDMDc(X2,dX2,U,Psi,dPsi,gedmdopt);
fprintf(sprintf('[%s] Finished gEDMD\n', datetime('now')));

save('result_gEDMD.mat');