
clear; close all; format compact; clc;

addpath('../../SOS Toolbox');
addpath('../../SOS Toolbox/multipoly/multipoly');
addpath('../../SOS Toolbox/sosopt/sosopt');
addpath('../../SOS Toolbox/sedumi-master/sedumi-master');

load('../sampling/samples.mat');
fprintf(sprintf('Number of samples (Case 1) = %d\n',size(X1,2)));
fprintf(sprintf('Number of samples (Case 2) = %d\n',size(X2,2)));

load('../gEDMD/result_gEDMD.mat');
Psi = Psi_poly;
gamma = 0;
x = s2p(x);
x1 = x(1);
x2 = x(2);
bx = s2p(bx);

ax = C_a_Psi'*Psi;
cx = C_c_Psi'*Psi;
sx = C_s_Psi'*Psi;

L2 = L2 - L1;

%% APPROXIMATION OF THE TERMS IN STABILITY EQUATIONS
% -----------------------------------------------------------------------
% Then, the following terms are approximated:
% term1 = \Delta \cdot (Fa) = (I - K1)*a(x)/dt                -- (2)
% term2 = \Delta \cdot (Gc) = (K1 - K2)*c(x)/dt               -- (3)
% term3 = \Delta \cdot (bFa) = (I - K1)*a(x)*b(x)/dt          -- (4)
% term4 = \Delta \cdot (nGc) = (K3 - K2)*c(x)*b(x)/dt         -- (5)
% -----------------------------------------------------------------------

ndic = length(Psi); % number of dictionary functions

% Approximation of the terms
if opts.OCP.deg_a(end) == 0
    C_a_Psi = subs(C_a_Psi,decvars(1),1);
    C_ab_Psi = subs(C_ab_Psi,decvars(1),1);
end

nx = length(x); % number of states
% L1 = (K1'-eye(ndic))/dt; % Koopman generator for Case 1
% L2 = (K2'-K1')/dt; % Koopman generator for Case 2
Div_F = 0; Div_G = 0;
F = C_x_Psi'*L1*Psi;
G = C_x_Psi'*L2*Psi;

for i1 = 1 : nx
    Div_Fi = diff(F(i1), x(i1));
    Div_F = Div_F + Div_Fi;
    Div_Gi = diff(G(i1), x(i1));
    Div_G = Div_G + Div_Gi;
end

Div_F.coefficient(find(abs(Div_F.coefficient) <= 1e-9)) = 0;
Div_G.coefficient(find(abs(Div_G.coefficient) <= 1e-9)) = 0;

term(1) = C_a_Psi'*L1*Psi + C_a_Psi'*Psi*Div_F;
term(2) = C_c_Psi'*L2*Psi + C_c_Psi'*Psi*Div_G;
term(3) = C_ab_Psi'*L1*Psi + C_ab_Psi'*Psi*Div_F;
term(4) = C_bc_Psi'*L2*Psi + C_bc_Psi'*Psi*Div_G;

% Formulate stability equation
sos_const1 = (1+alpha)*bx*sum(term(1:2)) - alpha*sum(term(3:4));
sos_const1.coefficient(find(abs(sos_const1.coefficient)<=1e-3)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if deg_a == 0 % if a is constant
    C_a_amon = mpvar('a',1,1);
else
    deg_amon = deg_a; % minimum degree : maximum degree of a(x)
    amon = monomials(x, deg_amon); % M_a(x): monomials of a(x)
    namon = length(amon); % number of monomials of a(x)
    C_a_amon = mpvar('a',namon,1); % C_a_amon
end

deg_cmon = deg_c; % minimum degree : maximum degree of c(x)
cmon = monomials(x, deg_cmon); % M_c(x): monomials for c(x)
ncmon = length(cmon); % number of elements in M_c(x)
C_c_cmon = mpvar('c',ncmon,1); % C_c_cmon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -----------------------------------------------------------------------
% Solve L1 OCP
% -----------------------------------------------------------------------
% Find maximum-degree monomial vector among a(x), c(x), and s(x)
% z_deg = max([floor(ax.maxdeg/2),floor(cx.maxdeg/2),floor(sx.maxdeg/2)])+1;
zx_deg = floor(Psi.maxdeg/2) + 1;
zx = monomials(x,0:zx_deg);

Mx_info = sprintf('Mx_%d_%d.mat',zx_deg,length(x));

if isfile(Mx_info)
    load(Mx_info,'Mx','Z','Xi_Z_Xi');
    disp('Mx and Xi_Z_Xi loaded from previous save data.');
else
    disp('Start calculating Mx and Xi_Z_Xi.');
    % Big PSD matrix, M(x) = Tx'*Z*Tx
    Tx = kron(zx,eye(2));
    tmp = monomials(x,0:zx_deg*2);
    [Xi_Z_Xi,Z,~] = sosdecvar('Z',tmp(1:size(Tx,1)));
    % Mx = Tx'*Z*Tx;
    
%     nZrow = size(Z,1);
%     nTxcol = size(Tx,2);
%     tmp = mpvar('tmp',nZrow,nTxcol); % initialization
%     for i1 = 1 : nZrow
%         for i2 = 1 : nTxcol
%             tmp(i1,i2) = Z(i1,:)*Tx(:,i2);
%         end
%     end
%     
%     Mx = mpvar('Mx',nTxcol,nTxcol);
%     for i1 = 1 : nTxcol
%         for i2 = 1 : nTxcol
%             Mx(i1,i2) = Tx(:,i1)'*tmp(:,i2);
%         end
%     end
    
    nZrow = size(Z,1);
    nTxrow = size(Tx,1);
    nTxcol = size(Tx,2);    
    tmp = mpvar('tmp',nZrow,nTxcol); % initialization
    for i1 = 1 : nZrow
        for i2 = 1 : nTxcol
            tmp2 = mpvar('tmp2',nTxrow,1);
            parfor i3 = 1 : nTxrow
                tmp2(i3) = Z(i1,i3)*Tx(i3,i2);
            end
            tmp(i1,i2) = ones(1,nTxrow)*tmp2;            
        end
      
    end
    clearvars tmp2
    
    Mx = mpvar('Mx',nTxcol,nTxcol);
    for i1 = 1 : nTxcol
        for i2 = 1 : nTxcol
            tmp3 = mpvar('tmp3',nTxrow,1);
            parfor i3 = 1 : nTxrow
                tmp3(i3) = Tx(i3,i1)'*tmp(i3,i2);
            end
            Mx(i1,i2) = ones(1,nTxrow)*tmp3;            
        end
      
    end
    
%     % quadratic cost function relaxation constraint
%     Xi = monomials(x,0:zx_deg*2);
%     Xi = Xi(1:size(Z,1));
%     
%     nZrow = size(Z,1);
%     nXirow = size(Tx,2);
%     
%     tmp2 = mpvar('tmp2',nZrow,1);
%     for i1 = 1 : nZrow
%         tmp2(i1) = Z(i1,:)*Xi;
%     end
%     Xi_Z_Xi = Xi'*tmp2;
    
    
    disp('Finish calculating Mx and Xi_Z_Xi and save the mat file.');
end

% Equality constraints
eq_const1 = sx-Mx(1);
eq_const2 = cx-Mx(2);
eq_const3 = ax*R^(-1)-Mx(4);
eq_const1.coefficient(find(abs(eq_const1.coefficient)<=0)) = 0;
eq_const2.coefficient(find(abs(eq_const2.coefficient)<=0)) = 0;
eq_const3.coefficient(find(abs(eq_const3.coefficient)<=0)) = 0;

sos_const2 = Xi_Z_Xi;
sos_const2.coefficient(find(abs(sos_const2.coefficient)<=0)) = 0;

% cost function
qx = x1^2 + x2^2;
d1 = p2s(qx)/p2s((bx)^alpha)*p2s(Psi(1:size(C_a_amon,1),1));
d2 = 1/p2s((bx)^alpha)*p2s(Psi(1:size(C_c_cmon,1),1));
rx1n = [-5, -0.1]; rx1p = [0.1, 5]; rx2n = [-5, -0.1]; rx2p = [0.1, 5];
d1i = []; d2i = [];
for i1 = 1 : length(d1)
    if hasSymType(d1(i1),'variable') == 1
        varnames = symvar(d1(i1));
        if length(varnames) == 2 
            tmp = integral2(matlabFunction(d1(i1)),rx1n(1),rx1n(2),rx2n(1),rx2n(2)) ...
                + integral2(matlabFunction(d1(i1)),rx1p(1),rx1p(2),rx2p(1),rx2p(2));
        elseif isequal(varnames,p2s(x1))
            tmp = integral2(matlabFunction(d1(i1)),rx1n(1),rx1n(2)) ...
                + integral2(matlabFunction(d1(i1)),rx1p(1),rx1p(2));
        elseif isequal(varnames,p2s(x2))
            tmp = integral2(matlabFunction(d1(i1)),rx2n(1),rx2n(2)) ...
                + integral2(matlabFunction(d1(i1)),rx2p(1),rx2p(2));
        end
    else
        tmp = double(d1(i1));
    end
    d1i = [d1i; tmp];
end
for i1 = 1 : length(d2)
    if hasSymType(d2(i1),'variable') == 1
        varnames = symvar(d2(i1));
        if length(varnames) == 2 
            tmp = integral2(matlabFunction(d2(i1)),rx1n(1),rx1n(2),rx2n(1),rx2n(2)) ...
                + integral2(matlabFunction(d2(i1)),rx1p(1),rx1p(2),rx2p(1),rx2p(2));
        else
            if isequal(varnames,p2s(x1))
                tmp = integral(matlabFunction(d2(i1)),rx1n(1),rx1n(2)) ...
                    + integral(matlabFunction(d2(i1)),rx1p(1),rx1p(2));
            else
                tmp = integral(matlabFunction(d2(i1)),rx2n(1),rx2n(2)) ...
                    + integral(matlabFunction(d2(i1)),rx2p(1),rx2p(2));
            end
        end
    else
        tmp = double(d2(i1));
    end
    d2i = [d2i; tmp];
end
d1i = double(d1i); d2i = double(d2i);
d1i = [d1i;zeros(ndic-size(C_a_amon,1),1)];
d2i = [d2i;zeros(ndic-size(C_c_cmon,1),1)];
obj = d1i'*C_a_Psi + d2i'*C_s_Psi;

% sos formulation
MM = 100;
MC = 2;
eps = 0;
pconstr = polyconstr;
pconstr(1) = sos_const1 >= gamma*ax*bx;
pconstr(2) = sos_const2 >= eps;
pconstr(3) = eq_const1 == eps;
pconstr(4) = eq_const2 == eps;
pconstr(5) = eq_const3 == eps;
if opts.OCP.deg_a(end) ~= 0
    pconstr(6) = ax >= eps;
end
%pconstr(7) = MM*ax - MC*sx >= eps;

% solve SOS problem
sosopts = sosoptions;
% sosopts.simplify = 'off';
% sosopts.form = 'kernel';
% sosopts.feastol = 1e-12;
[info,dopt,sossol] = sosopt(pconstr, x, obj, sosopts);

if info.feas == 1 % SOS is feasible
    disp('SOS solution was found.');
else % SOS is infeasible
    disp('SOS was infeasible.');
    return;
end

% solution for a(x), denoted by 'ax'
C_a_Psi_sol = subs(C_a_Psi, dopt);
ax_sol = C_a_Psi_sol'*Psi;
% ax.coefficient = round(ax.coefficient,6);

% solution for c(x), denoted by 'cx'
C_c_Psi_sol = subs(C_c_Psi, dopt);
cx_sol = C_c_Psi_sol'*Psi;
% cx.coefficient = round(cx.coefficient,6);

% solution for ab(x), denoted by 'cx'
C_ab_Psi_sol = subs(C_ab_Psi, dopt);
ab_sol = C_ab_Psi_sol'*Psi;

% solution for ab(x), denoted by 'cx'
C_bc_Psi_sol = subs(C_bc_Psi, dopt);
bc_sol = C_bc_Psi_sol'*Psi;

ZZ = subs(Z,dopt);


%% Test estimated controller
ys = -5;
ye = 5;
dy = 0.5;
y1s = -6;
y1e = 6;

tspan = [0 40];

% -----------------------------------------------------------------------
% Draw a vector field
% -----------------------------------------------------------------------
[y1i, y2i] = meshgrid(y1s:dy:y1e, y1s:dy:y1e);
yi = [y1i(:), y2i(:)];

ucont = zeros(size(y1i));
vcont = zeros(size(y2i));

for i1 = 1 : numel(y1i)
    xdot = feval(@(y) model_VDP_control(0,y,x,ax_sol,cx_sol), yi(i1,:)');
    ucont(i1) = xdot(1);
    vcont(i1) = xdot(2);
end


% -----------------------------------------------------------------------
% Simulate the controlled system from arbitrary initial points
% -----------------------------------------------------------------------
xinits = [-4,-4;4,-4;-4,4;4,4];
Ninit = size(xinits,1);
hbar = parfor_progressbar(numel(y1i), 'Simulating trajectories..');
parfor i1 = 1 : length(xinits)
    init = xinits(i1,:)';
    [tcont{i1}, ycont{i1}] = ode15s(@(t,y) model_VDP_control(t,y,x,ax_sol,cx_sol), tspan, init);
    hbar.iterate(1);
end
close(hbar);


figure;
hold on; 
quiver(y1i,y2i,ucont,vcont);

% Open loop trajectory plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nsim = 2000;
deltaT = 0.01;
f_open =  @(t, x, u)( [x(2,:) ; -x(1,:)+x(2,:).*(1-x(1,:).^2) + u] );
f_d = @(t,x,u) (x + deltaT*f_open(t, x, u));

x444 = [y1s, y1s, y1e, y1e, y1s]';  % for SOS automatica paper
y444 = [y1s, y1e, y1e, y1s, y1s];

cc1 = cell2mat(ycont(1,1));
rr2 = plot(cc1(:,1),cc1(:,2),'k','LineWidth',2);
cellfun(@(x) plot(x(:,1),x(:,2),'k','LineWidth',2), ycont); 
plot(x444, y444,'color','k', 'LineWidth', 2);

for j = 1:size(xinits,1)
    x_open = xinits(j,:)';
    u_value = [];    
    for i = 0:Nsim-1
        x_open = [x_open, f_d(0,x_open(:,end), 0)];    
    end
    rr1 = plot(x_open(1,:),x_open(2,:), '--r','LineWidth',2);
    plot(ycont{j}(1,1),ycont{j}(1,2),'o','MarkerSize',10,'markerfacecolor','g');
    plot(0,0,'o','MarkerSize',10,'markerfacecolor','m');
    
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hold off;
axis tight equal;
xlim([y1s y1e]); 
ylim([y1s y1e]);
ylabel('$x_2$','Interpreter','Latex','FontSize', 25);
xlabel('$x_1$','Interpreter','Latex','FontSize', 25);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 25, 'TickLabelInterpreter', 'Latex');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 25, 'TickLabelInterpreter', 'Latex');
lgr = legend([rr1,rr2],'$Open\;\;loop$', '$Closed\;\;loop$','Interpreter','Latex');
lgr.FontSize = 15;



%Cost_fn = (p2s(qx)*p2s(ax))/(p2s(bx)^alpha) + (abs(p2s(cx)))/(p2s(bx)^alpha);
% Cost_fn = (p2s(qx)) + (p2s(cx_sol)/p2s(ax_sol))^2;
% figure;
% Cost = zeros(size(ycont));
% for ii = 1:length(ycont)
%     F = cell2mat(ycont(1,ii));
%     for jj = 1:size(F,1)
%        temp11 = subs(Cost_fn,p2s(x),F(jj,:)');
%        Cost(1,ii) = Cost(1,ii) + temp11;
%        u11(1,jj) =  subs(p2s(cx_sol),p2s(x), F(jj,:)') / subs(p2s(ax_sol),p2s(x), F(jj,:)');
%     end
%     tt = 1:size(F,1);
%     %figure;
%     plot(tt,u11(1:size(F,1)),'LineWidth',1.5);
%     hold on;
%     fontsize = 20;
% ylabel('$k(x)$','Interpreter','Latex','FontSize', fontsize);
% xlabel('$t$','Interpreter','Latex','FontSize', fontsize);
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 20, 'TickLabelInterpreter', 'Latex');
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 20, 'TickLabelInterpreter', 'Latex');
% axis equal;
% %xlim([0 200]); ylim([-10 10]);
% end








