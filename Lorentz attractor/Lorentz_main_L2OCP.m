
clear; close all; format compact; clc;

addpath('../SOS Toolbox');
addpath('../SOS Toolbox/multipoly/multipoly');
addpath('../SOS Toolbox/sosopt/sosopt');
addpath('../SOS Toolbox/sedumi-master/sedumi-master');


% load sample data for Case 1 and 2, X_i and Y_i
load('./sampling/sample_Lorenz.mat');
fprintf(sprintf('Number of samples (Case 1) = %d\n',size(X1,2)));
fprintf(sprintf('Number of samples (Case 2) = %d\n',size(X2,2)));

rho = param.rho;
sigma = param.sigma;
beta = param.beta;

% -----------------------------------------------------------------------
% Now, denote Koopman operator estimated from the case (a) as K1, 
% (b) as K2. Estimate K1 and K2 using EDMD algorithm with monomials 
% of polynomial function as dictionary functions.
% -----------------------------------------------------------------------

% create polynomial variables x and u
pvar x1 x2 x3 real% state variables
pvar u % control
x = [x1; x2; x3]; % state variable vector

% model, dxdt = fx + gx*u

fx = [sigma*(x2-x1);
    x1*(rho-x3) - x2;
    x1*x2 - beta*x3;];
gx = [0; 1; 0];

% constant parameters, alpha and b(x)
alpha = 4;
gamma = 0;
% bmon = monomials(x, 0:1);
% bvar_mat = sprandsym(length(bmon), 2/50, 0.1, 1);
% bx = bmon' * bvar_mat * bmon;
% bx = 3*x1^2 + 2*x1*x2 + 2*x2^2;
bx = x'*x;

% finding b(x) as a constrol Lyapunov function of linearization
% [Vlin,Alin,Plin] = linstab(fx,x);
% bx = Vlin;

% finding b(x) as a control Lyapunov function of linearization using LQR
% Alin = double(subs(jacobian(fx,x),x,[0;0;0]));
% Q = [1 0.5 0.5;
%     0.5 1 0.5;
%     0.5 0.5 1];
% R = 3.5;
% N = [1.5; 0; 1];
% [~,Plin,~] = lqr(Alin,gx,Q,R,N);
% bx = x'*Plin*x;

% checking if b(x) >= 0
% if issos(bx) ~= 1
%     disp('b(x) is not a SOS polynomial');
%     return
% end

% -----------------------------------------------------------------------
% Create orthogonal dictionary functions, and coefficient matries of
% a(x), c(x), a(x)*b(x), and b(x)*c(x) in terms of Psi(x)
% -----------------------------------------------------------------------
% degree of monomials of a(x) and c(x), i.e., M_a(x), M_c(x)
deg_a = 0 : 1; % minimum degree : maximum degree of a(x)
deg_c = 1 : 3; % minimum degree : maximum degree of c(x)
deg_s = 1 : 4; % minimum degree : maximum degree of s(x)

% type and order of orthogonal polynomial basis
order = 7;
type = 'Monomial';
% type = 'Legendre';
% type = 'Hermite';

% in case of returnlegendre polynomial, set interval
opts.interval = [-5 5];

% options for numerical precisions of coefficients
% for example, some of constant coefficients can occur in the conversion
% from coefficient vectors in terms of monomials to coefficient vectors
% in terms of basis vectors, and due to numerical 
opts.round = 1e-6;

% options for tolerance of error
opts.errnrm = 1e-6;

% generates coefficient matrices and dictioanry functions
[C_s_Psi,C_a_Psi, C_c_Psi, C_ab_Psi, C_bc_Psi, C_x_Psi, Psi, polys, decvars,decvars_1] = ...
    FindOrthoPoly_1(x,bx,deg_a,deg_c,deg_s,order,type,opts);

% [C_a_Psi, C_c_Psi, C_ab_Psi, C_bc_Psi, C_x_Psi, Psi, polys, decvars] = ...
%     FindMonoPoly(x,bx,deg_a,deg_c,order,type,opts);

% estimate K1 and K2
edmdopt.method = 1;
edmdopt.batch = 10;

% estimate K1 (Case 1)
[K1, G1, A1, Xr1, Yr1, Psi_X1, Psi_Y1] = EDMD_poly(X1, Y1, x, Psi, edmdopt);
fprintf(sprintf('Estimated ||A1 - G1*K1|| = %f\n', norm(A1-G1*K1, 'fro')));

% estimate K1 (Case 1)
[K2, G2, A2, Xr2, Yr2, Psi_X2, Psi_Y2] = EDMD_poly(X2, Y2, x, Psi, edmdopt);
fprintf(sprintf('Estimated ||A2 - G2*K2|| = %f\n', norm(A2-G2*K2, 'fro')));


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
if deg_a == 0
    C_a_Psi = subs(C_a_Psi,decvars(1),1);
    C_ab_Psi = subs(C_ab_Psi,decvars(1),1);
end
a_x = C_a_Psi'*Psi;
cx = C_c_Psi'*Psi;
sx = C_s_Psi'*Psi;
nx = length(x); % number of states
L1 = (K1' - eye(ndic))/dt; % Koopman generator for Case 1
L2 = (K2'- K1')/dt; % Koopman generator for Case 2
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

Div_F_model = double(sum(diag(jacobian(fx,x))));

term(1) = C_a_Psi'*L1*Psi + C_a_Psi'*Psi*Div_F;
term(2) = C_c_Psi'*L2*Psi + C_c_Psi'*Psi*Div_G;
term(3) = C_ab_Psi'*L1*Psi + C_ab_Psi'*Psi*Div_F;
term(4) = C_bc_Psi'*L2*Psi + C_bc_Psi'*Psi*Div_G;

% Formulate stability equation
sos_poly = (1+alpha)*bx*sum(term(1:2)) - alpha*sum(term(3:4));
sos_poly.coefficient(find(abs(sos_poly.coefficient)<=1e-0)) = 0;

qx = x'*x;

d1 = (p2s(qx)/p2s((bx)^alpha))*p2s(Psi);
d2 = (1/p2s((bx)^alpha))*p2s(Psi);
rx1n = [-2.5, -0.1]; rx1p = [0.1, 2.5]; rx2n = [-2.5, -0.1]; rx2p = [0.1, 2.5]; rx3n = [-2.5, -0.1]; rx3p = [0.1, 2.5];
d1i = []; d2i = [];


parfor i1 = 1 : length(d1)
    if hasSymType(d1(i1),'variable') == 1
        varnames = symvar(d1(i1));
        if length(varnames) == 3 
            tmp = integral3(matlabFunction(d1(i1)),rx1n(1),rx1n(2),rx2n(1),rx2n(2),rx3n(1),rx3n(2)) ...
                + integral3(matlabFunction(d1(i1)),rx1p(1),rx1p(2),rx2p(1),rx2p(2),rx3p(1),rx3p(2));
        elseif isequal(varnames,p2s(x1))
            tmp = integral1(matlabFunction(d1(i1)),rx1n(1),rx1n(2)) ...
                + integral1(matlabFunction(d1(i1)),rx1p(1),rx1p(2));
        elseif isequal(varnames,p2s(x2))
            tmp = integral1(matlabFunction(d1(i1)),rx2n(1),rx2n(2)) ...
                + integral1(matlabFunction(d1(i1)),rx2p(1),rx2p(2));
        elseif isequal(varnames,p2s(x3))
            tmp = integral1(matlabFunction(d1(i1)),rx3n(1),rx3n(2)) ...
                + integral1(matlabFunction(d1(i1)),rx3p(1),rx3p(2));
        elseif length(varnames) == 2
            if isequal(varnames,p2s(x1x2))
                tmp = integral2(matlabFunction(d1(i1)),rx1n(1),rx1n(2),rx2n(1),rx2n(2)) ...
                    + integral2(matlabFunction(d1(i1)),rx1p(1),rx1p(2),rx2p(1),rx2p(2));
            elseif isequal(varnames,p2s(x2x3))
                tmp = integral2(matlabFunction(d1(i1)),rx2n(1),rx2n(2),rx3n(1),rx3n(2)) ...
                    + integral2(matlabFunction(d1(i1)),rx2p(1),rx2p(2),rx3p(1),rx3p(2));
            elseif isequal(varnames,p2s(x1x3))
                tmp = integral2(matlabFunction(d1(i1)),rx1n(1),rx1n(2),rx3n(1),rx3n(2)) ...
                    + integral2(matlabFunction(d1(i1)),rx1p(1),rx1p(2),rx3p(1),rx3p(2));
            end
        end
    else
        tmp = double(d1(i1));
    end
    d1i = [d1i; tmp];
 end 
for i1 = 1 : length(d2)
    if hasSymType(d2(i1),'variable') == 1
        varnames = symvar(d2(i1));
        if length(varnames) == 3 
            tmp1 = integral3(matlabFunction(d2(i1)),rx1n(1),rx1n(2),rx2n(1),rx2n(2),rx3n(1),rx3n(2)) ...
                + integral3(matlabFunction(d2(i1)),rx1p(1),rx1p(2),rx2p(1),rx2p(2),rx3p(1),rx3p(2));
        elseif isequal(varnames,p2s(x1))
            tmp1 = integral1(matlabFunction(d2(i1)),rx1n(1),rx1n(2)) ...
                + integral1(matlabFunction(d2(i1)),rx1p(1),rx1p(2));
        elseif isequal(varnames,p2s(x2))
            tmp1 = integral1(matlabFunction(d2(i1)),rx2n(1),rx2n(2)) ...
                + integral1(matlabFunction(d2(i1)),rx2p(1),rx2p(2));
        elseif isequal(varnames,p2s(x3))
            tmp1 = integral1(matlabFunction(d2(i1)),rx3n(1),rx3n(2)) ...
                + integral1(matlabFunction(d2(i1)),rx3p(1),rx3p(2));
        elseif length(varnames) == 2
            if isequal(varnames,p2s(x1x2))
                tmp1 = integral2(matlabFunction(d2(i1)),rx1n(1),rx1n(2),rx2n(1),rx2n(2)) ...
                    + integral2(matlabFunction(d2(i1)),rx1p(1),rx1p(2),rx2p(1),rx2p(2));
            elseif isequal(varnames,p2s(x2x3))
                tmp1 = integral2(matlabFunction(d2(i1)),rx2n(1),rx2n(2),rx3n(1),rx3n(2)) ...
                    + integral2(matlabFunction(d2(i1)),rx2p(1),rx2p(2),rx3p(1),rx3p(2));
            elseif isequal(varnames,p2s(x1x3))
                tmp1 = integral2(matlabFunction(d2(i1)),rx1n(1),rx1n(2),rx3n(1),rx3n(2)) ...
                    + integral2(matlabFunction(d2(i1)),rx1p(1),rx1p(2),rx3p(1),rx3p(2));
            end
        end
    else
        tmp1 = double(d2(i1));
    end
    d2i = [d2i; tmp1];
end   

obj = d1i'*C_a_Psi + d2i'*C_s_Psi;

disp('Cost function integration done.');

L1min = 'Pro';

if strcmp(L1min, 'yes')
    % -----------------------------------------------------------------------
    % Solve SOS with L1 norm minimization
    % -----------------------------------------------------------------------
    % formulate SOS objective function, ojb = c^T*t
    ndecvar = length(decvars);
    t = mpvar('t',ndecvar,1);
    c = ones(ndecvar,1);
    obj = c'*t;
    
    % formulate SOS constraints
    pconstr = polyconstr; % initialize polynomial constraint
    pconstr(1) = sos_poly >= gamma*a_x*bx; % SOS constraint
    for i1 = 1 : ndecvar
        pconstr(i1+1) = -decvars(i1) + t(i1) >= 0;
    end
    for i1 = 1 : ndecvar
        pconstr(i1+ndecvar+1) = decvars(i1) + t(i1) >= 0;
    end
    
    % solve SOS problem
    sosopts = sosoptions;
    % opts.simplify = 'on';
    % opts.form = 'kernel';
    % opts.feastol = 1e-9;
    [info,dopt,sossol] = sosopt(pconstr,x,obj,sosopts);
elseif strcmp(L1min, 'no')
    % -----------------------------------------------------------------------
    % Solve SOS as a feasibility problem
    % -----------------------------------------------------------------------
    % formulate SOS constraints
    pconstr = polyconstr; % initialize polynomial constraint
    pconstr(1) = sos_poly >= 0; % SOS constraint
    
    % solve SOS problem
    sosopts = sosoptions;
    % sosopts.simplify = 'on';
    % sosopts.form = 'kernel';
    % sosopts.feastol = 1e-12;
    [info,dopt,sossol] = sosopt(pconstr, x, sosopts);
elseif strcmp(L1min,'Pro')
    pconstr = polyconstr;
pconstr(1) = sos_poly >= gamma*a_x*bx; % SOS constraint
pconstr(2) = (sx - cx) >= 0;
pconstr(3) = (sx + cx) >= 0;

if deg_a(end) ~= 0
    pconstr(4) = a_x >= 0;
end

sosopts = sosoptions;
[info,dopt,sossol] = sosopt(pconstr,x,obj,sosopts);
    
end

if info.feas == 1 % SOS is feasible
    disp('SOS solution was found.');
else % SOS is infeasible
    disp('SOS was infeasible.');

    return;
end

% solution for a(x), denoted by 'ax'
C_a_Psi_sol = subs(C_a_Psi, dopt);
ax = C_a_Psi_sol'*Psi;
% ax.coefficient = round(ax.coefficient,6);

% solution for c(x), denoted by 'cx'
C_c_Psi_sol = subs(C_c_Psi, dopt);
cx = C_c_Psi_sol'*Psi;
% cx.coefficient = round(cx.coefficient,6);

% solution for s(x), denoted by 'sx'
C_s_Psi_sol = subs(C_s_Psi, dopt);
sx = C_s_Psi_sol'*Psi;
% sx.coefficient = round(sx.coefficient,6);

% solution for ab(x), denoted by 'cx'
C_ab_Psi_sol = subs(C_ab_Psi, dopt);
ab = C_ab_Psi_sol'*Psi;

% solution for ab(x), denoted by 'cx'
C_bc_Psi_sol = subs(C_bc_Psi, dopt);
bc = C_bc_Psi_sol'*Psi;


%% Test estimated controller
ys = opts.interval(1);
ye = opts.interval(2);
dy = 0.1;

tspan = [0 5];


xxs1 = -30;
xxe1 = 30;
xxs2 = -30;
xxe2 = 30;
xxs3 = 0;
xxe3 = 50;
xxd = 8;
[x1i, x2i, x3i] = meshgrid(xxs1:xxd:xxe1, xxs2:xxd:xxe2, xxs3:xxd:xxe3);
xxi = [x1i(:), x2i(:), x3i(:)];

ucont = zeros(size(x1i));
vcont = zeros(size(x2i));
wcont = zeros(size(x3i));

for i1 = 1 : numel(x1i)
    xdot = feval(@(y) model_Lorentz_control(0,y,x,ax,cx,param), xxi(i1,:)');
    ucont(i1) = xdot(1);
    vcont(i1) = xdot(2);
    wcont(i1) = xdot(3);
end

% -----------------------------------------------------------------------
% Simulate the controlled system from arbitrary initial points
% -----------------------------------------------------------------------
xst = -10; xend = 10; Nsample = 5;
clearvars xinits
xinits(:,1) = -xst + (xend-xst).*rand(Nsample,1);
xinits(:,2) = -xst + (xend-xst).*rand(Nsample,1);
xinits(:,3) = -xst + (xend-xst).*rand(Nsample,1);
xinits = [-16,16,16; -16,16,26; -16,16,36; 16,-16,16; 16,-16,26; 16,-16,36];
u = 0;


for i1 = 1 : size(xinits,1)
    init = xinits(i1,:)';
    [tcont{i1}, ycont{i1}] = ode15s(@(t,y) model_Lorentz_control(t,y,x,ax,cx,param), tspan, init);
    [tmodel{i1}, ymodel{i1}] = ode15s(@(t,x) model_Lorentz(t,x,u,param), tspan, init);
    
end


% plot 1: states vs time
fontsize = 25;
fontsize2 = 32;
linewidth = 1.5;
% figure(1); hold on;
% for i1 = 1 : size(xinits,1)
%         plot(tmodel{i1}, ymodel{i1}(:,1), '--r');
%         plot(tmodel{i1}, ymodel{i1}(:,2), '--g'); 
%         plot(tmodel{i1}, ymodel{i1}(:,3), '--b'); 
%         plot(tcont{i1}, ycont{i1}(:,1), 'r');
%         plot(tcont{i1}, ycont{i1}(:,2), 'g');
%         plot(tcont{i1}, ycont{i1}(:,3), 'b');
%         
%         xlabel('$t$','Interpreter','Latex','FontSize', fontsize);
%         ylabel('$x_{1\sim3}$','Interpreter','Latex','FontSize', fontsize);
%         xAX = get(gca,'XAxis');
%         set(xAX,'FontSize', fontsize, 'TickLabelInterpreter', 'Latex');
%         yAX = get(gca,'YAxis');
%         set(yAX,'FontSize', fontsize, 'TickLabelInterpreter', 'Latex');
%         
% %         if i1 == 1
% %             L = legend('$x_1$ (no control)','$x_2$ (no control)','$x_3$ (no control)','$x_1$ (controlled)','$x_2$ (controlled)','$x_3$ (controlled)');
% %             set(L,'Interpreter','Latex','FontSize',fontsize2);
% %         end
% end


% plot 2: 3 dimensional plot
figure(2); hold on; 
% quiver3(x1i,x2i,x3i,ucont,vcont,wcont,'LineWidth',1);
grid on;
for i1 = 1 : size(xinits,1)
    rr1 = plot3(ymodel{i1}(:,1),ymodel{i1}(:,2),ymodel{i1}(:,3), '--r','LineWidth',linewidth);
    rr2 = plot3(ycont{i1}(:,1),ycont{i1}(:,2),ycont{i1}(:,3), 'b','LineWidth',linewidth);
    ll1 = plot3(ycont{i1}(1,1),ycont{i1}(1,2),ycont{i1}(1,3),'o','MarkerSize',15,'markerfacecolor','g');
    ll2 = plot3(0,0,0,'o','MarkerSize',15,'markerfacecolor','m');
    
    xlabel('$x_1$','Interpreter','Latex','FontSize', fontsize);
    ylabel('$x_2$','Interpreter','Latex','FontSize', fontsize);
    zlabel('$x_3$','Interpreter','Latex','FontSize', fontsize);
    xAX = get(gca,'XAxis');
    set(xAX,'FontSize', fontsize, 'TickLabelInterpreter', 'Latex');
    yAX = get(gca,'YAxis');
    set(yAX,'FontSize', fontsize, 'TickLabelInterpreter', 'Latex');
    zAX = get(gca,'ZAxis');
    set(zAX,'FontSize', fontsize, 'TickLabelInterpreter', 'Latex');
    
%     if i1 == 1
%         L = legend('No control', 'Controlled');
%         set(L,'Interpreter','Latex','FontSize',fontsize2);
%         grid on;
%     end
lgr = legend([rr1, rr2], '$Uncontrolled \;\;trajectory$', '$Controlled \;\;trajectory$','Interpreter','Latex');
lgr.FontSize = 10;
end
axis equal;
hold off;


% figure(3);
% for ii = 1:length(ycont)
%     F = cell2mat(ycont(1,ii));
%     for jj = 1:size(F,1)
% %        temp11 = subs(Cost_fn,p2s(x),F(jj,:)');
% %        Cost(1,ii) = Cost(1,ii) + temp11;
%        u11(1,jj) =  subs(p2s(cx),p2s(x), F(jj,:)') / subs(p2s(ax),p2s(x), F(jj,:)');
%     end
%     tt = 1:size(F,1);
%     %figure;
%     plot(tt,u11(1:size(F,1)),'LineWidth',1.5);
%     hold on;
%     fontsize = 25;
% ylabel('$k(x)$','Interpreter','Latex','FontSize', fontsize);
% xlabel('$t$','Interpreter','Latex','FontSize', fontsize);
% xAX = get(gca,'XAxis');
% set(xAX,'FontSize', 25, 'TickLabelInterpreter', 'Latex');
% yAX = get(gca,'YAxis');
% set(yAX,'FontSize', 25, 'TickLabelInterpreter', 'Latex');
% xlim([0 200]);
% %ylim([-10 10]);
% axis equal;
% end
% hold off;
