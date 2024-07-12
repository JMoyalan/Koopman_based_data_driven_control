clc; clear all; close all;

load('../gEDMD/result_gEDMD.mat');

UU = [2;2];
dT = 0.01;
dT1 = 0.1;
dT2 = 0.2;
dT3 = 0.3;
ns = 1e3;
ns1 = 1e4;
U = zeros(2,ns);

U(:,1) = UU;
for  i = 1:ns-1
    F(1) = UU(2);
    F(2) = (1-UU(1)^2)*UU(2)- UU(1);
    U(1,i+1) = UU(1,1) + F(1)*dT;
    U(2,i+1) = UU(2,1) + F(2)*dT;
    UU = U(:,i+1);
end

Psi_X = zeros(size(L1,1),ns1);

Psi_X(:,1) = double(subs(p2s(Psi_poly),x, U(:,1)));

for j = 1:ns1-1
    Psi_dot = L1*Psi_X(:,j);
    Psi_X(:,j+1) = Psi_X(:,j) + Psi_dot*0.001;
%     Psi_X(:,j+1) = Psi_dot*0.01;
end
aa = zeros(2,size(L1,1));
for k = 1:ns
    a = U(:,k)*pinv(Psi_X(:,k));
    aa = aa + a;
end
aa = aa/ns;

Y = zeros(2,ns1);

for i = 1:ns1
%     Y(1,i) = Psi_X(2,i);
%     Y(2,i) = Psi_X(11,i); 
    Y(:,i) = aa*Psi_X(:,i);
end

% Y = zeros(2,ns);
% YY = [2;2];
% Y(:,1) = YY;
% for  i = 1:ns-1
%     FF(1) = YY(2);
%     FF(2) = (1-YY(1)^2)*YY(2)- YY(1);
%     Y(1,i+1) = YY(1,1) + FF(1)*dT1;
%     Y(2,i+1) = YY(2,1) + FF(2)*dT1;
%     YY = Y(:,i+1);
% end

% W = zeros(2,ns);
% WW = [2;2];
% W(:,1) = WW;
% for  i = 1:ns-1
%     FF(1) = WW(2);
%     FF(2) = (1-WW(1)^2)*WW(2)- WW(1);
%     W(1,i+1) = WW(1,1) + FF(1)*dT2;
%     W(2,i+1) = WW(2,1) + FF(2)*dT2;
%     WW = W(:,i+1);
% end
% 
% V = zeros(2,ns);
% VV = [2;2];
% V(:,1) = VV;
% for  i = 1:ns-1
%     FF(1) = VV(2);
%     FF(2) = (1-VV(1)^2)*VV(2)- VV(1);
%     V(1,i+1) = VV(1,1) + FF(1)*dT3;
%     V(2,i+1) = VV(2,1) + FF(2)*dT3;
%     VV = V(:,i+1);
% end



figure;
%quiver(y1i,y2i,ucont,vcont);
hold on;
rr1 = plot(U(1,:),U(2,:),'LineWidth',1.5);
rr2 = plot(Y(1,:),Y(2,:),'LineWidth',1.5);
% rr3 = plot(W(1,:),W(2,:),'LineWidth',1.5);
% rr4 = plot(V(1,:),V(2,:),'LineWidth',1.5);


ylabel('$x_2$','Interpreter','Latex','FontSize', 40);
xlabel('$x_1$','Interpreter','Latex','FontSize', 40);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 40, 'TickLabelInterpreter', 'Latex');
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 40, 'TickLabelInterpreter', 'Latex');
axis tight equal;
xlim([-6 6]); 
ylim([-5 5]); 
% lgr = legend([rr1, rr2, rr3, rr4], '$Uncontrolled \; dynamics$', '$Order \; 9$','$Order \; 7$','$Order \; 5$','Interpreter','Latex');
% lgr.FontSize = 15;

hold off;


