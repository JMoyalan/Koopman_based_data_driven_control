function dy = model_kmc_cont(t,x,Psi,L,C_x_Psi,cp)
% % norm normalization
% xn = (x - norms.Xmean')./norms.Xnorm';
% un = (u - norms.Umean')./norms.Unorm';
% 
% Psi_X = [Psi([xn;cp(:)]); un];
% dyn = C_x_Psi'*L*Psi_X;
% dy = dyn.*norms.dXnorm' + norms.dXmean';

% % range normalization
% xn = (x - norms.Xmin')./(norms.Xmax'-norms.Xmin');
% un = (u - norms.Umin')./(norms.Umax'-norms.Umin');
% 
% Psi_X = [Psi([xn;cp(:)]); un];
% dyn = C_x_Psi'*L*Psi_X;
% dy = dyn.*(norms.dXmax'-norms.dXmin') + norms.dXmin';

% no normalization
Psi_X = Psi([x;cp(:)]);
dy = C_x_Psi'*L*Psi_X;