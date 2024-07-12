function [C_s_Psi,C_a_Psi,C_c_Psi,C_ab_Psi,C_bc_Psi,C_x_Psi,Psi,Psi_poly,dPsi,c] = FindBasis(x,dx,opts)
% Polynomial basis
if opts.polybasis.enable == 1
%     [C_x_Psi_poly, Psi_poly, dPsi_poly] = FindPolyBasis(s2p(x),s2p(dx),opts.polybasis.order,opts.polybasis.type,opts.polybasis);
    [C_s_Psi,C_a_Psi, C_c_Psi, C_ab_Psi, C_bc_Psi, C_x_Psi_poly, Psi_poly, dPsi_poly] = ...
        FindPolyBasis(s2p(x),s2p(dx),s2p(opts.OCP.bx),opts.OCP.deg_a,opts.OCP.deg_c,opts.OCP.deg_s,opts.polybasis.order,opts.polybasis.type,opts.polybasis);
    Psi_poly_sim = sym({});
    dPsi_poly_sim = sym({});
    for i1 = 1 : length(Psi_poly)
        Psi_poly_sim(i1,1) = p2s(Psi_poly(i1));
        dPsi_poly_sim(i1,1) = p2s(dPsi_poly(i1));
    end 
    fprintf(sprintf('[%s] Constructed polynomial basis\n', datetime('now')));
else
    Psi_poly = [];
    dPsi_poly = [];
    C_x_Psi_poly = [];
end

% Sinusoidal basis
if opts.sinubasis.enable == 1
    [Psi_sinu, dPsi_sinu] = FindSinuBasis(x,dx,opts.sinubasis);
    C_x_Psi_sinu = zeros(length(Psi_sinu),length(x));
    fprintf(sprintf('[%s] Constructed sinusoidal function basis\n', datetime('now')));
else
    Psi_sinu = [];
    dPsi_sinu = [];
    C_x_Psi_sinu = [];
end

% TPS RBF
if opts.tpsrbf.enable == 1
    [Psi_TPSRBF, dPsi_TBFRBF, c] = FindTPSRBF(x,dx,opts.tpsrbf);
    C_x_Psi_TPSRBF = zeros(length(Psi_TPSRBF),length(x));
    fprintf(sprintf('[%s] Constructed thin plate spline radial basis function (TPS RBF) basis\n', datetime('now')));
else
    Psi_TPSRBF = [];
    dPsi_TPSRBF = [];
    C_x_Psi_TPSRBF = [];
    c = [];
end

% combine basis
if opts.polybasis.enable == 1
    if opts.tpsrbf.enable == 1
        C_x_Psi = [C_x_Psi_poly; ...
            C_x_Psi_sinu; ...
            C_x_Psi_TPSRBF];
        Psi = matlabFunction([Psi_poly_sim; ...
            Psi_sinu; ...
            Psi_TPSRBF], ...
            'Vars', {[x;c(:)]});
        dPsi = matlabFunction([dPsi_poly_sim; ...
            dPsi_sinu; ...
            dPsi_TBFRBF], ...
            'Vars', {[x;dx;c(:)]});
    else
        C_x_Psi = [C_x_Psi_poly; C_x_Psi_sinu];
        Psi = matlabFunction([Psi_poly_sim; Psi_sinu], 'Vars', {x});
        dPsi = matlabFunction([dPsi_poly_sim; dPsi_sinu], 'Vars', {[x;dx]});
    end
else
    if opts.tpsrbf.enable == 1
        C_x_Psi = [C_x_Psi_sinu; ...
            C_x_Psi_TPSRBF];
        Psi = matlabFunction([Psi_sinu; ...
            Psi_TPSRBF], ...
            'Vars', {[x;c(:)]});
        dPsi = matlabFunction([dPsi_sinu; ...
            dPsi_TBFRBF], ...
            'Vars', {[x;dx;c(:)]});
    else
        C_x_Psi = [C_x_Psi_sinu];
        Psi = matlabFunction([Psi_sinu], 'Vars', {x});
        dPsi = matlabFunction([dPsi_sinu], 'Vars', {[x;dx]});
    end
end