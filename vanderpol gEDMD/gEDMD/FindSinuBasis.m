function [Psi_sinu, dPsi_sinu] = FindSinuBasis(y,dy,opts)
ny = length(y);

% Sinusoidal basis
Psi_sinu = [];
for i1 = 1 : opts.n    
    Psi_sinu = [Psi_sinu; sin(i1*y); cos(i1*y)];
end

% Derivative of Sinusoidal basis
dPsi_sinu = jacobian(Psi_sinu,y)*dy;