function [Psi_TPSRBF, dPsi_TPSRBF, c] = FindTPSRBF(x,dx,opts)
% NOTE: Structure of a matrix of center points is
%       number of center points by number of states, and later,
%       in FindBasis.m, c is vectorized by the command c(:).
nc = opts.nc; % number of center points
nx = length(x); % number of states
syms c [nc,nx] real; % symbolic expression of center points

Psi_TPSRBF = [];
for i1 = 1 : nc
    r = norm(x-transpose(c(i1,:)));
    Psi_TPSRBF = [Psi_TPSRBF; r^2*log(r)];
end

% Derivative of Sinusoidal basis
dPsi_TPSRBF = jacobian(Psi_TPSRBF,x)*dx;