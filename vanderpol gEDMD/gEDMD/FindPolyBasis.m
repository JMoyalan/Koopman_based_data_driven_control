% ==================================================================
% CODE TO CREATE VECTORS OF COEFFICIENT OF POLYNOMIALS AND
% MONOMIALS, AND POLYNOMIAL BASIS VECTORS
%
% HISTORY:
% 3/14 Created [HC]
% 4/18 Added decision variables in the output [HC]
% 5/11 Made orthogonal polynomial generation more efficient [HC]
% ==================================================================
% function [C_x_Psi, Psi, dPsi] = FindPolyBasis(x,dx,order,type,opts)
function [C_s_Psi,C_a_Psi, C_c_Psi, C_ab_Psi, C_bc_Psi, C_x_Psi, Psi, dPsi] = FindPolyBasis(x,dx,bx,deg_a,deg_c,deg_s,order,type,opts)
%% Create orthogonal polynomials
% ------------------------------------------------------------------
% Convert polynomial variables to symbolic variables
% (This is because MATLAB function "hermiteH" cannot deal with
% polynomial variables. So, we create orthogonal polynomials
% using symbolic toolbox, then convert it back to polynomial vars.
% ------------------------------------------------------------------
x = p2s(x); % symbolic state variables

% ------------------------------------------------------------------
% create orthogonal polynomials H_i(x)
% ------------------------------------------------------------------
if strcmp(type, 'Hermite')
    for i1 = 1 : order+1
        H{i1} = hermiteH(i1-1,x)';
    end
    Hmat = cat(1, H{:});
elseif strcmp(type, 'Legendre')
    P = legendrepol(order,opts.interval(1),opts.interval(2));
    p = [];
    for i1 = 1 : order+1
        p = [p; transpose(x.^(order-i1+1))];
    end
    Hmat = P*p;
    for i1 = 1 : size(Hmat,1)
        H{i1} = Hmat(i1,:);
    end
end


%% create orthogonal basis vector (dictionary functions), Psi(x)
% such that Psi = [H_0(x1)*H_0(x2), H_1(x1)*H_0(x2), ...,
%                  H_0(x1)*H_1(x2), H_1(x1)*H_1(x2), ...,
%                  H_0(x1)*H_2(x2), H_1(x1)*H_2(x2), ...].

% ------------------------------------------------------------------
% convert symbolic variables back to polynomial variables
% ------------------------------------------------------------------
x = s2p(x); % polynomial state variables
nx = length(x); % number of states

% ------------------------------------------------------------------
% create Psi(x) consisting of orthogonal polynomials
% ------------------------------------------------------------------
if strcmp(type, 'Hermite') || strcmp(type, 'Legendre')
    Psi = ones((order+1)^nx,1);
    for i1 = 1 : nx
        Psi = Psi.* ...
            repmat(kron(Hmat(:,i1),ones((order+1)^(i1-1),1)),(order+1)^(nx-i1),1);
    end
    Psi = s2p(Psi);
elseif strcmp(type, 'Monomial')
    if opts.noconst == 1
        Psi = monomials(x, 1:order);
    else
        Psi = monomials(x, 0:order);
    end
end

% ------------------------------------------------------------------
% number of dictionary functions
% ------------------------------------------------------------------
ndic = length(Psi);

%% Create common monomial basis vector, M_x(x)
% : This is a monomial basis vector common to all polynomials
%   including a(x), c(x), and even Psi(x).
%   We create those polynomials in the basis of M_x(x), then
%   later convert them to the basis of Psi(x).

% ------------------------------------------------------------------
% create common monomial vector, M_x(x)
% : max degree of common monomial vector should be greater than 
%   or equal to the max degree of monomial vector for Psi(x).
% ------------------------------------------------------------------
deg_xmon = 0 : Psi.maxdeg; % common monomial vector should be 
xmon = monomials(x, deg_xmon); % M_x(x)
nxmon = length(xmon); % number of elements in M_x(x)

% ------------------------------------------------------------------
% Unit test: check max degrees
% ------------------------------------------------------------------
if Psi.maxdeg > xmon.maxdeg
    error('Max degree of dictionary function is larger than the max degree of common monomials. Increase the degree of common monomials and try again.');
end

%% Reorganize Psi(x) in terms of M_x(x). Find a coefficient matrix
% C_Psi_xmon such that Psi(x) = C_Psi_xmon'*M_x(x).
% ------------------------------------------------------------------
% Find monomial vector of Psi(x), "C_Psi_Psimon", coefficient vector
% of Psi(x), "Psimon", such that Psi(x) = C_Psi_Psimon'*Psimon.
% ------------------------------------------------------------------
nPsimon = size(Psi.degmat,1);
tmp = repmat(x',nPsimon,1).^Psi.degmat(:,end-nx+1:end);
Psimon = ones(nPsimon,1);
for i1 = 1 : size(tmp,2)
    Psimon = Psimon.*tmp(:,i1); % M_Psi(x)
end
C_Psi_Psimon = Psi.coefficient; % C_Psi_Psimon

% ------------------------------------------------------------------
% unit test: check if Psi(x) = C_Psi_Psimon'*M_Psi(x)
% ------------------------------------------------------------------
err = Psi - (C_Psi_Psimon'*Psimon);
if norm(full(err.coefficient)) >= opts.errnrm
    error('Unit test failed: Psi(x) ~= transpose(C_Psi_Psimon)*Psimon. Check your coefficient calculation and try again.');
end

% ------------------------------------------------------------------
% Create coefficient vectors for Psi(x) in terms of M_x(x)
% ------------------------------------------------------------------
% arrange coefficient matrix of M_Psi(x) in terms of M_x(x)
% such that M_Psi(x) = cmat_Psi_xmon'*M_x(x).
cmat_Psi_xmon = zeros(nxmon,nPsimon);
[~,iPsimon2xmon] = ismember(xmon.varname, Psimon.varname);
[ixmon_Psi,iPsimon] = ismember(xmon.degmat, Psimon.degmat(:,iPsimon2xmon), 'rows','legacy');
cmat_Psi_xmon(find(ixmon_Psi==1),:) = Psimon.coefficient(nonzeros(iPsimon),:);
% [ixmon_Psi,iPsimon] = ismember(Psimon.degmat, xmon.degmat, 'rows');
% cmat_Psi_xmon(find(ixmon_Psi==1),:) = xmon.coefficient(nonzeros(iPsimon),:);


% create a new coefficient vector of Psi(x) in terms of M_x(x),
% "C_Psi_xmon", such that Psi(x) = C_Psi_xmon'*xmon.
C_Psi_xmon = cmat_Psi_xmon*C_Psi_Psimon;

% ------------------------------------------------------------------
% unit test: check if Psi(x) = C_Psi_xmon'*M_x(x)
% ------------------------------------------------------------------
err = Psi - (C_Psi_xmon'*xmon);
if norm(full(err.coefficient)) >= opts.errnrm
    error('Unit test failed: Psi(x) ~= transpose(C_Psi_Psimon)*Psimon.');
end

%% Create polynomial, a(x)
if deg_a == 0 % if a is constant
    C_a_amon = mpvar('a',1,1);
    C_a_Psi = [C_a_amon; zeros(ndic-1,1)];
    ax = C_a_amon;
else % if a(x) is a polynomial function
    % ------------------------------------------------------------------
    % create monomial vector, M_a(x), and coefficient vector, C_a_amon
    % such that a(x) = C_a_amon'*M_a(x). Subsequently, construct a(x).
    % ------------------------------------------------------------------
    deg_amon = deg_a; % minimum degree : maximum degree of a(x)
    amon = monomials(x, deg_amon); % M_a(x): monomials of a(x)
    namon = length(amon); % number of monomials of a(x)
    C_a_amon = mpvar('a',namon,1); % C_a_amon
    ax = C_a_amon'*amon; % a(x)
    
    % ------------------------------------------------------------------
    % Unit test: check max degrees
    % ------------------------------------------------------------------
    % max degree of common monomials
    for i1 = 1 : nx
        maxdeg1 = max(amon.degmat(:,i1));
        maxdeg2 = max(xmon.degmat(:,i1));
        if maxdeg1 > maxdeg2
            error('Max degree of a(x) is larger than the max degree of common monomials. Increase the degree of common monomials and try again.');
        end
    end
    
    % max degree of Psi
    for i1 = 1 : nx
        maxdeg1 = max(amon.degmat(:,i1));
        maxdeg2 = max(Psi.degmat(:,i1));
        if maxdeg1 > maxdeg2
            error('Max degree of a(x) is larger than the max degree of Psi. Increase the order of orthogonal basis vector and try again.');
        end
    end
    
    % ------------------------------------------------------------------
    % create coefficient vectors for a(x) in terms of M_x(x)
    % ------------------------------------------------------------------
    % arrange coefficient matrix of M_a(x) in terms of M_x(x)
    % such that M_a(x) = cmat_ax_xmon'*M_x(x)
    cmat_ax_xmon = zeros(nxmon,namon);
    [ixmon_a,iamon] = ismember(xmon.degmat, amon.degmat, 'rows');
    cmat_ax_xmon(find(ixmon_a==1),:) = amon.coefficient(nonzeros(iamon),:);
    
    % ------------------------------------------------------------------
    % unit test: check if M_a(x) = cmat_ax_xmon'*M_x(x)
    % ------------------------------------------------------------------
    err = amon - cmat_ax_xmon'*xmon;
    if norm(full(err.coefficient)) >= opts.errnrm
        error('Unit test failed: M_a(x) ~= transpose(cmat_ax_xmon)*xmon. Check your coefficient calculation and try again.');
    end
    
    % ------------------------------------------------------------------
    % create a new coefficient vector of a(x) in terms of M_x(x),
    % "C_a_xmon", such that a(x) = C_a_xmon'*M_x(x).
    % ------------------------------------------------------------------
    C_a_xmon = cmat_ax_xmon*C_a_amon;
    
    % ------------------------------------------------------------------
    % unit test: check if a = C_a_xmon'*M_x(x).
    % ------------------------------------------------------------------
    err = ax - C_a_xmon'*xmon;
    if norm(full(err.coefficient)) >= opts.errnrm
        error('Unit test failed: a(x) ~= transpose(C_a_xmon)*M_x(x).');
    end
    
    % ------------------------------------------------------------------
    % create a new coefficient vector of a(x) in terms of Psi(x),
    % "C_a_Psi", such that a(x) = C_a_Psi'*Psi;
    % ------------------------------------------------------------------
    % programmer's note: we convert C_a_xmon to symbolic before
    % computing inversion and back to polynomial since polynomial toolbox
    % can't handle "\" inverse. Also, if calculation involves inversion,
    % polynomial toolbox does not provide good numerical precision.
    % use conversion to symbolic when it comes to inversion calculation.
    C_a_Psi = s2p( C_Psi_xmon\p2s(C_a_xmon) );
%     C_a_Psi = pinv(C_Psi_xmon'*C_Psi_xmon)*C_Psi_xmon'*C_a_xmon;
    C_a_Psi.coefficient(find(abs(C_a_Psi.coefficient)<=opts.round)) = 0;
    
    % ------------------------------------------------------------------
    % unit test: check if a(x) = C_a_amon'*amon = C_a_Psi'*Psi
    % ------------------------------------------------------------------
    err = ax - C_a_Psi'*Psi;
    if norm(full(err.coefficient)) >= opts.errnrm
        error('Unit test failed: a(x) ~= transpose(C_a_Psi)*Psi.');
    end
end

%% check max degree of polynomial, b(x)
% max degree of common monomials
for i1 = 1 : nx
    maxdeg1 = max(bx.degmat(:,i1));
    maxdeg2 = max(xmon.degmat(:,i1));
    if maxdeg1 > maxdeg2
        error('Max degree of b(x) is larger than the max degree of common monomials. Increase the degree of common monomials and try again.');
    end
end

% max degree of Psi
for i1 = 1 : nx
    maxdeg1 = max(bx.degmat(:,i1));
    maxdeg2 = max(Psi.degmat(:,i1));
    if maxdeg1 > maxdeg2
        error('Max degree of b(x) is larger than the max degree of Psi. Increase the order of orthogonal basis vector and try again.');
    end
end

%% Create polynomial, c(x)
% ------------------------------------------------------------------
% check if c(x) is not a constant
% ------------------------------------------------------------------
if deg_c == 0
    error('Degree of c(x) cannot be a constant.');
end

% ------------------------------------------------------------------
% create monomial vector of c(x), "M_c(x)", and coefficient vector
% of c(x), "C_c_cmon", such that c(x) = C_c_cmon'*M_c(x).
% ------------------------------------------------------------------
deg_cmon = deg_c; % minimum degree : maximum degree of c(x)
cmon = monomials(x, deg_cmon); % M_c(x): monomials for c(x)
ncmon = length(cmon); % number of elements in M_c(x)
C_c_cmon = mpvar('c',ncmon,1); % C_c_cmon
cx = C_c_cmon'*cmon; % c(x) = C_c_cmon'*M_c(x).

% ------------------------------------------------------------------
% check max degrees
% ------------------------------------------------------------------
% max degree of common monomials
for i1 = 1 : nx
    maxdeg1 = max(cmon.degmat(:,i1));
    maxdeg2 = max(xmon.degmat(:,i1));
    if maxdeg1 > maxdeg2
        error('Max degree of c(x) is larger than the max degree of common monomials. Increase the degree of common monomials and try again.');
    end
end

% max degree of Psi
for i1 = 1 : nx
    maxdeg1 = max(cmon.degmat(:,i1));
    maxdeg2 = max(Psi.degmat(:,i1));
    if maxdeg1 > maxdeg2
        error('Max degree of c(x) is larger than the max degree of Psi. Increase the order of orthogonal basis vector and try again.');
    end
end

% ------------------------------------------------------------------
% create coefficient vectors for c(x) in terms of M_x(x),
% "C_c_xmon", such that c(x) = C_c_xmon'*M_x(x)
% ------------------------------------------------------------------
cmat_cx_xmon = zeros(nxmon,ncmon);
[ixmon_c,icmon] = ismember(xmon.degmat, cmon.degmat, 'rows');
cmat_cx_xmon(find(ixmon_c==1),:) = cmon.coefficient(nonzeros(icmon),:);
C_c_xmon = cmat_cx_xmon*C_c_cmon; % C_c_xmon

% ------------------------------------------------------------------
% unit test: check if c(x) = C_c_cmon'*cmon = C_c_xmon'*xmon
% ------------------------------------------------------------------
err = cx - C_c_xmon'*xmon;
if norm(full(err.coefficient)) >= opts.errnrm
    error('Unit test failed: c(x) ~= transpose(C_c_xmon)*xmon.');
end

% ------------------------------------------------------------------
% create a new coefficient vector of c(x) in terms of Psi(x)
% such that c(x) = C_c_Psi'*Psi;
% ------------------------------------------------------------------
% programmer's note: we convert C_a_xmon to symbolic before 
% computing inversion and back to polynomial since polynomial toolbox
% can't handle "\" inverse. Also, if calculation involves inversion,
% polynomial toolbox does not provide good numerical precision.
% use conversion to symbolic when it comes to inversion calculation.
C_c_Psi = s2p( C_Psi_xmon\p2s(C_c_xmon) );
% C_c_Psi = (C_Psi_xmon'*C_Psi_xmon)*C_Psi_xmon'*C_c_xmon;
C_c_Psi.coefficient(find(abs(C_c_Psi.coefficient)<=opts.round)) = 0;

% ------------------------------------------------------------------
% unit test: check if c(x) = C_c_cmon'*cmon = C_c_Psi'*Psi
% ------------------------------------------------------------------
err = cx - C_c_Psi'*Psi;
if norm(full(err.coefficient)) >= opts.errnrm
    error('Unit test failed: c(x) ~= transpose(C_c_Psi)*Psi.');
end

%% Create polynomial, s(x)
% ------------------------------------------------------------------
% check if s(x) is not a constant
% ------------------------------------------------------------------
if deg_s == 0
    error('Degree of s(x) cannot be a constant.');
end

% ------------------------------------------------------------------
% create monomial vector of s(x), "M_s(x)", and coefficient vector
% of s(x), "C_s_smon", such that s(x) = C_s_smon'*M_s(x).
% ------------------------------------------------------------------
deg_smon = deg_s; % minimum degree : maximum degree of s(x)
smon = monomials(x, deg_smon); % M_s(x): monomials for s(x)
nsmon = length(smon); % number of elements in M_s(x)
C_s_smon = mpvar('s',nsmon,1); % C_s_smon
sx = C_s_smon'*smon; % s(x) = C_s_smon'*M_s(x).

% ------------------------------------------------------------------
% check max degrees
% ------------------------------------------------------------------
% max degree of common monomials
for i1 = 1 : nx
    maxdeg1 = max(smon.degmat(:,i1));
    maxdeg2 = max(xmon.degmat(:,i1));
    if maxdeg1 > maxdeg2
        error('Max degree of s(x) is larger than the max degree of common monomials. Increase the degree of common monomials and try again.');
    end
end

% max degree of Psi
for i1 = 1 : nx
    maxdeg1 = max(smon.degmat(:,i1));
    maxdeg2 = max(Psi.degmat(:,i1));
    if maxdeg1 > maxdeg2
        error('Max degree of s(x) is larger than the max degree of Psi. Increase the order of orthogonal basis vector and try again.');
    end
end

% ------------------------------------------------------------------
% create coefficient vectors for s(x) in terms of M_x(x),
% "C_s_xmon", such that s(x) = C_s_xmon'*M_x(x)
% ------------------------------------------------------------------
cmat_sx_xmon = zeros(nxmon,nsmon);
[ixmon_s,ismon] = ismember(xmon.degmat, smon.degmat, 'rows');
cmat_sx_xmon(find(ixmon_s==1),:) = smon.coefficient(nonzeros(ismon),:);
C_s_xmon = cmat_sx_xmon*C_s_smon; % C_s_xmon

% ------------------------------------------------------------------
% unit test: check if s(x) = C_s_cmon'*smon = C_s_xmon'*xmon
% ------------------------------------------------------------------
err = sx - C_s_xmon'*xmon;
if norm(full(err.coefficient)) >= opts.errnrm
    error('Unit test failed: s(x) ~= transpose(C_s_xmon)*xmon.');
end

% ------------------------------------------------------------------
% create a new coefficient vector of s(x) in terms of Psi(x)
% such that s(x) = C_s_Psi'*Psi;
% ------------------------------------------------------------------
% programmer's note: we convert C_a_xmon to symbolic before 
% computing inversion and back to polynomial since polynomial toolbox
% can't handle "\" inverse. Also, if calculation involves inversion,
% polynomial toolbox does not provide good numerical precision.
% use conversion to symbolic when it comes to inversion calculation.
C_s_Psi = s2p( C_Psi_xmon\p2s(C_s_xmon) );
% C_s_Psi = (C_Psi_xmon'*C_Psi_xmon)*C_Psi_xmon'*C_s_xmon;
C_s_Psi.coefficient(find(abs(C_s_Psi.coefficient)<=opts.round)) = 0;

% ------------------------------------------------------------------
% unit test: check if s(x) = C_s_smon'*smon = C_s_Psi'*Psi
% ------------------------------------------------------------------
err = sx - C_s_Psi'*Psi;
if norm(full(err.coefficient)) >= opts.errnrm
    error('Unit test failed: s(x) ~= transpose(C_s_Psi)*Psi.');
end

%% create polynomial, z1(x) = a(x)*b(x)
ab = ax*bx; % z1(x) = a(x)*b(x)

% ------------------------------------------------------------------
% Coefficient variable vector of z1(x)
% ------------------------------------------------------------------
coef_vars = [];
for i1 = 1 : length(ab.varname)-nx
    coef_vars = [coef_vars; s2p(str2sym(ab.varname{i1}))];
end

% ------------------------------------------------------------------
% Find monomial vector of z1(x), M_ab(x), and coefficient vector,
% C_ab_abmon such that C_ab_abmon'*M_ab(x).
% ------------------------------------------------------------------
nabmon = size(ab.degmat,1); % number of monomials in M_ab(x)
tmp = repmat(x',nabmon,1).^ab.degmat(:,end-nx+1:end);
abmon = ones(nabmon,1);
for i1 = 1 : size(tmp,2)
    abmon = abmon.*tmp(:,i1); % M_ab(x)
end
C_ab_abmon = ab.coefficient.*ab.degmat(:,1:end-nx)*coef_vars; % C_ab_abmon

% ------------------------------------------------------------------
% unit test: check if a(x)*b(x) = C_ab_abmon'*abmon
% ------------------------------------------------------------------
if double(ab - (C_ab_abmon'*abmon)) ~= 0
    error('Unit test failed: a(x)*b(x) ~= transpose(C_ab_abmon)*abmon.');
end

% ------------------------------------------------------------------
% unit test: check max degrees of z1(x)
% ------------------------------------------------------------------
% max degree of common monomials
for i1 = 1 : nx
    maxdeg1 = max(abmon.degmat(:,i1));
    maxdeg2 = max(xmon.degmat(:,i1));
    if maxdeg1 > maxdeg2
        error('Max degree of a(x)*b(x) is larger than the max degree of common monomials. Increase the degree of common monomials and try again.');
    end
end

% max degree of Psi
for i1 = 1 : nx
    maxdeg1 = max(abmon.degmat(:,i1));
    maxdeg2 = max(Psi.degmat(:,i1));
    if maxdeg1 > maxdeg2
        error('Max degree of a(x)*b(x) is larger than the max degree of Psi. Increase the order of orthogonal basis vector and try again.');
    end
end

% ------------------------------------------------------------------
% Create coefficient vectors for z1(x) in terms of M_x(x) such that
% z1(x) = a(x)*b(x) = C_ab_xmon'*M_x(x).
% ------------------------------------------------------------------
% arrange coefficient matrix of M_ab(x) in terms of M_x(x)
cmat_ab_xmon = zeros(nxmon,nabmon);
[ixmon_ab,iabmon] = ismember(xmon.degmat, abmon.degmat, 'rows');
cmat_ab_xmon(find(ixmon_ab==1),:) = abmon.coefficient(nonzeros(iabmon),:);

% create a new coefficient vector of z1(x) in terms of M_x(x)
% such that z1(x) = C_ab_xmon'*xmon.
C_ab_xmon = cmat_ab_xmon*C_ab_abmon;

% ------------------------------------------------------------------
% unit test: check if z1(x) = C_ab_abmon'*abmon = C_ab_xmon'*xmon
% ------------------------------------------------------------------
err = ab - C_ab_xmon'*xmon;
if norm(full(err.coefficient)) >= opts.errnrm
    error('Unit test failed: a(x)*b(x) ~= transpose(C_ab_xmon)*xmon.');
end

% ------------------------------------------------------------------
% create a new coefficient vector of z1(x) in terms of Psi(x)
% such that z1(x) = C_ab_Psi'*Psi;
% ------------------------------------------------------------------
C_ab_Psi = s2p( C_Psi_xmon\p2s(C_ab_xmon) );
% C_ab_Psi = pinv(C_Psi_xmon'*C_Psi_xmon)*C_Psi_xmon'*C_ab_xmon;
C_ab_Psi.coefficient(find(abs(C_ab_Psi.coefficient)<=opts.round)) = 0;

% ------------------------------------------------------------------
% unit test: check if z1(x) = C_ab_abmon'*abmon = C_an_Psi'*Psi
% ------------------------------------------------------------------
err = ab - C_ab_Psi'*Psi;
if norm(full(err.coefficient)) >= opts.errnrm
    error('Unit test failed: a(x)*b(x) ~= transpose(C_ab_Psi)*Psi.');
end

%% create polynomial, z2(x) = b(x)*c(x)
% ------------------------------------------------------------------
% create polynomial, z2(x) = b(x)*c(x)
% ------------------------------------------------------------------
bc = bx*cx;

% ------------------------------------------------------------------
% Coefficient variable vector of z2(x)
% ------------------------------------------------------------------
coef_vars = [];
for i1 = 1 : length(bc.varname)-nx
    coef_vars = [coef_vars; s2p(str2sym(bc.varname{i1}))];
end

% ------------------------------------------------------------------
% Find monomial vector of z2(x), bcmon
% ------------------------------------------------------------------
nbcmon = size(bc.degmat,1);
tmp = repmat(x',nbcmon,1).^bc.degmat(:,end-nx+1:end);
bcmon = ones(nbcmon,1);
for i1 = 1 : size(tmp,2)
    bcmon = bcmon.*tmp(:,i1);
end

% ------------------------------------------------------------------
% Find coefficients of z2(x) such that z2(x) = C_bc_bcmon'*bcmon
% ------------------------------------------------------------------
C_bc_bcmon = bc.coefficient.*bc.degmat(:,1:end-nx)*coef_vars;

% unit test: check if b(x)*c(x) = C_bc_bcmon'*bcmon;
if double(bc - (C_bc_bcmon'*bcmon)) ~= 0
    error('Unit test failed: b(x)*c(x) ~= transpose(C_bc_bcmon)*bcmon. Check your coefficient calculation and try again.');
end

% ------------------------------------------------------------------
% check max degrees
% ------------------------------------------------------------------
% max degree of common monomials
for i1 = 1 : nx
    maxdeg1 = max(bcmon.degmat(:,i1));
    maxdeg2 = max(xmon.degmat(:,i1));
    if maxdeg1 > maxdeg2
        error('Max degree of b(x)*c(x) is larger than the max degree of common monomials. Increase the degree of common monomials and try again.');
    end
end

% max degree of Psi
for i1 = 1 : nx
    maxdeg1 = max(bcmon.degmat(:,i1));
    maxdeg2 = max(Psi.degmat(:,i1));
    if maxdeg1 > maxdeg2
        error('Max degree of b(x)*c(x) is larger than the max degree of Psi. Increase the order of orthogonal basis vector and try again.');
    end
end

% ------------------------------------------------------------------
% Create coefficient vectors for z2(x) in terms of M_x(x) and Psi(x)
% ------------------------------------------------------------------
% initialize coefficient matrix of z2(x) in terms of M_x(x)
cmat_bc_xmon = zeros(nxmon,nbcmon);

% arrange coefficient matrix of M_bc(x) in terms of M_x(x)
[ixmon_bc,ibcmon] = ismember(xmon.degmat, bcmon.degmat, 'rows');
cmat_bc_xmon(find(ixmon_bc==1),:) = bcmon.coefficient(nonzeros(ibcmon),:);

% create a new coefficient vector of z2(x) in terms of M_x(x)
% such that z2(x) = C_bc_xmon'*xmon.
C_bc_xmon = cmat_bc_xmon*C_bc_bcmon;

% create a new coefficient vector of z2(x) in terms of Psi(x)
% such that z2(x) = C_bc_Psi'*Psi;
% C_bc_Psi = (C_Psi_xmon'*C_Psi_xmon)\C_Psi_xmon'*C_bc_xmon;
C_bc_Psi = s2p( C_Psi_xmon\p2s(C_bc_xmon) );
C_bc_Psi.coefficient(find(abs(C_bc_Psi.coefficient)<=opts.round)) = 0;

% ------------------------------------------------------------------
% unit test: check if z2(x) = C_bc_bcmon'*bcmon = C_bc_xmon'*xmon
% ------------------------------------------------------------------
err = bc - C_bc_xmon'*xmon;
if norm(full(err.coefficient)) >= opts.errnrm
    error('Unit test failed: b(x)*c(x) ~= transpose(C_bc_xmon)*xmon.');
end

% ------------------------------------------------------------------
% unit test 3: check if z2(x) = C_bc_bcmon'*abmon = C_bc_Psi'*Psi
% ------------------------------------------------------------------
err = bc - C_bc_Psi'*Psi;
if norm(full(err.coefficient)) >= opts.errnrm
    error('Unit test failed: b(x)*c(x) ~= transpose(C_bc_Psi)*Psi.');
end

%% Find coefficient vector, C_x_Psi, such that x = C_x_Psi'*Psi
% ------------------------------------------------------------------
% arrange coefficient matrix of a state vector x in terms of M_x(x)
% such that x = cmat_x_xmon'*M_x(x)
% ------------------------------------------------------------------
cmat_x_xmon = zeros(nxmon,nx);
[~,ix2xmon] = ismember(xmon.varname, x.varname);
[ixmon_x,ix] = ismember(xmon.degmat, x.degmat(:,ix2xmon), 'rows');
cmat_x_xmon(find(ixmon_x==1),:) = x.coefficient(nonzeros(ix),:);

% ------------------------------------------------------------------
% unit test: check if x = cmat_x_xmon'*M_x(x)
% ------------------------------------------------------------------
err = x - cmat_x_xmon'*xmon;
if norm(full(err.coefficient)) >= opts.errnrm
    error('Unit test failed: x ~= transpose(cmat_x_xmon)*xmon.');
end

% ------------------------------------------------------------------
% create a new coefficient vector of the state vector x 
% in terms of Psi(x) such that x = C_x_Psi'*Psi;
% ------------------------------------------------------------------
C_x_Psi = C_Psi_xmon\cmat_x_xmon;
C_x_Psi(find(C_x_Psi <= opts.round)) = 0;

% ------------------------------------------------------------------
% unit test: check if x = C_x_Psi'*Psi(x)
% ------------------------------------------------------------------
err = x - C_x_Psi'*Psi;
if norm(full(err.coefficient)) >= opts.errnrm
    error('Unit test failed: x ~= transpose(C_x_Psi)*Psi.');
end

%% Find derivative of basis
dPsi = jacobian(Psi,x)*dx;