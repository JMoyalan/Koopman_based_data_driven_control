function [K, G, A, Xr, Yr, Psi_X, Psi_Y] = EDMD_poly(X, Y, xvar, basis, opt)
if opt.method == 1 % using parallel cellfun using batch data
    % parallel batch preparation
    M = size(X, 2);
    nx = size(X, 1);
    dX = round(M/opt.batch);
    cols = [repmat(dX,1,opt.batch-1), M - (opt.batch-1)*dX];
    data_X = mat2cell(X, nx, cols);
    data_Y = mat2cell(Y, nx, cols);
    
    % calculate values of basis functions for each sub-data
    hbar = parfor_progressbar(opt.batch, 'Computing Koopman operator..');
    parfor i1 = 1 : opt.batch
        % sub data length
        Mi = size(data_X{i1}, 2);
        
        % convert each column of X and Y to cells
        Xcell = mat2cell(data_X{i1}, nx, ones(Mi,1));
        Ycell = mat2cell(data_Y{i1}, nx, ones(Mi,1));
        
        % calculate Psi(X) and Psi(Y)
        sub_Psi_X{i1} = cell2mat( ...
            cellfun(@(x) double(subs(basis, xvar, x)), ...
            Xcell, 'UniformOutput', false) );
        sub_Psi_Y{i1} = cell2mat( ...
            cellfun(@(x) double(subs(basis, xvar, x)), ...
            Ycell, 'UniformOutput', false) );
        
        % data matrix reordered in parallel computing
        sub_X{i1} = data_X{i1};
        sub_Y{i1} = data_Y{i1};
        
        hbar.iterate(1);
    end
    close(hbar);
    
    % concatenate all sub matrices
    Psi_X = cat(2, sub_Psi_X{:});
    Psi_Y = cat(2, sub_Psi_Y{:});
    
    % concaternate all sub data
    Xr = cat(2, sub_X{:});
    Yr = cat(2, sub_Y{:});
    
    % calculate Koopman operator matrix
    G = 1/M*( Psi_X*Psi_X' );
    A = 1/M*( Psi_X*Psi_Y' );
    K = pinv(G)*A;
    
else % using parallel for loop
    G = 0; A = 0; M = size(X,2);
    hbar = parfor_progressbar(M, 'Computing Koopman operator..');
    for i1 = 1 : M
        Psi_X(:,i1) = double(subs(basis,xvar,X(:,i1)));
        Psi_Y(:,i1) = double(subs(basis,xvar,Y(:,i1)));
        G = G + Psi_X(:,i1)*Psi_X(:,i1)'; % compact form of the lifted data points for X
        A = A + Psi_X(:,i1)*Psi_Y(:,i1)'; % compact form of the lifted data points for Y
        hbar.iterate(1);
    end
    close(hbar);
    
    % Estimate Koopman operator K
    G = 1/M*G;
    A = 1/M*A;
    K = pinv(G)*A;
    
    % data matrix
    Xr = X;
    Yr = Y;
end