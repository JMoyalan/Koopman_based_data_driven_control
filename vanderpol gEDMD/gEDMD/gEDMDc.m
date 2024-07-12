function [op, normPsi, normdPsi, cp] = gEDMDc(X, dX, U, Psi, dPsi, opt)
% find center points of RBFs if the basis includes RBFs
if opt.tpsrbf_enable == 1
    fprintf(sprintf('[%s] Started calculating RBF center points using k-means clustering\n', datetime('now')));
    nc = opt.tpsrbf_nc;
    nx = size(X,1);
    
    % Compute RBF center points
    if isempty(opt.tpsrbf_ctps)
        if strcmp(opt.type,'edmd')
            data = unique([X'; dX'],'row','stable');
        else
            data = X';
        end
        [~, cp] = kmeans(data, nc, 'MaxIter', 10000);
        
        % Remove overlapping center points
        cnt = 1;
        ndup = 0;
        max_iter = 100;
        while 1
            if strcmp(opt.type,'edmd')
                iovepts = double([ismember(cp,X','row'), ismember(cp,dX','row')]);
                idup = find(sum(iovepts, 2) >= 1);
            else
                iovepts = ismember(cp,X','row');
                idup = find(iovepts >= 1);
            end
            
            if isempty(idup)
                break;
            end
            
            cp(idup, :) = cp(idup, :) + 1e-1*randn(length(idup), nx);
            ndup = ndup + length(idup);
            fprintf(sprintf('[%s] Removed %d overlapping data point\n', datetime('now'), length(idup)));
            
            cnt = cnt + 1;
            if cnt > max_iter
                fprintf(sprintf('[%s] Stopped the loop of removing overlapping points since it reached the maximum %d iteration number\n', datetime('now'), max_iter));
                break;
            end
        end
        
        ctps = cp;
        
        fprintf(sprintf('[%s] Found RBF basis center points\n', datetime('now')));
    else
        ctps = opt.tpsrbf_ctps;
        fprintf(sprintf('[%s] Loaded existing RBF basis center points\n', datetime('now')));
    end    
else
    ctps = [];
end

if opt.grbf_enable == 1
    fprintf(sprintf('[%s] Started calculating RBF center points using k-means clustering\n', datetime('now')));
    nc = opt.grbf_nc;
    nx = size(X,1);
    
    % Compute RBF center points
    if isempty(opt.grbf_cgrbf)
        if strcmp(opt.type,'edmd')
            data = unique([X'; dX'],'row','stable');
        else
            data = X';
        end
        [~, cp] = kmeans(data, nc, 'MaxIter', 10000);
        
        % Remove overlapping center points
        cnt = 1;
        ndup = 0;
        max_iter = 100;
        while 1
            if strcmp(opt.type,'edmd')
                iovepts = double([ismember(cp,X','row'), ismember(cp,dX','row')]);
                idup = find(sum(iovepts, 2) >= 1);
            else
                iovepts = ismember(cp,X','row');
                idup = find(iovepts >= 1);
            end
            
            if isempty(idup)
                break;
            end
            
            cp(idup, :) = cp(idup, :) + 1e-1*randn(length(idup), nx);
            ndup = ndup + length(idup);
            fprintf(sprintf('[%s] Removed %d overlapping data point\n', datetime('now'), length(idup)));
            
            cnt = cnt + 1;
            if cnt > max_iter
                fprintf(sprintf('[%s] Stopped the loop of removing overlapping points since it reached the maximum %d iteration number\n', datetime('now'), max_iter));
                break;
            end
        end
        
        cgrbf = cp;
        
        fprintf(sprintf('[%s] Found RBF basis center points\n', datetime('now')));
    else
        cgrbf = opt.tpsrbf_cgrbf;
        fprintf(sprintf('[%s] Loaded existing RBF basis center points\n', datetime('now')));
    end
else
    cgrbf = [];
end

% combine center points
cp = [ctps(:); cgrbf(:)];

if opt.method == 1 % using parallel cellfun using batch data
    % parallel batch preparation
    M = size(X, 2);
    nx = size(X, 1);
    Xi = round(M/opt.batch);
    cols = [repmat(Xi,1,opt.batch-1), M - (opt.batch-1)*Xi];
    if strcmp(opt.type,'gedmd')
        if opt.tpsrbf_enable == 1
            cpmat = repmat(cp,1,M);
            ndX = nx+size(cpmat,1);
            ndZ = 2*nx+size(cpmat,1);
            data_X = mat2cell([X;cpmat], ndX, cols);
            data_Z = mat2cell([X;dX;cpmat], ndZ, cols);            
        else
            ndX = nx;
            ndZ = 2*nx;
            data_X = mat2cell(X, ndX, cols);
            data_Z = mat2cell([X;dX], ndZ, cols);
        end
    elseif strcmp(opt.type,'edmd')
        if opt.tpsrbf_enable == 1
            cpmat = repmat(cp,1,M);
            ndX = nx+size(cpmat,1);
            ndZ = nx+size(cpmat,1);
            data_X = mat2cell([X;cpmat], ndX, cols);
            data_Z = mat2cell([dX;cpmat], ndZ, cols);
        else
            ndX = nx;
            ndZ = nx;
            data_X = mat2cell(X, ndX, cols);
            data_Z = mat2cell(dX, ndZ, cols);
        end
    end
    
    % calculate values of basis functions for each sub-data
    hbar = parfor_progressbar(opt.batch, 'Computing Koopman operator..');
    parfor i1 = 1 : opt.batch
        % sub data length
        Mi = size(data_X{i1}, 2);
        
        % convert each column of X and Z to cells        
        Xcell = mat2cell(data_X{i1}, ndX, ones(Mi,1));
        Zcell = mat2cell(data_Z{i1}, ndZ, ones(Mi,1));        
        
        % calculate Psi(X) and dPsi(X)
        sub_Psi_X{i1} = cell2mat( cellfun(Psi, Xcell, 'UniformOutput', false) );
        sub_dPsi_X{i1} = cell2mat( cellfun(dPsi, Zcell, 'UniformOutput', false) );
        
        % data matrix reordered in parallel computing
        sub_X{i1} = data_X{i1};
        sub_Z{i1} = data_Z{i1};
        
        hbar.iterate(1);
    end
    close(hbar);
    
    % concatenate all sub matrices
    Psi_X = double(cat(2, sub_Psi_X{:}));
    dPsi_X = double(cat(2, sub_dPsi_X{:}));
    
    % augmented observables
    if ~isempty(U)
        Z_aug = [Psi_X; U];
    else
        Z_aug = Psi_X;
    end
    
    % concaternate all sub data
%     Xr = cat(2, sub_X{:});
%     Zr = cat(2, sub_Z{:});
%     dXr = Zr(nx+1:end,:);
    
    % normalize rows (correspoding to time) of Psi(X)
    if opt.normal.X == 1
        normPsi = vecnorm(Z_aug,1,2);
        fprintf(sprintf('[%s] Completed normalization of X\n', datetime('now')));
    else
        normPsi = ones(size(Z_aug,1),1);
    end
    Z_aug_n = Z_aug./normPsi;
    
    % normalize rows (correspoding to time) of dPsi(X)
    if opt.normal.dX == 1
        normdPsi = vecnorm(dPsi_X,1,2);
        normdPsi(find(normdPsi == 0)) = 1;
        fprintf(sprintf('[%s] Completed normalization of dX\n', datetime('now')));
    else
        normdPsi = ones(size(dPsi_X,1),1);
    end
    dPsi_X_n = dPsi_X./normdPsi;

    % calculate Koopman generator matrix (initial guess)
    Z_aug_n_T = Z_aug_n'; % transpose of normalized Psi(X)
    dPsi_X_n_T = dPsi_X_n'; % transpose of normalized dPsi(X)

    % solution 1
    G = 1/M*Z_aug_n*Z_aug_n_T;
    A = 1/M*Z_aug_n*dPsi_X_n_T;
    L = pinv(G)*A;
%     L = G\A;
    
    % solution 2
%     L = Psi_X_n_T\dPsi_X_n_T;

    fprintf(sprintf('[%s] Found initial solution with error = %.12f\n', datetime('now'), norm(dPsi_X_n_T - Z_aug_n_T*L,'fro')));
    
    % iterative sparsification
    if opt.sparse == 1
        lambda = 0.05;
        for k = 1 : 10
            smallinds = (abs(L)<lambda);    % find small coefficients
            L(smallinds)=0;                 % and threshold
            for ind = 1:nx                   % n is state dimension
                biginds = ~smallinds(:,ind);
                % Regress dynamics onto remaining terms to find sparse Xi
                L(biginds,ind) = Psi_X_n_T(:,biginds)\dPsi_X_n_T(:,ind);
            end
            fprintf(sprintf('[%s] Found sparse solution at %d iteration with error = %.12f\n', datetime('now'), k, norm(dPsi_X_n_T - Psi_X_n_T*L,'fro')));
        end
    end
    
    % sparse solution
    op = L'./normPsi';
    fprintf(sprintf('[%s] Found final solution with error = %.12f\n', datetime('now'), norm(dPsi_X./normdPsi - op*Z_aug,'fro')));
    
    op = op.*normdPsi; % reverse normalization of dPsi(x)
    
    % Draw sparsity map
    figure;
    spy(op);
    title('Sparsity map of the estimated Koopman generator');
    
else % using parallel for loop
     % This needs to be modified to add normalization procedure !!!
     % This needs to be modified to reflect the changes in method 1 !!!
    G = 0; A = 0; M = size(X,2);
    hbar = parfor_progressbar(M, 'Computing Koopman operator..');
    for i1 = 1 : M
        % Lifting data points
        if strcmp(opt.type,'gedmd')
            if opt.tpsrbf_enable == 1
                Psi_X = Psi([X(:,i1); cp]);
                Psi_X = [Psi_X; U(:,i1)];
                dPsi_X = dPsi([X(:,i1); dX(:,i1); cp]);
            else
                Psi_X = Psi(X(:,i1));
                Psi_X = [Psi_X; U(:,i1)];
                dPsi_X = dPsi([X(:,i1); dX(:,i1)]);
            end
        elseif strcmp(opt.type,'edmd')
            if opt.tpsrbf_enable == 1
                Psi_X = Psi([X(:,i1); cp]);
                dPsi_X = dPsi([dX(:,i1); cp]);
            else
                Psi_X = Psi(X(:,i1));
                dPsi_X = dPsi(dX(:,i1));
            end
        end
        G{i1} = 1/M*Psi_X*Psi_X';
        A{i1} = 1/M*dPsi_X*Psi_X';
        hbar.iterate(1);
    end
    close(hbar);
    
    % compute G and A
    G = sum(cat(3, G{:}),3);
    A = sum(cat(3, A{:}),3);
    
    % Estimate Koopman operator K
    op = A*pinv(G);
    
    fprintf(sprintf('[%s] Found solution with error = %.12f\n', datetime('now'), norm(dPsi_X - op*Psi_X,'fro')));
    
    % data matrix
    Xr = X;
    dXr = dX;
    
    % norms
    normPsi = [];
    normdPsi = [];
end