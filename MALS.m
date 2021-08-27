function A = MALS(X, coef)
% X:                cell array, 1 by view_num, each array is d_v by num

%    v = max(size(X, 2), size(X, 1));
		v = length(X);
    num = size(X{1}, 2);
    NITER = 30;
    
    %% ==================== Initialization =====================
	rho = 1.2;
    miu = 1.25;
    threshold = 10^-7;
    
    E = cell(v, 1);
    R = cell(v, 1);
    Q = cell(v, 1);
    Z1 = cell(v, 1);
    Z2 = cell(v, 1);
    Wv = ones(v, 1) / v;
   
    for i = 1 : v
       E{i} = zeros(size(X{i}));
       Z1{i} = zeros(size(X{i}));
       
       R{i} = zeros(num);
       Q{i} = zeros(num);
       Z2{i} = zeros(num);       
    end    
    
    lambQ = coef.lambQ; % for Tr(QLQ^t); 
    lambA = coef.lambA; % for ||A||^2;
    lambR = 1;
    lambE = 1; % for ||E||-norm;
    lambWv = 2;
    k = coef.k;
    islocal = coef.islocal;
    
    A = zeros(num);
    
    %% ===================== Normalization =====================
    % pre-caculate X^tX
    XTX = cell(v, 1);
    
    for i = 1 :v
        X{i} = NormalizeFea(X{i}, 0); % unit norm
        XTX{i} = X{i}' * X{i}; 
    end
    
    
    %% ======================= Initilization ====================
    %initialize weighted_distX
    SUM = zeros(num);
    distX_initial = zeros(num, num, v);
    for i = 1:v
        distX_initial(:, :, i) =  L2_distance_1(X{i}, X{i});                  %initialize X
        SUM = SUM + distX_initial(:, :, i);
    end
    distX = 1 / v * SUM;
    [distXs, idx] = sort(distX,2);

    %initialize S
    S = zeros(num);
%     rr = zeros(num,1);
    for i = 1:num
        di = distXs(i, 2 : k + 2);
%         rr(i) = 0.5 * (k * di(k + 1) - sum(di(1 : k)));
        id = idx(i, 2 : k + 2);
        S(i, id) = (di(k + 1) - di) / ( k * di(k + 1) - sum(di(1 : k)) + eps);               %initialize S
        
    end
    
    
    %% ======================= Iteration =======================
    for iter = 1 : NITER
%         fprintf('Starting the %d-th iteration...\n', iter);
        
        % construct LA
        D = diag(sum(A));
        L = D - A;
        
        for iv = 1 : v
            % update E_v
            temp1 = X{iv} - X{iv} * R{iv} + Z1{iv} / miu;
            temp2 = lambE / miu; % lamdba2 = 1;
            % |E|^1-norm
            E{i} = max(0, temp1 - temp2) + min(0, temp1 + temp2);
            clear temp1 temp2;       

            % update R_v
            R{iv} = solve_nuclear_norm(X{iv}, R{iv}, Z1{iv}, Z2{iv}, miu, lambR, E{iv}, Q{iv});

            % update Q_v
            Q{iv} = (miu * R{iv} - Z2{iv}) / (4 * lambQ * Wv(iv) * L + miu * eye(num));         
            
            % update multipliers
            Z1{iv} = Z1{iv} + miu * (X{iv} - X{iv} * R{iv} - E{iv});
            Z2{iv} = Z2{iv} + miu * (Q{iv} - R{iv});
                     
        end        
        
        % update A
        % update weighted distQ
        SUM = zeros(num);
        distQ = zeros(num, num, v);
        for i = 1 : v
            distQ(:, :, i) =  (Wv(i)^lambWv) * L2_distance_1(Q{i}, Q{i});
            SUM = SUM + distQ(:, :, i);
        end    
        
        A = zeros(num);
        for i = 1:num
            if islocal == 1
                idxa0 = idx(i, 2 : k + 1);
            else
                idxa0 = 1:num;
            end;
            dqi = SUM(i, idxa0);
            ad = -(0.5 * lambQ * dqi) ./ (2 * lambA); % lamb3 or alpha
            A(i,idxa0) = EProjSimplex_new(ad);
        end;
        A = (A + A') / 2;
        
        % update W_v
        for i = 1 : v
            Wv(i) = 0.5 / sqrt(sum(sum(distQ(:, :, i) .* A)));      
        end

        diff1 = 0;
        diff2 = 0;
        for i = 1 : v
            diff1 = diff1 + max(sum(abs(X{i} - X{i} * R{i} - E{i}), 1));
            diff2 = diff2 + max(sum(abs(Q{i} - R{i}), 1));
        end
        
        if diff2 <= threshold && diff1 <= threshold
            fprintf('Algorithm has reached convergence at %d-th iteration', iter);
            break;
        end        

        miu = miu * rho;
    end

end