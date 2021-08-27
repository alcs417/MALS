function [idx_eg, idx_rc] = MALS_wrapper(X, lambda_R, lambda_A)

    coef.lambR = 1;
    coef.lambQ = lambda_R;
    coef.lambA = lambda_A; 
    coef.lambE = 1;
    coef.k = 9;       % k nearest neighbors
    coef.islocal = 1;

    A = MALS(X, coef);
    NUMC = 2 : 15;
    [K1, ~, K12, ~] = Estimate_Number_of_Clusters_given_graph(A, NUMC);
    fprintf('The number of clusters estimated by eigen gap : %d\n', K1);
    idx_eg = computeLabels(A, K1);
    fprintf('The number of clusters estimated by rotation cost : %d\n', K12);
    idx_rc = computeLabels(A, K12);

end

function idx = computeLabels(T, k)
    Z = (abs(T)+ abs(T')) / 2;
    idx = clu_ncut(Z, k);
end
