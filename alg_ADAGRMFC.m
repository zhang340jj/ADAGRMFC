function y3=alg_ADAGRMFC(Y,Sd,St,cv_setting,nr_fold,left_out,use_WKNKN,K,eta,use_W_matrix)
%alg_grmf predicts DTIs based on the ADA-GRMFC algorithm described in the following paper: 
% Junjun Zhang and Minzhu Xie
% (2022) Graph regularized non-negative matrix factorization with prior knowledge consistency 
%constraint for drug-target interaction prediction
% INPUT:
%  Y:           interaction matrix
%  Sd:          pairwise drug similarities matrix
%  St:          pairwise target similarities matrix
%  cv_setting:  cross validation setting ('cv_d'or 'cv_t')
%  nr_fold:     number of folds in cross validation experiment
%  left_out:    if cv_setting=='cv_d' --> left_out is 'drug' indices that are left out
%               if cv_setting=='cv_t' --> left_out is 'target' indices that are left out
%              
%
% OUTPUT:
%  y3:  prediction matrix

    % get best parameters
    [k,lambda_l,lambda_d,lambda_t,p,num_iter] = alg_ADAGRMFC_parest(Y,Sd,St,cv_setting,nr_fold,left_out,use_WKNKN,K,eta,use_W_matrix);
    fprintf('k%g\t\t%g\t%g\t%g\t\t',k,lambda_l,lambda_d,lambda_t);

    % preprocessing Y
    if use_WKNKN
        Y = preprocess_WKNKN(Y,Sd,St,K,eta);
    end

    % preprocessing Sd & St
    Sd = preprocess_PNN(Sd,p);
    St = preprocess_PNN(St,p);

    % Laplacian Matrices
    Dd = diag(sum(Sd));
    Dt = diag(sum(St));
    Ld = Dd - Sd;
    Ld = (Dd^(-0.5))*Ld*(Dd^(-0.5));
    Lt = Dt - St;
    Lt = (Dt^(-0.5))*Lt*(Dt^(-0.5));

    % initialize A & B
    [A,B] = initializer(Y,k);


    % if weight matrix W is not being used (i.e. classifier='grmf')
    if ~use_W_matrix
        [A,B] = alg_ADAGRMFC_predict(Y,A,B,Ld,Lt,lambda_l,lambda_d,lambda_t,num_iter);

    % else if weight matrix W is being used (i.e. classifier='wgrmf')
    else
        test_ind = get_test_indices(Y,cv_setting,left_out);
        W = ones(size(Y));
        W(test_ind) = 0;
        [A,B] = alg_ADAGRMFC_predict(Y,A,B,Ld,Lt,lambda_l,lambda_d,lambda_t,num_iter,W);
    end

    % compute prediction matrix
    y3 = A*B';

    %--------------------------------------------------------------------

end