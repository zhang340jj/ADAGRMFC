function [A,B]=alg_srcmf_predict(Y,A,B,Sd,St,lambda_l,lambda_d,lambda_t,num_iter,W)
%alg_cmf_predict is a helper function of CMF that executes the alternating
%least squares method to obtain the solution. 
    
    K = size(A,2);
    lambda_d_Sd = lambda_d*Sd;
    lambda_t_St = lambda_t*St;
    lambda_l_eye_K = lambda_l*eye(K);
   
    
    AA = initializer(Sd,K);
    BB  = initializer(St,K);
    
    
    if nargin < 10
        AtA = A'*A;
        BtB = B'*B;
        for z=1:num_iter
            A = (Y*B + lambda_d*AA)  / (BtB +( lambda_l+lambda_d)*eye(K));
            AtA = A'*A;
            B = (Y'*A + lambda_t*BB) / (AtA +( lambda_l+lambda_t)*eye(K));
            BtB = B'*B;
            
             %update AA
        AA = (2*Sd*AA+A)/(2*(AA')*AA+eye(K));
        
        %update BB
        BB = (2*St*BB+B)/(2*(BB')*BB+eye(K));
        end
        
    else
        H = W .* Y;
        for z=1:num_iter
            B_old = B;
            HB_plus_lambda_d_Sd_A_old = H*B + lambda_d*AA;
            lambda_l_eye_k_plus_lambda_d_A_oldt_A_old = ( lambda_l+lambda_d)*eye(K);
            for a=1:size(A,1)
                A(a,:) = HB_plus_lambda_d_Sd_A_old(a,:) / (lambda_l_eye_k_plus_lambda_d_A_oldt_A_old + B'*diag(W(a,:))*B);
            end
            A_old = A;
            HtA_plus_lambda_t_St_B_old = H'*A + lambda_t*BB;
            lambda_l_eye_k_plus_lambda_t_B_oldt_B_old = ( lambda_l+lambda_t)*eye(K);
            for b=1:size(B,1)
                B(b,:) = HtA_plus_lambda_t_St_B_old(b,:) /(lambda_l_eye_k_plus_lambda_t_B_oldt_B_old + A'*diag(W(:,b))*A);
            end

                         %update AA
        AA = (2*Sd*AA+A)/(2*(AA')*AA+eye(K));
        
        %update BB
        BB = (2*St*BB+B)/(2*(BB')*BB+eye(K));
            
%             % for readability...
%             A_old = A;
%             lambda_d_A_oldt_A_old = lambda_d*(A_old'*A_old);
%             for a=1:size(A,1)
%                 A(a,:) = (H(a,:)*B + lambda_d_Sd(a,:)*A_old) / (B'*B + lambda_l_eye_k + lambda_d_A_oldt_A_old);
%             end
%             B_old = B;
%             lambda_t_B_oldt_B_old = lambda_t*(B_old'*B_old);
%             for b=1:size(B,1)
%                 B(b,:) = (H(:,b)'*A + lambda_t_St(b,:)*B_old) / (A'*A + lambda_l_eye_k + lambda_t_B_oldt_B_old);
%             end
        end
    end
    
end