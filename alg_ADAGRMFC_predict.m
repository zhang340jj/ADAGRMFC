function [A,B]=alg_ADAGRMFC_predict(Y,A,B,Ld,Lt,lambda_l,lambda_d,lambda_t,num_iter,W)
%alg_ADAGRMFC_predict is a helper function of ADAGRMFC that executes the
%alternating  direction algorithm to obtain the solution.
    
    K = size(A,2);
    lambda_d_Ld = lambda_d*Ld;
    lambda_t_Lt = lambda_t*Lt;
    lambda_l_eye_K = lambda_l*eye(K);
 
 if nargin < 10
    
    opts.tol = 1e-5;
    opts.maxit = 500;
    opts.print = 1;
    
    
    U1 = zeros(size(A));
    V1 = zeros(size(B));
    Lam = zeros(size(A));
    pai = zeros(size(B));
%判断投影是否为满秩
    Pomega=Y-A*B';

    m1=445;
    n1=664;

    beta1=0.01;
    af=0.4;
    gamma=1.618;



L = length(Pomega);
 Zfull = (L/(m1*n1) > 0.2 ) || m1*n1 < 5e5;

  
for iter = 1:num_iter
      
    % updating variables A, U1 and Lam
gX = lambda_l*U1 - Lam-lambda_d_Ld*A;
if Zfull    
       A=(Y*B+gX)/(B'*B + lambda_l_eye_K);
else       
        A=(B+Y*(B*B')+gX)/(B'*B + lambda_l_eye_K);
    end

        U1= max(0, A);
           
    % updating variables B, V1 and pai
    
%
    if Zfull

     B =(Y'*A+lambda_l*V1 - pai -lambda_t_Lt*B) / (A'*A + lambda_l_eye_K);
    else

     B = (Y'*A+Y'*(A*A')+lambda_l*V1 - pai -lambda_t_Lt*B) / (A'*A + lambda_l_eye_K);
    end
    
    V1= max(0, B + pai/beta1);    
    pai = pai + gamma*beta1*(B-V1);
    
    % updating variable Y
    if Zfull   
        Y = A*B'; 
        Res = Y-A*B'; 
    else
        Y=A*B'+Res;
    end 

   


  
   
        end
    end
    
end










