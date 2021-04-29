function  [k, prop, ym, yvar, beta0,lkhd,step] = ep_EM_bic(y, beta, i)
% EM algorithm for ep Mixture Model with mulitivarite form through BIC
% 
% y -- data n*1
% mmax -- maximum no. of mixtures
% 
% k -- no. of mixtures
% prop -- mixture proportions 
% ym -- mean vectors
% yv -- sigma2
% beta -- 1*k
%
%
k = i;
[n, d] = size(y);
 beta0=beta';

% K-means Initialization
[prop0, ymean0, var0] = gmm_kmean(y, i);
yvar0=var0(:);
% 
if (sum(prop0==0) > 0) 
   ind = (prop0>0);
   k = sum(ind);
   prop0= prop0(ind);
   ymean0 = ymean0(ind, :);
   yvar0 = yvar0(ind);
   beta0=beta0(ind);
end;
%  ymean0=zeros(k,1)+randn(k,1);
yvar0=yvar0.*gamma(0.5./beta0)./(2.^(1./beta0).*gamma(1.5./beta0));
 
epsilon = 1e-7;
step = 0;
delta = 1;
                %lambda = 1;
Dof = 1 + d + d*(d+1)/2; 
lkhd0 = 0;
%lkhd0_p = -10000;
options = optimoptions('fsolve','Display','none');


while (delta >= epsilon)&(step < 1000)
step=step+1;
    
    %E-step
                
 for j=1:k
     if (yvar0(j)==0)
         yvar0(j)=1e-10;
     end
     pdf_est(:,j)=eppdf1(y,ymean0(j),yvar0(j),beta0(j));
    
 end

prob0 = repmat(prop0, n, 1).*pdf_est;
    prob1 = sum(prob0, 2);
    idx=find(prob1==0);
     prob1(idx)=1e-10;
    h_est = prob0./repmat(prob1, 1, k);
     

    %M-step
    % maximize mean and variance
    xx=repmat(y',k,1);
%     ymean0=fsolve(@(m) diag(((xx-repmat(m,1,n)).^repmat(2*beta0-ones(k,1),1,n))*h_est),ymean0);
%     ymean0=fsolve(@(m) beta0.*diag(((xx-repmat(m,1,n)).^repmat(2,k,n).^repmat(beta0,1,n)./(xx-repmat(m,1,n)) )*h_est),ymean0,options);
    ymean0=fsolve(@(m) beta0.*diag(((xx-repmat(m,1,n)).^repmat(2,k,n).^repmat(beta0-0.5,1,n))*h_est),ymean0,options);
%     minfun=@(m) abs(diag(((xx-repmat(m,1,n)).^repmat(2,k,n).^repmat(beta0,1,n)./(xx-repmat(m,1,n)) )*h_est)-zeros(k,1));
%     [ymean0,fval]=fmincon(minfun,ymean0);
%     ymean0=fsolve(@(m) diag(abs(xx-repmat(m,1,n)).^repmat(p0-ones(k,1),1,n)*h_est*sign(xx-repmat(m,1,n))),ymean0);
%     ymean0=fsolve(@(m) ((beta0.*diag((xx-repmat(ymean0,1,n)).^repmat(2,k,n).^repmat(beta0,1,n)*h_est)./sum(h_est,1)').^(1./beta0)).*beta0.*diag(((xx-repmat(m,1,n)).^repmat(2,k,n).^repmat(beta0,1,n)./(xx-repmat(m,1,n)) )*h_est),ymean0);

    yvar0=(beta0.*diag((xx-repmat(ymean0,1,n)).^repmat(2,k,n).^repmat(beta0,1,n)*h_est)./sum(h_est,1)').^(1./beta0);
%     for i=1:k
%         eta0(i)=sum(h_est(:,i))/(p0(i)*(abs(y'-repmat(ymean0(i),1,n))).^p0(i)*h_est(:,i));
%     end
  
   
    
    % maximize prop
    
    prop1 = sum(h_est, 1)/n;
    
    prop1 = prop1/sum(prop1);
    
         
    

       prop0 = prop1;
    
     %  lkhd0 = lkhd;
     %  lkhd0_p = lkhd_p;
       
   
%     ymean0
%     yvar0
%     prop0
prob0 = repmat(prop0, n, 1).*pdf_est;
prob1 = sum(prob0, 2);
lkhd1 = 2*sum(log(prob1))-k*Dof*log(n);% BIC
delta = abs(lkhd0-lkhd1)  ;
lkhd0=lkhd1;
likeli(step)=lkhd0;
like(step)=sum(log(prob1));

end;


prob0 = repmat(prop0, n, 1).*pdf_est;
prob1 = sum(prob0, 2);
lkhd = 2*sum(log(prob1))-k*Dof*log(n);% BIC


prop = prop0;
ym = ymean0;
yvar= yvar0;
