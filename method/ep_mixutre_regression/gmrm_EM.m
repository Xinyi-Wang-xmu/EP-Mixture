function  [k, prop, yb, yv, lkhd] = gmrm_EM(y,x, lambda, mmax);

% EM algorithm for Gaussian Mixture Model 
% 
% y -- data
% mmax -- maximum no. of mixtures
% 
% k -- no. of mixtures
% prop -- mixture proportions 
% ym -- mean vectors
% yv -- covariance matrices
%
%

k = mmax;
[n, p] = size(x);

% K-means Initialization
[prop0, ymean0, yvar0] = gmm_kmean(y, mmax);
if (sum(isnan(yvar0))>0)
    [prop0, ymean0, yvar0] = gmm_kmean(y, mmax);
end

yvar0=yvar0(:);

yb0=ymean0./mean(x);% k*p

if (sum(prop0==0) > 0) 
   ind = (prop0>0);
   k = sum(ind);
   prop0= prop0(ind);
   yb0 = yb0(ind, :);
   yvar0 = yvar0(ind);
end;



epsilon = 1e-4;
step = 0;
delta = 1;
                %lambda = 1;
Dof = 2 + p ; 
%lkhd0 = -10000;
%lkhd0_p = -10000;

while (delta >= epsilon)&(step < 200)
    
    step = step+1;
    %E-step
                %for i=1:k;
                    %for j=1:n;
                        %pdf_est(j, i) = mvnpdf(y(j, :)', ymean0(i, :)', yvar0(:, :, i));
                    %end;
                %end;
    for i=1:k;
        yc = y - x*yb0(i,:)';
        pdf_est(:, i) = normpdf(yc,0,sqrt(yvar0(i)));    
    end;
     
    prob0 = repmat(prop0, n, 1).*pdf_est;
    prob1 = sum(prob0, 2);
     idx=find(prob1==0);
     prob0(idx,:)=1/k;
%      if (length(idx)==n)
%          prob0=1/k*ones(n,k);
%      else
%           prob1(idx,:)=1e-10;
%      end
         prob1 = sum(prob0, 2);
    h_est = prob0./repmat(prob1, 1, k);

    %%% LKHD difference    
    % lkhd = sum(log(prob1));
    % lkhd_p = lkhd - lambda* Dof * sum(log( (epsilon + prop0)/epsilon ));
    % disp(lkhd - lkhd0);
    % disp(lkhd_p - lkhd0_p);
   
    
    %M-step
    % maximize mean and variance
    ymean0 = h_est'*y./repmat(sum(h_est, 1)', 1, p);
    
    for i=1:k;
        
        yb0(i,:)=inv(x'*diag(h_est(:,i))*x)*x'*diag(h_est(:,i))*y;
%         yb0(i,:)=fsolve(@(m) x'*((y-x*m').*h_est(:,i)),yb0(i,:));
     
        yvar0(i)=(h_est(:,i)'*((y-x*yb0(i,:)').^2))/sum(h_est(:,i));

        % yvar0(:, :, i) = (y - x*yb0(i,:)')'*diag(h_est(:,i))*(y - x*yb0(i,:)')/sum(h_est(:,i));
    end;
    
    % maximize prop
    if step==1
        LL=k;
    else
        LL=sum(prop1<lambda/sqrt(n));
    end
    
    prop1 = max(0, (sum(h_est, 1) - lambda*Dof)/(n - k*lambda * Dof)) ;% lambda here is n*lambda in paper?
   
 %   prop1 = max(0, (sum(h_est, 1))/(n - k*lambda * Dof));
 %   prop1(prop1<epsilon) = 0;
    prop1 = prop1/sum(prop1); 
    
    delta = sum(abs(prop1-prop0)); 
         
    if (sum(prop1==0) > 0) 
       ind = (prop1 >0);
       yb0 = yb0(ind, :);
       yvar0 = yvar0(ind);
       pdf_est = pdf_est(:, ind);
       prop0 = prop1(ind);
        
       k = sum(prop1 > 0);
       delta = 1;
       
     %  lkhd0 = -10000;
     %  lkhd0_p = -10000;
     %  disp('=====');
    else
       prop0 = prop1;
    
     %  lkhd0 = lkhd;
     %  lkhd0_p = lkhd_p;
       
    end;
    
end;


prob0 = repmat(prop0, n, 1).*pdf_est;
prob1 = sum(prob0, 2);
lkhd = 2*sum(log(prob1))-k*Dof*log(n);% BIC


prop = prop0;
yb = yb0;
yv = yvar0;

 
    