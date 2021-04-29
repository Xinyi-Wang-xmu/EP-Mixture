function  [k, prop, yb, yvar, beta0,lkhd,step] = epr_EM(y,x, beta,lambda, mmax)
% EM algorithm for ep Mixture Model with mulitivarite form
% 
% y -- dependent data n*1
% x -- independent data n*d
% mmax -- maximum no. of mixtures
% k -- no. of mixtures
% prop -- mixture proportions 
% yb -- coefficient vectors
% yv -- sigma2
% beta -- 1*k
%
%
% global y x n k d h_est beta0

% handles.y=y;
% handles.x=x; 
k = mmax;
[n, d] = size(x);
 beta0=beta';
% handles.n=n;
% handles.k=k;
% handles.d=d;
% handles.beta0=beta0;

% K-means Initialization
[prop0, ymean0, var0] = gmm_kmean(y, mmax);
% yvar0=var0(:);
% yvar0=ones(1,k);
% yb0=ymean0./mean(x); % k*d
% 
% uniform random initialization
yvar0=unifrnd(0,10,1)*ones(k,1);
% yvar0=unifrnd(0,10,k,1);
yb0=unifrnd(-10,10,k,d);

if (sum(prop0==0) > 0) 
   ind = (prop0>0);
   k = sum(ind);
   prop0= prop0(ind);
   yb0 = yb0(ind, :);
   yvar0 = yvar0(ind);
   beta0=beta0(ind);
end;
%  ymean0=zeros(k,1)+randn(k,1);
% yvar0=yvar0.*gamma(0.5./beta0)./(2.^(1./beta0).*gamma(1.5./beta0));
 
epsilon = 1e-4;
step = 0;
delta = 1;
                %lambda = 1;
Dof = 2 + d ; 
lkhd0 = 0;
%lkhd0_p = -10000;
options = optimoptions('fsolve','Display','none');


while (delta >= epsilon)&(step < 1000)
step=step+1;
    
    %E-step
                
 for j=1:k
     
     pdf_est(:,j)=eppdf2(y,x,yb0(j,:)',yvar0(j),beta0(j));
    
 end

prob0 = repmat(prop0, n, 1).*pdf_est;
    prob1 = sum(prob0, 2);
    idx=find(prob1==0);
    prob0(idx,:)=1e-10*ones(length(idx),k);
%     for i=1:length(idx)
%      prob0(idx(i),:)=1e-10*ones(1,k);
%     end
    prob1 = sum(prob0, 2);
    h_est = prob0./repmat(prob1, 1, k); %n*k
%     handles.h_est=h_est;
 

    %M-step
%     yb0=fzero(@findb, yb0')
%     yb0=fsolve(@findb, yb0', options);
%     yb0=yb0';
    for j=1:k
        try
        yb0(j,:)=fsolve(@(m) beta0(j)*x'*((y-x*m').^2.^(beta0(j)-1).*(y-x*m').*h_est(:,j)),yb0(j,:),options);
        end
        yvar0(j)=((beta0(j)*h_est(:,j)'*((y-x*yb0(j,:)').^2.^beta0(j)))/sum(h_est(:,j)))^(1/beta0(j));
    end
%     [a,b]=size(h_est);
%     if b==1
%         label=ones(1,n);
%     else
%         [M,label]=max(h_est');
%     end
%         xb=[];
% %         if (length(label)==1)
% %             label=label*ones(1,n);
% %         end
%         for no=1:n
%             xb(no,:)=x(no,:)*yb0(label(no),:)';
%         end
%         
%         yvar0=var(y-xb)*gamma(0.5./beta0)./(2.^(1./beta0).*gamma(1.5./beta0));
%          for j=1:k
%              ind=(label==j);
%              yvar0(j)=var(y(ind)-x(ind,:)*yb0(j,:)')*gamma(0.5/beta0(j))/(2^(1/beta0(j))*gamma(1.5/beta0(j)));
%          end

%        yvar0(j)=(beta0(j)*h_est(:,j)'*((y-x*yb0(j,:)').^2.^beta(j)))/sum(h_est(:,j))^(1/beta0(j));

%         error=0;
%         b0=yb0(j,:);
%         v0=yvar0(j);
%         while(error>=epsilon)
%         try
%         b1=fsolve(@(m) beta0(j)*x'*((y-x*m').^2.^(beta0(j)-1).*(y-x*m').*h_est(:,j))/v0,b0,options);
% %           syms m
% %           yb0(j,:)=vpasolve(beta0(j)*x'*((y-x*m').^2.^(beta0(j)-1).*(y-x*m').*h_est(:,j))==0,m)
%         end
%         v1=(beta0(j)*h_est(:,j)'*((y-x*b1').^2.^beta(j)))/sum(h_est(:,j))^(1/beta0(j));
%         error=abs(b1-b0)+abs(v1-v0);
%         b0=b1;
%         v0=v1;
%         end
%     end
    
%        yb0=fsolve(@(m) beta0.*((repmat(y,1,k)-x*yb0').^2.^repmat((beta0-0.5)',n,1).*h_est)'*x,yb0,options);
%        yvar0=(beta0.*diag(((repmat(y,1,k)-x*yb0')'.^repmat(2,k,n).^repmat(beta0,1,n)*h_est)./sum(h_est,1)')).^(1./beta0);

    
  
    
    % maximize prop
    if step==1
        LL=k;
    else
        LL=sum(prop1<lambda/sqrt(n));
    end
    prop1 = max(0, (sum(h_est, 1) - n*lambda*Dof)/(n - k*lambda*n* Dof)) ;
    
    prop1 = prop1/sum(prop1);
    
    delta = sum(abs(prop1-prop0)); 
         
    if (sum(prop1==0) > 0) 
       ind = (prop1 >0);
       yb0 = yb0(ind, :);
       yvar0 = yvar0(ind);
       pdf_est = pdf_est(:, ind);
       prop0 = prop1(ind);
       beta0=beta0(ind); 
       k = sum(prop1 > 0);
       delta = 1;
    else
       prop0 = prop1;
         
    end;
         
%     ymean0
%     yvar0
%     prop0
prob0 = repmat(prop0, n, 1).*pdf_est;
prob1 = sum(prob0, 2);
lkhd1 = 2*sum(log(prob1))-k*Dof*log(n);% BIC
lkhd0=lkhd1;
likeli(step)=lkhd0;
like(step)=sum(log(prob1))-n*lambda*Dof*sum(log(prop0));

end;


prob0 = repmat(prop0, n, 1).*pdf_est;
prob1 = sum(prob0, 2);
lkhd = 2*sum(log(prob1))-k*Dof*log(n);% BIC


prop = prop0;
yb = yb0;
yvar= yvar0;