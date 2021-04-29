function [model,R,p,echo] = maximizationModel(x,R,eta,mu,p,C)

n = length(x);
k = size(R,2);

% update Pi
Tmp1 = 1/(1-2*k*C)*(sum(R)/(n)-2*C);
Pi = max(zeros(k,1),Tmp1')';

if length(find(Pi==0))~=length(Pi)
    Pi = Pi/sum(Pi);
end

ind = find(Pi==0); len_ind = length(ind);

if k > len_ind
    % delete component
    ind = find(Pi==0);
    k = k - length(ind);
    Pi(ind) = [];
    mu(ind)=[];
    R(:,ind) = [];
    Suma=sum(R,2);% 处理sum(R,2)元素为0，出现NaN的情况
    a=find(Suma==0);
    Suma(a)=1e-6;
    R = R./repmat(Suma,1,k);
    eta(ind) = [];
    p(ind) = [];
    
    N = sum(R);
   
    for l=1:k
      fun=@(mu) sum(p(l)*R(:,l)'.*((abs(x-mu*ones(1,n))).^(p(l)-1)).*sign(x-mu*ones(1,n))); %??
      mu(l)=fsolve(fun,mean(x));
    end
    % update lambda
    for l = 1:k
        tmp = sum(R(:,l)'.*((abs(x-mu(l)*ones(1,n))).^(p(l))));
        eta(l) = N(l)/(p(l)*tmp);
    end
    
    model.eta = eta;
    model.Pi = Pi;
    model.mu=mu;
    echo = 0;
else
    model.eta = eta;
    model.Pi = Pi;
    model.mu=mu;
    echo = 1;
end