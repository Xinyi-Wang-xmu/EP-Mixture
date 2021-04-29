function [model,R,p] = NmaximizationModel(x,R,p,C)

n = length(x);
k = size(R,2);

% update Pi
Tmp1 = 1/(1-2*k*C)*(sum(R)/(n)-2*C);
Pi = max(zeros(k,1),Tmp1')';
N = sum(R);
for l=1:k
   fun=@(mu) sum(p(l)*R(:,l)'.*((abs(x-mu*ones(1,n))).^(p(l)-1)).*sign(x-mu*ones(1,n))); %??
   mu(l)=fsolve(fun,mean(x));
end

for l = 1:k
        tmp = sum(R(:,l)'.*((abs(x-mu(l)*ones(1,n))).^(p(l))));
        eta(l) = N(l)/(p(l)*tmp);
end


    model.mu=mu;
    model.eta = eta;
    model.Pi = Pi;
  