function [R, llh, llh_BIC] = Nexpectation(x, model, p, C)
eta = model.eta;  
Pi = model.Pi;
mu=model.mu;
n = length(x);
k = length(Pi);
logRho = zeros(n,k);
epsilon = 1e-6;

for l = 1:k
    logRho(:,l) = Nlogeppdf(x,mu(l),eta(l),p(l));
end
logRho = bsxfun(@plus,logRho,log(Pi));
T = logsumexp(logRho,2); %log(sum(Rho))
llh_BIC = sum(T); 
llh = sum(T)-n*C*2*(sum(log(epsilon+Pi))-k*log(epsilon)); % uncomplete_penalized_loglikelihood conditional??
% llh = sum(T)/n - n*C*(sum(log(epsilon+Pi))-k*log(epsilon)); % loglikelihood %%divide n is not correct
logR = bsxfun(@minus,logRho,T); %R=Rho/sum(Rho)
R = exp(logR);