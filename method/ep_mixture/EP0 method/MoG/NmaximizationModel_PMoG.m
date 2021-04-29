function [model,R,k] = NmaximizationModel_PMoG(X,R,lambda)

[d,n] = size(X);
k = size(R,2);

% nk = sum(R,1);

% mu = bsxfun(@times, X*R, 1./nk);

Tmp1 = 1/(1-k*lambda*2)*(sum(R)/n-lambda*2);
w = max(zeros(1,k),Tmp1);
nk = sum(R,1);
for i=1:k
    mu(i)=sum(R(:,i)'.*X)/nk(i);
end

sqrtR = sqrt(R); 
    for i = 1:k
        Xo = bsxfun(@minus,X,mu(i));
        Xo = bsxfun(@times,Xo,sqrtR(:,i)');
        Sigma(i) = (Xo*Xo')/nk(i);
        Sigma(i) = Sigma(i)+(1e-6); % add a prior for numerical stability
    end

    model.mu = mu; model.Sigma = Sigma; model.eta = 1./(2*Sigma); model.Pi = w; 
