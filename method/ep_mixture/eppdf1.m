function pro=eppdf1(x,mu,sigma2,beta)
% x is n*p mu is 1*p pr is p*n
%sigma2=sigma^2
[n,p]=size(x);
k=2*beta/(gamma(1/2/beta)*2^(1+1/2/beta)*sigma2^(1/2));
% pr=k*(det(sigma2))^(-0.5)*exp(-0.5*(x-repmat(mu,n,1)).^beta);
pro=k*exp(-0.5/sigma2^beta*(x-repmat(mu,n,1)).^2.^beta);
end
