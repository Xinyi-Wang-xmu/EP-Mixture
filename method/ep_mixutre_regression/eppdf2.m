function pro=eppdf2(y,x,b,sigma2,beta)
% yis n*1,x is n*p, b is p*1
%sigma2=sigma^2
[n,p]=size(x);
k=2*beta/(gamma(1/2/beta)*2^(1+1/2/beta)*sigma2^(1/2));
% pr=k*(det(sigma2))^(-0.5)*exp(-0.5*(x-repmat(mu,n,1)).^beta);
pro=k*exp(-0.5/(sigma2^beta)*(y-x*b).^2.^beta);
end