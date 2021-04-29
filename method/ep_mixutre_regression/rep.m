function y = rep(n,mu,sigma2,beta)
% random number generator of the exponential power distribution
% input: n -  the number of the random number
%        mu - location parameter
%        sigma2 - scale parameter
%        beta - shape parameter

if sigma2<=0
    error('sigma2 must be positive');
end

 if beta==1
     y = normrnd(mu, sqrt(sigma2), 1,n);

 else
            qg = gamrnd(0.5/beta,2,1,n);
            z = qg.^(0.5/beta);
            c1 = zeros(1,n); c2 = zeros(1,n);
            tmp = rand(1,n);
            c1(find(tmp<=0.5)) = -z(find(tmp<=0.5));
            c2(find(tmp>0.5)) = z(find(tmp>0.5));
            z = c1 + c2;
            y = mu + z*sqrt(sigma2);%mean of y is mu
 end
 
end
% [nElems, centers] = hist(y, 500);
% bar(centers, nElems/n);axis([-10,60,0,0.01]);
