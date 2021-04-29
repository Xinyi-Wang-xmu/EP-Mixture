function y = Nlogeppdf(x, mu,lambda, p)
n=length(x);
    for i=1:n
    m=x(i)-mu;
    y(i,:) = log(p) + log(lambda)/p - log(2) - log(gamma(1/p))-lambda*(abs(m)).^(p);%log pdf pf EP
end
