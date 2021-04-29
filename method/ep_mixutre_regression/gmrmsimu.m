function []=gmrmsimu(id,rept,n,k0,b,prop,bet,ysigma,mmax);
% b k0*d
% prop 1*k0
% ysigma 1*k0
% bet 1*1

d=length(b)-1;
beta0=bet*ones(1,k0);
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    x=normrnd(0,1,n,d);
    x=[ones(n,1),x];
    k0=length(prop);
    y = [];
    xx=x;
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,0, ysigma(m), beta0(m))'+xx(1:ni,:)*b(m,:)';
    y = [y; y1];
    xx=xx((ni+1):length(xx),:);
end;
ni=n-length(y);
y1 = rep(ni,0, ysigma(k0),beta0(k0))'+xx(1:ni,:)*b(k0,:)';
y = [y; y1];

[k_max, prop_max, yb_max, yvar_max,lambda_max,lkhd_max] = gmrm_EstSel(y,x, mmax);

k0_g(i)=k_max;
propo_g(i,1:k_max)=prop_max;
yb_g(1:k_max,1:(d+1),i)=yb_max;
sigma_g(1:k_max,i)=yvar_max;
BIC_g(i)=lkhd_max;
% likelii(i,1:step)=likeli;
% likei(i,1:step)=like;

end
save(['simu_gmrm_',num2str(id)], 'k0_g', 'propo_g', 'yb_g', 'sigma_g', 'BIC_g')
