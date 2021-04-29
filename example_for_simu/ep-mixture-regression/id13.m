clear all;
clc;
rept=100;
n=100;
b=[1.9,-2.2,0.8,-1.3;2.4,2.5,-2.4,0.3];
d=3;
ysigma=[2,1];
prop=[0.5,0.5];
beta0=1.5*[1,1];
mmax=10;
beta=1.5*ones(1,mmax);
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

lambda=sqrt(log(length(y)))*logspace(-4,1,100);


[k_max, prop_max, yb_max, yvar_max,beta_max,lambda_max,lkhd_max,step_max] = epr_EstSel(y,x, beta, lambda,mmax);

k0_e(i)=k_max;
propo_e(i,1:k_max)=prop_max;
yb_e(1:k_max,1:(d+1),i)=yb_max;
sigma_e(1:k_max,i)=yvar_max;
BIC_e(i)=lkhd_max;
% likelii(i,1:step)=likeli;
% likei(i,1:step)=like;
stepi_e(i)=step_max;

end
save simu_epmr_1_1 k0_e propo_e yb_e sigma_e BIC_e stepi_e
%%
clear all;
clc;
rept=100;
n=100;
b=[1.9,-2.2,0.8,-1.3;2.4,2.5,-2.4,0.3];
d=3;
ysigma=[2,1];
prop=[0.5,0.5];
beta0=1.5*[1,1];
mmax=10;
beta=1.5*ones(1,mmax);
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
save simu_gmrm_1_1 k0_g propo_g yb_g sigma_g BIC_g 
%%