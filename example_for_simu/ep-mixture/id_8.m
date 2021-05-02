clear all;
clc;
rept=100;
n=300;
ymean=[-3,1,3];
ysigma=[0.25,4,0.125];
prop=[0.4,0.3,0.3];
beta0=1.5*[1,1,1];
beta=1.5*ones(1,10);
mmax=10;
id=8;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,ymean(m), ysigma(m), beta0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rep(ni,ymean(k0), ysigma(k0),beta0(k0));
y = [y, y1];
y=y';
lambda=sqrt(log(length(y)))*logspace(-4,-3,10);

[k_maxe, prop_maxe, ymean_maxe, yvar_maxe,beta_maxe,lambda_maxe,lkhd_maxe,step_maxe]=ep_EstSel1(y, beta, lambda, mmax);
k0_ep(i)=k_maxe;
propo_ep(i,1:k_maxe)=prop_maxe;
yme_ep(1:k_maxe,i)=ymean_maxe;
sigma_ep(1:k_maxe,i)=yvar_maxe;
lambda0_ep(i)=lambda_maxe;
BIC_ep(i)=lkhd_maxe;
% likelii(i,1:step)=likeli;
% likei(i,1:step)=like;
stepi_ep(i)=step_maxe;
end

save (['Simu3_epmupenalty_knn_normsimu_',num2str(n),num2str(rept),num2str(id)], 'k0_ep', 'propo_ep', 'yme_ep', 'sigma_ep', 'lambda0_ep', 'BIC_ep', 'stepi_ep')
%%
clear all;
clc;
rept=100;
n=600;
ymean=[-3,1,3];
ysigma=[0.25,4,0.125];
prop=[0.4,0.3,0.3];
beta0=1.5*[1,1,1];
beta=1.5*ones(1,10);
mmax=10;
id=8;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,ymean(m), ysigma(m), beta0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rep(ni,ymean(k0), ysigma(k0),beta0(k0));
y = [y, y1];
y=y';
lambda=sqrt(log(length(y)))*logspace(-4,-3,10);

[k_maxe, prop_maxe, ymean_maxe, yvar_maxe,beta_maxe,lambda_maxe,lkhd_maxe,step_maxe]=ep_EstSel1(y, beta, lambda, mmax);
k0_ep(i)=k_maxe;
propo_ep(i,1:k_maxe)=prop_maxe;
yme_ep(1:k_maxe,i)=ymean_maxe;
sigma_ep(1:k_maxe,i)=yvar_maxe;
lambda0_ep(i)=lambda_maxe;
BIC_ep(i)=lkhd_maxe;
% likelii(i,1:step)=likeli;
% likei(i,1:step)=like;
stepi_ep(i)=step_maxe;
end

save (['Simu3_epmupenalty_knn_normsimu_',num2str(n),num2str(rept),num2str(id)], 'k0_ep', 'propo_ep', 'yme_ep', 'sigma_ep', 'lambda0_ep', 'BIC_ep', 'stepi_ep')
%%
clear all;
clc;
rept=100;
n=1000;
ymean=[-3,1,3];
ysigma=[0.25,4,0.125];
prop=[0.4,0.3,0.3];
beta0=1.5*[1,1,1];
beta=1.5*ones(1,10);
mmax=10;
id=8;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,ymean(m), ysigma(m), beta0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rep(ni,ymean(k0), ysigma(k0),beta0(k0));
y = [y, y1];
y=y';
lambda=sqrt(log(length(y)))*logspace(-4,-3,10);

[k_maxe, prop_maxe, ymean_maxe, yvar_maxe,beta_maxe,lambda_maxe,lkhd_maxe,step_maxe]=ep_EstSel1(y, beta, lambda, mmax);
k0_ep(i)=k_maxe;
propo_ep(i,1:k_maxe)=prop_maxe;
yme_ep(1:k_maxe,i)=ymean_maxe;
sigma_ep(1:k_maxe,i)=yvar_maxe;
lambda0_ep(i)=lambda_maxe;
BIC_ep(i)=lkhd_maxe;
% likelii(i,1:step)=likeli;
% likei(i,1:step)=like;
stepi_ep(i)=step_maxe;
end

save (['Simu3_epmupenalty_knn_normsimu_',num2str(n),num2str(rept),num2str(id)], 'k0_ep', 'propo_ep', 'yme_ep', 'sigma_ep', 'lambda0_ep', 'BIC_ep', 'stepi_ep')
%%
clear all;
clc;
rept=100;
n=300;
ymean=[-3,1,3];
ysigma=[0.25,4,0.125];
prop=[0.4,0.3,0.3];
beta0=1.5*[1,1,1];
beta=1.5*ones(1,10);
mmax=10;
id=8;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,ymean(m), ysigma(m), beta0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rep(ni,ymean(k0), ysigma(k0),beta0(k0));
y = [y, y1];
y=y';
[ BestModel, minBIC,numComponents] = gmmSel(y, 10);
mu_g(i,1:numComponents)=BestModel.mu';
var_g(i,1:numComponents)=BestModel.Sigma;
Pi_g(i,1:numComponents)=BestModel.ComponentProportion;
BIC_g(i)=-minBIC;
yy_g(i,:)=y;
kk_g(i)=numComponents;
end

save(['Simu3_gmm_',num2str(n),num2str(rept),num2str(id)], 'mu_g', 'var_g', 'Pi_g', 'BIC_g', 'yy_g', 'kk_g') 

%%
clear all;
clc;
rept=100;
n=600;
ymean=[-3,1,3];
ysigma=[0.25,4,0.125];
prop=[0.4,0.3,0.3];
beta0=1.5*[1,1,1];
beta=1.5*ones(1,10);
mmax=10;
id=8;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,ymean(m), ysigma(m), beta0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rep(ni,ymean(k0), ysigma(k0),beta0(k0));
y = [y, y1];
y=y';
[ BestModel, minBIC,numComponents] = gmmSel(y, 10);
mu_g(i,1:numComponents)=BestModel.mu';
var_g(i,1:numComponents)=BestModel.Sigma;
Pi_g(i,1:numComponents)=BestModel.ComponentProportion;
BIC_g(i)=-minBIC;
yy_g(i,:)=y;
kk_g(i)=numComponents;
end

save(['Simu3_gmm_',num2str(n),num2str(rept),num2str(id)], 'mu_g', 'var_g', 'Pi_g', 'BIC_g', 'yy_g', 'kk_g') 

%%
clear all;
clc;
rept=100;
n=1000;
ymean=[-3,1,3];
ysigma=[0.25,4,0.125];
prop=[0.4,0.3,0.3];
beta0=1.5*[1,1,1];
beta=1.5*ones(1,10);
mmax=10;
id=8;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,ymean(m), ysigma(m), beta0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rep(ni,ymean(k0), ysigma(k0),beta0(k0));
y = [y, y1];
y=y';
[ BestModel, minBIC,numComponents] = gmmSel(y, 10);
mu_g(i,1:numComponents)=BestModel.mu';
var_g(i,1:numComponents)=BestModel.Sigma;
Pi_g(i,1:numComponents)=BestModel.ComponentProportion;
BIC_g(i)=-minBIC;
yy_g(i,:)=y;
kk_g(i)=numComponents;
end

save(['Simu3_gmm_',num2str(n),num2str(rept),num2str(id)], 'mu_g', 'var_g', 'Pi_g', 'BIC_g', 'yy_g', 'kk_g') 

%%
clear all;
clc;
rept=100;
n=300;
ymean=[-3,1,3];
ysigma=[0.25,4,0.125];
prop=[0.4,0.3,0.3];
beta0=1.5*[1,1,1];
beta=1.5*ones(1,10);
mmax=10;
id=8;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,ymean(m), ysigma(m), beta0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rep(ni,ymean(k0), ysigma(k0),beta0(k0));
y = [y, y1];
y=y';
[k, pro, ymeann, yvar,lambda,lkhd] = gmm_EstSel(y, 10);
mu_gp(i,1:k)=ymeann';
var_gp(i,1:k)=yvar;
Pi_gp(i,1:k)=pro;
BIC_gp(i)=lkhd;
yy_gp(i,:)=y;
kk_gp(i)=k;
llambda_gp(i)=lambda;

end

save (['Simu3_gmmpenalty_',num2str(n),num2str(rept),num2str(id)], 'mu_gp', 'var_gp', 'Pi_gp', 'BIC_gp', 'yy_gp', 'kk_gp', 'llambda_gp')

%%
clear all;
clc;
rept=100;
n=600;
ymean=[-3,1,3];
ysigma=[0.25,4,0.125];
prop=[0.4,0.3,0.3];
beta0=1.5*[1,1,1];
beta=1.5*ones(1,10);
mmax=10;
id=8;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,ymean(m), ysigma(m), beta0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rep(ni,ymean(k0), ysigma(k0),beta0(k0));
y = [y, y1];
y=y';
[k, pro, ymeann, yvar,lambda,lkhd] = gmm_EstSel(y, 10);
mu_gp(i,1:k)=ymeann';
var_gp(i,1:k)=yvar;
Pi_gp(i,1:k)=pro;
BIC_gp(i)=lkhd;
yy_gp(i,:)=y;
kk_gp(i)=k;
llambda_gp(i)=lambda;

end

save (['Simu3_gmmpenalty_',num2str(n),num2str(rept),num2str(id)], 'mu_gp', 'var_gp', 'Pi_gp', 'BIC_gp', 'yy_gp', 'kk_gp', 'llambda_gp')

%%
clear all;
clc;
rept=100;
n=1000;
ymean=[-3,1,3];
ysigma=[0.25,4,0.125];
prop=[0.4,0.3,0.3];
beta0=1.5*[1,1,1];
beta=1.5*ones(1,10);
mmax=10;
id=8;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,ymean(m), ysigma(m), beta0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rep(ni,ymean(k0), ysigma(k0),beta0(k0));
y = [y, y1];
y=y';
[k, pro, ymeann, yvar,lambda,lkhd] = gmm_EstSel(y, 10);
mu_gp(i,1:k)=ymeann';
var_gp(i,1:k)=yvar;
Pi_gp(i,1:k)=pro;
BIC_gp(i)=lkhd;
yy_gp(i,:)=y;
kk_gp(i)=k;
llambda_gp(i)=lambda;

end

save (['Simu3_gmmpenalty_',num2str(n),num2str(rept),num2str(id)], 'mu_gp', 'var_gp', 'Pi_gp', 'BIC_gp', 'yy_gp', 'kk_gp', 'llambda_gp')

%%
clear all;
clc;
rept=100;
n=300;
ymean=[-3,1,3];
ysigma=[0.25,4,0.125];
prop=[0.4,0.3,0.3];
beta0=1.5*[1,1,1];
beta=1.5*ones(1,10);
mmax=10;
id=8;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,ymean(m), ysigma(m), beta0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rep(ni,ymean(k0), ysigma(k0),beta0(k0));
y = [y, y1];
y=y';
[k_max, prop_max, ymean_max, yvar_max,beta_max,lkhd_max,step_max] = ep_EstSel_bic(y, beta, mmax);

k0_e(i)=k_max;
propo_e(i,1:k_max)=prop_max;
yme_e(1:k_max,i)=ymean_max;
sigma_e(1:k_max,i)=yvar_max;
BIC_e(i)=lkhd_max;
% likelii(i,1:step)=likeli;
% likei(i,1:step)=like;
stepi_e(i)=step_max;

end

save (['Simu3_epmu_knn_normsimu_',num2str(n),num2str(rept),num2str(id)], 'k0_e', 'propo_e', 'yme_e', 'sigma_e', 'BIC_e', 'stepi_e')

%%
clear all;
clc;
rept=100;
n=600;
ymean=[-3,1,3];
ysigma=[0.25,4,0.125];
prop=[0.4,0.3,0.3];
beta0=1.5*[1,1,1];
beta=1.5*ones(1,10);
mmax=10;
id=8;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,ymean(m), ysigma(m), beta0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rep(ni,ymean(k0), ysigma(k0),beta0(k0));
y = [y, y1];
y=y';
[k_max, prop_max, ymean_max, yvar_max,beta_max,lkhd_max,step_max] = ep_EstSel_bic(y, beta, mmax);

k0_e(i)=k_max;
propo_e(i,1:k_max)=prop_max;
yme_e(1:k_max,i)=ymean_max;
sigma_e(1:k_max,i)=yvar_max;
BIC_e(i)=lkhd_max;
% likelii(i,1:step)=likeli;
% likei(i,1:step)=like;
stepi_e(i)=step_max;

end

save (['Simu3_epmu_knn_normsimu_',num2str(n),num2str(rept),num2str(id)], 'k0_e', 'propo_e', 'yme_e', 'sigma_e', 'BIC_e', 'stepi_e')

%%
clear all;
clc;
rept=100;
n=1000;
ymean=[-3,1,3];
ysigma=[0.25,4,0.125];
prop=[0.4,0.3,0.3];
beta0=1.5*[1,1,1];
beta=1.5*ones(1,10);
mmax=10;
id=8;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:rept
    i
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rep(ni,ymean(m), ysigma(m), beta0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rep(ni,ymean(k0), ysigma(k0),beta0(k0));
y = [y, y1];
y=y';
[k_max, prop_max, ymean_max, yvar_max,beta_max,lkhd_max,step_max] = ep_EstSel_bic(y, beta, mmax);

k0_e(i)=k_max;
propo_e(i,1:k_max)=prop_max;
yme_e(1:k_max,i)=ymean_max;
sigma_e(1:k_max,i)=yvar_max;
BIC_e(i)=lkhd_max;
% likelii(i,1:step)=likeli;
% likei(i,1:step)=like;
stepi_e(i)=step_max;

end

save (['Simu3_epmu_knn_normsimu_',num2str(n),num2str(rept),num2str(id)], 'k0_e', 'propo_e', 'yme_e', 'sigma_e', 'BIC_e', 'stepi_e')
