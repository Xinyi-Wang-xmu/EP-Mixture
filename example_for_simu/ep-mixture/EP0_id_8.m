% EP0 method for model 8 (as an example)
tic;
clear all;clc;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:100
    i
    ymean=[-3,1,3];
    ysigma=[0.25,4,0.125];
    prop=[0.4,0.3,0.3];
    n=300;
    sigmap=[0.4368,1.7472,0.3089];%[(2*0.25^1.5/3)^(1/3),,]
    p0=[3,3,3];
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rnormp(ni,ymean(m),sigmap(m),p0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rnormp(ni,ymean(k0),sigmap(k0),p0(k0));
y = [y, y1];
CC_PMoEP = sqrt(log(length(y)))*logspace(-4,-2,50);
len_lam_PMoEP = length(CC_PMoEP); 
BIC_max=-10^40;
BIC_back = [];
    
for l = 1:len_lam_PMoEP 

    % parameter setting 
    param_PMoEP.maxiter = 100;      % the maximum iteration number
    param_PMoEP.k = 10;              % the initialized number of mixture components
    param_PMoEP.display = 0;        % the display setting
    param_PMoEP.tol = 1.0e-10;      % the tolerance
    p=3*ones(1,10);
    C= CC_PMoEP(l);
    disp(['lambda is: ',num2str(C)]);

    % call main function 
    [label,model,llh,llh_BIC,p] = EM_PMoEP(y,param_PMoEP,p,C);
    hat_K = length(model.Pi);
    N=length(y);
    BIC_PMoEP = 2*llh_BIC - 3*hat_K*log(N);
    BIC_PMo(l) = 2*llh_BIC - 3*hat_K*log(N);

    BIC_back = [BIC_back BIC_PMoEP];
   
   if BIC_max<BIC_PMoEP
       p_max=p;
       model_max=model;
       BIC_max=BIC_PMoEP;
       lambda_max=CC_PMoEP(l);
   end
   
end

k=length(model_max.Pi);

kk(i)=k;

end

save Simu3_caoeppenalty_3001008 kk

toc;

%%
tic;
clear all;clc;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:100
    i
    ymean=[-3,1,3];
    ysigma=[0.25,4,0.125];
    prop=[0.4,0.3,0.3];
    n=600;
    sigmap=[0.4368,1.7472,0.3089];%(2*0.25^1.5/3)^(1/3)
    p0=[3,3,3];
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rnormp(ni,ymean(m),sigmap(m),p0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rnormp(ni,ymean(k0),sigmap(k0),p0(k0));
y = [y, y1];
CC_PMoEP = sqrt(log(length(y)))*logspace(-4,-2,50);
len_lam_PMoEP = length(CC_PMoEP); 
BIC_max=-10^40;
BIC_back = [];

for l = 1:len_lam_PMoEP 

    % parameter setting 
    param_PMoEP.maxiter = 100;      % the maximum iteration number
    param_PMoEP.k = 10;              % the initialized number of mixture components
    param_PMoEP.display = 0;        % the display setting
    param_PMoEP.tol = 1.0e-10;      % the tolerance
    p=3*ones(1,10);
    C= CC_PMoEP(l);
    disp(['lambda is: ',num2str(C)]);

    % call main function 
    [label,model,llh,llh_BIC,p] = EM_PMoEP(y,param_PMoEP,p,C);
    hat_K = length(model.Pi);
    N=length(y);
    BIC_PMoEP = 2*llh_BIC - 3*hat_K*log(N);
    BIC_PMo(l) = 2*llh_BIC - 3*hat_K*log(N);

    BIC_back = [BIC_back BIC_PMoEP];
   
   if BIC_max<BIC_PMoEP
       p_max=p;
       model_max=model;
       BIC_max=BIC_PMoEP;
       lambda_max=CC_PMoEP(l);
   end
   
end
k=length(model_max.Pi);

kk(i)=k;

end

save Simu3_caoeppenalty_6001008 kk

toc;

%%
tic;
clear all;clc;
addpath(genpath(pwd));
s = RandStream('mt19937ar','Seed', 5489);
RandStream.setGlobalStream(s);
for i=1:100
    i
    ymean=[-3,1,3];
    ysigma=[0.25,4,0.125];
    prop=[0.4,0.3,0.3];
    n=1000;
    sigmap=[0.4368,1.7472,0.3089];%(2*0.25^1.5/3)^(1/3)
    p0=[3,3,3];
    k0=length(prop);
    y = [];
for m = 1:(k0-1);
    ni = floor(n*prop(m));
    y1 = rnormp(ni,ymean(m),sigmap(m),p0(m));
    y = [y, y1];
end;
ni=n-length(y);
y1 = rnormp(ni,ymean(k0),sigmap(k0),p0(k0));
y = [y, y1];
CC_PMoEP = sqrt(log(length(y)))*logspace(-4,-2,50);
len_lam_PMoEP = length(CC_PMoEP); 
BIC_max=-10^40;
BIC_back = [];

for l = 1:len_lam_PMoEP 

    % parameter setting 
    param_PMoEP.maxiter = 100;      % the maximum iteration number
    param_PMoEP.k = 10;              % the initialized number of mixture components
    param_PMoEP.display = 0;        % the display setting
    param_PMoEP.tol = 1.0e-10;      % the tolerance
    p=3*ones(1,10);
    C= CC_PMoEP(l);
    disp(['lambda is: ',num2str(C)]);

    % call main function 
    [label,model,llh,llh_BIC,p] = EM_PMoEP(y,param_PMoEP,p,C);
    hat_K = length(model.Pi);
    N=length(y);
    BIC_PMoEP = 2*llh_BIC - 3*hat_K*log(N);
    BIC_PMo(l) = 2*llh_BIC - 3*hat_K*log(N);

    BIC_back = [BIC_back BIC_PMoEP];
   
   if BIC_max<BIC_PMoEP
       p_max=p;
       model_max=model;
       BIC_max=BIC_PMoEP;
       lambda_max=CC_PMoEP(l);
   end
   
end

k=length(model_max.Pi);
kk(i)=k;

end

save Simu3_caoeppenalty_10001008 kk

toc;