function [label,model,llh,llh_BIC,p] = EM_PMoEP(X,param,p,lambda)
% Description: EM algorithm for fitting PMoEP model.
% USAGE: [label,model,TW,OutU,OutV,llh,llh_BIC,p] = EM_PMoEP(InW,InX,r,param,p,lambda)
%Input:
%   X: d x n input data matrix
%   param:
%      --param.maxiter: maximal iteration number
%      --param.k: the number of mixture components
%      --param.display: display the iterative process
%      --param.tol: the thresholding for stop
%   p: the candidates components
%   lambda: the tuning parameter

%Output:
%   label: the labels of the noises
%   model: model.eta, the precisions of the different EPs
%          model.Pi,the mixing coefficients
%   W: d x n weighted matrix
%   OutU: the final factorized matrix U
%   OutV: the final factorized matrix V
%   llh:  the log likelihood
%   llh_BIC:  the log likelihood used in BIC criterion
%   p: the selected components

% Author: Xiangyong Cao(caoxiangyong45@gmail.com)
% If you have any quesion, please contact Xiangyong Cao

% Reference paper: 
% @InProceedings{Cao_2015_ICCV,
% author = {Cao, Xiangyong and Chen, Yang and Zhao, Qian and Meng, Deyu and Wang, Yao and Wang, Dong and Xu, Zongben},
% title = {Low-Rank Matrix Factorization Under General Mixture Noise Distributions},
% journal = {The IEEE International Conference on Computer Vision (ICCV)},
% month = {June},
% year = {2015}
% }

% version 1.0, date: 12-24-2015

% initialization
[d,n] = size(X);
if (~isfield(param,'maxiter'))
    maxiter = 100;
else
    maxiter = param.maxiter;
end


if (~isfield(param,'k'))
    k = 3;
else
    k = param.k;
end

if (~isfield(param,'display'))
    display = 0;
else
    display = param.display;
end

if (~isfield(param,'tol'))
    tol = 1e-7;
else
    tol = param.tol;
end

param.method = 2;
%Initialize the parameters of MoEP
if length(p)~=length(find(p==2))
    R = initialization_PMoEP(X,k,'random'); % initialize label posterior probability??
    [~,label(1,:)] = max(R,[],2);% R每行最大值的位置？？
    R = R(:,unique(label));
    %eta=5*ones(1,k)+randn(1,k);
    eta = 10*rand(1,k); % eta = 10*ones(1,k);
    nk = sum(R,1);
    Pi=1/k*ones(1,k)+randn(1,k);
    %Pi =  nk/size(R,1);
    %mu=2*mean(X) * sort(rand(k,1));
    model.eta = eta;
    model.Pi = Pi;
    %model.mu=mu;
else
    R = initialization_PMoG(X,k);
    [~,label(1,:)] = max(R,[],2);
    R = R(:,unique(label));
    model.Sigma = rand(1,k);
    nk = sum(R,1);
    model.Pi = nk/size(R,1)+randn(1,k);
    %model.mu=2*mean(X) * sort(rand(k,1));
    model.mu = zeros(1,k);
end

% llh = -inf(1,maxiter);
converged = false;


t = 1;

%Initialized E Step R(t)
if length(p)~=length(find(p==2))
    [R, llh(t),llh_BIC] = expectation(X,model, p, lambda);
else
    [R, llh(t),llh_BIC] = expectation_PMoG(X,model,lambda);
end

while ~converged && t < maxiter
   
    t = t+1;
   
    % M Step pi,eta(t), not t+1
    if length(p)~=length(find(p==2))
        %[model,~,p,echo] = maximizationModel(X,R,model.eta,model.mu,p,lambda);
        [model,~,p,echo] = maximizationModel(X,R,model.eta,p,lambda);
    else
        [model,R,echo,k] = maximizationModel_PMoG(X,R,lambda);
        p = 2*ones(1,k);
    end
    
    if echo==1 %all k component has pi=0
        llh=llh(t-1);
        llh_BIC=llh; label=[]; 
        return;
    end
    
    % E Step update R 
    if length(p)~=length(find(p==2))
        [R, llh(t),~] = expectation(X,model, p, lambda);
        L1 = llh(t);
    else
        Sigma = 1./(2*model.eta);
        model.Sigma = Sigma;
        [R, llh(t),llh_BIC] = expectation_PMoG(X,model,lambda);
        L1 = llh(t);
    end
    
    [~,label] = max(R,[],2);
    u = unique(label);   % non-empty components
    if size(R,2) ~= size(u,2)
        R = R(:,u);   % remove empty components
        R = R./repmat(sum(R,2),1,size(R,2));
        if length(p)~=length(find(p==2))
            model.eta = model.eta(u);
            %model.mu=model.mu(u);
            model.Pi = model.Pi(u);
            model.Pi = model.Pi/sum(model.Pi);
            p = p(u);
        else
            model.mu = model.mu(u);
            model.Sigma = model.Sigma(u);
            model.Pi = model.Pi(u);
        end
    end
    converged = abs(llh(t)-llh(t-1)) < tol;
    k = length(u);
    if mod(t,10)==0
    disp(['Iteration ',num2str(t),': there are ',num2str(k),...
        ' EP components and the corresponding p are ',num2str(p)])
    end
end

[~,label] = max(R,[],2);

if ~display
    disp(['The likelihood in this step is ',num2str(L1),';']);
    if length(p)~=length(find(p==2))
        disp(['There are ',num2str(k),' hylaplace noises mixed in data']);
        disp(['The selected p is ',num2str(p)]);
        disp(['with precisions ',num2str(model.eta)]);
        disp(['with weights ',num2str(model.Pi),'.']);
    else
        disp(['There are ',num2str(k),' Gaussian noises mixed in data']);
        disp(['with means ',num2str(model.mu)]);
        disp(['with variances ',num2str(model.Sigma)]);
        disp(['with weights ',num2str(model.Pi),'.']);
    end
    disp(['L2 RMSE is ', num2str(sqrt(mean(X.^2)))]);
    disp(['L1 RMSE is ', num2str(mean(abs(X)))]);
end
if converged
    fprintf('Converged in %d steps.\n',t-1);
else
    fprintf('Not converged in %d steps.\n',maxiter);
end
