function [ BestModel, minBIC,numComponents] = gmmSel(y, mmax);

% Gaussian Normal Mixture through BIC
% Estimation and Selection 
%
% y -- data n*p
% mmax -- max no. of normal mixtures
%
% k -- estimated no. of normal mixtures
% prop -- estimated proportions
% ymean -- estimated normal means
% yvar -- estimated normal variance matrices

% EM estimation
BIC=[];
GMModels=[];

for i=1:mmax
    
    %i
   %GMModels{i}=fitgmdist(y,i);
   %GMModels{i}=fitgmdist(y,i,'Start','plus');
   GMModels{i}=fitgmdist(y,i,'RegularizationValue',0.01);
   BIC(i)=GMModels{i}.BIC;
end
[minBIC, numComponents]=min(BIC);%BIC=-loglikelihood+...????????
BestModel=GMModels{numComponents}