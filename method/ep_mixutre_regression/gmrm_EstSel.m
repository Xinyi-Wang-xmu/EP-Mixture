function [k_max, prop_max, yb_max, yvar_max,lambda_max,lkhd_max] = gmrm_EstSel(y,x, mmax);

% Gaussian Normal Mixture 
% Estimation and Selection 
%
% y -- data
% mmax -- max no. of normal mixtures
%
% k -- estimated no. of normal mixtures
% prop -- estimated proportions
% ymean -- estimated normal means
% yvar -- estimated normal variance matrices

% EM estimation

lambda=sqrt(log(length(y)))*logspace(-1,1,40);

%lambda=2;
lkhd_max=-10^40;
lkhd_back = [];
k_max=0;


for i=1:length(lambda)
    
    i;
    
   [k, prop, yb, yvar,lkhd] = gmrm_EM(y,x, lambda(i), mmax);
   lkhd_back = [lkhd_back lkhd];
   
   if lkhd_max<lkhd
       k_max=k;
       prop_max=prop;
       yb_max=yb;
       yvar_max=yvar;
       lkhd_max=lkhd;
       lambda_max=lambda(i);
   end
end 
pause(0);