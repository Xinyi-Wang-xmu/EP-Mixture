function [k_max, prop_max, yb_max, yvar_max,beta_max,lambda_max,lkhd_max,step_max] = epr_EstSel(y, x,beta, lambda, mmax);

% EP Mixture Regression
% Estimation and Selection 
%
% y -- dependent data n*1
% x -- independent data n*d
% mmax -- max no. of normal mixtures
% k -- estimated no. of normal mixtures
% prop -- estimated proportions
% yb -- estimated coefficients
% yvar -- estimated sigma2

% EM estimation

% lambda=sqrt(log(length(y)))*logspace(-1,1,40)/length(y);

%lambda=2;
lkhd_max=-10^40;
lkhd_back = [];
k_max=0;


for i=1:length(lambda)
    
%     i
    
   [k, prop, yb, yvar,betab,lkhd,step] = epr_EM(y, x, beta, lambda(i), mmax);
   lkhd_back = [lkhd_back lkhd];
   
   if lkhd_max<lkhd
       k_max=k;
       prop_max=prop;
       yb_max=yb;
       yvar_max=yvar;
       beta_max=betab;
       lkhd_max=lkhd;
       lambda_max=lambda(i);
       step_max=step;
   end
end 
pause(0);
          