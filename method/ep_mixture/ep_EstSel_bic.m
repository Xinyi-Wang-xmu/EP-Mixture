function [k_max, prop_max, ymean_max, yvar_max,beta_max,lkhd_max,step_max] = ep_EstSel_bic(y, beta, mmax);

% EP Mixture through bic
% Estimation and Selection 
%
% y -- data n*1
% mmax -- max no. of normal mixtures
%
% k -- estimated no. of normal mixtures
% prop -- estimated proportions
% ymean -- estimated means
% yvar -- estimated sigma2

% EM estimation

% lambda=sqrt(log(length(y)))*logspace(-1,1,40)/length(y);

%lambda=2;
lkhd_max=-10^40;
lkhd_back = [];
k_max=0;


for i=1:mmax
    
    %i
    beta0=beta(1:i);
   [k, prop, ymean, yvar,betab,lkhd,step] = ep_EM_bic(y, beta0, i);
   lkhd_back = [lkhd_back lkhd];
   
   if lkhd_max<lkhd
       k_max=k;
       prop_max=prop;
       ymean_max=ymean;
       yvar_max=yvar;
       beta_max=betab;
       lkhd_max=lkhd;
       step_max=step;
   end
end 
pause(0);
           