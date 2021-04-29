function [prop, ymean, yvar] = gmm_kmean(y, k);

% K-means clustering
% 
% y    - (n x p) p-dimensional input data
% k    - no. of clusters
%
% prop -- proportions of each clusters
% ymean -- estimated mean vectors
% yvar -- estimated covariance matrices
%
% Tao Huang, June 8 2009, 
%
%  covariance ??? equal ???

[n, p] = size(y);

epsilon = 1e-4;   

ind = randperm(n);
ymean0 = y(ind(1:k), :);

step = 0;
rss0 = 100;

dist_mat = sqdist(y, ymean0);           % n x k matrix
  
% update class
[ymin, ycls] = min(dist_mat, [], 2);
rss1 = mean(ymin);
delta = rss0 - rss1;



while (delta>=epsilon)&(step <100)
      
  step = step+1;

  % update mean vector
  for i=1:k
      if sum(ycls==i)>0
         ymean0(i, :) = mean(y(ycls==i, :));
      end
  end;
  
  % update class 
  dist_mat = sqdist(y, ymean0);
  [ymin, ycls] = min(dist_mat, [], 2);
  
  % update rss
  rss0 = rss1;
  rss1 = mean(ymin);
  delta = rss0 - rss1;
  
end

ymean = ymean0;

for i=1:k
    if sum(ycls==i)<3
        prop(1,i)=0;
    else
         prop(1, i) = mean(ycls==i);
    yvar(:, :, i) = cov(y(ycls==i, :)); 
    end 
end

prop = prop/sum(prop);
