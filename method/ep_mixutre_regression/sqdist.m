function dist = sqdist(y, ymean);

[n, p] = size(y);
[k, p] = size(ymean);

for i = 1:k;
    dist(:, i) = sum((y - repmat(ymean(i,:), n, 1)).^2, 2);
end;
