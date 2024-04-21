function [Om] = make_Om_random_scan(m,n,r)

% ランダムスキャン
t0 = rand(m*n,1);
[~, indx] = sort(t0);

num_nonzero = floor(m*n*r);

Om = zeros(m,n);
Om( indx(1:num_nonzero)) = 1;

