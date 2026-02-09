function [rho, num, den] = fast_pearsoncor(x,y)
% [rho, num, den] = fast_pearsoncor(x,y)
% rho = num./den

[nt nv] = size(x);

dx = x - repmat(mean(x),[nt 1]);
dy = y - repmat(mean(y),[nt 1]);

num = sum(dx.*dy);
den = sqrt(sum(dx.*dx).*sum(dy.*dy));
rho = num./den;
ind = find(den==0);
rho(ind) = 0;

return;

