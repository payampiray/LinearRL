function se = serr(x,dim)
if nargin<2, dim=1; end;
s = std(x,[],dim);
n = size(x,dim);
se = s./sqrt(n);