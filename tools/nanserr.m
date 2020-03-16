function se = nanserr(x,dim)
if nargin<2, dim=1; end;
s = nanstd(x,[],dim);
n = sum(~isnan(x),dim);
se = s./sqrt(n);