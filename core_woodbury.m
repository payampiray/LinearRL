function [U, z, D, L] = core_woodbury(L0,D0,P,q)
% P0 = core_griding(5);
% P0(end,:) = 0; P0(end,end)=1;
% q = .1*ones(size(P0,1),1);
% L0 = (diag(exp(q))-P0);
% D0 = L0^-1;
% ds = 5;
% P = P0;
% P(ds,:) = 0; P(ds,ds) = 1;

L = diag(exp(q))-P;

ds = find(sum(abs(L-L0),2)~=0);
e = zeros(size(P,1),1); e(ds) = 1;

d = L(ds,:) - L0(ds,:);
terms = diag(P)==1;
p = P(:,terms); p(terms,:)=0;
m0 = D0*e;

z0 = D0*p;

alpha = 1/(1+d*m0);
zc = z0 - alpha*m0*(d*z0);

z = zc*exp(-q(terms));
z(terms) = exp(-q(terms));

G = P*z;
zg = z'./G;
U = P.*zg;

D = D0 - 1/(1+d*m0)*m0*d*D0;

% DD = (diag(exp(q))-P)^-1;
% [uu, zz] = core_lmdp(P,q);
end