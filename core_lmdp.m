function [U, z, MNN] = core_lmdp(P,c,MNN)

terminals = diag(P)==1;

if nargin<3
L = diag(exp(c)) - P;
L(terminals,:) = [];
L(:,terminals) = [];
MNN = (L^-1);
end

Pnt = P(~terminals,terminals);
Qt = exp(-c(terminals));

znt = MNN*Pnt*Qt;
z = zeros(size(c));
z(~terminals) = znt;
z(terminals) = exp(-c(terminals));

G = P*z;

zg = z'./G;
U = P.*zg;
end
