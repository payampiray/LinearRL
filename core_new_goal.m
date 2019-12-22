function [U, z] = core_new_goal(P0,M0,gs,q)
% plans towards a new goal using woodbury
% gs: the goal state 

P = P0;
P(gs,:) = 0; P(gs,gs) = 1;
e = zeros(size(P,1),1); e(gs) = 1;


dp = P0(e==1,:) - P(e==1,:);
m0 = M0(:,e==1);
p = P(:,e==1); p(e==1)=0;
z0 = M0*p;

alpha = (dp*z0)/(1+dp*m0);
z2 = z0(~e) - alpha*m0(~e);

z = nan(size(P,1),1);
z(~e) = z2*exp(-q(e==1));
z(e==1) = exp(-q(e==1));

G = P*z;

zg = z'./G;
U = P.*zg;
end