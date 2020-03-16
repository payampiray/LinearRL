function [U, expv] = core_new_goal(T0,M0,gs,q)
% plans towards a new goal using woodbury
% gs: the goal state 

T = T0;
T(gs,:) = 0; T(gs,gs) = 1;
e = zeros(size(T,1),1); e(gs) = 1;

dp = T0(e==1,:) - T(e==1,:);
m0 = M0(:,e==1);
p = T(:,e==1); p(e==1)=0;
z0 = M0*p;

alpha = (dp*z0)/(1+dp*m0);
expv2 = z0(~e) - alpha*m0(~e);

expv = nan(size(T,1),1);
expv(~e) = expv2*exp(-q(e==1));
expv(e==1) = exp(-q(e==1));

G = T*expv;
expvg = expv'./G;
U = T.*expvg;
end