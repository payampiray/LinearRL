function [pii, expv, M] = core_lrl(T,c,M)
% implenetation of the linear RL model
% 
% [pii, expv, M] = core_lrl(T,c)
% T is the transition probability under the default policy
% c is the cost vector (i.e. negative reward) across all states
% pii is the optimized decision policy
% expv = exp(v), where v is the value function
% M is the default representation (DR)
% 
% [pii, expv, M] = core_lrl(T,c,M)
% This one relies on a given DR matrix (M) as the input.
% 
% See Piray and Daw (2019), "A common model explaining flexible decision
% making, grid fields and cognitive control", biorxiv.
% ------------------------

% reward vector across all states
r = -c;

% terminal states
terminals = diag(T)==1;

% computing M (if not given)
if nargin<3
    L = diag(exp(c)) - T;
    L(terminals,:) = [];
    L(:,terminals) = [];
    M = (L^-1);
end

% P = T_{NT}
P = T(~terminals,terminals);
expr = exp(r(terminals));

expv_N = M*P*expr;
expv = zeros(size(r));
expv(~terminals) = expv_N;
expv(terminals) = exp(r(terminals));

% A matrix formulation of equation 6 of manuscript
G = T*expv;
zg = expv'./G;
pii = T.*zg;
end
