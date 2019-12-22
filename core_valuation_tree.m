function soft_pi = core_valuation_tree(P,q,n,depth,beta)

% n = 7;
% k = 2;
% [P, A, xy, paths] = core_treeing(n,k);
% smax = paths(end,1:end-1);
% q = rand(size(P,1),1);
% depth = 3;

%-----------------------------
if nargin<4
    depth = n;
end
if nargin<5
    beta = 1;
end

availables = P>0;
terminals = diag(P)==1;

ns = size(P,1);
V = zeros(ns,1);
V(terminals) = q(terminals);

for i=n:-1:1
    snt = (2^(i-1)):(2^i-1);
    for s=snt
        nexts = P(s,:)>0;
        if i>depth
            V(s) = q(s) + mean(V(nexts));
        else
            p = P(s,nexts)';
            V(s) = q(s) + min(V(nexts));
        end
    end
end


soft_pi = zeros(ns,ns);
policy = zeros(ns,ns);
for s=1:ns
    nexts = find(availables(s,:));
    p = P(s,nexts)';
    [~,a] = min(p.*V(nexts));
    policy(s, nexts(a) ) = 1;
    
    v = p.*V(nexts);
    q = exp( -beta*v);
    q = q/sum(q);
    soft_pi(s,nexts) = q;
end
value = -V;
end