function [policy, value, pii] = core_value_iteration(P,q,beta)
if nargin<3
    beta = 1;
end

ns = size(P,1);

availables = P>0;
terminals = diag(P)==1;
snt = find(~terminals)';

V = rand(ns,1);    
V(terminals) = q(terminals);
tolv = .01; i = 0; dv = 1;
while dv > tolv
    Vpre = V; i = i+1;

    for s=snt
        nexts = (availables(s,:));
        p = P(s,nexts)';
        V(s) = q(s) + min(V(nexts));
    end
    dv = max(abs(V-Vpre));
%     dvv(i) = dv;
end

pii = zeros(ns,ns);
policy = zeros(ns,ns);
for s=1:ns
    nexts = find(availables(s,:));
    p = P(s,nexts)';
    [~,a] = min(p.*V(nexts));
    policy(s, nexts(a) ) = 1;
    
    v = p.*V(nexts);
    q = exp( -beta*v);
    q = q/sum(q);
    pii(s,nexts) = q;
end
value = -V;

end
