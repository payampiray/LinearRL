function [P, k, L, D] = core_train_default(s,P,cost,alpha,L,D,learn_p)
if nargin<5
    learn_p = 0;
    L = diag(exp(cost))-P;
    D = L^-1;
end

ss = s;
ok = 1;
k = 0;
t = diag(P)==1;
M = D;
M(t,:) = [];
M(:,t) = [];    
while ok    
    k = k+1;
    U = core_lrl(P,cost,M);
    
    u = U(s,U(s,:)>0);
    [a, r]=randp(u);
    suc = find(U(s,:));
    next = suc(a);
    
    if learn_p
        P = learn_P(s,next,P,alpha);
        [~, ~, D, L] = core_woodbury(L,D,P,cost);
        M = D;
        M(t,:) = [];
        M(:,t) = [];        
    end
    
    s = next;
    ss = [ss s];
    if P(s,s)==1
        ok = 0;
    end
end
end


function P = learn_P(cur,next,P,alpha)
ps = P(cur,:);
ps(next) = ps(next) + alpha*(1-ps(next));
ps = ps/sum(ps);

P(cur,:) = ps;
end

function [a, r]=randp(u)
    csu = [0 cumsum(u,2)];
    r = rand;
    a = find(r>csu,1,'last');
end
