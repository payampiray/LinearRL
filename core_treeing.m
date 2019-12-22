function [P, A, xy, paths] = core_treeing(N,K)
% N = 3;
% K = 2;

states = [1 1+(1:K)];
% p{1} = zeros(1,length(states));
% p{1}(1+(1:K)) = 1;
children(1,:) = 1+(1:K);
levels = 1;
parents = nan;

x = 0;
y = 1;
s = 1;
for n = 2:N
    parent = n-1;
    parent_node = find(levels==parent);
    for i = 1:length(parent_node)
        
        dx = -1:(2/(K-1)):1;
        for k=1:K        
            s = s+1;
            child = states(end)+(1:K);
            states = [states child];            
            levels(s) = n;

    %         p{s} = zeros(1,length(states));
    %         p{s}(children) = 1;

            children(s,:) = child;
            parents(s,:) = parent_node(i);

            y(s,:) = n;
            x(s,:) = x(parent_node(i))+dx(k);
        end
    end
end
ns = length(states);
cs = length(levels);
levels = [levels (n+1)*ones(1,ns-cs)];
children = [children; nan(ns-cs,K)];

P = zeros(ns,ns);
for s=1: length(children)
    if ~isnan(children(s,:))
        P(s,children(s,:)) = 1/K;     
    else
        P(s,s) = 1;
    end
end


x = nan(ns,1); x(1) = 0;
y = nan(ns,1); y(1) = N+1;
ypre = y(1);
for n=2:(N+1)
    sn = find(levels==n);
    x0 = ceil(-length(sn)/2);
    xs = (x0):(-x0);
    if mod(length(sn),2)==0, xs(xs==0) = []; end
    
    y(sn) = ypre - 1; ypre = ypre - 1;
    x(sn) = xs;
end

A = P;
A = ((A + A')>0) + 0.0;
xy = [x y];


terminals = find(diag(P)==1);
paths = nan(length(terminals),N+1);
for i=1:length(terminals)
    s = terminals(i);
    st= nan(1,N+1);
    st(end) = s;
    for n=N:-1:1
        [s,~] = find(children==s);    
        st(n) = s;
    end
    paths(i,:) = st;
end
end