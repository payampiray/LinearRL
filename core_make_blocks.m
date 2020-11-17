
function [lij, P, c] = core_make_blocks(lij,P,c,blocked,terminals,goal)
ns = size(P,1);
c = c*ones(ns,1);
c(terminals) = goal;

for i=1:length(terminals)
    P(terminals(i),:) = 0; P(terminals(i),terminals(i)) = 1;
end

P(blocked,:) = [];
P(:,blocked) = [];
A = (P>0)+0.0;
D = diag(sum(A,2));
P = D^-1*A;
c(blocked) = [];
lij(blocked,:) = [];


end