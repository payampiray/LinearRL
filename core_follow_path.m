function [costs, path, intens, U1] = core_follow_path(P,U,q,s0)

termianls = find(diag(P)==1);

s = s0;
path = s0;
costs = q(s0);
intens = [];
U1 = zeros(size(U));
while ~any(termianls==s)
    [us,next] = max(U(s,:));
    s = next;
    u = us(1);
    
    U1(s,next) = 1;    
    costs = [costs q(s)];
    path = [path s];
    intens = [intens u];
end
intens = [intens 1];

end