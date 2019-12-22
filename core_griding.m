function [P, xy, lij] = core_griding(I,J)
if nargin<2, J = I; end

ns = I*J;
P = zeros(ns,ns);
A = zeros(ns,ns);
xy = zeros(ns,2);
for i=1:I
    for j=1:J
        k = sub2ind([I J],i,j);
        IJ = [i-1 j;
              i j-1;
%               i j;
              i j+1;
              i+1 j];
        As = [1;2;3;4];        
        next = nansub2ind([I J],IJ);
        nonan = ~isnan(next);
        next = next(nonan);        
        As = As(nonan);
        P(k, next) = 1/length(next);
        A(k,next) = As;
        
        xy(k,:) = [j i];
    end
end

idx = (1:ns)';
x = xy(:,1);
y = xy(:,2);

lij = [idx y x]; % label, row, column
end

function neigbs = nansub2ind(n,IJ)
neigbs = nan(1,size(IJ,1));
for k=1:size(IJ,1)
    inan = any(IJ(k,:)<1 | IJ(k,:)>n);
    if ~inan
        neigbs(k) = sub2ind([n n],IJ(k,1),IJ(k,2));
    end
    
end

end