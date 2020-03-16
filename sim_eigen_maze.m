function h = sim_eigen_maze(plt_nr,plt_nc,plt_np)
do_plot = 1;

eignums = [15 20 32];
A = run(eignums);

if ~do_plot
    return;
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 1;
    plt_nc = 3;
    plt_np = 1:3;
    fsiz = [0.3536    0.6907    0.6    0.2204];    
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end
%----------------------

for i= 1:length(A)
    h(i) = subplot(plt_nr,plt_nc,plt_np(i));
    imagesc(A{i});
    set(gca,'xtick',[],'ytick',[],'box','on');    
end

end

function A = run(eignums)
n = 50;
I = n;
J = n;

[T, xy] = core_griding(I,J);
c= .1*ones(size(T,1),1);

D = (diag(exp(c))-T)^-1;

i = max(eignums);
[U,~] = eigs(D,i);
idx = sub2ind([n n],xy(:,1),xy(:,2));

A = cell(1,length(eignums));
for i=1:length(eignums)
    Ap = zeros(n,n);
    Ap(idx) = U(:,eignums(i));
    A{i} = Ap;
end

end
