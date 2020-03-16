function h = sim_border_cells(plt_nr,plt_nc,plt_np)
do_plot = 1;

eignums = [1 10 6 12];
A = run(eignums);

if ~do_plot
    return;
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 2;
    plt_nc = 2;
    plt_np = 1:4;
    fsiz = [0.3536    0.6907    0.4    0.4];    
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end
%----------------------

for i= 1:length(A)
    h(i) = subplot(plt_nr,plt_nc,plt_np(i));
    imagesc(A{i}); hold on;
    set(gca,'xtick',[],'ytick',[],'box','on');
end

end

function As = run(eignums)
n = 20; type = 1;

pipedir = def('pipedir');
fdir = fullfile(pipedir,mfilename); makedir(fdir);
fname = fullfile(fdir,sprintf('sim%d_n%d.mat',type,n));

[T,xy] = core_griding(n);
c = .1*ones(size(T,1),1);

sb1 = [0:n:(n^2-n)]+n-1;
sb2 = [0:n:(n^2-n)]+n-2;

P0 = T;
A = (T>0)+0;
for i=1:length(sb1)
    s1 = sb1(i);
    s2 = sb2(i);
    A(s1,s2) = 0;
end
D = sum(A,2);
T = diag(D.^-1)*A;


%----
L0 = diag(exp(c))-P0;
L  = diag(exp(c))-T;

ds = find(sum(abs(L - L0),2)~=0);
D0 = L0^-1;

m0 = D0(:,ds);
d = L(ds,:) - L0(ds,:);
alpha = (eye(length(ds))+d*m0)^-1;

if ~exist(fname,'file')
    [UD0,ED0]=eig(D0);
    UD0 = UD0(:,1:50);
    ED0 = ED0(:,1:50);
    save(fname,'UD0','ED0','D0');
else
    f = load(fname);
    UD0 = f.UD0;
end

% UD0 = real(UD0);

Dx =   m0*alpha*d*UD0(:,eignums);

idx = sub2ind([n n],xy(:,1),xy(:,2));
As = cell(1,length(eignums));
for i=1:length(eignums)
    u = Dx(:,i);
    if mean(u>0)<.5, u = -u; end
    
    A = zeros(n,n);
    A(idx) = u;    
    As{i} = A;
        
end
end
