function fig4
% grid fields
def('addpath');

fsiz = [0.3526    0.5259    .45    .9];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 3;
nc = 3;

h(4) = bar_derdikman(nr,nc,4);
h(7:9) = sim_eigen_maze(nr,nc,[7 8 9]);

for i = [1 5]
h(i) = subplot(nr,nc,i);
set(h(i),'visible','off');
end

h([2 3 6 8 9]) = [];

% --------
fs = def('fs');
fn = def('fn');
fsy = def('fsy');
alf = def('alf');
fsA = def('fsA');
xsA = -.05 + def('xsA');
ysA = def('ysA');
abc = def('abc');
bw  = .15;
cols = def('col');
cols = cols([3 2 1],:);

for i= 2:length(h)
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

ys1 = .5;
text(xsA,ys1,abc(1),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(1));
end

function h = bar_derdikman(plt_nr,plt_nc,plt_np)


mx = [0.008659793814433062, 0.49360824742268034, 0.578762886597938];
ex = [0.018762886597938167, 0.5239175257731958, 0.604742268041237] - mx;
labels = {'Hairpin','Virtual Hairpin','Open field 2'};

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 1;
    plt_nc = 1;
    plt_np = 1;
    fsiz = [0.2    0.4    0.1    0.4];
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end


fsy = def('fsy');
bw = .2;

h = subplot(plt_nr,plt_nc,plt_np(1));
hp = errorbarKxN(mx,ex,labels,struct('barwidth',bw,'colmap',[1 1 1]*.5));
hax=get(gca,'XAxis');
set(hax,'fontsize',fsy,'TickLabelRotation', 20);
set(hp,'ylim',[0 .7]);
ylabel(sprintf('Correlation with\n open field 1'),'fontsize',fsy);
end