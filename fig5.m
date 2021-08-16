function fig5
% border cells
def('addpath');

fsiz = [0.3526    0.5259    .35    0.25];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 2;
nc = 4;
h(5:8) = sim_border_cells(nr,nc,[3 4 7 8]);


for i = [1 2]
h(i) = subplot(nr,nc,i);
set(h(i),'visible','off');
end

h([2:4 6:8]) = [];

% --------
fn = def('fn');
fsA = def('fsA');
ysA = def('ysA');
abc = def('abc');
cols = def('col');

xsA = [-.7 -0.2];
for i= 1:length(h)
    text(xsA(i),ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

end
