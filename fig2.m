function fig2

def('addpath');

fsiz = [0.3526    0.5259    .3    0.2630*3];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 3;
nc = 2;

sim_DR_tree(nr,nc,2);
h(3:4) = sim_DR_maze(nr,nc,3:4);
h(5:6) = sim_approx_navigation(nr,nc,5:6);

for i = 1
h(i) = subplot(nr,nc,i);
set(h(i),'visible','off');
end
h([2 6]) = [];

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

for i= 1:length(h)
    set(h(i),'fontsize',fs,'fontname',fn);
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

end

