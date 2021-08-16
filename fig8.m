function fig8
def('addpath');

fsiz = [0.3526    0.5259    .35    0.25*3];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 3;
nc = 2;

h(3:4) = sim_control_cost(nr,nc,3:4);
h(5:6) = sim_control_Krebs(nr,nc,5:6);

for i = [1 2]
h(i) = subplot(nr,nc,i);
set(h(i),'visible','off');
end
h([3 4]) = [];

% --------
fs = def('fs');
fn = def('fn');
fsA = def('fsA');
abc = def('abc');

xsA = -.05;
ysA = [.7 .7 1.1 1.1];
for i= 1:length(h)
    set(h(i),'fontsize',fs,'fontname',fn);
    text(xsA,ysA(i),abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

end