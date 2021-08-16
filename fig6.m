function fig6

% fsiz = [0.3526    0.5259    .35    0.45];
fsiz = [0.3526    0.5259    .55    0.2];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 1;
nc = 3;

h(2) = sim_2step(nr,nc,2);

for i = [1 3]
h(i) = subplot(nr,nc,i);
set(h(i),'visible','off');
end

% --------
fs = def('fs');
fn = def('fn');
fsy = def('fsy');
alf = def('alf');
fsA = def('fsA');
xsA = -.05 + def('xsA');
ysA = 1.03;def('ysA');
abc = def('abc');
bw  = .15;
cols = def('col');
cols = cols([3 2 1],:);

for i= 1:length(h)
    set(h((i)),'fontsize',fs,'fontname',fn);
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

end
