function fig7

fsiz = [0.3526    0.5259    .45    0.25*3];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 3;
nc = 3;

h(2:3) = sim_train_default(nr,nc,2:3);
h(4:9) = sim_habit(nr,nc,4:9);

for i = 1
h(i) = subplot(nr,nc,i);
set(h(i),'visible','off');
end
% h(3) = [];

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

for i= 1:length(h)
    set(h(i),'fontsize',fs,'fontname',fn);
    ht(i) = text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

hx2 = ht(2);
set(hx2,'units','inch');
pos2 = get(hx2,'position');

hx = ht(end-1);
set(hx,'units','inch')
pos = get(ht(end-1),'position');

pos2(1) = pos(1);
set(hx2,'position',pos2);
end
