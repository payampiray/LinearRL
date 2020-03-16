function fig5
% border cells

fsiz = [0.3526    0.5259    .25    0.2630*3];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 4;
nc = 2;
h(5:8) = sim_border_cells(nr,nc,5:8);


for i = [1 2]
h(i) = subplot(nr,nc,i);
set(h(i),'visible','off');
end

h([2:4 6:8]) = [];

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
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

end
