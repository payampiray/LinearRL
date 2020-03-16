function fig3
def('addpath');

fsiz = [0.3526    0.5259    .6    0.2630*3];

close all;
figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);

nr = 3;
nc = 3;

h(2:3) = sim_tolman_latent(nr,nc,[2 3]);
h(6) = sim_policy_revaluation(nr,nc,6); 
h(7:9) = sim_detour(nr,nc,[7 8 9]);

for i = [1 4]
h(i) = subplot(nr,nc,i);
set(h(i),'visible','off');
end

h(5) = [];

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
    set(h((i)),'fontsize',fs,'fontname',fn);
    text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
end

end
