function h = sim_tolman_latent(plt_nr,plt_nc,plt_np)
do_plot = 1;

[mx_tll, tll_labels, TB, AA] = run_tolman_latent;


if ~do_plot
    return;
end
%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 1;
    plt_nc = 2;
    plt_np = [1 2];
    fsiz = [0.1375    0.6907    0.45    0.2139];
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end
%----------------------
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

h(1) = subplot(plt_nr,plt_nc,plt_np(1));
imagesc(TB,'AlphaData',AA);
set(gca,'box','on','ytick',[],'xtick',[]);

h(2) = subplot(plt_nr,plt_nc,plt_np(2));
errorbarKxN(mx_tll',0*mx_tll',tll_labels,struct('colmap',cols,'barwidth',bw));
alpha(alf);
ylabel('Probability','fontsize',fsy);
legend(tll_labels,'fontsize',fsy,'location','north','box','off');
hax=get(gca,'XAxis');
set(hax,'fontsize',fsy);

end

function [mx, tll_labels, TB2, AA, cols, alf] = run_tolman_latent
I = 9;
J = 9;
c0 = 1;
reward = [-5 -5];
blocked = [2:9 11:18 20:27 29:36 47:54 56:63 65:72 74:81];
terminals = [1 73];
cp = 37;

[P0, ~, lij] = core_griding(I,J);
% plot_grids(P0,[],I,J);

[lij0, P, c] = make_blocks(lij,P0,c0,blocked,terminals,reward);
cp = find(lij0(:,1) == cp);
U1 = core_lrl(P,c);

config.add_arrow = 0;
config.add_labels = 0;

figure;
[TB, AA] = plot_grids(U1,config,I,J,[],[],U1,lij0); close;

TB2 = nan([size(TB),3]);
TB2(:,:,1) = TB;
% TB2(:,:,2) = TB*.3;
% TB2(:,:,3) = TB*.3;
aa = (TB==0).*(AA==.5); aa = aa==1;
[it,jt] = find(aa);

cols = def('col');cols(1,:) = [];
alf = .8;
TB2(it(1),jt(1),:) =  cols(1,:);
TB2(it(2),jt(2),:) =  cols(2,:);

AA(AA==1) = 0;
AA(aa) = alf;
AA(AA==.5) = .1;

pcp1 = U1(cp,:);
pcp1 = pcp1(pcp1>0);
%-----

reward(1) = -reward(1);
c(diag(P)==1) = reward;
U2 = core_lrl(P,c);
pcp2 = U2(cp,:);
pcp2 = pcp2(pcp2>0);

mx = [pcp1([1 3]); pcp2([1 3])];
tll_labels = {'Learning','Test'};
end

function [lij, T, c] = make_blocks(lij,T,c,blocked,terminals,goal)
ns = size(T,1);
c = c*ones(ns,1);
c(terminals) = goal;

for i=1:length(terminals)
    T(terminals(i),:) = 0; T(terminals(i),terminals(i)) = 1;
end

T(blocked,:) = [];
T(:,blocked) = [];
A = (T>0)+0.0;
d = diag(sum(A,2));
T = d^-1*A;
c(blocked) = [];
lij(blocked,:) = [];


end