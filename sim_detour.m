function h = sim_detour(plt_nr,plt_nc,plt_np)
do_plot = 1;

[P, cost, barrierrs, cp, scp, xb, yb, TB, AA, xt, yt] = detour;
M = (diag(exp(cost))-P)^-1;

[U1] = core_lrl(P,cost);
p1 = U1(cp,scp);

[U2] = blocking(P,barrierrs,M,cost);

p2 = U2(cp,scp);

mx = [p1; p2]';
directions = {'left','straight','right'};
labels = {'Training','Test'};

if ~do_plot
    return;
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 1;
    plt_nc = 3;
    plt_np = [1 2 3];
    fsiz = [0.1375    0.6907    0.6    0.2139];
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
title(labels{1},'fontsize',fsy);
% text(xt,yt,'Food','fontsize',fs,'HorizontalAlignment','center');

h(2) = subplot(plt_nr,plt_nc,plt_np(2));
imagesc(TB,'AlphaData',AA);
set(gca,'box','on','ytick',[],'xtick',[]);
hold on;
plot(xb,yb,'linewidth',6,'color','k');
title(labels{2},'fontsize',fsy);
% text(xt,yt,'Food','fontsize',fs,'HorizontalAlignment','center');

h(3) = subplot(plt_nr,plt_nc,plt_np(3));
errorbarKxN(mx,0*mx,labels,struct('barwidth',bw,'colmap',cols));
alpha(alf);
set(h,'fontname',fn);
ylabel('Probability','fontsize',fsy);
hax=get(gca,'XAxis');
set(hax,'fontsize',fsy);
legend(directions,'location','north','fontsize',fsy,'box','off');

end

function [U, z, P] = blocking(P0,ss,M0,q)
P = P0;
A = P>0;
A(ss(1),ss(2))=0;
D = sum(A,2);
P = diag(D.^-1)*A;

st = ss(1);
dp = P0(st,:) - P(st,:);
e  = zeros(size(P0,1),1); e(st)=1;

terms = diag(P)==1;
p = P(:,terms); p(terms)=0;
m0 = M0*e;

z0 = M0*p;

alpha = (dp*z0)/(1+dp*m0);
zNT = z0(~terms) - alpha*m0(~terms);

z = nan(size(P,1),1);
z(~terms) = zNT*exp(-q(terms));
z(terms) = exp(-q(terms));

G = P*z;
zg = z'./G;
U = P.*zg;

end

function [P, cost, barriers, cp, scp, xb, yb, TB2, AA, jt, it, scp0] = detour
I = 9;
J = 9;
c0 = 1;
reward = -5;
barriers0 = [33 32];
terminals = [28];
cp = 34; % critical point
scp0 = [25 33 43];


bb = [];
n = 9;
for i = [2 3 5 6 7 8]
    b = [1 2 4 5 6 8 9]+(i-1)*n;
    bb = [bb b];
end
blocked = bb;
blocked = [1 2 8 9 73 74 80 81 blocked];


[P0, ~, lij] = core_griding(I,J);
% plot_grids(P0,[],I,J); close;

[lij0, P, cost] = make_blocks(lij,P0,c0,blocked,terminals,reward);
cp = find(lij0(:,1) == cp);
scp = nan(size(scp0));
for i=1:length(scp0)
    scp(i) = find(lij0(:,1) == scp0(i));
end
% terminals = find(lij0(:,1) == terminals);
U1 = core_lrl(P,cost);

config.add_arrow = 0;
config.add_labels = 0;

figure;
[TB, AA] = plot_grids(U1,config,I,J,[],[],U1,lij0); close;

TB2 = nan([size(TB),3]);
TB2(:,:,1) = TB;
% TB2(:,:,2) = TB*.3;
% TB2(:,:,3) = TB*.3;
aa = (TB==0).*(AA==.5); aa = aa==1;

[it,jt] = ind2sub([I,J],terminals);
[ib,jb] = ind2sub([I,J],barriers0);

cols = def('col');
col  = cols(3,:);

alf = .8;
TB2(it(1),jt(1),:) =  col;

AA = zeros(size(TB));
AA(TB==1) = .1;
AA(terminals) = alf;

xb = ib-1.5+[0 0];
yb = J - jb +.5;

barriers = nan(size(barriers0));
for i=1:length(barriers0)
    barriers(i) = find(lij0(:,1) == barriers0(i));
end

end

function [lij, P, c] = make_blocks(lij,P,c,blocked,terminals,goal)
ns = size(P,1);
c = c*ones(ns,1);
c(terminals) = goal;

for i=1:length(terminals)
    P(terminals(i),:) = 0; P(terminals(i),terminals(i)) = 1;
end

P(blocked,:) = [];
P(:,blocked) = [];
A = (P>0)+0.0;
D = diag(sum(A,2));
P = D^-1*A;
c(blocked) = [];
lij(blocked,:) = [];


end