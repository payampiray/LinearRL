function h = sim_train_default(plt_nr,plt_nc,plt_np)
do_plot = 1;

[u,p,lambdas,labels,xt]=run;

if ~do_plot
    return;
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 1;
    plt_nc = 2;
    plt_np = [1 2];
    fsiz = [0.2    0.4    0.4    0.25];
    figure; set(gcf,'units','normalized'); set(gcf,'position',fsiz);
end

col = def('col');
fs = def('fs');
fn = def('fn');
fsy = def('fsy');
alf = def('alf');
fsA = def('fsA');
xsA = -.05 + def('xsA');
ysA = def('ysA');
abc = def('abc');
bw = .25;
col = [col(3,:); col(2,:)];

mxx = {p,u};
locs = {'northwest','southwest'};
for j=1:2
    mx = [p(:,j) u(:,j)];
    mx = mxx{j};
    
    h(j) = subplot(plt_nr,plt_nc,plt_np(j));
    for i=1:size(mx,2)
        plot(lambdas,mx(:,i),'-x','color',[col(i,:) 1],'linewidth',2); hold on;
    end
    ylim([0,1.05]);    
    xlim([0,1.1*max(lambdas)]);
    xtick = [0 lambdas];
    set(gca,'xtick',xtick);
    
    if j==2
        plot([0,1.1*max(lambdas)],.5*[1 1],'k');
    end
    xlabel('control cost parameter (\lambda)','Interpreter','tex','fontsize',fsy);
    ylabel('Probability of A','fontsize',fsy);
    legend(xt,'location',locs{j},'fontsize',fsy,'box','off');
    title(labels{j},'fontsize',fsy);
end

% for j=1:2
%     mx = [p(:,j) u(:,j)];
%     mx(mx==0.5) = 0.5001;
%     h(j) = subplot(plt_nr,plt_nc,plt_np(j));
%     errorbarKxN(mx',0*mx',lambdas,struct('colmap',col,'barwidth',bw,'basevalue',0));
%     alpha(alf);
%     set(h(j),'fontname',fn,'fontsize',fs);
%     ylabel('Probability','fontsize',fsy);
%     ylim([0 1.05]);
%     % hax=get(gca,'XAxis');
%     % set(hax,'fontsize',fsy);
%     xlabel('control cost parameter (\lambda)','Interpreter','tex','fontsize',fsy);
%     if j==2
%         hlg = legend(labels,'location','northeast','fontsize',fsy,'box','off','AutoUpdate','off');
%         x = get(gca,'xlim');
%         plot(x,.5*[1 1],'k');
%     end
%     title(xt{j});
% end

end

function [u,p,lambdas,labels,xtitles]=run
n = 2;
k = 2;

ntrain = 1000;
alpha = 0.01;

lambdas = .5:.5:4;

Pb = [0 .5 .5;0 1 0;0 0 1];


ns = size(Pb,1);
rbs = zeros(ns,1);
rbs(2) = 5;

rs = zeros(ns,1);
rs(3) = 5;

labels = {'Default policy','Decision policy'};
xtitles = {'no training','over-trained'};

pipedir = def('pipedir'); 
makedir(fullfile(pipedir,mfilename));
fname = fullfile(pipedir,mfilename,sprintf('sim.mat'));
ok = exist(fname,'file');

if ~ok
    for i=1:length(lambdas)        
        rng(0);

        lambda = lambdas(i);
        P0 = Pb;
        
        P =run_base(P0,-rbs,lambda,ntrain,alpha);
        p(i,:) = [P0(1,3) P(1,3)];

        q = -rs;
        [U0] = core_lrl(P0,q,[],lambda);
        [c0(i,:), path0(i,:)] = core_follow_path(P0,U0,q,1);
        u(i,1) = U0(1,3);

        [U] = core_lrl(P,q,[],lambda);
        [c(i,:), path(i,:)] = core_follow_path(P,U,q,1);    
        u(i,2) = U(1,3);        
    end
    save(fname,'u','p','lambdas','labels','xtitles');
end

x = load(fname);
u = x.u;
p = x.p;
lambdas = x.lambdas;
labels = x.labels;
xtitles = x.xtitles;

end

function [P,D0,D]=run_base(P,q,lambda,nsim,alpha)

q(q==0) = 10^-10;
q = q/lambda;

L = diag(exp(q))-P;
D = L^-1;

% terminals = abs(diag(P)-1)<eps;
% D(1:size(D,1):numel(D)) = nan;
% M = L(~terminals,~terminals)^-1;
% D(~terminals,~terminals) = M;
% D(terminals,terminals) = diag(nan(sum(terminals),1));

D0 = D;
for i=1:nsim
    [P, ~, L, D] = core_train_default(1,P,q,alpha,L,D,1);
end

terms = diag(P)==1;
D0 (terms,:) = [];
D0 (:,terms) = [];
D (terms,:) = [];
D (:,terms) = [];
end
