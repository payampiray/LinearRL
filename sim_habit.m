function h = sim_habit(plt_nr,plt_nc,plt_np)
do_plot = 1;

[siz,P0,lij0,U1_trained,U2_no_train,U2_trained,paths1,paths2,minpath1,minpath2,path2_no_train,path2_trained] = habit_problem;

if ~do_plot
    return;
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 2;
    plt_nc = 3;
    plt_np = 1:6;
    fsiz = [0.2    0.4    0.45    0.45];
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

% config = struct('pos0',.3,'arrow_shift0',.3,'arrow_length0',.6,'colmap','summer','alpha',.4,'col_terminal',[.2 .4 1],'colorful',0);
config = struct('pos0',.3,'arrow_shift0',.3,'arrow_length0',.6,'colmap','summer','alpha',.4,'col_terminal',[.2 .4 1],'colorful',1);
config.labels = repmat('',size(lij0,1),1);
config.labels = cell(size(lij0,1),1);
config.labels{1} = 'H';
config.add_labels = 1;
config.str_label = '%s';

h(1)=subplot(plt_nr,plt_nc,plt_np(1));
plot_grids(P0,config,siz,siz,[],[],lij0); %close;

config.col_terminal = cols(1,:);
h(2)=subplot(plt_nr,plt_nc,plt_np(2));
plot_grids(U1_trained,config,siz,siz,[],[],lij0); %close;

mx1 = mean(paths1);
ex1 = serr(paths1);
mx2 = mean(paths2);
ex2 = serr(paths2);
bv1 = minpath1;
bv2 = minpath2;
labels = {'no training','over-trained'};

h(3) = subplot(plt_nr,plt_nc,plt_np(3));
errorbarKxN(mx1,ex1,labels,struct('barwidth',bw,'colmap',cols,'basevalue',bv1));
alpha(gca,alf);
set(h,'fontname',fn);
ylabel('Number of steps','fontsize',fsy);
hax=get(gca,'XAxis');
set(hax,'fontsize',fsy);
title('Same room training','fontsize',fsy);
% ----------------
config.add_arrow = 1;
config.path = path2_no_train;
h(4) = subplot(plt_nr,plt_nc,plt_np(4));
plot_grids(U2_no_train,config,siz,siz,[],[],lij0); %close;
title('no training','fontsize',fsy);

config.path = path2_trained;
h(5) = subplot(plt_nr,plt_nc,plt_np(5));
plot_grids(U2_trained,config,siz,siz,[],[],lij0); %close;
title('over-trained','fontsize',fsy);

h(6) = subplot(plt_nr,plt_nc,plt_np(6));
errorbarKxN(mx2,ex2,labels,struct('barwidth',bw,'colmap',cols,'basevalue',bv2));
alpha(gca,alf);
set(h(6),'fontname',fn);
ylabel('Number of steps','fontsize',fsy);
hax=get(gca,'XAxis');
set(hax,'fontsize',fsy);
title('Different room training','fontsize',fsy);



% % for i= 1:length(h)
% %     text(xsA,ysA,abc(i),'fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(i));
% % end
end

function [siz,P0,lij0,U1_trained,U2_no_train,U2_trained,paths1,paths2,minpath1,minpath2,path2_no_train,path2_trained] = habit_problem
pipedir = def('pipedir');
fdir = fullfile(pipedir,mfilename); makedir(fdir);
fname = fullfile(fdir,sprintf('sim.mat'));
if exist(fname,'file')
    load(fname);
    return;
end

[P, cost_step, reward, scp, lij0] = define;
siz = sqrt(max(lij0(:,1)));
alpha_P = .01;
nsim = 1000;
nsim_train = 1000;

P0 = P;
cost = cost_step*ones(size(P,1),1); cost(diag(P)==1)=reward;
L0 = diag(exp(cost))-P0;
D0 = L0^-1;

L = L0;
D = D0;
s = 1;
rng(0);
for i=1:nsim_train
    [P, ~, L, D] = train_episode(s,P,cost,alpha_P,L,D,1);
end

[P1_no_train, cost] = new_terminal(P0,scp(2),cost_step,reward);
[P1_trained] = new_terminal(P,scp(2),cost_step,reward);
for i=1:nsim    
    [~, paths1(i,1)] = train_episode(s,P1_no_train,cost,alpha_P);
    [~, paths1(i,2)] = train_episode(s,P1_trained,cost,alpha_P);
end
% U1_no_train = core_lrl(P1_no_train,cost);
U1_trained = core_lrl(P1_trained,cost);
% [~, path1_no_train] = core_follow_path(P1_no_train,U1_no_train,cost,s);
[~, path1_trained] = core_follow_path(P1_trained,U1_trained,cost,s);

% goal in a new room
[P2_no_train, cost] = new_terminal(P0,scp(3),cost_step,reward);
[P2_trained] = new_terminal(P,scp(3),cost_step,reward);
for i=1:nsim
    [~, paths2(i,1)] = train_episode(s,P2_no_train,cost,alpha_P);
    [~, paths2(i,2)] = train_episode(s,P2_trained,cost,alpha_P);
end
U2_no_train = core_lrl(P2_no_train,cost);
U2_trained = core_lrl(P2_trained,cost);
[~, path2_no_train] = core_follow_path(P2_no_train,U2_no_train,cost,s);
[~, path2_trained] = core_follow_path(P2_trained,U2_trained,cost,s);

minpath1 = length(path1_trained)-1;
minpath2 = length(path2_no_train)-1;

save(fname);

end

% ---------------------------------------------------
function [P, cost, L, D] = new_terminal(P,terminal,cost_step,reward)
t = find(diag(P)==1);
neighbours = find(P(:,t));
P(t,t) = 0;
P(t,neighbours) = 1/length(neighbours);

P(terminal,:) = 0;
P(terminal,terminal) = 1;

cost = cost_step*ones(size(P,1),1); cost(terminal)=reward;

L = diag(exp(cost))-P;
D = L^-1;
end

function [P, k, L, D] = train_episode(s,P,cost,alpha,L,D,learn_p)
if nargin<5
    learn_p = 0;
    L = diag(exp(cost))-P;
    D = L^-1;
end

ss = s;
ok = 1;
k = 0;
t = diag(P)==1;
M = D;
M(t,:) = [];
M(:,t) = [];    
while ok    
    k = k+1;
    U = core_lrl(P,cost,M);
    
    u = U(s,U(s,:)>0);
    [a, r]=randp(u);
    suc = find(U(s,:));
    next = suc(a);
    
    if learn_p
        P = learn_P(s,next,P,alpha);
        [~, ~, D, L] = core_woodbury(L,D,P,cost);
        M = D;
        M(t,:) = [];
        M(:,t) = [];        
    end
    
    s = next;
    ss = [ss s];
    if P(s,s)==1
        ok = 0;
    end
end
end

function P = learn_P(cur,next,P,alpha)
ps = P(cur,:);
ps(next) = ps(next) + alpha*(1-ps(next));
ps = ps/sum(ps);

P(cur,:) = ps;
end

function [a, r]=randp(u)
    csu = [0 cumsum(u,2)];
    r = rand;
    a = find(r>csu,1,'last');
end
% ---------------------------------------------------

function [P, c0, reward, terminals, lij, room] = define
close all;
I = 11;
J = 11;

c0 = .5;
reward = -1;
terminals0 = [71 115 75];

blocked = [56 57 59:63 65 66, 6 17 39 50 72 83 105 116];

[P0, ~, lij] = core_griding(I,J);
% plot_grids(P0,[],I,J); 
% close;

[lij, P] = core_make_blocks(lij,P0,c0,blocked,terminals0(1),reward);
terminals = nan(size(terminals0));
for i=1:length(terminals0)
    terminals(i) = find(lij(:,1) == terminals0(i));
end

room = nan(size(lij,1),1);
for i=1:size(lij,1)
    if lij(i,2)<=5 && lij(i,3)<=5
        room(i) = 1;
    elseif lij(i,2)<=5 && lij(i,3)>6
        room(i) = 2;
    elseif lij(i,2)>6 && lij(i,3)<=5
        room(i) = 3;
    elseif lij(i,2)>6 && lij(i,3)>6
        room(i) = 4;
    end
end

config.add_arrow = 0;
config.add_labels = 1;
config.str_label = '%d';

config.labels = lij(:,1);
% config.labels = (1:length(P))';

figure;
[TB, AA] = plot_grids(P,config,I,J,[],[],P,lij); 
close;

end