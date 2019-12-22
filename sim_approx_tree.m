function h = sim_approx_tree(plt_nr,plt_nc,plt_np)
do_plot = 1;

[P, ~, xy, paths] = core_treeing(7,2);

n = 7;
pipedir = def('pipedir');
fname = fullfile(pipedir,'sim_approx_tree',sprintf('sim_n%d.mat',n));

if ~exist(fname,'file')
    nsim = 100;
    error_cost = nan(nsim,n);
    error_rank = nan(nsim,n);
    for i=1:nsim
        [e_cost, e_rank, labels] = run_approximation(n);
        error_cost(i,:) = e_cost;
        error_rank(i,:) = e_rank;
        fprintf('sim %03d\n',i);
    end
    sims = struct('n',n,'error_cost',error_cost,'error_rank',error_rank,'labels',{labels}); %#ok<NASGU>
    save(fname,'-struct','sims');
end
sims = load(fname); %sims = fsim.sims;
error_cost = sims.error_cost;
error_rank = sims.error_rank;
labels = sims.labels;
nsim = size(error_cost,1);

for i=1:length(labels)-1
    labels{i} = sprintf('D%d',i);
end
labels{end} = 'LRL';

mcost = mean(error_cost);
mrank = mean(error_rank);
ecost = serr(error_cost);
erank = serr(error_rank);

mx1 = mcost;
mx2 = mrank;
ex1 = ecost;
ex2 = erank;

if ~do_plot
    return;
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 1;
    plt_nc = 1;
    plt_np = 1;
    fsiz = [0.3536    0.6907    0.2    0.2204];    
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
bw  = .3;
cols = def('col'); cols(1,:) = [];
col  = [.5 .5 .5];

h(1) = subplot(plt_nr,plt_nc,plt_np(1));
errorbarKxN(mx1,ex1,labels,struct('colmap',col,'barwidth',bw));
ylabel('Cost','fontsize',fsy);
xlabel('Models','fontsize',fsy);
alpha(alf);
hax=get(gca,'XAxis');
set(hax,'fontsize',fsy);
end

function [error_cost, error_path, labels] = run_approximation(n)

k = 2;
[P, ~, xy, paths] = core_treeing(n,k);

ns = length(P);

q = rand(ns,1);
q = 10*q;

% optimal path
for i=1:size(paths,1)
    for j=1:size(paths,2)
        costs(i,j) = q(paths(i,j));
    end
end
[costsorted,idx]=sort(sum(costs,2));

[U] = core_lmdp(P,q);
[cost_lmdp, path_lmdp] = core_follow_path(P,U,q,1);

beta = 10;
% [pi1] = core_value_iteration(P,q,beta);
% [cost_val, path_val] = core_follow_path(P,pi1,q,1);

for i=1:(n-1)
    pi = core_valuation_tree(P,q,n,i,beta);
    [cost_val, path_val] = core_follow_path(P,pi,q,1);
    costpar(i) = sum(cost_val);
    
    labels{i} = sprintf('Depth %d',i);    
end
% labels{n} = sprintf('Full MB');


cost = [costpar sum(cost_lmdp)];
goodness = nan(size(cost));
for i=1:length(cost)
    [j] = find(cost(i)==costsorted);
    if length(j)>1, j = j(1); end
    goodness(i) = j;
end

error_cost = cost - min(costsorted);
error_path = goodness-1;
labels = [labels 'linear RL'];
end
