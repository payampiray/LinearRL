function h = sim_approx_tree(plt_nr,plt_nc,plt_np)
do_plot = 1;

[P, ~, xy, paths] = core_treeing(7,2); %#ok<ASGLU>

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
sims = load(fname);
error_cost = sims.error_cost;
error_rank = sims.error_rank;
labels = sims.labels;

for i=1:length(labels)-1
    labels{i} = sprintf('D%d',i);
end
labels{end} = 'LRL';

mcost = mean(error_cost);
mrank = mean(error_rank);
ecost = serr(error_cost);
erank = serr(error_rank);

x1 = error_cost;
mx1 = mcost;
ex1 = ecost;

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
fsy = def('fsy');
alf = def('alf');

h(1) = subplot(plt_nr,plt_nc,plt_np(1));
raincloud1xN(x1,mx1,ex1,labels);
ylim([0,12]);

ylabel('Cost','fontsize',fsy);
xlabel('Models','fontsize',fsy);
alpha(alf);
hax=get(gca,'XAxis');
set(hax,'fontsize',fsy);
end

function [error_cost, error_path, labels] = run_approximation(n)

k = 2;
[T, ~, ~, paths] = core_treeing(n,k);

ns = length(T);

q = rand(ns,1);
q = 10*q;

% optimal path
for i=1:size(paths,1)
    for j=1:size(paths,2)
        costs(i,j) = q(paths(i,j));
    end
end
[costsorted]=sort(sum(costs,2));

[U] = core_lrl(T,q);
[cost_lrl] = core_follow_path(T,U,q,1);

beta = 10;
costpar = nan(1,n-1);
labels = cell(1,n-1);
for i=1:(n-1)
    pi = core_valuation_tree(T,q,n,i,beta);
    [cost_val, ~] = core_follow_path(T,pi,q,1);
    costpar(i) = sum(cost_val);
    
    labels{i} = sprintf('Depth %d',i);    
end

cost = [costpar sum(cost_lrl)];
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
