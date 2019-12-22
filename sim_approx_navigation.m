function h = sim_approx_navigation(plt_nr,plt_nc,plt_np)
do_plot = 1;

n = 10;
[maze, co, cl] = run_navigation(n);

if ~do_plot
    return;
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 1;
    plt_nc = 2;
    plt_np = [1 2];
    fsiz = [0.3536    0.6907    0.3    0.2204];    
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

mclr = 'k';
mlw = 1;

subplot(plt_nr,plt_nc,plt_np(1));
show_maze(maze, mlw, mclr);
% title('10x10 maze','fontsize',fsy,'fontname',fn,'fontweight','normal')
h(1) = gca;

%-------------------------------------
mm = max(max([co cl]-1));
mstep = 10;
mm = ceil(mm/mstep)*mstep;
dd = 0:mstep:mm;

h(2) = subplot(plt_nr,plt_nc,plt_np(2));

plot(co-1,cl-1,'k','linewidth',2);
xlabel('Optimal exhaustive search','fontsize',fsy,'fontname',fn);
ylabel('Linear RL approximation','fontsize',fsy,'fontname',fn);
set(gca,'fontsize',fs,'fontname',fn);
xlim(dd([1 end]));
ylim(dd([1 end]));
set(gca,'xtick',dd);
set(gca,'ytick',dd);
end

function [maze, co, cl, dd, to, tl] = run_navigation(n)
nsim = 'all'; snsim = 'all';
q = 1;

bigpipedir = def('bigpipedir');
pipedir = def('pipedir');
fdir = fullfile(bigpipedir,'run_navigation'); makedir(fdir);
fmaze = fullfile(fdir,sprintf('n%d_maze.mat',n));
fname0 = fullfile(fdir,sprintf('n%d_M0.mat',n));

fdir = fullfile(pipedir,'run_navigation'); makedir(fdir);
fname = fullfile(fdir,sprintf('n%d_%s.mat',n,snsim));
ftname= fullfile(fdir,sprintf('n%d_%s_temp.mat',n,snsim));

if ~exist(fmaze,'file')
    rng(0);
    [id, rr, cc, ptr_left, ptr_up, ptr_right, ptr_down] = core_maze(n,n,'v'); %#ok<ASGLU>
    save(fmaze,'id','rr','cc','ptr_left', 'ptr_up', 'ptr_right', 'ptr_down');
end

fm = load(fmaze);
maze = [fm.rr, fm.cc, fm.ptr_left, fm.ptr_up, fm.ptr_right, fm.ptr_down];

P = zeros(size(maze,1),size(maze,1));
for i=1:size(P,1)
    d = maze(i,3:end);
    di = abs(d(d<0));
    P(i,di) = 1/length(di);
end

P0 = P;
q = q*ones(size(P0,1),1);

np = size(P0,1);

if ~exist(fname,'file')
    if ~exist(fname0,'file')
        M0 = (diag(exp(q))-P0)^-1;
        save(fname0,'M0');
    else
        M0 = load(fname0); M0 = M0.M0;
    end

    
    if ~ischar(nsim)
        for i=1:nsim
            ss = randperm(np);
            states(i,1:2) = ss(1:2);    
        end
    else
        states = nchoosek(1:np,2);
%         states = [ones(np-1,1) (2:np)'];
        nsim = size(states,1);
    end
    
    for i=1:nsim
        ss = states(i,:);

        s0 = ss(1);
        st = ss(2);
        P = P0;
        P(st,:) = 0; P(st,st) = 1;

        tic;
        [policy, value, pi] = core_value_iteration(P,q);    
        to(i,1) = toc;
        [c] = core_follow_path(P,policy,q,s0);
        co(i,1) = length(c);        

        tic;
        U = core_new_goal(P0,M0,st,q);
        tl(i,1) = toc;
        [c, path, intens, U1] = core_follow_path(P,U,q,s0);
        cl(i,1) = length(c);


        if mod(i,50)==0
            fprintf('sim %04d\n',i);
            save(ftname,'states','co','cl','to','tl');
        end
    end
    save(fname,'states','co','cl','to','tl');
end

ff = load(fname);
co = ff.co;
cl = ff.cl;

to = ff.to;
tl = ff.tl;
dd = 1:2:(2*n);

end

function show_maze(maze, lw, colr)
rr = maze(:,1);
cc = maze(:,2);
% ptr_left = maze(:,3);
% ptr_up = maze(:,4); 
ptr_right = maze(:,5); 
ptr_down = maze(:,6);

row = max(rr);
col = max(cc);

line([.5,col+.5],[.5,.5],'linewidth',lw,'color',colr) % draw top border
line([.5,col+.5],[row+.5,row+.5],'linewidth',lw,'color',colr) % draw bottom border
line([.5,.5],[0.5,row+.5],'linewidth',lw,'color',colr) % draw left border
line([col+.5,col+.5],[.5,row+.5],'linewidth',lw,'color',colr)  % draw right border

for ii=1:length(ptr_right)
    if ptr_right(ii)>0 % right passage blocked
        line([cc(ii)+.5,cc(ii)+.5],[rr(ii)-.5,rr(ii)+.5],'linewidth',lw,'color',colr);
        hold on
    end
    if ptr_down(ii)>0 % down passage blocked
        line([cc(ii)-.5,cc(ii)+.5],[rr(ii)+.5,rr(ii)+.5],'linewidth',lw,'color',colr);
        hold on
    end
    
end
axis equal
axis([.5,col+.5,.5,row+.5])
% axis off
set(gca,'YDir','reverse')
set(gca,'ytick',[],'xtick',[])
end
