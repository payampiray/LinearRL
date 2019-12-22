function h = sim_control_cost(plt_nr,plt_nc,plt_np)
do_plot = 1;

P = [.5;.2];
u = .01:.01:.99;
for i=1:size(P,1)
    p = P(i,:);
    Elogu = log(u).*u + log(1-u).*(1-u);
    Elogp = log(p).*u + log(1-p).*(1-u);
    cc(i,:) = Elogu - Elogp;
end
P = [P 1-P];
lbls = {'A','B'};

if ~do_plot
    return;
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 1;
    plt_nc = 2;
    plt_np = 1:2;
    fsiz = [0.3536    0.6907    0.4    0.2204];    
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
bw  = .15;
cols = def('col');

ccmax = max(max(cc));
for i=1:plt_nc
    h(i) = subplot(plt_nr,plt_nc,plt_np(i) );
    plot(u,cc(i,:),'linewidth',2,'color','k');
%     set(h,'xtick',[0 p(1) 1]);
    ylim([0 ccmax]);
    set(h,'fontsize',fs,'fontname',fn);
    xlabel('Probability of choosing A');
    if i==1        
        ylabel(sprintf('Control\ncost'),'fontsize',fsy);    
    end
end

end