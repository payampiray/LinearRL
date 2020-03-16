function sim_DR_tree(plt_nr,plt_nc,plt_np)
do_plot = 1;

M = run_DR;


if ~do_plot
    return;
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 1;
    plt_nc = 1;
    plt_np = 1;
    fsiz = [0.3536    0.6907    0.13    0.2204];    
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

subplot(plt_nr,plt_nc,plt_np(1));
hm = heatmap(M/1000);
set(hm,'ColorbarVisible','off','fontsize',fsy);
title('DR');
end

function M = run_DR
N = 3; K = 2;
[T] = core_treeing(N,K);
c = ones(size(T,1),1);
c(8) = -.5;
c(9) = -.5;

[~,~,M] = core_lrl(T,c);

end
