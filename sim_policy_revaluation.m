function h = sim_policy_revaluation(plt_nr,plt_nc,plt_np)
do_plot = 1;

[mx, labels, action] = run;

if ~do_plot
    return;
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 1;
    plt_nc = 1;
    plt_np = 1;
    fsiz = [0.1375    0.6907    0.25    0.2139];
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
cols = cols([2 3 1],:);

subplot(plt_nr,plt_nc,plt_np);
h = errorbarKxN(mx',0*mx',labels,struct('colmap',cols,'barwidth',bw));
ylabel('Probability','fontsize',fsy);
alpha(alf);
legend(action,'fontsize',fs,'location','north','box','off');
hax=get(gca,'XAxis');
set(hax,'fontsize',fsy);

end

function [mx, labels, leglabels] = run

P = zeros(6,6);
P(1,2:3) = .5;
P(2,4:5) = .5;
P(3,5:6) = .5;
P(4:6,4:6) = eye(3);

lambda = 10;

c1 = -[0;0;0;0;15;30]/lambda;
[U1,~,MNN] = core_lrl(P,c1);

c2 = c1;
c2(4) = -45/lambda;
[U2] = core_lrl(P,c2,MNN);

mx = [U1(1,2:3); U2(1,2:3)];

labels = {'Training','Test'};
leglabels = {'2','3'};
end
