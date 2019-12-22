function h = sim_bandit(plt_nr,plt_nc,plt_np)
do_plot = 1;

[mx, u, leglabels, mo, uo] = run;

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
bw  = .15;
bw  = .15;
cols = def('col'); cols(1,:) = []; cols = cols([2 1],:);
cols(3,:) = [0 .2 .8];

h(1) = subplot(plt_nr,plt_nc,plt_np(1));
plot(u,mx(:,1),'linewidth',2,'color',[cols(1,:) alf]); hold on;
plot(u,mx(:,2),'linewidth',2,'color',[cols(2,:) alf]);
plot(u,mx(:,3),'linewidth',2,'color','k');
plot([0,1],[mo mo],'--','linewidth',1,'color',[0 0 0 alf]);
plot([uo,uo],[mo mo],'o','linewidth',1,'color',[0 0 0 alf],'markersize',10,'MarkerFaceColor','k');
% ht = text(uo,mo*.7,'$\pi^*_A=0.73$','fontsize',fsy+2,'fontweight','bold','Interpreter','latex','HorizontalAlignment','center');
text(uo,mo*.7,sprintf('$v^*_A=%0.2f$',mo),'fontsize',fsy+2,'fontweight','bold','Interpreter','latex','HorizontalAlignment','center');
% ht = text(uo,mo*.8,'\pi^*=0.73','fontsize',fsy,'fontweight','bold','Interpreter','tex','HorizontalAlignment','center');
% ylim([0 1]);
ylabel('Value','fontsize',fsy);
xlabel('Probability of choosing A','fontsize',fsy);
legend(leglabels,'location','southeast','fontsize',fsy,'box','off');
end

function [mx, u, leglabels, mo, uo] = run

p = [.5 .5];
c = -[1 0];

u  = (0.001:.001:.999)';
uu = [u 1-u];
Er = -uu*c';
kl = sum(uu.*log(bsxfun(@times,uu,1./p)),2);

mx = [Er kl Er-kl];
% mx(mx==0) = .005;

leglabels = {'Expected reward','Control cost',['Expected reward ' char(8211) ' Control cost']};

uo = exp(-c).*p;
uo = uo./sum(uo);
Ero = -uo*c';
klo = sum(uo.*log(bsxfun(@times,uo,1./p)),2);
mo = Ero-klo;
uo = uo(1);
end
