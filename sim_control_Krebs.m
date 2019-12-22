function h = sim_control_Krebs(plt_nr,plt_nc,plt_np)
% reward decreases error rate, Krebs et al. 2010 cognition
do_plot = 1;

% [esim, mdata, sdata, labels] = run;
[mx, ex, titlelabel, ylbl, labels] = run;

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
cols = def('col');


col(1,:) = [56 189 165]/255;
col(2,:) = [245 104 80]/255;
%----------------------


for i=1:2    
    h(i) = subplot(plt_nr,plt_nc,plt_np(i) );    
    errorbarKxN(mx{i},ex{i},{''},struct('barwidth',bw,'colmap',col));
    set(gca,'fontsize',fs);
    ylabel(ylbl{i},'fontsize',fsy);
    xlabel('Trial type','fontsize',fsy);
    title(titlelabel{i},'fontsize',fsy);
    hax=get(gca,'XAxis');
    set(hax,'fontsize',fsy);

    if i==1
        legend(labels,'fontsize',fsy,'location','northwest','box','off');
    end
    if i==2
        ylim([0 1]);
    end
end

end

function [mx, ex, titlelabel, ylbl, labels] = run


vc = [2 0];
p = [.8 .2];

lambda = 1;
for i=1:length(vc)
    value = [0 vc(i)];
    ui = exp(value/lambda).*p;
    ui = ui/sum(ui);
    u(i,:) = ui;
    v(i,:) = value;
end
esim = u(:,1);


% Krebs et al data, Table 1
N = 20;
m1 = [5.7 8.8 10.5 6.7];
e1 = [(2.6) (5.1) (4.3) (3.4)]/sqrt(N);

m2 = [10.2 19.7 15.6 12.5];
e2 = [(5.1) (10.4) (8.0) (7.5)]/sqrt(N);
labels = {'Potential reward','No reward'};

sdata = mean([e1; e2],2);
mdata = mean([m1; m2],2);

mx = {mdata,esim};
ex = {sdata,0*esim};
titlelabel = {'Data','Model'};
ylbl = {'Error (%)','Error probability'};

end