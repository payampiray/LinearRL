function h = sim_PIT_Corbit(plt_nr,plt_nc,plt_np)

do_plot = 1;

[mx, ex, ylbl, ttl, labels]=run;

if ~do_plot
    return;
end

%--------------------------------------------------------------------------
if nargin<1
    close all;    
    plt_nr = 2;
    plt_nc = 2;
    plt_np = 1:4;
    fsiz = [0.3536    0.6907    0.3    0.5];    
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

%----------------------

h = nan(1,4);
for i=1:4
    h(i) = subplot(plt_nr,plt_nc,plt_np(i) );    
    errorbarKxN(mx{i}',ex{i}',cell(1,2),struct('barwidth',bw));
    set(h(i),'fontsize',fs);
    ylabel(ylbl{i},'fontsize',fsy);
    title(ttl{i},'fontsize',fsy);    

    if i==2 || i==4
        legend(labels,'fontsize',fsy,'box','off');
        ry = max(get(h(i-1),'ylim'))/max(mx{i-1});
        ylim([0 max(ry)*max(mx{i})]);       
    end
    xlabel('Stimulus','fontsize',fsy);
end
text(xsA-.1,.5,'Hunger','fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(1),'horizontalalignment','center','Rotation',90);
text(xsA-.1,.5,'Satiety','fontsize',fsA,'Unit','normalized','fontname',fn,'parent',h(3),'horizontalalignment','center','Rotation',90);

end

function [mx, ex, ylbl, ttl, labels]=run
% 
% % without Pavlovian training
% % state 1: S1,
% % state 2: O1
% % state 3: O2
% % state 4: O3

% after Pavlovian training
pp = .5;
ns = 4;
resprate = 1;

P = zeros(ns,ns);
P(1,:) = [0 pp (1-pp)/(ns-2)*ones(1,ns-2)];
P(2:ns,2:ns) = eye(ns-1);
c = zeros(ns,1);
c(2:(ns-1)) = -5;
U1 = core_lmdp(P,c);

same1 = U1(1,2)*resprate;
diff1 = U1(1,3)*resprate;
msim1 = [same1 diff1];

% with Pavlovian training under satiation
% resprate = 20;

P = zeros(ns,ns);
P(1,:) = [0 pp (1-pp)/(ns-2)*ones(1,ns-2)];
P(2:ns,2:ns) = eye(ns-1);
c = zeros(ns,1);
U2 = core_lmdp(P,c);

same2 = U2(1,2)*resprate;
diff2 = U2(1,3)*resprate;
msim2 = [same2 diff2];

% extracted using this tool from Corbit et al. paper
% https://apps.automeris.io/wpd/
Dist0 =  385.62339181286563; % ==10
Dist1 =  873.8947368421052;
Dist2 =  964.0594774701386;
Dist3 =  568.726837823757;
Dist4 =  665.8260063337149;

dm_same1 = Dist1/Dist0*10;
de_same1 = (Dist2-Dist1)/Dist0*10;
dm_diff1 = Dist3/Dist0*10;
de_diff1 = (Dist4-Dist3)/Dist0*10;
mdata1 = [dm_same1 dm_diff1];
edata1 = [de_same1 de_diff1];
labels = {'same','different'};

Dist1 = 400.6770491661981;
Dist2 = 448.3067667558829;
Dist3 = 208.55383330219934;
Dist4 = 257.83144400558587;

dm_same2 = Dist1/Dist0*10;
de_same2 = (Dist2-Dist1)/Dist0*10;
dm_diff2 = Dist3/Dist0*10;
de_diff2 = (Dist4-Dist3)/Dist0*10;
mdata2 = [dm_same2 dm_diff2];
edata2 = [de_same2 de_diff2];
mx = {mdata1, msim1, mdata2, msim2};
ex = {edata1, msim1*0, edata2, 0*msim2};
ylbl = {'Lever presses','Probability','Lever presses','Probability'};
ttl = {'Data','Model','Data','Model'};
end
