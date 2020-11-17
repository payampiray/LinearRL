function h = sim_2step(plt_nr,plt_nc,plt_np)
do_plot = 1;

[mx, labels, action] = run;

if ~do_plot
    return;
end

%--------------------------------------------------------------------------
if nargin<1
%     close all;    
    plt_nr = 1;
    plt_nc = 1;
    plt_np = 1;
    fsiz = [0.1375    0.6907    0.2    0.25];
    fsiz = [0.2    0.4    0.2    0.25];    
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
cols = cols([2 3],:);

subplot(plt_nr,plt_nc,plt_np);
h = errorbarKxN(mx',0*mx',labels,struct('colmap',cols,'basevalue',.5,'barwidth',bw));
ylabel('stay probability','fontsize',fsy);
alpha(alf);
legend(action,'fontsize',fs,'location','north','box','off');
hax=get(gca,'XAxis');
set(hax,'fontsize',fsy);
ylim([.5 .89]);
end

function [x, rownames, colnames]=run
P = zeros(7,7);
P(1,2:3) = [.5 .5];
P(2,4:5) = .5;
P(3,6:7) = .5;
P(4:7,4:7) = eye(4);


r0 = 0;
rb = r0;
r(:,1) = [zeros(3,1);r0; +.25;r0;r0];
r(:,2) = [zeros(3,1);r0; -.25; r0;r0];

T = [.7 .3;.3 .7];
Tname = {'common', 'rare';'rare', 'common'};
pstay = .75;

Ti = T^-1;


p = [];
lbls = cell(2,0);
for i=1:size(r,2)
    [U] = core_lrl(P,-r(:,i));
    u = U(1,2:3)';    
    p0 = Ti*u;    
    
    for j=1:2
        ps = zeros(2,1);
        ps(j) = pstay;
        ps(3-j) = 1-pstay;
        pii = p0.*ps;
        pi(j) = pii(j)/sum(pii);
    end
    
    s = find(r(:,i)~=r0);
    if s<6
        Tnamei = Tname(1,:);
    else
        Tnamei = Tname(2,:);
    end
    if r(5,i)<rb
        cn = repmat({'rewarded'},1,2);
    else
        cn = repmat({'unrewarded'},1,2);
    end
    lbl = [cn;Tnamei];
    p = [p pi];    
    lbls = [lbls lbl ]; 
end

x = [p(1:2);p(3:4)];
rownames = {'rewarded','unrewarded'};
colnames = {'common','rare'};

end
