function h = sim_DR_maze(plt_nr,plt_nc,plt_np)
do_plot = 1;

[TB1,U1,path1] = run_DR_maze(1);

[TB2,U2,path2,AA,imp_states,labels,xi,yi,x,y] = run_DR_maze(2);

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

h(1) = subplot(plt_nr,plt_nc,plt_np(1));
add_arrow(TB1,AA,imp_states,labels,xi,yi,U1,x,y,path1);
set(gca,'xtick',[],'ytick',[],'box','on');


h(2) = subplot(plt_nr,plt_nc,plt_np(2));
add_arrow(TB2,AA,imp_states,labels,xi,yi,U2,x,y,path2);
set(gca,'xtick',[],'ytick',[],'box','on');
end

function [TB2,Um,path,AA,imp_states,labels,xi,yi,x,y] = run_DR_maze(goal)
n = 21;
I = n;
J = n;
[P, ~, lij] = core_griding(I,J);
terms =  [n n*n];
start = (n-1)*n/2+1;
c= 0.001*ones(size(P,1),1);
for i=1:2 %:length(terms)
    P(terms(i),:) = 0;
    P(terms(i),terms(i)) = 1;
end
c(terms) = 0;
c(terms(goal)) = -5;

U = core_lmdp(P,c);


[~, path] = core_follow_path(P,U,c,start);

Um = zeros(size(U));
for i=1:size(U,1)
   [~,ix] = max(U(i,:)); 
   Um(i,ix) = 1;
end

TB = ones(I,J);
AA = TB*.5;
imp_states = [start terms];

TB2 = nan([size(TB),3]);
TB2(:,:,1) = TB;
[xi,yi] = ind2sub(size(TB),imp_states);

cols = def('col');
cols(1,:) = [1 .7 0];
alf = .6;

for k=1:length(imp_states)
    TB2(xi(k),yi(k),:) =  cols(2,:);
end

AA(AA==1) = 0;
AA(imp_states) = alf;
AA(AA==.5) = 0;

labels = {'H','F','W'};

x = lij(:,3);
y = lij(:,2);

end

function add_arrow(TB2,AA,imp_states,labels,xi,yi,U,x,y,path)
imagesc(TB2,'AlphaData',AA);
set(gca,'box','off','ytick',[],'xtick',[]);

for k=1:length(imp_states)
    text(yi(k),xi(k),labels{k},'HorizontalAlignment','center',...
        'fontsize',12,'fontweight','bold');
end

hold on;

col = [1 0 0];
linewidth = 4;

pos0 = .68;
arrow_shift0  = .1;
arrow_length0 = .8;

gcfunit = get(gcf,'units');
set(gcf,'units','normalized');
axis equal
pos = get(gcf,'position');
set(gcf,'units',gcfunit);

arrow_shift = arrow_shift0*pos(3)/pos0;
arrow_length = arrow_length0*pos(3)/pos0;

states = 1:size(U,1);
if ~isempty(path)
    states = path;
end

for i=2:(length(states)-1)
    s = states(i);
    nexts = find(U(s,:)>0);
    if ~isempty(nexts) %&& ~terminals(s)
        u = U(s,nexts);
        
        u = .5;
        intens = u;

        for k=1:length(nexts)
            x1 = x(s);
            y1 = y(s);

            x2 = x(nexts(k));
            y2 = y(nexts(k));

            if x1 == x2 && y2>y1
                marker = 'v';
                y1 = y1 + arrow_shift;
                y2 = y2 - arrow_length;

            elseif x1 == x2 && y2<y1
                marker = '^';
                y1 = y1 - arrow_shift;
                y2 = y2 + arrow_length;

            elseif y1 == y2 && x2>x1
                marker = '>';
                x1 = x1 + arrow_shift;
                x2 = x2 - arrow_length;

            elseif y1 == y2 && x2<x1
                marker = '<';
                x1 = x1 - arrow_shift;
                x2 = x2 + arrow_length;
            else
                fprintf('something is wrong');
            end
            plot([x1 x2],[y1 y2],'linewidth',linewidth*intens(k),'color',[col u(k)]);
            a = scatter(x2,y2,'Marker',marker,'sizedata',20,'MarkerFaceColor','r','MarkerEdgeColor','r',...
                'MarkerFaceAlpha',u(k),'MarkerEdgeAlpha',u(k)); 
        end
    end
end

end