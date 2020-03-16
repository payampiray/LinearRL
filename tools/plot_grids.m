function [TB, AA] = plot_grids(U,config,I,J,IW,JW,P,lij)
if nargin<2, config = struct('add_arrow',0); end
if nargin<3, I = sqrt(size(U,1)); end
if nargin<4, J = sqrt(size(U,2)); end
if nargin<5, IW = []; end
if nargin<6, JW = []; end
if nargin<7, P = U; end
if nargin<8, lij = []; end
   
%--------------------------------------------------------------------------
if isempty(lij)
    l = (1:(I*J))';
    [i,j] = ind2sub([I,J],l);
    lij = [l i j];
end

if isempty(config)
    config.add_labels = 1;
end
% if nargin>6, config.add_arrow=1; end
% if ~isempty(config.path), config.add_arrow = 1; end

inp = inputParser;
inp.addParameter('add_arrow',nargin>6);
inp.addParameter('add_labels',1);
inp.addParameter('labels',[]);
inp.addParameter('str_label','%0.1f');
inp.addParameter('path',[]);
inp.parse(config);
config = inp.Results;
%--------------------------------------------------------------------------
% plot properties
fs = 12;
col = [1 0 0];
alf = .5;
linewidth = 4;

% these parameters determine the shape of arrow
pos0 = .68;
arrow_shift0  = .1;
arrow_length0 = .8;

try %#ok<TRYNC>
%     fs = get(gca,'fontsize');
end

%--------------------------------------------------------------------------
% A is the same as lij(:,1), but in a matrix format
ind = sub2ind([I,J],lij(:,2),lij(:,3));
A= nan(I,J);
A(ind) = lij(:,1);
% heatmap(A);


% manages barriers
if size(lij)<(I*J)
    idx_blked = false(I*J,1);
    lij_blked = nan(I*J,3);    
    l0 = (1:(I*J))';    
    for k=1:length(l0)
        index = lij(:,1)==l0(k);
        if sum(index)==0
            idx_blked(k) = 1;
        elseif sum(index)==1
            lij_blked(k,:) = lij(index,:);
        elseif sum(index)>1
            error('!');
        end
    end
else
    lij_blked = lij;
end
% 
% l = (1:(I*J))';
% [i,j] = ind2sub([I,J],l);
% lij_temp = [l i j];
% A= nan(I,J);
% A(ind) = lij_temp(:,1);
% figure;
% heatmap(A);

%--------------------------------------------------------------------------

% adjust U and P for blocked ones by adding zero row and columns
U0 = U; clear U;
P0 = P; clear P;
sb = isnan(lij_blked(:,1));
U = zeros(length(sb),length(sb));
P = zeros(length(sb),length(sb));
% k = 0;
% for i=1:size(P,1)
%     if ~isnan(lij_blked(i))
%         k = k+1;
%         P(i,~sb) = P0(k,:);
%         U(i,~sb) = U0(k,:);
%     end
% end
U(~sb,~sb) = U0;
P(~sb,~sb) = P0;

% create a binary matrix determining terminals and blocked ones (TB)
terms = (diag(P)==1);
terminals = false(size(P,1),1);
terminals(terms) = 1;
TB = ones(I,J);
TB(terminals) = 0;

TB(isnan(A)) = 0;

% create a transparency matrix: alf for all states but the blocked ones.
AA = alf*ones(size(TB));
AA(isnan(A)) = 1;

% plot A and add TB as an image with AA as the alpha value
h = heatmap(A);
set(h,'colormap',[1 1 1],'ColorbarVisible','off',...
    'YDisplayLabels',repmat({''},I,1),'XDisplayLabels',repmat({''},J,1) );
% fs = get(gca,'fontsize');
imagesc(TB,'alphadata',AA);
colormap(gca,'winter');


% add grid
[X,Y]=meshgrid(1:(J+1),1:(I+1));
hold on;
plot(X-.5,Y-.5,'k');
plot(X'-.5,Y'-.5,'k'); axis off
set(gca,'fontsize',fs);

% add state labels (note: do not label the blocked ones)
if ~isempty(config.labels)
    % if labels are given
    labels = config.labels;
    config.add_labels = 1; 

    L = nan + lij_blked(:,1);
    L(lij(:,1)) = labels;
    
else
    % labels are state numbers
    L = lij_blked(:,1);
    config.str_label = '%d';
end

x = lij_blked(:,2);
y = lij_blked(:,3);
if config.add_labels   
    for k=1:numel(L)
        if ~isnan(L(k))
            text(y(k),x(k),sprintf(config.str_label,L(k)),...
                'HorizontalAlignment','center','fontsize',fs);            
        end
    end
end

% add walls. note: make sure that doors are open!
n = I*J;
XJW = X(:,JW+1);
YJW = Y(:,JW+1);
AJW1 = A(:,JW);
AJW2 = A(:,JW+1);
for w=1:size(XJW,2)
    nans = isnan(AJW1(:,w)) | isnan(AJW2(:,w));
    doors = false(size(nans));
    doors(~nans) = P(sub2ind([n,n],AJW1(~nans,w),AJW2(~nans,w)));
    for i=1:I
        if ~doors(i)
            plot(XJW([i i+1],w)-.5,YJW([i i+1],w)-.5,'k','linewidth',3)
        end
    end
end

XIW = X(IW+1,:);
YIW = Y(IW+1,:);
AIW1 = A(IW,:);
AIW2 = A(IW+1,:);
for w=1:size(XIW,1)
    nans = isnan(AIW1(w,:)) | isnan(AIW2(w,:));
    doors = false(size(nans));
    doors(~nans) = P(sub2ind([n,n],AIW1(w,~nans),AIW2(w,~nans)));
    for j=1:J
        if ~doors(j)
            plot(XIW(w,[j j+1])-.5,YIW(w,[j j+1])-.5,'k','linewidth',3)
        end
    end
end

% add U-dependent arrows
if ~config.add_arrow
    return;
end

gcfunit = get(gcf,'units');
set(gcf,'units','normalized');
pos = get(gcf,'position');
set(gcf,'units',gcfunit);


x = lij_blked(:,3);
y = lij_blked(:,2);
path = config.path;

arrow_shift = arrow_shift0*pos(3)/pos0;
arrow_length = arrow_length0*pos(3)/pos0;

states = 1:size(U,1);
if ~isempty(path)
    states = path;
end

for i=1:length(states)
    s = states(i);
    nexts = find(U(s,:)>0);
    if ~isempty(nexts) && ~terminals(s)
        u = U(s,nexts);

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
            end
            plot([x1 x2],[y1 y2],'linewidth',linewidth*intens(k),'color',[col u(k)]);
            scatter(x2,y2,'Marker',marker,'MarkerFaceColor','r','MarkerEdgeColor','r',...
                'MarkerFaceAlpha',u(k),'MarkerEdgeAlpha',u(k)); 
        end
    end
end

end