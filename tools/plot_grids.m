function [TB, AA] = plot_grids(P,config,I,J,IW,JW,lij)
if nargin<2, config = struct('add_arrow',0); end
if nargin<3, I = sqrt(size(P,1)); end
if nargin<4, J = sqrt(size(P,2)); end
if nargin<5, IW = []; end
if nargin<6, JW = []; end
if nargin<7, lij = []; end

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

inp.addParameter('alpha',.5);
inp.addParameter('font_size',14);
inp.addParameter('colmap','pink');
inp.addParameter('colorful',0);
inp.addParameter('col_terminal',[1 .2 .2]);


% these are used for the path
inp.addParameter('col_path',[1 0 0]);
inp.addParameter('line_width',4);
inp.addParameter('pos0',.3);
inp.addParameter('arrow_shift0',.2);
inp.addParameter('arrow_length0',.8);

inp.parse(config);
config = inp.Results;
%--------------------------------------------------------------------------
% plot properties
fs = config.font_size;
col = config.col_path;
alf = config.alpha;
linewidth = config.line_width;
colmap_name = config.colmap;
colorful = config.colorful;
col_terminal = config.col_terminal;

% these parameters determine the shape of arrow
pos0 = config.pos0;
arrow_shift0  = config.arrow_shift0;
arrow_length0 = config.arrow_length0;

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

%--------------------------------------------------------------------------

% adjust U and P for blocked ones by adding zero row and columns
P0 = P; clear P;
sb = isnan(lij_blked(:,1));
P = zeros(length(sb),length(sb));
P(~sb,~sb) = P0;

% create a binary matrix determining terminals and blocked ones (TB)
terms = (diag(P)==1);
terminals = false(size(P,1),1);
terminals(terms) = 1;
TB = ones(I,J);
TB(terminals) = 0;
TB(isnan(A)) = 0;

% plot A and add TB as an image with AA as the alpha value
% h = heatmap(A);
% set(h,'colormap',[1 1 1],'ColorbarVisible','off',...
%     'YDisplayLabels',repmat({''},I,1),'XDisplayLabels',repmat({''},J,1) );


% create a transparency matrix: alf for all states but the blocked ones.
if colorful
    AA = zeros(size(TB));    
    AA(isnan(A)) = 1;
    TBC = nan([I,J,3]);
    colorful_sqs = find(terminals);
    col_colorfuls = repmat(col_terminal,length(colorful_sqs),1);

    for i=1:size(colorful_sqs)
        [it,jt] = ind2sub([I,J],colorful_sqs(i));
        TBC(it,jt,:) = col_colorfuls;
        AA(it,jt) = alf;
    end
    imagesc(TBC,'alphadata',AA);
else

    AA = alf*ones(size(TB));    
    AA(isnan(A)) = 1;
    imagesc(TB,'alphadata',AA);
    colormap(gca,colmap_name);
end

% add grid
[X,Y]=meshgrid(1:(J+1),1:(I+1));
hold on;
plot(X-.5,Y-.5,'k');
plot(X'-.5,Y'-.5,'k');
% axis off;
set(gca,'fontsize',fs,'Xtick',[],'Ytick',[]);

% add state labels (note: do not label the blocked ones)
if ~isempty(config.labels)
    % if labels are given
    labels = config.labels;
    config.add_labels = 1;
else
    % labels are state numbers
    labels = lij(:,1);
    config.str_label = '%d';
end

if ~iscell(labels)
    labels0 = labels;
    labels = num2cell(labels);
    labels(isnan(labels0)) = cell(sum(isnan(labels0)),1);
end
L = cell(size(lij_blked(:,1)));
L(lij(:,1)) = labels;

x = lij_blked(:,2);
y = lij_blked(:,3);
if config.add_labels   
    for k=1:numel(L)
        if ~isempty(L{k})
            text(y(k),x(k),sprintf(config.str_label,L{k}),...
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

gcfunit = get(gca,'units');
set(gcf,'units','normalized');
pos = get(gca,'position');
set(gca,'units',gcfunit);


x = lij_blked(:,3);
y = lij_blked(:,2);
path = config.path;

arrow_shift = arrow_shift0*pos(3)/pos0;
arrow_length = arrow_length0*pos(3)/pos0;

if ~isempty(path)
    path0 = path;    
    path = lij(path0,1);


    for i=1:(length(path)-1)
        s = path(i);
        nexts = path(i+1);
        if ~isempty(nexts) && ~terminals(s)
    %         u = U(s,nexts);
            u = 1;
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

end