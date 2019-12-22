function [h, hb] = errorbarKxN(mx,ex,facnames,config)
% examples
% mx = [ [2;2.5;2.8]  [.05;.25;.9] ];
% ex = .5*rand(size(mx));
% facnames = {'x1','x2'};
% legnames = {'leg1','leg2','leg3'};
% figure;errorbarKxN(mx,ex,facnames);
% another example
% mx = [ [2;2.5;]  [.05;.25;] ];
% ex = .5*rand(size(mx));
% facnames = {'x1','x2'};
% legnames = {'leg1','leg2'};

%-------------------------------------------
if nargin<3, facnames = []; end
if nargin<4
    config = struct('legnames',[]);
end

p = inputParser;
p.addParameter('legnames',[]);
p.addParameter('colmap',[]);
p.addParameter('basevalue',0);
p.addParameter('barwidth',[]);
p.parse(config);
config    = p.Results;

legnames = config.legnames;
colmap = config.colmap;
basevalue = config.basevalue;
barwidth = config.barwidth;
%-------------------------------------------

[K,N] = size(mx);
% if(K>3), error('Number of rows must be less than 3'); end

if size(ex,1)==K
    el    = mx-ex;
    eh    = mx+ex;
elseif size(ex,1)==(2*K)
    el    = ex(1:K,:);
    eh    = ex(K+(1:K),:);
else
    error('!');
end

% if size(ex,2)~=size(mx,2), ex=ex'; end
% if size(ex,2)~=size(mx,2), error('ex is not matched with mx'); end

kk = 0:(1/K):1; kk(end)=[];
dx = 2;

if isempty(colmap)
if K==1
    colmap = [.5 .5 .5];
else
    colmap = repmat( (0:(1/(K-1)):1)',1,3);
% colmap = repmat( ((1/(K)):(1/(K)):1)',1,3);
% colmap = repmat( ((1/(K+1)):(1/(K+1)):(K/(K+1)))',1,3);
end
end

if isempty(barwidth)
    wb = 1/(K+.5);
else
    wb = barwidth;
end
isleg = ~isempty(legnames);
a     = nan(1,N);
% figure;
% axs = zeros(N,2);
hb = zeros(N,K);
for i=1:N
    ax = -median(kk) + kk + +dx*(i-1);
    axs(i,:) = ax;
    a(i) = median(ax);
    for k=1:K
        hb(i,k) = bar(ax(k),mx(k,i),wb,'FaceColor',colmap(k,:),'EdgeColor','k','linewidth',1,'basevalue',basevalue);
        hold on;
    end
    if isleg, hleg = legend(legnames); set(hleg,'AutoUpdate','off'); isleg = 0; end
    for k=1:K        
        plot([ax(k);ax(k)],[el(k,i);eh(k,i)],'-','color','k','linewidth',2);
    end
end

if isempty(facnames)
    facnames = repmat({''},1,N);
end
set(gca,'xtick',a);
set(gca,'xticklabel',facnames);

if K>1
dd = axs(1,2)-axs(1,1);
xlims = [axs(1,1)-dd axs(end,end)+dd];
set(gca,'xlim',xlims);
end

% % set y-axis
% yrng = get(gca,'ytick');
% ystp = diff(yrng); 
% ystp = ystp(1);
% ymax = yrng(end)+.95*ystp;
% ymin = yrng(2)-0.95*ystp;
% set(gca,'ylim',[ymin;ymax]);

% set axes propertis
set(gca,'box','off');
set(gca,'ticklength', [0 0]);

% return axes handle
h = gca;

end