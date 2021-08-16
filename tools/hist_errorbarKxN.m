function [h, hb] = hist_errorbarKxN(x,y,mx,ex,facnames,config)


%-------------------------------------------
if nargin<4, facnames = []; end
if nargin<5
    config = struct('legnames',[]);
end

p = inputParser;
p.addParameter('legnames',[]);
p.addParameter('colmap',[]);
p.addParameter('basevalue',0);
p.addParameter('barwidth',[]);
p.addParameter('hist_dist',.5);
p.addParameter('hist_scale',1);
p.parse(config);
config    = p.Results;

legnames = config.legnames;
colmap = config.colmap;
basevalue = config.basevalue;
barwidth = config.barwidth;
hist_dist = config.hist_dist;
hist_scale = config.hist_scale;
%-------------------------------------------

[K,N] = size(mx);

if ~iscell(y)
    l = 0;
    z = cell(K,N);    
    for k=1:K
        for n=1:N
            l = l+1;
            z{k,n} = y(:,l);
        end
    end
    y = z;
    
%     l = 0;
%     u = cell(K,N);    
%     for k=1:K
%         for n=1:N
%             l = l+1;
%             u{k,n} = x(:,l);
%         end
%     end
%     x = u;    
end


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
        
        sgn = 2*i-3;        
        ah = ax(k)+sgn*hist_dist + [-y{k,i}/2 y{k,i}/2]*hist_scale;
        xx = [x x];
        for l=1:length(ah)
            plot(ah(l,:),xx(l,:),'color',[1 1 1]*0,'linewidth',1);
        end
        
    end
    if isleg, hleg = legend(legnames); set(hleg,'AutoUpdate','off'); isleg = 0; end
    for k=1:K        
        plot([ax(k);ax(k)],[el(k,i);eh(k,i)],'-','color','k','linewidth',2);
    end
end

yl = get(gca,'ylim');
yl(1) = basevalue;
set(gca,'ylim',yl);

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

% set axes propertis
set(gca,'box','off');
set(gca,'ticklength', [0 0]);

% return axes handle
h = gca;

end