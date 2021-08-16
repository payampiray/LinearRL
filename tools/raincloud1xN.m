function h = raincloud1xN(y, mdy, edy, labels)
dx =2;
for i=1:size(y,2)
%     face_color = face_color(1,:);
    ax = (i-1)*dx;    
    h(i,:) = draw(ax, y(:,i), mdy(i), edy(i), [1 1 1]*.5, .6);
    hold on;
end

xlim([-1 ax+1]);
set(gca,'xticklabel',labels);

end

function h = draw(ax, y, mdy, edy, face_color, alf)

if nargin<4
face_color = [0.7384    0.0515    0.1973];
end
if nargin<5
alf = .4;
end
isright = 1;


% patch
n_points = 100;
patch_alpha = alf;
box_mid = ax;
n_extra_y = 2;

% bar
bar_mid = .05;
bar_width = .3; % width
bar_linewidth = 1; % line width
rect_alpha = alf;

%------------------------------------
isright(isright==0) = -1;

N = length(y);
maxy=max(y);
miny=min(y);

extra_y = n_extra_y*(maxy-miny)/(n_points);
binranges = linspace(miny-extra_y,maxy+extra_y,n_points);
[dens,dens_pos] = ksdensity(y,binranges,'function','pdf','kernel','normal');
area=sum(dens)*2*(binranges(2)-binranges(1));

%Normalized area is a third of the max available area
max_area= maxy-miny;
dens=dens*max_area/(3*area);
% plot(-dens,1:npoints);

% if left
if ~isright
    xpatch=[box_mid-dens(1:end-1) ; ...
            box_mid+zeros(1,length(dens)-1) ; ...
            box_mid+zeros(1,length(dens)-1);  ...
            box_mid-dens(2:end)];
end

% if right    
if isright
    xpatch=[box_mid+zeros(1,length(dens)-1) ; ...
            box_mid+dens(1:end-1) ; ...
            box_mid+dens(2:end) ; ...
            box_mid+zeros(1,length(dens)-1)];    
end

ypatch=[dens_pos(1:end-1) ; ...
    dens_pos(1:end-1) ; ...
    dens_pos(2:end) ; ...
    dens_pos(2:end)];

minxp = min(min(xpatch,[],2));
maxxp = max(max(xpatch,[],2));
rxp = ceil((maxxp - minxp)*10)/10;

hpatch = patch(xpatch,ypatch,[1 1 1],'EdgeAlpha',1,'FaceColor',face_color,'EdgeColor','none','FaceAlpha',patch_alpha);    
hpatchl = line(box_mid +isright*[0 dens 0 0], [dens_pos(1) dens_pos dens_pos(end) dens_pos(1)],'color',face_color);
% xl1 = rxp*[-1 1]+box_mid;
hold on;

bar_mid = box_mid + -isright*(bar_mid + bar_width);


hrect = rectangle('position',[bar_mid-bar_width/2 0 bar_width mdy],'FaceColor',[face_color rect_alpha],'linewidth',bar_linewidth);
hrectl = plot([bar_mid bar_mid],mdy+[-edy edy],'-','color','k','linewidth',2);

h = [hpatch hpatchl hrect hrectl];
end