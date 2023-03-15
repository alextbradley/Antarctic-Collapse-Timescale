% Make figure 1 of the ms, showing:
% (a) constructed temperature profile along Pine Island
% (b) constructed temperature profile along Ross
% (c) Evolution of crevasse depth in idealized setup in cold and warm
% conditions
% (d) Evoluton of temperature profiles in these situations
% (e) 'Collapse' time as a function of initial thickness for different values
% of melt rate.
%
% ATB (aleey@bas.ac.uk), 15/02/2023. MIT Licence

%% Prelimiaries and figure setup
clear
addpath('../functions');
fig = figure(1);clf
positions = [0.05, 0.12, 0.24, 0.38;
             0.05, 0.58, 0.24, 0.38;
             0.41, 0.12, 0.24,  0.84; 
             0.71, 0.12,  0.24, 0.84];
sz        = size(positions);
fig.Position(3:4) = [1200,360];
fs = 13; %fontsize
for i = 1:sz(1)
    ax(i) = subplot('Position', positions(i,:)); hold on
    ax(i).FontSize = fs;
end

%% make (a), (b)
load('../gendata/figure1/figure1-data.mat')
cmap = cmocean('thermal',100);
cmap = cmap(1:end-10,:);

for i = 1:2
    %load data in and surf colours
    T = cell2mat(Tflowline(i));
    xx    = cell2mat(xflowline(i));
    z     = cell2mat(zflowline(i));
    surf(ax(i), xx/1e3,z, T, 'linestyle', 'none'); view(ax(i), [0,90])
 
    %tidy plot
    ax(i).YTick = -1000:500:0;
    ax(i).YLim = [-1000, 200];
    ax(i).YLabel.String = 'depth (m)';
    clim([-20,-2])
    colormap(ax(i),cmap);
    grid(ax(i), 'on');
    box(ax(i), 'on');
    

end
ax(1).XLabel.String = 'distance along flowline (km)';

% add the colorbar
c = colorbar(ax(2));
c.Label.String = 'temp (C)';
c.Label.FontSize = fs+2;
c.Position(1) = 0.305;
c.Position(2) = 0.12;
c.Position(3) = 0.0075;
c.Position(4) = 0.84;
c.Limits = [-22,2];

%% make the insets showing the flowline
if 0
for i = 1:2
    figure(i+1);clf; 
    xs = insetss(i).xs;
    ys = insetss(i).ys;
    us = insetss(i).us;
   
    us(~isnan(us)) = 1; %set a uniform colour
    %pl = imagesc(xs,ys, us');
    %set(pl, 'AlphaData', ~isnan(us'));
    contourf(xs,ys, (us'), [1,1]);
    hold on
    flc = cell2mat(flc_all(i));
    plot(flc(:,1), flc(:,2), 'w', 'linewidth', 2.75)
    colormap(0.3*[1,1,1])
    axis equal
    axs = gca;
    axs.Visible = 'off';
end
end

%% make plot showing crevasse depth as a function of time
colmap = [0, 33, 153; 153,0,63]/255;
ss = figc_data;
 for i =1:2
     plot(ax(3), [ss(i,:).t], [ss(i,:).dimless_crev_depth], 'linewidth', 1.8, 'color',colmap(i,:));
     box(ax(3), 'on');

 end
 ax(3).XLabel.String = 'time (yrs)';
 ax(3).XLabel.Interpreter = 'none';
 ax(3).XTick = 0:100:400;
 ax(3).XLim = [0,400];
 ax(3).YLim = [0,1];
 ax(3).YTick = 0:0.2:1;
 ax(3).YLabel.String = '$d_{c}/H$';
 ax(3).YLabel.Interpreter = 'latex';



% %% add inset showing evolution of temperature profile
% figure(2); clf;
% for i = 1:2
%     subplot(1,2,i); hold on
%     tt = [ss(i,:).t];
%     for it = 1:20:length(tt)
%         Tp = [ss(i,it).temp_profile];
%         plot(Tp-273.15, linspace(0,1,length(Tp)), 'linewidth', 1.5);
%         %pause
%         
%     end
% end

%%

axin = axes;
axin.Position = [0.545, 0.175, 0.1, 0.34];

axin.XTick = [-25, 0];
axin.YTick = [0,1];
Ts = -22;
Tb = -2;
axin.XLim = [-25,0];
box(axin, 'on');
hold(axin, 'on')
nl = 4; %number of lines
for i = 1:2
    tt = [ss(i,:).t];
    
    idxp = linspace(1,length(tt), nl); %time indices to plot
    idxp = round(idxp);

    %make a colormap
    if i == 1
        T = [180/255,180/255, 218/255 ; %first colour
        colmap(i,:)];
    else
        T = [218/255,180/255, 180/255 ; %first colour
        colmap(i,:)];
    end
    x = [0 1];
    map = interp1(x,T,linspace(0,1,nl));

    for it = 1:length(idxp)
        Tp = [ss(i,idxp(it)).temp_profile];
        zz = [ss(i,idxp(it)).zz];

        p = plot(axin, Tp-273.15, zz, 'linewidth', 1.5);
        palpha =  it/length(idxp);
        p.Color = [map(it,:)]; %set the color nad alpha

        %put the point on the main line
        scatter( ax(3),[ss(i,idxp(it)).t], [ss(i,idxp(it)).dimless_crev_depth],'filled', 'MarkerFaceColor' , map(it,:), 'MarkerEdgeColor', 'k', ...
            'LineWidth',1);
    end
   
end
axin.FontSize = 13;

%% part (d)

%get data from array
mm = figd_data.m;
HH = figd_data.h;
ct = figd_data.ct;

%setup colormap
colmap = flipud(cmocean('matter',length(mm)+2));
colmapcts = cmocean('matter');
colmapcts = colmapcts(40:end,:);

%plot curves sequentialy
for im = 1:length(mm)-1
    yy = (ct(:,1)./ct(:,im) - 1)*100;
    plot(ax(4),HH, yy, 'linewidth', 1.5, 'color', colmap(im+2,:));
    %drawnow; pause
end

ax(4).YLim = [0, 500];
ax(4).YTick = 0:100:500;
ax(4).XLim = [200,1000];
ax(4).XTick = 0:200:1000;

ax(4).XLabel.String = 'current ice thickness (m)';
ax(4).YLabel.String = '$\tau(h,m)/\tau_\infty(h)$';
ax(4).YLabel.Interpreter = 'latex';
box(ax(4), 'on')
%set(ax(4), 'YScale', 'log');

%% sort out colorbar
axghost = axes();
axghost.Visible = 'off';

nmd = log10(mm);
nm = nmd(end) - nmd(1); 
cc = colorbar; clabels = 10.^(nmd(1):nmd(end)); clim(log10(clabels([1 end]))); set(cc,'Ticks',log10(clabels),'TickLabels',clabels); 
cc.FontSize = c.FontSize;
cc.Label.String = 'melt rate (m/yr)';
cc.Label.FontSize = fs+2;
cc.Colormap = colmapcts;
cc.Position(3) = c.Position(3);
cc.Position(2) = c.Position(2);
cc.Position(4) = c.Position(4);
cc.Position(1) = 0.955;
cc.Label.Position(1) = 3;
%% tidy
for i = 1:sz(1)
    ax(i).FontSize = fs;
    ax(i).XLabel.FontSize = fs+2;
    ax(i).YLabel.FontSize = fs+2;
    ax(i).FontName = 'Arial';
end