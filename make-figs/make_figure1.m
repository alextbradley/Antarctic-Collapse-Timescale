% Make figure 1 of the ms, showing:
% (a) constructed temperature profile along Pine Island
% (b) constructed temperature profile along Ross
% (c) Evolution of crevasse depth in idealized setup in cold and warm
% conditions
% (c)inset Evoluton of temperature profiles in these situations
% (d) 'Collapse' time as a function of initial thickness for different values
% of melt rate.
%
% ATB (aleey@bas.ac.uk), 15/02/2023. MIT Licence

%% Prelimiaries and figure setup
clear
addpath('../functions');
fig = figure(1);clf
positions = [0.05, 0.12, 0.27, 0.38;
             0.05, 0.58, 0.27, 0.38;
             0.38, 0.12, 0.27,  0.84; 
             0.71, 0.12,  0.24, 0.84];
sz        = size(positions);
fig.Position(3:4) = [1200,360];
fs = 13; %fontsize
for i = 1:sz(1)
    ax(i) = subplot('Position', positions(i,:)); hold on
    ax(i).FontSize = fs;
end

%% Make panels (a), (b) showing temperature profiles 
load('../gendata/figure1/figure1-data.mat')
cmap = cmocean('thermal',100);
cmap = cmap(1:end-10,:);

for i = 1:2
    %load data in and surf colours
    T = cell2mat(Tflowline(i));
    xx    = cell2mat(xflowline(i));
    z     = cell2mat(zflowline(i));
    if i == 1 %Ross
        nsy = 30;
        nsx = 2; %how finely to sampel in x and y
    else
        nsy = 1;
        nsx = 1; %how finely to sampel in x and y
    end

    surf(ax(i), xx(1:nsy:end,1:nsx:end)/1e3,z(1:nsy:end,1:nsx:end), T(1:nsy:end,1:nsx:end), 'linestyle', 'none'); view(ax(i), [0,90])
    % have to plut only 1 in 20 values to make size manageable
 
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
ax(1).XLim = [0,1000];
ax(2).XLim = [0,60];
% add the colorbar
c = colorbar(ax(2));
c.Label.String = 'temp (C)';
% c.Label.FontSize = fs+2;
% c.Position(1) = 0.305;
% c.Position(2) = 0.12;
% c.Position(3) = 0.0075;
% c.Position(4) = 0.84;
c.Limits = [-22,2];

c.Location = 'south';
c.Position(4) = 0.02;
c.Position(1) = 0.2;
c.Position(3) = 0.1;

%% Make panel (a), (b) insets showing the flowlines
% for i = 1:2
%     figure(i+1);clf; 
%     xs = insetss(i).xs;
%     ys = insetss(i).ys;
%     us = insetss(i).us;
%    
%     us(~isnan(us)) = 1; %set a uniform colour
%     %pl = imagesc(xs,ys, us');
%     %set(pl, 'AlphaData', ~isnan(us'));
%     contourf(xs,ys, (us'), [1,1]);
%     hold on
%     flc = cell2mat(flc_all(i));
%     plot(flc(:,1), flc(:,2), 'w', 'linewidth', 2.75)
%     colormap(0.3*[1,1,1])
%     axis equal
%     axs = gca;
%     axs.Visible = 'off';
% end

%% Make panel (c), showing crevasse depth as a function of time for both ice shelves
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

%% Make panel (c) inset, showing temperature profiles through the thinning experiment

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

%% Make panel (d): enhancement of collapse time for different initial ice thickness and melt rates

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
    yy = (ct(:,1)./ct(:,im) - 1);
    plot(ax(4),HH, yy, 'linewidth', 1.5, 'color', colmap(im+2,:));
    %drawnow; pause
end

ax(4).YLim = [0, 5];
ax(4).YTick = 0:1:5;
ax(4).XLim = [200,1000];
ax(4).XTick = 0:200:1000;

ax(4).XLabel.String = 'initial ice thickness (m)';
ax(4).YLabel.String = '$(\tau -  \tau_\infty)/\tau_\infty$';
ax(4).YLabel.Interpreter = 'latex';
box(ax(4), 'on')
%set(ax(4), 'YScale', 'log');

% sort out colorbar
axghost = axes();
axghost.Visible = 'off';

nmd = log10(mm);
nm = nmd(end) - nmd(1); 
cc = colorbar; clabels = 10.^(nmd(1):nmd(end)); clim(log10(clabels([1 end]))); set(cc,'Ticks',log10(clabels),'TickLabels',clabels); 
cc.FontSize = c.FontSize;
cc.Label.String = 'melt rate (m/yr)';
cc.Label.FontSize = fs+2;
cc.Colormap = colmapcts;
cc.Position = [0.9050    0.5400    0.0080    0.4000];
cc.Label.Position(1) = 3;

% Tidy panels
for i = 1:sz(1)
    ax(i).FontSize = fs;
    ax(i).XLabel.FontSize = fs+2;
    ax(i).YLabel.FontSize = fs+2;
    ax(i).FontName = 'Arial';
end