% Make figure 1 of the ms, showing:
% (a) temperature profiles are the start of the simulation and throughout
% in the warm and cold conditions
% (b) evolution of crevasse depth in idealized setup in cold and warm
% conditions
% (c) 'Collapse' time as a function of initial thickness for different values
% of melt rate.
%
% ATB (aleey@bas.ac.uk), 15/02/2023. MIT Licence

%% Prelimiaries and figure setup
clear
addpath('../functions');
fig = figure(1);clf
positions = [0.05, 0.12, 0.27, 0.84;
             0.38, 0.12, 0.27,  0.84; 
             0.71, 0.12,  0.24, 0.84];
sz        = size(positions);
fig.Position(3:4) = [1200,360];
fs = 13; %fontsize
for i = 1:sz(1)
    ax(i) = subplot('Position', positions(i,:)); hold on; box on
    ax(i).FontSize = fs;
end

load('../gendata/figure1/figure1-data.mat')
ss = figc_data; %rename for some reason???

%% Make panel (b), showing crevasse depth as a function of time for both ice shelves
colmap = [0, 33, 153; 153,0,63]/255;
colmap = [0, 33, 153; 250,194,5]/255;

ss = figc_data;
 for i =1:2
     plot(ax(2), [ss(i,:).t], [ss(i,:).dimless_crev_depth], 'linewidth', 1.8, 'color',colmap(i,:));
     box(ax(2), 'on');

 end
 ax(2).XLabel.String = 'time (yrs)';
 ax(2).XLabel.Interpreter = 'none';
 ax(2).XTick = 0:100:400;
 ax(2).XLim = [0,400];
 ax(2).YLim = [0,1];
 ax(2).YTick = 0:0.2:1;
 ax(2).YLabel.String = '$d_{c}/H$';
 ax(2).YLabel.Interpreter = 'latex';


%% Make panel (a) showing temperature profiles 
colmap = [0, 33, 153; 153,0,63]/255;
colmap = [0, 33, 153; 250,194,5]/255;

Ts = -22;
Tb = -2;
nl = 5; %number of lines
for i = 1:2
    tt = [ss(i,:).t];
    
    idxp = linspace(1,length(tt), nl); %time indices to plot
    idxp = round(idxp);

    %make a colormap
    if i == 1
        T = [140/255,140/255, 220/255 ; %first colour
        colmap(i,:)];
        T = flipud(T);
    else
        T = [220/255,140/255, 140/255 ; %first colour
        colmap(i,:)];
        T = flipud(T);
        
    end
    x = [0 1];
    map = interp1(x,T,linspace(0,1,nl));

    for it = 1:length(idxp)
        Tp = [ss(i,idxp(it)).temp_profile];
        zz = [ss(i,idxp(it)).zz];

        p = plot(ax(1), Tp-273.15, zz, 'linewidth', 1.5);
        palpha =  it/length(idxp);
        p.Color = [map(it,:)]; %set the color nad alpha

        %put the point on the main line
        scatter( ax(2),[ss(i,idxp(it)).t], [ss(i,idxp(it)).dimless_crev_depth],'filled', 'MarkerFaceColor' , map(it,:), 'MarkerEdgeColor', 'k', ...
            'LineWidth',1);
    end

    % add colorbar
    axnew = axes(); axnew.Visible = 'off';
    cin1 = colorbar(axnew);
    shg
    if i == 1
    cin1.Position = [0.08, 0.15, 0.01, 0.3];
    else
        cin1.Position = [0.07, 0.15, 0.01, 0.3];

    end

    cin1.TickLabels ='';
    colormap(cin1, map)
   
end
ax(1).XLabel.String = 'temperature (C)';
ax(1).YLabel.String = "dimensionless depth";
ax(1).YTick = [0:0.2:1];



%% Make panel (c): enhancement of collapse time for different initial ice thickness and melt rates

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
    plot(ax(3),HH, yy, 'linewidth', 1.5, 'color', colmap(im+2,:));
    %drawnow; pause
end

ax(3).YLim = [0, 5];
ax(3).YTick = 0:1:5;
ax(3).XLim = [200,1000];
ax(3).XTick = 0:200:1000;

ax(3).XLabel.String = 'initial ice thickness (m)';
ax(3).YLabel.String = '$(\tau -  \tau_\infty)/\tau_\infty$';
ax(3).YLabel.Interpreter = 'latex';
box(ax(3), 'on')
%set(ax(4), 'YScale', 'log');

% sort out colorbar
axghost = axes();
axghost.Visible = 'off';

nmd = log10(mm);
nm = nmd(end) - nmd(1); 
cc = colorbar; clabels = 10.^(nmd(1):nmd(end)); clim(log10(clabels([1 end]))); set(cc,'Ticks',log10(clabels),'TickLabels',clabels); 
cc.FontSize =fs;
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