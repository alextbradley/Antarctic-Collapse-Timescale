% Make supplementary figure C, showing a comparison between the analytic
% temperature solution of Sergienko et al., and our approximation for the
% cases shown in figure 1.
% 
% 
%% Prelimiaries and figure setup
clear
addpath('../functions');
fig = figure(1);clf
 fig.Renderer = 'painters'; %help with saving
fig.Position(3:4) = [1200, 540];
for i = 1:6
    ax(i) = subplot(3,2,i); hold on; box on;
    ax(i).FontSize = 14; 
end


%% Panels (a), (b): temperatures from Sergienko et al.
load('../gendata/figure1/figure1-data.mat')
cmap = cmocean('thermal',100);
cmap = cmap(1:end-10,:);

for i = 1:2
    %load data in and surf colours
    T = cell2mat(Tflowline_act(i));
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
%ax(1).XLabel.String = 'distance along flowline (km)';
ax(1).XLim = [0,1000];
ax(2).XLim = [0,60];
% add the colorbar
c = colorbar(ax(1));
% c.Position(1) = 0.92;
% c.Position(3) = 0.006;
c.Location = 'south';
c.Label.String = 'temperature (C)';
c.Ticks = [-22, -12, -2];
c.Limits = [-22, -2];
c.Position(1) = 0.35;
c.Position(3) = 0.1;
%% Panels (c), (d): our approximation to the temperature

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

    surf(ax(i+2), xx(1:nsy:end,1:nsx:end)/1e3,z(1:nsy:end,1:nsx:end), T(1:nsy:end,1:nsx:end), 'linestyle', 'none'); view(ax(i), [0,90])
    % have to plut only 1 in 20 values to make size manageable
 
    %tidy plot
    ax(i+2).YTick = -1000:500:0;
    ax(i+2).YLim = [-1000, 200];
    ax(i+2).YLabel.String = 'depth (m)';
    clim([-20,-2])
    colormap(ax(i+2),cmap);
    grid(ax(i+2), 'on');
    box(ax(i+2), 'on');
    

end
ax(3).XLim = [0,1000];
ax(4).XLim = [0,60];

%% Panels (e), (f): difference between the two profiles
diff_cmap = cmocean('balance');
diff_cmap = diff_cmap(30:end-30, :);
for i = 1:2
    %load data in and surf colours
    T = cell2mat(Tflowline_act(i)) - cell2mat(Tflowline(i));
    
    % work out the abs mean diff
    TT = abs(T(~isnan(T)));
    txt(i) = mean(mean(TT));

    xx    = cell2mat(xflowline(i));
    z     = cell2mat(zflowline(i));
    if i == 1 %Ross
        nsy = 30;
        nsx = 2; %how finely to sampel in x and y
    else
        nsy = 1;
        nsx = 1; %how finely to sampel in x and y
    end


    surf(ax(i+4), xx(1:nsy:end,1:nsx:end)/1e3,z(1:nsy:end,1:nsx:end), T(1:nsy:end,1:nsx:end), 'linestyle', 'none'); view(ax(i), [0,90])
    %tidy plot
    ax(i+4).YTick = -1000:500:0;
    ax(i+4).YLim = [-1000, 200];
    ax(i+4).YLabel.String = 'depth (m)';
    clim([-4,4])
    colormap(ax(i+4),diff_cmap);
    grid(ax(i+4), 'on');
    box(ax(i+4), 'on');

end
ax(5).XLabel.String = 'distance along flowline (km)';
ax(6).XLabel.String = 'distance along flowline (km)';
ax(5).XLim = [0,1000];
ax(6).XLim = [0,60];
cd = colorbar(ax(6));
cd.Location = 'south';
cd.Label.String = 'temperature difference (C)';

cd.Position(1) = c.Position(1);
cd.Position(3) = c.Position(3);
cd.Position(4) = c.Position(4);

%% tidy stuff
for i = 1:3
    ax(2*i - 1).XLim = [-2, 1000]; %so we can still see the box
    ax(2*i - 1).XLim = [0, 1000]; 
    ax(2*i).XLim = [-0.2, 60]; %so we can still see the box
    ax(2*i).XLim = [0, 60];
end

for i = 1:2
    ax(2*i-1).XTickLabel = {};
    ax(2*i-1).XTick = ax(5). XTick;
end

for i = 1:2
    ax(2*i).XTickLabel = {};
    ax(2*i).XTick = ax(6). XTick;
end

for i = 1:2
    ax(i).Position(2) = 0.8; 
end

% put panels c and d just above e and f
gap = 0.02; 
for i = 3:4
    ax(i).Position(2) = ax(5).Position(2) + ax(5).Position(4) + gap;
end

%put panels a and b just above c and d
for i = 1:2
    ax(i).Position(2) = ax(3).Position(2) + ax(3).Position(4) + gap;
end

%adjust the colormap 
c.Position(2) = ax(2).Position(2);
c.Position(4) = 0.015;
c.Position(2) = 0.59;
cd.Position(4) = 0.015;
cd.Ticks = [-4,0,4];


% move panels left a bit
for i = 1:3
    ax(2*i).Position(1) =0.48;
end

for i = 1:3
    ax(2*i-1).Position(1) = 0.07;
end

c.Position(1) = 0.29;
cd.Position(1) = 0.29;

% add the rmse labels
text(ax(5), 300, -850, ['MAE: ',  sprintf('%.2f', txt(1)), ' C'], 'FontSize', 14)
text(ax(6), 23, -850, ['MAE: ', sprintf('%.2f', txt(2)), ' C'], 'FontSize', 14)