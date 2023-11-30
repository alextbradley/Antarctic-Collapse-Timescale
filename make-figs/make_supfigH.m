% Make supplementary figure "H" of the manuscript, showing (a) basal melt
% rates, (b) thinning rates, (c) strain rates and (d) ice thickness from
% bedmachine on ice shelves
%
%
%% Preliminaries
%data = load('../data/ice_sheet_data.mat', 'xx','yy', 'H', 'm', 'eflow', 'dMdt');


% Notes to self: use (a) to show that basal melt rates are indeed higher in
% Amundsen sea sector and where thermal forcing is high, (c) use to show
% that negative strain rates are fairly limited (we don't need to quantify
% this)
%% Prep Figure
figure(1); clf;
w = 0.3;
positions = [0     0.66  w w;
             0.5   0.66  w w;
             0     0.33  w w;
             0.5   0.33  w w;
             0     0.01  w w;
             0.5   0.01  w w];

positions = positions([5,6,3,4,1,2],:); %shift order

for i = 1:6
    ax(i) = subplot('Position', positions(i,:));
    %ax(i) = subplot(3,2,i);
    ax(i).Visible = 'off';
    hold on
    axis equal
    box on
end

% Basal melt rates
mm = data.m;
%mm = rot90(mm,3);
pl = imagesc(ax(1),data.xx, data.yy, mm);
set(pl, 'AlphaData', ~isnan(mm));
set(ax(1), 'YDir', 'reverse')
c(1) = colorbar(ax(1));
ax(1).CLim = [-5, 5];
cmap1 = cmocean('balance'); cmap1 = cmap1(20:end-20,:);
colormap(ax(1), cmap1);
ax(1).Visible = 'off';
c(1).Position(3) = 0.008;
c(1).Position(4) = 0.25;
c(1).Label.String = 'basal melt rate (m/yr)';
c(1).FontSize = 13;
%contour(ax(1), data.xx, data.yy, gl, [0.5, 0.5], 'k');

% Basal melt rates
mmpos = mm;
mmpos(mm <0.01) = 0.01;

pl = imagesc(ax(2),data.xx, data.yy, mmpos);
set(pl, 'AlphaData', ~isnan(mmpos));
set(ax(2), 'YDir', 'reverse')
c(2) = colorbar(ax(2));
ax(2).CLim = [-5, 5];
colormap(ax(2), cmap1);
ax(2).Visible = 'off';
c(2).Position(3) = 0.008;
c(2).Position(4) = 0.25;
c(2).Label.String = 'basal melt rate (m/yr)';
c(2).FontSize = 13;
%contour(ax(1), data.xx, data.yy, gl, [0.5, 0.5], 'k');


% strain rates
% Basal melt rates
eflow = data.eflow;
pl = imagesc(ax(3),data.xx, data.yy, eflow);
set(pl, 'AlphaData', ~isnan(eflow));
set(ax(3), 'YDir', 'reverse')
c(3) = colorbar(ax(3));
ax(3).CLim = 1e-2*[-2, 2];
cmap3 = cmocean('curl'); cmap3 = cmap3(20:end-20,:);
colormap(ax(3), cmap3);
ax(3).Visible = 'off';
c(3).Position(3) = 0.008;
c(3).Position(4) = 0.25;
c(3).Label.String = 'strain rate (1/yr)';
c(3).FontSize = 13;
% add the grounding line
hold on
gl = ~isnan(data.H);f
%p = contour(ax(3), data.xx, data.yy, gl, [0.5, 0.5], 'k', 'linewidth', 0.25);


%% Ice thickness
h= data.H;
pl = imagesc(ax(4),data.xx, data.yy, h);
set(pl, 'AlphaData', ~isnan(h));
set(ax(4), 'YDir', 'reverse')
c(4) = colorbar(ax(4));
ax(4).CLim = [0,1000];
cmap4 = cmocean('deep'); cmap4 = cmap4(20:end-20,:);
colormap(ax(4), cmap4);
ax(4).Visible = 'off';
c(4).Position(3) = 0.008;
c(4).Position(4) = 0.25;
c(4).Label.String = 'ice thickness (m)';
c(4).FontSize = 13;
% add the grounding line
hold on


%% Ice thinning rate
thinrate= -data.dMdt;
pl = imagesc(ax(5),data.xx, data.yy, thinrate);
set(pl, 'AlphaData', ~isnan(thinrate));
set(ax(5), 'YDir', 'reverse')
c(5) = colorbar(ax(5));
ax(5).CLim = [-3,3];
cmap5 = cmocean('balance'); cmap5 = cmap5(20:end-20,:);
colormap(ax(5), cmap5);
ax(5).Visible = 'off';
c(5).Position(3) = 0.008;
c(5).Position(4) = 0.25;
c(5).Label.String = 'thinning rate (m/yr)';
c(5).FontSize = 13;

%% Plot the mean thinning rate
jdir = dir('../data/ice-shelves/all-shelves/*.mat'); %where shelf files are stored
thinrate_used = nan(size(thinrate));

for i = 1:length(jdir)
    shelf = jdir(i).name;
    shelf = strrep(shelf,'.mat',''); %strip the .mat at the end

    fname = strcat('../data/ice-shelves/all-shelves/' ,shelf, '.mat');
    g = load(fname);
    thinrate_used(g.IN & ~isnan(thinrate)) = -g.mean_dhdt;
end

pl = imagesc(ax(6),data.xx, data.yy, thinrate_used);
set(pl, 'AlphaData', ~isnan(thinrate_used));
set(ax(6), 'YDir', 'reverse')
c(6) = colorbar(ax(6));
ax(6).CLim = [-3,3];
cmap6 = cmocean('balance'); cmap6 = cmap6(20:end-20,:);
colormap(ax(6), cmap6);
ax(6).Visible = 'off';
c(6).Position(3) = 0.008;
c(6).Position(4) = 0.25;
c(6).Label.String = 'thinning rate (m/yr)';
c(6).FontSize = 13;


%%
for i = [1,2,3,4,5,6]
    c(i).Position(1) = c(i).Position(1) + 0.1;
    c(i).Position(3) = 0.005;
    c(i).Position(4) = 0.15;

end
fig = gcf; 
fig.Position(3:4) = [900, 600];
