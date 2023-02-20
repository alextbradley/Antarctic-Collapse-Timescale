% Plot the collapse time for a specific ice shelf using the scattered
% interpolant (advective approximation)
addpath('../functions/')
%addpath('profiling/functions-oldgercrev')

% parallel tasks
% delete(gcp('nocreate'))
% local_cluster = parcluster('local')
% parpool(local_cluster)
%% Load the ice shelf data
f  = load('../data/ice_sheet_data.mat');
shelf_name= 'PineIsland';
shelf = shelf_name;
fname = strcat('../data/ice-shelves/all-shelves/' ,shelf, '.mat');
g = load(fname);

%% create restricted co-ordinates
[rId, cId] = find(g.IN) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
step = 1; %step in grid cell. Resolution becomes 500m*step
xs = f.xx(xminidx:step:xmaxidx);
ys = f.yy(yminidx:step:ymaxidx); %restricted co-ordinates
hs = f.H(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
ms = f.m(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
dhdts = f.dMdtadj(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
strains = f.eflow(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
Bs = f.B(xminidx:step:xmaxidx, yminidx:step:ymaxidx);


dhdts = -abs(dhdts);
ms(ms < 1e-1) = 1e-1; %set a minimum melt value
%% determine tags
% tags:
%    1: outside domain (no thickness data)
%    2: other missing data (no melt, strain rate, dhdt)
%    3: negative strain rate
%    4: positive thinning rate (currently stable)
%    5: positive thinning rate (currently unstable)
%    6: other (good data)

tags = 6*ones(size(ms));

tags(dhdts > 0) = 4;
tags(strains < 0) = 3;
tags(isnan(ms) | isnan(dhdts) | isnan(strains)) = 2;
tags(isnan(hs)) = 1;

%
figure(1); clf; %contourf(ys,xs,flipud(tags), 20, 'linestyle', 'none')
imagesc(ys,xs,flipud(tags))
c = colorbar;
caxis([1,6])
colormap(lines(6))
c.Ticks = linspace(2-(7/12),5+(7/12),6);
c.TickLabels = {'1','2','3','4','5','6'};
set(gca, 'YDir', 'normal');
%% compute the collapse time for those applicable
collapse_time_square = nan(size(ms));
tic

%constant parameters
Tb    = -2 + 273.15;     %basal temperature (kelvin)
Ts    = -20;
Ts    = Ts + 273.15;    %surface temp
B0    = 1.928;  %viscosity constant
rhoi  = 918.0;  %ice density
g     = 9.81;   %gravitational acceleration
kappa = 36;     %ice diffusivity
n     = 3;      %glen flow coeff
frac_tough = 150*1e3;
ghf = 40; %geothermal heat flux

parfor ix =  1:length(xs)
%for ix = 1:length(xs)
    %variable parameters
    row_H = hs(ix,:);
    row_dhdt = dhdts(ix,:);
    row_epsxx = strains(ix,:);
    row_tags = tags(ix,:);
    row_mdot = ms(ix,:);


    %get the collapse time on this row
    collapse_time_row = get_collapse_time_row(row_H, row_dhdt, row_epsxx, row_mdot, row_tags, ...
                                                Tb, Ts, B0, rhoi, g, kappa, n, frac_tough, ghf);

    collapse_time_square(ix,:) = collapse_time_row;
    ix;

end
toc

cc = collapse_time_square;
cc = cc(~isnan(cc));
mean(cc)
median(cc)

%% make the plot
figure(1); clf;  
ax(1) = gca; 
collapse_time_square(collapse_time_square == 0) =0.1;
cmap1 = cmocean('matter',100);
%contourf(flipud(log10(collapse_time_square)), 50, 'linestyle', 'none')
pl = imagesc((log10(collapse_time_square)));
set(pl, 'AlphaData', ~isnan(collapse_time_square)); 
axis equal
%
c = colorbar;
colormap(ax(1), cmap1);
cmin = 1; cmax = 4;

clabels = 10.^(cmin:1:cmax); caxis(log10(clabels([1 end]))); set(c,'Ticks',log10(clabels),'TickLabels',clabels); 
c.TickLabels = {'<10^1', '10^2','10^3', '>10^4'};
c.FontSize = 14;
c.FontName = 'GillSans';
c.Label.String = 'time to crevasse (yrs)';
c.Label.FontSize = 16;
ax(1).XTick = [];
ax(1).YTick = [];
ax(1).Visible = 'off';


%% add the missing data
ax(2) = axes();
tags_missing = nan(size(tags));
tags_missing(tags == 2) = 1;
pl = imagesc(tags_missing);
set(pl, 'AlphaData', ~isnan(tags_missing)); 
colormap(ax(2), 0.6* [1,1,1])
axis equal


% % add the negative strain rate as maximum
ax(3) = axes();
tags_negative_strain = nan(size(tags));
tags_negative_strain(tags == 3) = 1;
pl = imagesc( tags_negative_strain);
set(pl, 'AlphaData', ~isnan(tags_negative_strain)); 
colormap(ax(3), cmap1(end,:))
colormap(ax(3), [0,0,1])

axis equal


% add the positive thinning rate as maximum
ax(4) = axes();
tags_pos_thin = nan(size(tags));
tags_pos_thin(tags == 4) = 1;
pl = imagesc(tags_pos_thin);
set(pl, 'AlphaData', ~isnan(tags_pos_thin)); 
colormap(ax(4), cmap1(end,:))
axis equal


% add positive thinning rate, currently unstable as minimum 
ax(5) = axes();
tags_pos_thin = nan(size(tags));
tags_pos_thin(tags == 5) = 1;
pl = imagesc( tags_pos_thin);
set(pl, 'AlphaData', ~isnan(tags_pos_thin)); 
colormap(ax(5), cmap1(1,:))
axis equal


% % add zero collapse time as a minimum
% ax(6) = axes();
% zero_collapse = nan(size(collapse_time_square));
% zero_collapse(collapse_time_square == 0) = 1;
% pl = imagesc(zero_collapse);
% set(pl, 'AlphaData', ~isnan(zero_collapse)); 
% colormap(ax(6), cmap1(1,:))
% axis equal


for i = 1:5
    ax(i).Position = ax(1).Position;
    ax(i).Visible = 'off';

end