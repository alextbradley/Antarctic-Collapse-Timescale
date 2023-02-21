% Make figure 3 of the ms, showing the collapse time in (a) using the raw
% data and (b) using the absolute value of the thinning rate.
%
% ATB (aleey@bas.ac.uk), 17/02/2023. MIT Licence.
%
% (figure made as two individual figures)
%% Preliminaries
%clear
addpath('../functions');
fs = 14; %fontsize
cmap = cmocean('matter', 100);

%% Part (a)
figure(1);clf
ff = load('figure3-data.mat');
ct = ff.collapse_time; %collapse time data
tags = ff.tags;

% 1: collapse time
ax = gca;
pl = imagesc(ax, log10((ct))); 
set(pl, 'AlphaData', ~isnan(ct));
set(ax, "YDir", 'reverse')
axis equal
ax.Visible = 'off';
%pause

% 2: missing data
axlayer2 = axes();
tags_missing = nan(size(tags));
tags_missing(tags == 2) = 1;
pl = imagesc(tags_missing);
set(pl, 'AlphaData', ~isnan(tags_missing)); 
colormap(axlayer2, 0.6* [1,1,1])
axlayer2.Visible = 'off';
set(axlayer2, "YDir", 'reverse')

axis equal
axlayer2.Position = ax.Position;
%pause

% 3: negative strain rate
axlayer3 = axes();
tags_negative_strain =  nan(size(tags));
tags_negative_strain(tags == 3) = 1;
pl = imagesc( tags_negative_strain);
set(pl, 'AlphaData', ~isnan(tags_negative_strain)); 
colormap(axlayer3, cmap(end,:))
%colormap(ax(3),  0.3* [1,1,1])
axlayer3.Visible = 'off';
set(axlayer3, "YDir", 'reverse')
axis equal
axlayer3.Position = ax.Position;

% 4: negative dhdt
axlayer4 = axes();
tags_negative_strain =  nan(size(tags));
tags_negative_strain(tags == 4) = 1;
pl = imagesc( tags_negative_strain);
set(pl, 'AlphaData', ~isnan(tags_negative_strain)); 
colormap(axlayer4, cmap(end,:))
axlayer4.Visible = 'off';
set(axlayer4, "YDir", 'reverse')
axis equal
axlayer4.Position = ax.Position;

colormap(ax, cmap)
linkaxes([ax, axlayer2, axlayer3, axlayer4])
fig = gcf;
fig.Position(3:4) = [655,507];

%% Part (b) -- with approximation to thinning rate
figure(2);clf
clear ax axlayer2 axlayer3 axlayer4
ff = load('figure3-data.mat');
ct = ff.collapse_time; %collapse time data
tags = ff.tags;

% 1: collapse time
ax = gca;
pl = imagesc(ax, log10((ct))); 
c = colorbar(ax);c.Position(1) = 0.9;
c.Limits = [1,4];
set(pl, 'AlphaData', ~isnan(ct));
set(ax, "YDir", 'reverse')
colormap(ax, cmap)
axis equal
ax.Visible = 'off';
%pause

% 2: missing data
axlayer2 = axes();
tags_missing = nan(size(tags));
tags_missing(tags == 2) = 1;
pl = imagesc(tags_missing);
set(pl, 'AlphaData', ~isnan(tags_missing)); 
colormap(axlayer2, 0.6* [1,1,1])
axlayer2.Visible = 'off';
set(axlayer2, "YDir", 'reverse')
ct(tags_missing == 1)=  nan;

axis equal
axlayer2.Position = ax.Position;
%pause

% 3: negative strain rate
axlayer3 = axes();
tags_negative_strain = nan(size(tags));
tags_negative_strain(tags == 3) = 1;
pl = imagesc( tags_negative_strain);
set(pl, 'AlphaData', ~isnan(tags_negative_strain)); 
colormap(axlayer3, cmap(end,:))
%colormap(ax(3),  0.3* [1,1,1])
axlayer3.Visible = 'off';
set(axlayer3, "YDir", 'reverse')
axis equal
axlayer3.Position = ax.Position;
ct(tags_negative_strain == 1)=  nan;

linkaxes([ax, axlayer2, axlayer3])


fig = gcf;
fig.Position(3:4) = [655,507];
c.Ticks = 1:4;
c.FontSize = fs;
c.TickLabels = {'10^1', '10^2', '10^3', '10^4'};

%% Make figure 3c: violin plots of timescales
figure(3); clf; hold on; box on

% load the data from figure 2
fig2data = load('fig2_out.mat');

% restrict only to those shelves we're keeping
lk            = fig2data.ll(shelf_type~=0);
hk            = fig2data.h_ave(shelf_type~=0);
ct_avek       = fig2data.ct_ave(shelf_type~=0);
shelf_typek   = fig2data.shelf_type(shelf_type~=0);
shelf_countsk = fig2data.shelf_counts(shelf_type~=0);
shelf_namesk  = fig2data.shelf_names(shelf_type~=0);
bar_coldatak  = fig2data.shelf_cols(shelf_type~=0,:);


% sort the data in increasing mean
[~,I] = sort(ct_avek);
lks = lk(I);
shelf_namesks = shelf_namesk(I);
shelf_countsks = shelf_countsk(I);
bar_coldataks = bar_coldatak(I,:);
shelf_typeks = shelf_typek(I);
ct_aveks = ct_avek(I); 


w = 0.4; %width of the dist
for i = 1:length(shelf_typeks)

    %get the data as array
    counts = cell2mat(shelf_countsks(i));
    counts(counts < 10) = 10; %make the plotting a bit nicer

    %remove very long timescale pts
    counts = counts(counts < 5*1e3);

    %fit a kde to it
    kde = fitdist(counts,'kernel');

    %evaluate it
    x = logspace(1,4);
    y = pdf(kde,x);

    %scale the pdf
    y = y * w / max(y); 

    % create fill array
    cline = i; %centreline of distribution
    xf = [cline + y, flip(cline - y)];
    yf = [x,flip(x)];

    %fill the data
    fill(xf, yf, bar_coldataks(i,:), 'linewidth', 1)

    % add mean as a point
    plot(cline,mean(kde),'ko', 'markerfacecolor', 'k', 'markersize', 5)
    %or maybe the point of 90 mass or so?
%     [countsort,~] = sort(counts);
%     frac = 0.75;
%     idx = round(length(countsort)*frac); %index frac way along counts
%     plot(cline, countsort(idx),'ko', 'markerfacecolor', 'k', 'markersize', 5)
    


end

set(gca, 'YScale', 'log')

fig = gcf;
fig.Position(3:4) = [1000, 400];
ax3 = gca; ax3.FontSize = 14;

ax3.XLim = [0,length(shelf_namesk)+1];
ax3.XTick = 1:length(shelf_countsks);
ax3.XTickLabel = shelf_namesks;

