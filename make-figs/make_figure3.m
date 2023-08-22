% Make figure 3 of the ms, showing the (a) collapse time  and  (b) violin
% plots of the collapse time, showing the mean for both LEFM and NYE
% theory. 
%
% ATB (aleey@bas.ac.uk), 17/02/2023. MIT Licence.
%
% (figure made as three individual figures). NB: need to run figure2 first
% to get the statistics for (b).
%% Preliminaries
%clear
addpath('../functions');
fs = 14; %fontsize
cmap = flipud(cmocean('matter', 100));
%cmap = cmap(1:end-20,:);
%close all

% colormap a la Reese at al 2017
T = [0,   0,   0          %// dark
     249, 235, 81         %// yellow
     229, 122, 53        %// red  
     178, 40, 41           %dark red
     255, 255, 255        %// black
     255, 255, 255]./255; %// black again  -> note that this means values between 161 and 255 will be indistinguishable

x = [0
     floor(255/4)*1
     floor(255/4)*2
     floor(255/4)*3
     254
     255];

map = interp1(x/255,T,linspace(0,1,255));
%cmap = map;

%% Part (a)
fig1 = figure(1);clf;
fig1.Position(3:4) = [900, 900];
ff = load('figure3-data.mat');
ct = ff.collapse_time; %collapse time data
tags = ff.tags;
inf_col = [1,1,1]; %infinite time colour
misdat_col = 0.6* [1,1,1];


% 1: collapse time
ax = gca;
pl = imagesc(ax, log10((ct))); 
set(pl, 'AlphaData', ~isnan(ct));
set(ax, "YDir", 'reverse')
axis equal
ax.Visible = 'off';
clear clim; clim([0,4])
%pause
c = colorbar(ax);
c.Position(1) = 0.9; %move out of the way
c.Ticks = 0:4;
c.FontSize = 14;
c.TickLabels = {'10^0', '10^1', '10^2', '10^3', '10^4'};
c.Position(4) = 0.4;
c.Position(2) = 0.4;
c.Position(3) = 0.01;
%colorbar(ax);


% 2: missing data
axlayer2 = axes();
tags_missing = nan(size(tags));
tags_missing(tags == 2) = 1;
pl = imagesc(tags_missing);
set(pl, 'AlphaData', ~isnan(tags_missing)); 
colormap(axlayer2, misdat_col)
axlayer2.Visible = 'off';
set(axlayer2, "YDir", 'reverse')

axis equal
axlayer2.Position = ax.Position;
%pause

% 3: infinite collapse time
axlayer3 = axes();
infcollapse =  nan(size(ct));
infcollapse(isinf(ct)) = 1;
infcollapse(~isinf(ct)) = nan;
pl = imagesc( infcollapse);
set(pl, 'AlphaData', ~isnan(infcollapse)); 
colormap(axlayer3, inf_col)

axlayer3.Visible = 'off';
set(axlayer3, "YDir", 'reverse')
axis equal
axlayer3.Position = ax.Position;


%outline of shelves
axlayer5 = axes();
isshelf = ~isnan(f.H');

contour(isshelf', [0.5, 0.5],'k', 'linecolor', 0.3*[1,1,1] , 'linewidth', 0.75);
set(axlayer5, 'YDir', 'reverse')
axlayer5.Visible = 'off';
axis equal


linkaxes([ax, axlayer2, axlayer3, axlayer5])
fig = gcf;
%fig.Position(3:4) = [655,507];

colormap(ax, cmap);



%% Make panel b: violin plots of timescales
figure(2); clf; hold on; box on

% load the data from figure 2
fig2data = load('fig2_out_.mat');

% restrict only to those shelves we're keeping
lk            = fig2data.ll(shelf_type~=0);
hk            = fig2data.h_ave(shelf_type~=0);
ct_avek       = fig2data.ct_ave(shelf_type~=0);
ct_ave_Nye_k  = fig2data.ct_ave_Nye(shelf_type~=0);
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
ct_ave_Nye_ks = ct_ave_Nye_k(I); 

w = 0.3; %width of the dist
for i = 1:length(shelf_typeks)

    %get the data as array
    counts = cell2mat(shelf_countsks(i));
    counts(counts == 0) = 1; %make the plotting a bit nicer
    %counts = counts(counts < 5e4); %remove very long (this is done by the
    %make figure2 script)

    %fit a kde to it
    if ~isempty(counts)
    kde = fitdist(counts,'kernel');

    %evaluate it
    x = logspace(0,5, 5e2);
    y = pdf(kde,x);

    %scale the pdf
    y = y * w / max(y); 

    % create fill array
    cline = i; %centreline of distribution
    xf = [cline + y, flip(cline - y)];
    yf = [x,flip(x)];

    %fill the data
    fill(xf, yf, bar_coldataks(i,:), 'linewidth', 1, 'EdgeColor', 0.2*[1,1,1])

    % add mean as a point
    plot(cline,ct_aveks(i),'ko', 'markerfacecolor', 0*[1,1,1], 'markersize', 6,'LineWidth',1.1 )

    % add the Nye result
    plot(cline,ct_ave_Nye_ks(i),'ko','marker','o' ,'markerfacecolor',  1*[1,1,1], 'markersize', 6, 'LineWidth',1.1 )

    end
end

set(gca, 'YScale', 'log')

fig = gcf;
fig.Position(3:4) = [1200, 350];
ax3 = gca; ax3.FontSize = 14;

shelf_namesks(shelf_namesks ==  "PopeSmithKohler") = "PSK";
shelf_namesks(shelf_namesks ==  "SwinburneSulzbergerNickerson") = "SSN";
ax3.XLim = [0,length(shelf_namesk)+1];
ax3.YLim = [10^0,4*10^4];
ax3.XTick = 1:length(shelf_countsks);
ax3.XTickLabel = shelf_namesks;
ax3.XTickLabelRotation = 45;
ax3.YLabel.String = 'collapse timescale';
