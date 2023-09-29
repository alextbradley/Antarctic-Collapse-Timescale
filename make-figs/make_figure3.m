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
f = load('../../data/ice-sheet-data.mat', 'tf_max');
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

ocean_cmap = cmocean('haline', 100);
ocean_cmap = flipud(ocean_cmap(50:end,:));



%% Part (a)
fig1 = figure(1);clf;
fig1.Position(3:4) = [900, 900];
ff = load('figure3-data.mat');
ct = ff.collapse_time; %collapse time data
tags = ff.tags;
inf_col = [1,1,1]; %infinite time colour
misdat_col = 0.6* [1,1,1];

% Background: thermal forcing
axlayer0 = gca; 
tfmax = f.tf_max;
pl0 = imagesc(axlayer0, f.tf_max);
set(pl0, 'AlphaData', ~isnan(tfmax));
colormap(axlayer0, ocean_cmap)
axlayer0.Visible = 'off';
set(axlayer0, "YDir", 'reverse')
axis equal
clear clim; clim([0,4])
c0 = colorbar(axlayer0);
c0.Position(1) = 0.9; %move out of the way
c0.Position(4) = 0.3;
c0.Position(2) = 0.5;
c0.Position(3) = 0.01;
c0.Ticks = 0:4;
c0.FontSize = 12;

% 1: collapse time
%ax = gca;
ax = axes();
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
c.FontSize = 12;
c.TickLabels = {'10^0', '10^1', '10^2', '10^3', '>10^4'};
c.Position(4) = 0.3;
c.Position(2) = 0.1;
c.Position(3) = 0.01;
%colorbar(ax);
axlayer0.Position = ax.Position;
c0.TickLabels = {"0","1","2","3",">4"};


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


linkaxes([axlayer0, ax, axlayer2, axlayer3, axlayer5])
fig = gcf;
%fig.Position(3:4) = [655,507];

colormap(ax, cmap);

%% print the coverage of this map (i.e. what percenatage of points are not missing data)
misdata = zeros(size(tags));
misdata(tags == 2) = 1;

ntot = prod(size(tags)) - sum(sum(tags == 1));  %total number of ice shelf points (tag = 1 corresponds to no thickness data)


per_miss = sum(sum(tags == 2))/ntot * 100;

fprintf('Percentage coverage of data is %.3f percent \n', 100- per_miss)

%% Make panel b: violin plots of timescales
figure(2); clf; hold on; box on

% load the data from figure 2
fig2data = load('fig2_out_.mat');

% restrict only to those shelves we're keeping
shelf_type(17) = 0; %remove the pine island subsection
lk            = fig2data.ll(shelf_type~=0);
hk            = fig2data.h_ave(shelf_type~=0);
ct_avek       = fig2data.ct_ave(shelf_type~=0);
ct_ave_Nye_k  = fig2data.ct_ave_Nye(shelf_type~=0);
shelf_typek   = fig2data.shelf_type(shelf_type~=0);
shelf_countsk = fig2data.shelf_counts(shelf_type~=0);
shelf_namesk  = fig2data.shelf_names(shelf_type~=0);
bar_coldatak  = fig2data.shelf_cols(shelf_type~=0,:);


% put back ground levels on
levs = [80, 180, 280, 980];
for i = 1:4
plot([0,30], levs(i)*[1,1], 'linewidth', 1, 'Color',[80,62,183]/255)
end

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
fig.Position(3:4) = [1200, 240];
ax3 = gca; ax3.FontSize = 14;
ax3.YTick = logspace(0,4,5);

% Adjust the names
shelf_namesks(shelf_namesks ==  "Thwaites") = "THW";
shelf_namesks(shelf_namesks ==  "Crosson") = "CRO";
shelf_namesks(shelf_namesks ==  "PopeSmithKohler") = "PSK";
shelf_namesks(shelf_namesks ==  "Dotson") = "DOT";
shelf_namesks(shelf_namesks ==  "PineIslandFast") = "PIF";
shelf_namesks(shelf_namesks ==  "PineIsland") = "PIG";
shelf_namesks(shelf_namesks ==  "Larsen") = "LAR";
shelf_namesks(shelf_namesks ==  "Wilkins") = "WIL";
shelf_namesks(shelf_namesks ==  "Brunt") = "BRU";
shelf_namesks(shelf_namesks ==  "George6") = "GVI";
shelf_namesks(shelf_namesks ==  "Getz") = "GET";
shelf_namesks(shelf_namesks ==  "Shackleton") = "SHA";
shelf_namesks(shelf_namesks ==  "West") = "WES";
shelf_namesks(shelf_namesks ==  "Cook") = "COO";
shelf_namesks(shelf_namesks ==  "Nansen") = "NAN";
shelf_namesks(shelf_namesks ==  "Cosgrove") = "COS";
shelf_namesks(shelf_namesks ==  "Borchgrevink") = "BOR";
shelf_namesks(shelf_namesks ==  "Abbot") = "ABB";
shelf_namesks(shelf_namesks ==  "SwinburneSulzbergerNickerson") = "SSN";
shelf_namesks(shelf_namesks ==  "FimbulJelbart") = "FIM";
shelf_namesks(shelf_namesks ==  "KingBaudoin") = "BAU";
shelf_namesks(shelf_namesks ==  "TottenMoscow") = "TOT";
shelf_namesks(shelf_namesks ==  "RiiserLarsen") = "RII";
shelf_namesks(shelf_namesks ==  "Amery") = "AME";
shelf_namesks(shelf_namesks ==  "Ronne") = "RON";
shelf_namesks(shelf_namesks ==  "Ross") = "ROSS ";
shelf_namesks(shelf_namesks ==  "Filchner") = "FIL";

ax3.XLim = [0,length(shelf_namesk)+1];
ax3.YLim = [10^0,4*10^4];
ax3.XTick = 1:length(shelf_countsks);
ax3.XTickLabel = shelf_namesks;
ax3.XTickLabelRotation = 45;
ax3.YLabel.String = 'collapse timescale';

%% Make circles showing collapse time 
figure(3); clf; hold on
count = 1;
levs  = logspace(0,4,length(cmap)); %these are the levels of the colourmap in log
for i = 1:27 
    subplot(3,9,count);
    
    diam =2*  sqrt(ct_aveks(i));

    %work out what the colour should be
    [~,idx] = min(abs(ct_aveks(i) - levs) ); %get index of nearest colourmap
    colpt = cmap(idx,:);
    
    
    plot(0,0,'o', 'markerfacecolor', colpt, 'markersize', diam, 'MarkerEdgeColor', 'k', 'LineWidth',1)
    count = count + 1;
    title(shelf_namesks(i))
    ax = gca; ax.Visible = 'off';
    hold on
    text(0,0, shelf_namesks(i))
end
