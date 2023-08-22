% Make figure 2 of the ms:
% (a) contour plot of l = kappa/(mH), with ice shelves labelled
% (b) locations of ice shelves in (l, strain rate,  dhdt) space.
%
% This code also produces figures 3c, showing distributions and mean
% of the collapse timescale for individual shelves
%
% 17/02/23, ATB (aleey@bas.ac.uk), MIT licence
%% Preliminaries
clear
addpath('../functions');
saveout = 1; %flag to specify saving output. (Need output from this script to run figure3 and 4 scripts)
f = load('../data/ice_sheet_data.mat');
fig = figure(1);clf
positions = [0.05, 0.12, 0.39, 0.83;
             0.58, 0.12, 0.38, 0.83];
sz        = size(positions);
fig.Position(3:4) = [1040,430];
fs = 14; %fontsize
for i = 1:2
    ax(i) = subplot('Position', positions(i,:)); hold on
    ax(i).FontSize = fs;
end

%colourmap
cmap = cmocean('matter', 100);
cmap = flipud(cmap(1:end-60,:));
colormap(ax(1) ,cmap)

cmap = flipud(cmocean('balance', 100)); cmap = cmap(20:end-20,:);
%cmap = flipud(cmocean('matter', 100)); %cmap = cmap(10:end,:);
colormap(ax(1), cmap);


warm_col = [203,0,63]/255; %warm shelves colour
cold_col = [0, 63, 203]/255; %cold shelves colour
%% Make panel a

%
% Make the contour map
%
kappai = 36; %thermal diffusivity
H = linspace(50,1000,100); %ice thicknesses
mdot = logspace(-1.5,log10(20),100);
[HH,mm] = meshgrid(H,mdot);
ell = kappai ./ HH ./ mm;

% contour plot
subplot('Position', positions(1,:))
%contourf(H, mdot, log10(ell), 100, 'linestyle', 'none');
pl= imagesc(H,mdot, log10(ell)); 
%set(pl, 'AlphaData', 0.9*ones(size(ell)))
set(ax(1), 'YDir', 'normal');
set(ax(1), 'YScale', 'log');

%labels
xlabel('ice thickness (m)', 'interpreter', 'none')
ylabel('melt rate (m/yr)', 'interpreter', 'none')

%colourbar
cc = colorbar;
cl = [-2.7,0];
clim = cl;
clabels = 10.^(ceil(min(clim):1:max(clim))); set(cc,'Ticks',log10(clabels),'TickLabels',clabels);
ax(1).CLim = cl;
cc.Label.Interpreter = 'latex';
cc.Label.String = '$\ell$';
cc.Label.FontSize = 15;
cc.Position(1) = 0.45;

% sort out limits
xl = [50,1000];
yl = 10.^[min(log10(mdot)),max(log10(mdot))]; %x and y lims of plot
ax(1).XTick = (0:200:1000);
ax(1).XLim = [xl(1)-2, xl(2)];
ax(1).YLim = [yl(1), yl(2)+0.1];
box(ax(1), 'on')

%
% add contours
%
hold on
levs = [-1.1, -1.5];
for i = 1:length(levs)
contour(H, mdot, log10(ell), levs(i)*[1,1], 'LineStyle','--', 'Color', 0.5*[1,1,1], 'LineWidth', 2);
end


%
%%%%%% add shelf points %%%%%%%
%

% init storage
jdir = dir('../data/ice-shelves/major/*.mat'); %where shelf files are stored

%we'll also work out shelf by shelf collapse data here
fig3data = load('figure3-data.mat');
ct = fig3data.collapse_time; %collapse time data
ct_Nye = fig3data.collapse_time_Nye;
tags = fig3data.tags;

%initialize storage
m_ave = zeros(1,length(jdir));
h_ave = zeros(1,length(jdir));
epsxx_ave = zeros(1,length(jdir));
thinrate_ave = zeros(1,length(jdir));


percov = zeros(1,length(jdir)); %store the % coverage of input data
collapse_coverage  =zeros(1,length(jdir)); %percentage of points with non-nan or non-inf output
shelf_type = zeros(1,length(jdir)); %flags for keeping the shelf or not (0), cold (1), warm (2)

shelf_counts = cell(1,length(jdir)); %for storing individual shelf counts for collapse time

ct_ave = zeros(1,length(jdir)); %mean values of the collapse time for LEFM
ct_ave_Nye = zeros(1,length(jdir)); %mean values of the collapse time for Nye
ll = zeros(1,length(jdir)); %store the lengthscales l
shelf_cols = zeros(length(jdir), 3);%store plot colours
area_shelf = zeros(1,length(jdir)); %store shelf areas


% loop over shelves
shelf_names = strings;
for i = 1:length(jdir)
    %load data for this shelf
    shelf = jdir(i).name;
    shelf = strrep(shelf,'.mat',''); %strip the .mat at the end
    shelf_names(i) = shelf;
    fname = strcat('../data/ice-shelves/major/' ,shelf, '.mat');
    g = load(fname);

    % restrict data to those points in shelf
    m_shelf = f.m; 
    m_shelf = m_shelf(g.IN); %only the points in this shelf
    m_shelf(m_shelf < 0) = 0; %add hoc adjustment
    h_shelf = f.H; 
    h_shelf = h_shelf(g.IN);
    epsxx_shelf = f.eflow; 
    epsxx_shelf = epsxx_shelf(g.IN);

    %restrict to where we have data
    idx =  (~isnan(h_shelf)) &  (~isnan(m_shelf))  & (~isnan(epsxx_shelf)); %points where we have data for melt, thichkness and strain rate
    h_shelf = h_shelf(idx);
    m_shelf = m_shelf(idx); %arrays with points in particular shelf with both melt and thicknes
    epsxx_shelf = epsxx_shelf(idx);
    epsxx_shelf = epsxx_shelf(epsxx_shelf > 0); % keep only positive strain rates

    %compute the shelf area as a sanity check?
    area_shelf(i) = length(h_shelf(~isnan(h_shelf))); %size in km^3 (grid size = 1e3 * 1e3)

    %get means of 
    m_ave(i) = median((m_shelf));
    h_ave(i) = median((h_shelf));
    epsxx_ave(i) = median(epsxx_shelf);
    thinrate_ave(i) = g.mean_dhdt;  %set the thinning rate to pre-determined value

    % compute the collapse time for this shelf 
    ct_shelf = ct(g.IN); %only the points in this shelf
    ct_keep = (ct_shelf(~isnan(ct_shelf))); %all non nan and non infpoints in shelf
    ct_keep = ct_keep(~isinf(ct_keep));
    ct_keep = ct_keep(ct_keep < 5e4); %remove very long points

    shelf_counts{i} = ct_keep;
    if ~isempty(ct_keep)
        pd = fitdist(ct_keep,'kernel');
        ct_ave(i) = mean(pd);
    end

    %compute collapse coverage
    h = f.H;
    hs = h(g.IN);
    collapse_coverage(i) = sum(sum(~isnan(ct_shelf)))/sum(sum(~isnan(hs))) * 100; %number of points with non nan collapse time

    
    %compute mean for Nye
    ct_Nye_shelf = ct_Nye(g.IN);
    ct_Nye_keep = (ct_Nye_shelf(~isnan(ct_Nye_shelf))); %all non nan points in shelf
    ct_Nye_keep = ct_Nye_keep(~isinf(ct_Nye_keep));
    ct_Nye_keep  = ct_Nye_keep(ct_Nye_keep < 5e4);

    if ~isempty(ct_Nye_keep)
        pd = fitdist(ct_Nye_keep,'kernel');
        ct_ave_Nye(i) = mean(pd);
    end
    

    % plot point and add name
    shelf_keep  = sum(idx); %number of points with data for thickness, melt, strain
    percov(i) = shelf_keep/area_shelf(i) * 100;
    threshold_keep = 0;
    l = kappai / h_ave(i) / m_ave(i);
    ll(i) = l;
    if (percov(i) > threshold_keep) || shelf == "Thwaites" %only plot if > threshold% coverage
        
        
        if l > 10^(max(levs)) %cold shelves
            shelf_type(i) = 1;
            ptcol = cold_col;
            t(i) = text(h_ave(i)-81,m_ave(i),shelf, 'FontSize', fs, 'FontName', 'Arial');
        else
            shelf_type(i) = 2;
            ptcol = warm_col;
            
            t(i) = text(ax(1),h_ave(i)+11,m_ave(i),shelf, 'FontSize', fs, 'FontName', 'Arial');
        end

        plot(ax(1), h_ave(i),m_ave(i), 'o', 'MarkerFaceColor', 0.8*[1,1,1], 'MarkerEdgeColor', 'k', 'MarkerSize', 6.5);

    else
        h_ave(i) = nan; m_ave(i) = nan; %force no plot
    end
    % work out what the colour should be 
    fr = (log10(l) - min(cl))/diff(cl); %fraction of way along colourmap
    fr = min(fr,1);
    fr = max(fr,0);
 
    if ~isnan(fr)
    shelf_cols(i,:) = cmap(round(fr*length(cmap)),:);
    end
    
         
end %end loop over shelves


%% Make panel b

ell_ave = kappai./ h_ave ./m_ave; %mean of the lengthscales
idx = shelf_type > 0;
ls = sum(idx);

%plot setup
ms = 50; %markersize
ax(2).XLabel.String = '$\ell$'; ax(2).XLabel.Interpreter = 'latex';
ax(2).YLabel.String= 'strain rate (1/yr)';
ax(2).ZLabel.String = 'inverse thinning rate (yr/m)';
set(ax(2), 'XScale', 'log')
set(ax(2), 'YScale', 'log')
set(ax(2), 'ZScale', 'log')
ax(2).View = [45,36];
grid(ax(2), 'on')


ax(2).XLim = [min(ell_ave(idx))*0.9,max(ell_ave(idx)) *1.1];            %ell limits
ax(2).YLim = [min(epsxx_ave(idx))*0.9,max(epsxx_ave(idx))*1.1 ];        %strain rate limits
ax(2).ZLim = [1./max(thinrate_ave(idx)) * 0.9,1./min(thinrate_ave(idx))*1.1 ]; %inverse thinning rate limits

 hold(ax(2), 'on');
 box(ax(2), 'on')
lw = 1;

%project onto the l axis
scatter3(ax(2),min(ax(2).XLim)*ones(1,ls), epsxx_ave(idx),1./thinrate_ave(idx),ms, shelf_cols(idx,:), 'Filled', 'MarkerEdgeColor','k', 'linewidth', lw);

% project onto the epsxx axis
scatter3(ax(2), ell_ave(idx), max(ax(2).YLim)*ones(1,ls),1./thinrate_ave(idx),ms, shelf_cols(idx,:), 'Filled', 'MarkerEdgeColor','k', 'linewidth', lw);

%project onto the thinning rate axis
scatter3(ax(2), ell_ave(idx),  epsxx_ave(idx) ,min(ax(2).ZLim)*ones(1,ls),ms,shelf_cols(idx,:), 'Filled', 'MarkerEdgeColor','k', 'linewidth', lw);
shg


%% Final tidying
for i = 1:2
    ax(i).FontSize = fs;
    ax(i).XLabel.FontSize = fs+2;
    ax(i).YLabel.FontSize = fs+2;
    ax(i).FontName = 'Arial';
end

    ax(2).ZLabel.FontSize = fs+2;



%% Save data for use in figures 3 and 4
if saveout
    save('fig2_out_.mat', 'shelf_names', 'shelf_type', 'shelf_cols','ll', 'h_ave', 'ct_ave',"shelf_counts", "ct_ave_Nye", "m_ave", 'thinrate_ave', "epsxx_ave")
end


