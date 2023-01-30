% Plot the collapse time for a specific ice shelf using the scattered
% interpolant (advective approximation)

%% Load the ice shelf data and scatter on here
f  = load('../data/ice_sheet_data.mat');
shelf_name= 'PineIsland';
shelf = shelf_name;
fname = strcat('../data/ice-shelves/' ,shelf, '.mat');
%% Load the scattered interpolant and plot
collapse_time = load('collapse_time_scattered_interpolant');
collapse_time = collapse_time.f;

m = logspace(-1,3);
H = linspace(100,1000);
[mm,HH] = meshgrid(m,H);
clf;contourf(m,H, collapse_time(mm,HH), 50, 'linestyle', 'none')
set(gca, 'XScale', 'log');
c =colorbar; 
c.Position(1) = 0.9;
ax = gca;
ax.Position(3) = 0.7;
%%

% Add the shelf data
g = load(fname);
m_iceshelf = f.m;
m_iceshelf = m_iceshelf(g.IN); %only the points in this shelf
h_iceshelf = f.H;
h_iceshelf = h_iceshelf(g.IN);
idx =  (~isnan(h_iceshelf)) &  (~isnan(m_iceshelf) & m_iceshelf > 1e-6); %points where we have point thickness and melt rate > 0
%idx =  (~isnan(h_iceshelf)) &  (~isnan(m_iceshelf)); %points where we have point thickness and melt rate > 0
h_iceshelf = h_iceshelf(idx);
m_iceshelf = m_iceshelf(idx); %arrays with points in particular shelf with both melt and thicknes

hold on
axnew = axes;
[binCounts, xbin, ybin] = histcounts2(log10(m_iceshelf),(h_iceshelf),[100,100]);
% Normalize bin counts to 0:1
binCountsNorm = (binCounts - min(binCounts(:))) ./ (max(binCounts(:)) - min(binCounts(:)));
% Plot the results *
h= histogram2('XBinEdges',10.^(xbin),'YBinEdges',ybin,'BinCounts',binCountsNorm, ...
    'DisplayStyle','tile','ShowEmptyBins','off'); % or you may what "off"
% Add color bar and make sure the color ranges from 0:1
%c2 = colorbar();
clim([0,1])
set(gca, 'XScale', 'log');
axnew.XLim = ax.XLim;
axnew.YLim = ax.YLim;
   colormap(gca, (cmocean('matter')));
   colorbar;
   axnew.Position = ax.Position;
axnew.Visible = 'off';
shg
%% Create a contour plot of collapse timescale

% create restricted co-ordinates
[rId, cId] = find(in) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
xs = f.xx(xminidx:xmaxidx);
ys = f.yy(yminidx:ymaxidx); %restricted co-ordinates
hs = f.H(xminidx:xmaxidx, yminidx:ymaxidx);
ms = f.m(xminidx:xmaxidx, yminidx:ymaxidx);
dhdts = f.dhdtadj(xminidx:xmaxidx, yminidx:ymaxidx);
ms(ms < 1e-1) = 1e-1;
% loop over grid points
collapse_time_square = nan(size(ms));
for ix = 1:length(xs)
    for iy = 1:length(ys)
        if ~isnan(ms(ix,iy)) && ~isnan(hs(ix,iy)) && ~isnan(dhdts(ix,iy)) 
              collapse_time_square(ix,iy) = collapse_time(ms(ix,iy),hs(ix,iy))/(-dhdts(ix,iy));  

        end
    end
end

collapse_time_square(collapse_time_square < 0) = 1e3; %really I should plot this in another colour
colormap(parula)
figure(1); clf; subplot(2,2,1); imagesc(ms); colorbar; title('melt rate'); axis equal
subplot(2,2,2); imagesc(hs); colorbar; title('ice thickness'); axis equal
subplot(2,2,3); imagesc(dhdts); colorbar; title('dhdt'); axis equal
subplot(2,2,4); pl = imagesc(log10(collapse_time_square)); c = colorbar; ; axis equal
clim([0,2]); c.Ticks = [0,1,2]; c.TickLabels = {'10^0', '10^1', '10^2'};
ax = gca; colormap(ax, cmocean('matter'))
set(pl, 'AlphaData', ~isnan(collapse_time_square));

% to do: add the missing data but is in shelf as a colour (e.g. grey) and
% points which are always stable in red
