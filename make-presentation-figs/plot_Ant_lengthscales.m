% Plot (a) the lengthscale \ell = Kappa / (mdot * H) as a function of mdot
% and H and (b) plots of ice shelf locations on these maps.

%% Parameters
kappai = 36; %thermal diffusivity
H = linspace(1,1200,1e2); %ice thicknesses
mdot = logspace(-2,2,1e2);
[HH,mm] = meshgrid(H,mdot);
ell = kappai ./ HH ./ mm;

xl = [0,1200];
yl = 10.^[-2,2]; %x and y lims of plot
% Make plot
figure(1); clf;
contourf(H, mdot, log10(ell), 100, 'linestyle', 'none');
xlabel('$H$', 'interpreter', 'latex')
ylabel('$\dot{m}$', 'interpreter', 'latex')
cc = colorbar;
ax = gca;
ax.FontSize = 16;
clabels = 10.^(-3:1); clim(log10(clabels([1 end]))); set(cc,'Ticks',log10(clabels),'TickLabels',clabels);
cc.Label.Interpreter = 'latex';
cc.Label.String = '$\ell$';
cc.Label.FontSize = 20;
ax = gca; ax.Position(3) = 0.65;
ax.FontSize = 14;
%cc.Position(1) = cc.Position(1)-(1e-5);
cmap = cmocean('ice', 100);
cmap = cmap(10:end-10,:);
colormap(cmap)
set(gca, 'YDir', 'normal');
set(gca, 'YScale', 'log');
hold on
%contour(H, mdot, log10(ell), -3:2, 'Color', 0.6*[1,1,1], 'linewidth', 1)
ax =gca;
ax.FontSize = 16;
ax.FontName = 'GillSans';
xticks(0:200:1200)
xlim(xl);
ylim(yl);
%contour(H, mdot, log10(ell), [-1,-1], 'k'); %

%% Make the same plot with the shelves on
figure(2); clf;

shelf_names = ["Ross", "Filchner", "Amery", "Ronne", "PineIsland", "PopeSmithKohler","Thwaites", "Getz","Larsen"];
f  = load('../data/ice_sheet_data.mat');
for i = 1:9
    subplot(3,3,i); hold on; box on
    %add the background
    contourf(H, mdot, log10(ell), 100, 'linestyle', 'none');
    xlabel('');
    ylabel('');
    set(gca, 'YDir', 'normal');
    set(gca, 'YScale', 'log');
    ax(i) = gca;
    
    colormap(ax(i), cmap);
    ax(i).XLim = xl;
    ax(i).YLim = yl;
    ax(i).YTick = 10.^(-2:2:2);
    ax(i).XTick = 0:400:1200;
    ax(i).FontSize = 15;
    ax(i).FontName = 'GillSans';
    clim([-3,1])
   
    shelf = shelf_names(i);
    fname = strcat('../data/ice-shelves/' ,shelf, '.mat');
    %title(shelf, 'FontName', 'GillSans');
    g = load(fname);
    m_iceshelf = f.m;
    m_iceshelf = m_iceshelf(g.IN); %only the points in this shelf
    h_iceshelf = f.H;
    h_iceshelf = h_iceshelf(g.IN);
    idx =  (~isnan(h_iceshelf)) &  (~isnan(m_iceshelf) & m_iceshelf > 1e-6); %points where we have point thickness and melt rate > 0
    %idx =  (~isnan(h_iceshelf)) &  (~isnan(m_iceshelf)); %points where we have point thickness and melt rate > 0
    h_iceshelf = h_iceshelf(idx);
    m_iceshelf = m_iceshelf(idx); %arrays with points in particular shelf with both melt and thicknes
   
    ax2 = axes;
    box on
    [binCounts, xbin, ybin] = histcounts2((h_iceshelf),log10(m_iceshelf),[100,100]);
    % Normalize bin counts to 0:1
    binCountsNorm = (binCounts - min(binCounts(:))) ./ (max(binCounts(:)) - min(binCounts(:)));
    % Plot the results *
    h= histogram2('XBinEdges',xbin,'YBinEdges',10.^(ybin),'BinCounts',binCountsNorm, ...
        'DisplayStyle','tile','ShowEmptyBins','off'); % or you may what "off"
    % Add color bar and make sure the color ranges from 0:1
    %c2 = colorbar();
    clim([0,1])
    grid off
    ax2.Visible = 'off';
    set(gca, 'YScale', 'log')
    colormap(gca, (cmocean('matter')));
    xlim(ax(i).XLim)
    ylim(ax(i).YLim)
    ax2.Position = ax(i).Position;
    ax2.XLim = [xl(1), xl(2)];
    hold on
    contour(H, mdot, log10(ell), [-1,-1], 'Color',0.05*[1,1,1], 'linewidth', 1.5) %add 0.1 contour (remember log!)
%     drawnow
%     pause
end %end loop over shelves
fig = gcf;
fig.Position(3:4) = [1255, 760];


