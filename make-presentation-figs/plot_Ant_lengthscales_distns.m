% Contour plot the lengthscale \ell = kappa / mdot H as a function of mdot
% and H with locations of ice shelves added.
%f  = load('../data/ice_sheet_data.mat');
%%
figure(1); clf;
% make the background image
kappai = 36; %thermal diffusivity
H = linspace(1,1200,1e2); %ice thicknesses
mdot = logspace(-1,2,1e2);
[HH,mm] = meshgrid(H,mdot);
ell = kappai ./ HH ./ mm;

xl = [0,1200];
yl = 10.^[min(log10(mdot)),max(log10(mdot))]; %x and y lims of plot
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
cmap = cmap(20:end-20,:);
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

%%
shelf_names = ["Ross", "Filchner", "Amery", "Ronne", "PineIsland", "PopeSmithKohler","Thwaites", "Getz","Larsen"];
mmean = zeros(1,length(shelf_names));
hmean = zeros(1,length(shelf_names));
mstd = zeros(1,length(shelf_names));
hstd = zeros(1,length(shelf_names));
for i = 1:9
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
    mmean(i) = median((m_iceshelf));
    hmean(i) = median((h_iceshelf));
    plot(hmean(i),mmean(i), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    if shelf_names(i) == "PopeSmithKohler"
        t(i) = text(hmean(i)-240,mmean(i),shelf_names(i), 'FontSize', 18, 'FontName', 'GillSans');

    else
        t(i) = text(hmean(i)+10,mmean(i),shelf_names(i), 'FontSize', 18, 'FontName', 'GillSans');
    end

    % add boxes around
%     mstd(i) = std(m_iceshelf);
%     hstd(i) = std(h_iceshelf);
%     xf = [hmean(i) - hstd(i), hmean(i) + hstd(i),  hmean(i) + hstd(i), hmean(i) - hstd(i)];
%     yf = [mmean(i) - mstd(i), mmean(i) - mstd(i),  mmean(i) + mstd(i), mmean(i) + mstd(i)];
%     yf = max(yf,0.1);
%     fill(xf, yf, 0.5*[1,1,1], 'FaceAlpha',0.3, 'EdgeAlpha',0.7);
     
end %end loop over shelves
fig = gcf;
fig.Position(3:4) = [1034, 639];
%% add contours
hold on
levs = [-1.2, -1.6];
for i = 1:length(levs)
contour(H, mdot, log10(ell), levs(i)*[1,1], 'LineStyle','--', 'Color', 0.7*[1,1,1], 'LineWidth', 2);
end
shg



