

% Flowline info:
% AM01: flowline 538, km 362
% AM02: flowline 1276, km 235
% AM03: flowline 805, km 310
% AM04: flowline 538, km 304
% AM05: flowline 538, km 240
% AM06: flowline 1276, km 173

flowlines    = [538, 1276, 805, 538, 538, 1276];
borehole_kms = [362, 235, 310, 304, 240, 173];
labels       = ["AM01","AM02","AM03","AM04","AM05","AM06"];

colors       = [143, 183, 92;
    79, 115, 175;
    179, 70, 70;
    118, 85, 147;
    79, 162, 187;
    228, 170, 61]/255;

Wang_data = load('../gendata/supfigureF/Wang2022Fig6_profiles.mat');
Wang_data = Wang_data.Wang2022data;
%% Get the data
ds = 1e3; %1km resolution
fname = '../data/ice-shelves/all-shelves/Amery.mat'; %shelf name
%data = load('../data/ice_sheet_data.mat');
ss = get_flowline_data(fname, data, ds);

% restrict to the grid
xs = data.yy(xminidx:xmaxidx);
ys = data.xx(yminidx:ymaxidx);
hs = data.H(xminidx:xmaxidx, yminidx:ymaxidx);
vs = data.v(xminidx:xmaxidx, yminidx:ymaxidx);
us = data.u(xminidx:xmaxidx, yminidx:ymaxidx);

%% Plot flowline locations on outline of shelf
figure(1); clf; hold on; axis equal

%outline of the shelf
isshelf = ~isnan(us);
pl = contour(xs, ys, isshelf', [0.5, 0.5], 'k', 'linewidth', 1.5); %outline of shelf

% put the flowlines on
for iF = 1:length(flowlines)
    flowline_idx = flowlines(iF);
    plot(ss(flowline_idx).flowline(:,1), ss(flowline_idx).flowline(:,2), 'k--', 'linewidth', 1.65)
end

% mark the boreholes
ms = 10; %markersize for boreholes
dx = 5e3; %x shift of labels
dy = 5e3; %y shift of labels
for iB = 1:length(borehole_kms)
    plot(ss(flowlines(iB)).flowline(borehole_kms(iB),1), ss(flowlines(iB)).flowline(borehole_kms(iB),2), 'ko', 'markersize', ms, 'markerfacecolor', colors(iB,:), 'MarkerEdgeColor',colors(iB,:))
    text(ss(flowlines(iB)).flowline(borehole_kms(iB),1) + dx,  ss(flowlines(iB)).flowline(borehole_kms(iB),2) + dy, labels(iB), 'FontSize',14)
end

% tidy
ax = gca;
ax.XLim = [1365000, 1640000];
ax.Visible = 'off';

%% Plot the modelled and observed temperatures at boreholes
figure(2); clf;
other_model_color = 0.7*[1,1,1];
ghf = 48;
Ts  = -22;
dzeta = 1e-3; %vertical dimensionless spacing
zeta = dzeta:dzeta:(1-dzeta);


for iF = 1:length(flowlines)
    s = ss(flowlines(iF)); %get the flowline
    flc = s.flowline;
    x = flc- flc(1,:);
    flc = s.flowline;
    x = sqrt(x(:,1).^2 + x(:,2).^2); %along flowline distance
    S = 34.6 *ones(size(x)); %salinity along flowline


    % construct temperature along flowline
    TgfF = get_grounding_line_temp(ghf, Ts, s.h(1));
    [T,z,xx] = get_flowline_temp(x, zeta,s.h,s.speed,S, abs(s.melt), TgfF); %thsi gets the 'actual' Sergienko profile, but we use it for the co-ordinates and melt rate later
    T_exp = zeros(size(T));
    for ix = 1:length(x)
        mm = max(s.melt(ix), .1);
        l  = kappa / mm / s.h(ix);
        anonT = @(z) ((Tb -  TgfF(z))*exp(-z/l)) ;
        TgfF = get_grounding_line_temp(ghf, Ts, s.h(ix));
        T_exp(:, ix) = TgfF(zeta) + (Tb - TgfF(zeta)).*exp(- zeta/l);
    end

        % add other modelled results first
         axs(iF) = subplot(2,3,iF); hold on; box on;
    plot(Wang_data(iF).T_BMB_ISMIP6, Wang_data(iF).depth_BMB_ISMIP6, 'Color',other_model_color, 'linewidth', 1.5)
    plot(Wang_data(iF).T_BMB_ROMS, Wang_data(iF).depth_BMB_ROMS, 'Color',other_model_color, 'linewidth', 1.5)
    plot(Wang_data(iF).T_BMB_CAL, Wang_data(iF).depth_BMB_CAL, 'Color',other_model_color, 'linewidth', 1.5)
    plot(Wang_data(iF).T_BMB_CAL2, Wang_data(iF).depth_BMB_CAL2, 'Color',other_model_color, 'linewidth', 1.5)


   
    plot(T_exp(:, borehole_kms(iF)), zeta*s.h(borehole_kms(iF)) - s.h(borehole_kms(iF)), '--', 'linewidth', 2, 'color', colors(iF,:), 'HandleVisibility', 'off');

    % add the obs
    plot(Wang_data(iF).T, Wang_data(iF).depth, 'Color',colors(iF,:), 'linewidth',2)
    %legendinfo{iF} = labels(iF);


    %tidy
    title(labels(iF));
    axs(iF).FontSize = 13;
    axs(iF).XLabel.String = 'temp (C)';
    axs(iF).YLabel.String = 'depth (m)';
    grid(axs(iF), 'on');
    axs(iF).YTick = flip(0:-100:(min(Wang_data(iF).depth)));
    axs(iF).YLim = [min(Wang_data(iF).depth), 0];
    %axs(iF).XLim = [-25,0];
    %axs(iF).XTick = -25:5:0;

end


