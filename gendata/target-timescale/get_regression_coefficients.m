function [s_epsxx, s_dhdt] = get_regression_coefficients(shelf,make_plot)
% return regression coefficients: how much do mean strain rate and mean
% dhdt increase with a unit increase in melting?



f  = load('../../data/ice_sheet_data.mat');
fname = strcat('../../data/ice-shelves/all-shelves/' ,shelf, '.mat');
g = load(fname);

%% create restricted co-ordinates
[rId, cId] = find(g.IN) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns

ms = f.m(xminidx:xmaxidx, yminidx:ymaxidx);
dhdts = f.dMdtadj(xminidx:xmaxidx, yminidx:ymaxidx);
strains = f.eflow(xminidx:xmaxidx, yminidx:ymaxidx);

% ms = f.m;
% dhdts = f.dMdtadj;
% strains = f.eflow; %uncomment for circum-Antarctic values

%% regress dhdt and strain with melt
idx = (~isnan(ms)) & (~isnan(dhdts)) & (~isnan(strains));
m = ms(idx); 
dhdt = dhdts(idx);
epsxx = strains(idx);
pdhdt = polyfit(m,dhdt,1);
pepsxx = polyfit(m,epsxx,1);


if make_plot
    figure(100); clf;
    ax(1) = subplot(1,2,1); scatter(m, dhdt); xlabel('melt'); ylabel('dhdt');
    hold on
    plot(ax(1).XLim, polyval(pdhdt, ax(1).XLim), 'k--', 'linewidth', 2)    
    ax(2) = subplot(1,2,2); scatter(m, epsxx); xlabel('melt'); ylabel('strain rate');
    hold on
    plot(ax(2).XLim, polyval(pepsxx, ax(2).XLim), 'k--', 'linewidth', 2)
end

s_epsxx = pepsxx(1) %expect > 0
s_dhdt = pdhdt(1)   %expect < 0
