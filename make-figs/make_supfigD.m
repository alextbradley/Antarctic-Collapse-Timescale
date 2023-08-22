% Make supplementary figure D of the ms: scatter plots of strain rate and
% thinning rate versus melt rate
%
% ATB (aleey@bas.ac.uk), 21/04/23. MIT licence.
%
%
figure(1); clf; 
%
% load data from figure 2 
%
fig2data = load("fig2_out_.mat");
idx = fig2data.shelf_type >0;
ms     = fig2data.m_ave(idx);
negdhdts  = fig2data.thinrate_ave(idx);
epsxxs = fig2data.epsxx_ave(idx);

%
% scatter these
% 
scatsize = 40;
scatcol  = [0, 33, 163]/255;
lw       = 1;
ax(1) = subplot(1,2,1); 
scatter(ms, negdhdts,  scatsize, shelf_cols(idx,:), 'filled','MarkerEdgeColor','k', 'linewidth', lw);

ax(2) = subplot(1,2,2); 
scatter(ms, epsxxs, scatsize, shelf_cols(idx,:),'filled','MarkerEdgeColor','k', 'linewidth', lw);

for i = 1:2
    ax(i).FontSize = 15;
    box(ax(i), 'on')
    hold(ax(i), 'on')
    ax(i).XLabel.String = 'melt rate (m/yr)';
    
end
 
ax(1).YLabel.String = 'thinning rate (m/yr)';
ax(2).YLabel.String = 'strain rate (1/yr)';


% 
% do the regression
%
lm_dhdt = fitlm(ms,dhdts);
lm_epsxx = fitlm(ms,epsxxs);

%get coeffs
p_dhdts = table2array(lm_dhdt.Coefficients(:,1));
p_epsxx = table2array(lm_epsxx.Coefficients(:,1));

%
% plot regression
%
ml = [0, 15];
plot(ax(2), ml, p_epsxx(1)+p_epsxx(2)*ml, 'k', 'linewidth', 1.5 );
plot(ax(1), ml, p_dhdts(1)+p_dhdts(2)*ml, 'k', 'linewidth', 1.5 );

%
% tidy
%
ax(2).YLim = [0, 0.015];
ax(1).YLim = [0, 8];
ax(1).YTick = 0:2:8;