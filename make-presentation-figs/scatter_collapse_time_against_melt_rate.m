% scatter the collapse time against melt rate
data = load('../make-figs/fig2_out_.mat');
fs = 14;

figure(1);clf; hold on; box on; grid on;
for i = 1:length(data.ave_crevasse_time)
    plot(data.ave_melt_rate(i), data.ave_crevasse_time(i), 'ko', 'markersize', 7, 'markerfacecolor', 'k');
    t(i) = text(data.ave_melt_rate(i)*1.05, data.ave_crevasse_time(i)*1.01,data.shelf_names_adj(i), 'FontSize', fs, 'FontName', 'Arial');

end

ax= gca;
ax.FontSize = fs;
set(ax, 'XScale', 'log');
set(ax, 'YScale', 'log');

% add the linear fit
lm = fitlm(log10(data.ave_melt_rate),log10(data.ave_crevasse_time));
p = table2array(lm_thinning_rate.Coefficients(:,1));

xp = ax.XLim;
plot(xp, 10^p(1)*xp.^p(2), 'k--', 'linewidth', 1.5);

ylabel('Collapse timescale (years)' );
xlabel('Average melt rate (m/yr)')
shg