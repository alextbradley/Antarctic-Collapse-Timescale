% Plot the Peclet number around Antarctica

addpath('../data/');
addpath('../functions')
f = load('../data/ice_sheet_data.mat', 'H', 'm', 'u', 'v');
speed = sqrt((f.u).^2 + (f.v).^2);
kappai = 36; %ice thermal diffusivity
L      = 50*1e3; %lengthscale (50km for WAIS, 500km for 
Peclet = f.H.^2 .* speed ./L ./kappai;
%%
figure(1); clf; 
surf(log10(flipud(Peclet)), 'linestyle', 'none'); view([0,90]);
xticks([])
yticks([]);
c = colorbar;
clim([-2,2]);
colormap(cmocean('balance'))
ax = gca;
ax.Visible = 'off';
c.Ticks = -2:2;
c.TickLabels = {'10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}'};
c.FontSize = 20;
c.FontName = 'GillSans';