% Make supplementary figure "A" of the manuscript, showing the Peclet
% number for (a) a lengthscale of 500km and (b) a lengthscale of 50km.
%
%% Preliminaries
clear
addpath('../functions');
f = load('supfigAdata.mat', 'H', 'vx', 'vy');


%% Compute the Peclet number
kappa = 36; %ice thermal diffusivity 
speed = sqrt((f.vx).^2 + (f.vy).^2);
Pe_50 = (f.H).^2 .* speed ./ (50*1e3) ./ kappa;
Pe_500 = (f.H).^2 .* speed ./ (500*1e3) ./ kappa;
Pe_500 = rot90(Pe_500);
Pe_50 = rot90(Pe_50);
%% Make the plot
figure(1); clf; 
ax(1) =subplot(1,2,1); 

pl = imagesc(Pe_500); 
axis equal
set(pl, 'AlphaData', ~isnan(log10(Pe_500)));
hold on

ax(2) = subplot(1,2,2);

axis equal
subplot(1,2,2); pl = imagesc(log10(Pe_50));
set(pl, 'AlphaData', ~isnan(Pe_50));

for i = 1:2
    ax(i).CLim = [-3,3];

    colormap(ax(i), cmocean('balance'));
    ax(i).Visible = 'off';
end
ax(1).Position(1) = 0.2;
c = colorbar(ax(2));
c.Position(1) = 0.93;
c.Position(3) = 0.01;
c.FontSize = 14;
c.Label.String = '$Pe = H^2 |\mathbf{u}|/(L\kappa)$';
c.Label.Interpreter = 'latex';
c.TickLabels = {'10^{-3}','10^{-2}','10^{-1}','10^{0}', '10^{1}', "10^{2}", "10^{3}"};