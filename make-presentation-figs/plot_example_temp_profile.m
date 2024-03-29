% Make a plot of an example temperature profile for a warm (PIG) ice shelf

%% Prelimnaries
%clear
addpath('../functions/')
%% Data info
%data= load('../data/ice_sheet_data.mat');
ds = 1e3; %flowline spacing
fname = '../data/ice-shelves/Ross.mat';
ss = get_flowline_data(fname, data, ds);
%nf = [50]; %flowline number (for pig, nf = 50 suggested, for Ross, nf =24)
nf = 24;
%% plot the ice velocity and flowline co-ordinates along flowline
f = load(fname);
in = f.IN;
[rId, cId] = find(in) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
xs = data.xx(xminidx:xmaxidx);
ys = data.yy(yminidx:ymaxidx); %restricted co-ordinates
vs = data.v(xminidx:xmaxidx, yminidx:ymaxidx);
us = data.u(xminidx:xmaxidx, yminidx:ymaxidx);
speeds   = sqrt(us.^2 + vs.^2);                 %ice velocities on restricted grid

%%
% plot speeds
%
figure(1); clf;
%contourf(xs, ys, log10(speeds'), 100, 'linestyle', 'none');
pl = imagesc(xs,ys, log10(speeds'));
set(pl, 'AlphaData', ~isnan(speeds'));
set(gca, 'YDir', 'normal');
clim([1,4])
axis equal
c = colorbar;
c.FontSize = 15;
c.FontName = 'GillSans';
c.Position(1) =0.85;
colormap((cmocean('matter', 100)))
colormap(0.2*[1,1,1])
ax = gca; ax.Visible = 'off';
c.Ticks = 1:4;
c.TickLabels = {'10^1','10^2','10^3','10^4'};

% add the flowline
hold on
flc = ss(nf).flowline;
hold on
plot(flc(:,1), flc(:,2), 'w', 'linewidth', 1.75)

%% Construct and plot the temperature profile
%get distance along flowline
x = flc- flc(1,:);
x = sqrt(x(:,1).^2 + x(:,2).^2);
dzeta = 1e-4; %vertical dimensionless spacing
zeta = dzeta:dzeta:(1-dzeta);
S = 34.6 *ones(size(x));
TgfF = @(z) -2- 18*z;% assume a linear profile between -2 at base and -20 at surface
Tgf = TgfF(zeta);


Ts  = ss(nf).temp(1) - 273.15;
ghf = 40; 
TgfF = get_grounding_line_temp(ghf, Ts, ss(nf).h(1));


[T,z,xx] = get_flowline_temp(x, zeta,ss(nf).h,ss(nf).speed,S, abs(ss(nf).melt), TgfF);

% Contour plot the temp
figure(2); clf
surf(xx/1e3,z, T, 'linestyle', 'none'); view([0,90])
colorbar
colormap(cmocean('thermal'))

ax = gca;
ax.FontSize = 15;
yticks(-800:200:200);
ylim([-800, 200])
grid off;
box on;
xlabel('distance along flowline (km)', 'FontSize', 16);
ylabel('depth (m)', 'FontSize', 16);
ax.FontName = 'GillSans';
ax.FontSize = 15;
c = colorbar; 
c.Label.String = 'temp (C)';
c.Label.FontSize = 16;
clim([-20,-2])
title('sergienko et al. flowline temp')
%% Repeat this with the approximate
% This profile assumes that there is no vertical advection, i.e. the
% grounding line profile advects horizontally. 
T_exp = zeros(size(T));


kappa = 36;
for ix = 1:length(x)
    Ts = -20;
    Tb = -2;
    mm = max(s.melt(ix), 0.1); 
    l  = kappa / mm / s.h(ix);

    ll(ix) = l;
    anonT = @(z) ((Tb -  TgfF(z))*exp(-z/l)) ;
    Tstar = -2; %basal temp 
   
    T_exp(:, ix) = TgfF(zeta) + (Tstar - TgfF(zeta)).*exp(- zeta/l);

end

figure(3); clf
surf(xx/1e3,z, T_exp, 'linestyle', 'none'); view([0,90])
colorbar
colormap(cmocean('thermal'))

ax = gca;
ax.FontSize = 15;
yticks(-800:200:200);
ylim([-800, 200])
grid off;
box on;
xlabel('distance along flowline (km)', 'FontSize', 16);
ylabel('depth (m)', 'FontSize', 16);
ax.FontName = 'GillSans';
ax.FontSize = 15;
c = colorbar; 
clim([-20,-2])
c.Label.String = 'temp (C)';
c.Label.FontSize = 16;
title('approximation to flowline temp')