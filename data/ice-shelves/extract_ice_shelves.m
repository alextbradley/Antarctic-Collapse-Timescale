% Generate indices and lengthscales for some ice shelves

fname = 'Filchner'; %name of the ice shelf
%f = load('../ice_sheet_data.mat');
%%
%f  = load('../ice_sheet_data.mat')
velocs = sqrt((f.u).^2 + (f.v).^2);
figure(1); clf; h = imagesc(f.xx, f.yy, velocs);
set(h, 'AlphaData', ~isnan(f.u)); caxis([0, 4000])
c = colorbar;
c.Label.String = 'velocity (m/yr)';
axis equal
title('zoom and click around ice shelf')
zoom on;

waitfor(gcf, 'CurrentCharacter', char(13)) %what until we hit return
zoom reset
zoom off; % to escape the zoom mode

[x,y] = ginput();
hold on; plot(x,y,'r');

title('click on ice front and grounding line')
[lx,ly] = ginput(2);
plot(lx,ly, 'ko-', 'markerfacecolor', 'k')
ax = gca;
l = sqrt(diff(lx).^2 + diff(ly).^2); %lengthscale of ice shelf
%%
[xx,yy] = meshgrid(f.xx, f.yy);
IN = inpolygon(xx, yy, x,y);
%%
% make a plot of the new region
% uu = f.u;
% uu(~IN) = nan;
% figure(2); h = imagesc(f.xx, f.yy,uu); set(h, 'AlphaData', ~isnan(uu));
% ax2 = gca; ax2.XLim = ax1
%% save 
fname = strcat(fname, '.mat');
save(fname, 'IN', 'l')

pause(3);
close(gcf)
