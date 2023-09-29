% Generate ice shelf boundaries based on Park et al., 2023 data

fname = 'Park2023_Larsen.mat';

melt = park2023(3,1).melt;
melt = squeeze(melt(:,:,165)); %the 2015 melt
melt(melt== 0) = nan; 
figure(1); clf; pl = imagesc(melt); set(pl, 'AlphaData', ~isnan(melt))

[x,y] = ginput();
hold on; plot(x,y,'ro-', 'markerfacecolor', 'r')

s = size(melt);
[xx,yy] = meshgrid(1:s(1),1:s(2));
IN = inpolygon(xx, yy, x,y);

melt_red = melt; melt_red(~IN) = nan;
figure(2);clf; 
pl = imagesc(melt_red); set(pl, 'AlphaData', ~isnan(melt_red))

save(fname, 'IN');
