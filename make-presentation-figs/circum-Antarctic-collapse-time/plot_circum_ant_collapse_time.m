% Make a plot of the circum-Antarctic collapse time, assembled from mat
% files listed in this directory.

addpath('../../functions/');

collapse_time = nan(size(f.H)); 
jdir = dir('*.mat');
for i = 1:length(jdir)
fname = jdir(i).name;
cdata = load(fname);
iidx = cdata.idx(:,cdata.nrow);  %get the indices this file corresponds to
if length(iidx(iidx ~= 0)) < length(iidx)
    iidx = iidx(iidx ~= 0);
end
collapse_time(iidx) = cdata.collapse_time_row;
end
%

% get data tags
tags = nan(size(f.m));
tags = 6*ones(size(f.m));
tags(f.dhdtadj > 0) = 4;
tags(f.eflow < 0) = 3;
tags(isnan(f.m) | isnan(f.dhdtadj) | isnan(f.eflow)) = 2;
tags(isnan(f.H)) = 1;

%% make the plot:
cmap = cmocean('matter', 100);
% 1: collapse time
figure(2); clf; 
ax(1) = gca; pl = imagesc(log10(collapse_time)); 
c = colorbar;
c.FontName = 'GillSans';
clim([1,4])
set(pl, 'AlphaData', ~isnan(collapse_time));
axis equal
colormap(ax(1), cmap)


% 2: missing data
ax(2) = axes();
tags_missing = nan(size(tags));
tags_missing(tags == 2) = 1;
pl = imagesc(tags_missing);
set(pl, 'AlphaData', ~isnan(tags_missing)); 
colormap(ax(2), 0.6* [1,1,1])
ax(2).Visible = 'off';
axis equal


% 3: negative strain rate
ax(3) = axes();
tags_negative_strain = nan(size(tags));
tags_negative_strain(tags == 3) = 1;
pl = imagesc( tags_negative_strain);
set(pl, 'AlphaData', ~isnan(tags_negative_strain)); 
colormap(ax(3), cmap(end,:))
%colormap(ax(3),  0.3* [1,1,1])
axis equal

% % 4: positive dhdt 
% ax(4) = axes();
% tags_thicken = nan(size(tags));
% tags_thicken(tags == 4) = 1;
% tags_thicken((tags_thicken == 1) & (collapse_time < 10)) = nan; %set anything with short collapse time to nan 
% pl = imagesc( tags_thicken);
% set(pl, 'AlphaData', ~isnan(tags_thicken)); 
% colormap(ax(4), cmap(end,:))
% axis equal

for i = 1:4
    ax(i).Position = ax(1).Position;
    ax(i).Visible = 'off';

end


linkaxes([ax])



