% how much of the sheet do we have in the shelves?

jdir = dir('./all-shelves/*.mat'); %where shelf files are stored

% loop over shelves
h_shelf = nan(size(f.H)); %we'll store shelf data in here
for i = 1:length(jdir)
    shelf = jdir(i).name;
    shelf = strrep(shelf,'.mat',''); %strip the .mat at the end
    shelf_names(i) = shelf;
    
    fname = strcat('./all-shelves/' ,shelf, '.mat');
    
    g = load(fname);
    h_thisshelf = f.H; h_thisshelf = h_thisshelf(g.IN);
    h_shelf(g.IN) = h_thisshelf;

end

% plot side by side with all shelves
clf;ax(1) =  subplot(1,3,1);
pl = imagesc(h_shelf);
set(pl, 'AlphaData', ~isnan(h_shelf));
title('shelf reconstruction')
axis equal

ax(2) = subplot(1,3,2);
pl = imagesc(f.H);
set(pl, 'AlphaData', ~isnan(f.H));
title('data')
axis equal

ax(3) = subplot(1,3,3);
idx = ~isnan(f.H) & isnan(h_shelf); %shows where the missing shelves are
pl = imagesc(idx);
set(pl, 'AlphaData', idx);
title('missing pts')
axis equal
colormap(ax(3), 0*[1,1,1])
linkaxes([ax(1),ax(2), ax(3)])


% say what percentage coverage we have in the shelves
percovsh = sum(sum(~isnan(h_shelf)))/ sum(sum(~isnan(f.H))) * 100