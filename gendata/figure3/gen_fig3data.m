% Generate data file for figure 3 of the ms, containing:
% crevasse_time     : LEFM crevasse time
% crevasse_time_Nye : Nye crevasse time
% tags              : specifies type of data point.

% assemble LEFM results
jdir = dir('./LEFM/step_1/*.mat');
collapse_time = nan(4401,5301);
for i = 1:length(jdir)
    data = load(strcat('./LEFM/step_1/', jdir(i).name)); %load data for this shelf
    
    % get the in points 
    shelf_info = load(strcat('../../data/ice-shelves/all-shelves/', data.shelf_name, '.mat'));
    
    %fill in the matrix 
    idx = ~isnan(data.collapse_time) & shelf_info.IN;


    collapse_time(idx) = data.collapse_time(idx);
end

save('figure3-data.mat', 'collapse_time');
%%


% 
%% Repeat for Nye
jdir = dir('./Nye/step_1/*.mat');
collapse_time_Nye = nan(4401,5301);
for i = 1:length(jdir)
    data = load(strcat('./Nye/step_1/', jdir(i).name)); %load data for this shelf
    
    %fill in the matrix
    idx = ~isnan(data.collapse_time);
    collapse_time_Nye(idx) = data.collapse_time(idx);
end

save('figure3-data.mat', 'collapse_time_Nye', '-append');

%%
% figure(1); clf; p = imagesc(log10(collapse_time_Nye./collapse_time));
% set(p, 'AlphaData', ~isnan(collapse_time));
% colormap(cmocean('balance')); clim([-1,1])
% colorbar

%% Get the tags
f  = load('../../data/ice_sheet_data.mat');
tags = 6*ones(size(f.m));

tags(f.dMdtadj > 0) = 4;
tags(f.eflow < 0) = 3;
tags(isnan(f.m) | isnan(f.dMdt) | isnan(f.eflow)) = 2;
tags(isnan(f.H)) = 1;

save('figure3-data.mat', 'tags', '-append');


%% plot the tags
figure(1); clf; 
imagesc(flipud(tags))
c = colorbar;
clim([1,6])
colormap(lines(6))
c.Ticks = linspace(2-(7/12),5+(7/12),6);
c.TickLabels = {'1','2','3','4','5','6'};
set(gca, 'YDir', 'normal');

