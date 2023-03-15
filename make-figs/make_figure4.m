% Make figure 4 of the ms, showing something to do with melt rate changes
% required.

% ATB (aleey@bas.ac.uk), 21/02/23. MIT licence.

%%
addpath('../functions');
f = load('../data/ice_sheet_data.mat', 'H', 'm'); %get the ice sheet thickness data (use this for shelf mask)


target_time = [77,177,277]; %77 corresponds to 2100, 177 to 2200 etc.
lt = length(target_time);
step = 10;        %grid resolution in km

shelf_info_folder = "../data/ice-shelves/all-shelves/"; %where to find co-ordinate files
dM_folder = strcat('../gendata/figure4/step_', num2str(step), '/target_time_', num2str(target_time(1)), '/'); %where to find dM results (need for size of dir here)
jdir = dir(strcat(dM_folder, '*.mat'));
ls = length(jdir);

%setup storage
dM = zeros(lt,ls);
shelf_names = strings(1,ls);
total_flux = zeros(1,ls);


dM_full = cell(lt,1);  %store arrays of full antarctic with shelves set to shelf wide dM
curr_mean_melt = nan(size(f.H));

for it = 1:lt
    dM_folder = strcat('../gendata/figure4/step_', num2str(step), '/target_time_', num2str(target_time(it)), '/'); %where to find dM results
    dM_full_this_targettime = nan(size(f.H)); 

    for is = 1:ls% loop over shelves

        shelf = jdir(is).name;
        data = load(strcat(dM_folder, shelf));
        dM(it,is) = data.dM;
        shelf_names(is) = strrep(shelf, '.mat', '');

        % get the shelf info
        shelf_info = load(strcat(shelf_info_folder, shelf));
        idx = ~isnan(f.H) & shelf_info.IN; %points in a shelf and within the shelf region

        % fill in the array
        dM_full_this_targettime(idx) = dM(it,is);

        % get the current mean melt rate in the shelf
        idxm = ~isnan(f.m) & shelf_info.IN; %points within the shelf region which have a melt rate
        curr_mean_melt(idx) = mean(mean(f.m(idxm))); %set points in the shelf to the mean melt rate

        total_flux(is) = sum(sum(f.m(idxm)))*1e3*1e3/1e9; %total flux from this region (in gt/yr)


    end %end loop over shelves
    dM_full{it} = dM_full_this_targettime;
end %end loop over target times

%% Make a plot
% figure(1); clf;
% dM_full_copy = dM_full; dM_full_copy(dM_full_copy == 0) = 1e-6; %set minimum to 1e-6 so that we can plot on log scale if we want
% %p = imagesc(log10(dM_full_copy));clim([-1,3]); %colormap(flipud(cmocean('amp')));
% p = imagesc(log10(dM_full_copy./curr_mean_melt)); colormap('parula');clim([-1,2])
% 
% 
% set(p, 'AlphaData', ~isnan(dM_full));
% 
% colorbar
% ax = gca; ax.FontSize = 15;


%%
figure(1); clf; plot(1:3,dM, 'linewidth', 2)
set(gca,'YScale','log')
ylim([0,2])
ylim(10.^[0,2])
xticks(1:3);
ax = gca; ax.FontSize = 15;
xticklabels = {'2100', '2200', '2300'};
ax.XTickLabel = xticklabels;
grid on

%% Make a bar chart
bar(dM'); 
ax = gca;
ax.XTick = 1:26;
ax.XTickLabel = shelf_names;
ax.XTickLabelRotation = 45;
ax.FontSize = 14;
set(ax,'YScale', 'log')
shg
ylabel('melt increase (m/yr)');
legend('2100', '2200', '2300')
