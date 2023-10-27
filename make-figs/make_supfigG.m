% Make supplementary figure "G", showing the basal melt rate in the Park et
% al. simulations for selected Antarctic ice shelves
%
% NB: this code takes several minutes to run (looping over time points in
% data is expensive). 
%
% 05/10/2023, ATB (aleey@bas.ac.uk), MIT Licence
%% Preliminaries

%park2023 = load('../gendata/figure4/park2023.mat'); % get the melt rates
addpath('../functions')


%% Get the data
shelves = ["Ross", "RonneFilchner", "Larsen", "Amery"];
nt = 300; %number of time points
melt_scenarios = nan(4, 3,10,nt); %4 shelves, 3 scenarios, 10 ensemble members, nt time points


for ish = 1:4 %for each shelf
    % load data for this shelf
    fname = strcat("Park2023_", shelves(ish), ".mat"); %shelf outlines for Park2023
    indata = load(fname);

    for isc = 1:3 %loop over scenarios
        for ie = 1:10 %loop over ensemble members
            for it = 1:nt %loop over time points

                grid_area = park2023(isc,ie).darea(:,:,it); %grid area (data in lat/lon)
                melt = park2023(isc,ie).melt(:,:,it);

                idx = (indata.IN) & (melt ~= 0) & (~isnan(melt)); %points to keep
                melt_keep = melt(idx);
                grid_keep = grid_area(idx);

                melt_total = sum(melt_keep .* grid_keep);
                area_total = sum(grid_keep);

                melt_scenarios(ish, isc,ie,it) = melt_total/area_total; %mean melt rate on shelf

            end
        end
    end
end

%% Make the plot
figure(1); clf;
colmap = [0,1,1; %cyan: SSP1-1.9
    1, 0,1; %magenta SSP2-4.5
    1,0,0]; %red SSP5-8.5

labels = ["Ross", "Ronne-Filchner", "Larsen", "Amery"];

colmap = linspecer(3);
colmap = colmap([3,1,2],:); %put green as ssp1, blue as ssp2 and red as ssp 3
ylims = [10, 2, 1.5, 8];
for ish = 1:4
    ax(ish) = subplot(2,2,ish); hold on; box on;
    for isc = 1:3
        melt_mean = squeeze(mean(melt_scenarios(ish, isc,:,:),3, 'omitnan')); %mean over 10 scenarios
        melt_std  = squeeze(std(melt_scenarios(ish, isc,:,:),1, 3, 'omitnan')); %std over 10 scenarios
        if isc < 3 %SSP1 and SSP2 are from 2015, SSP5 is from 1950
            tt = 2015:(2500-(486 - length(melt_mean)));
        else
            tt = 1850:2149;
        end
        tf = [tt, flip(tt)];
        mf = [(melt_mean - melt_std); flip(melt_mean+melt_std)];
        fill(tf, mf, colmap(isc,:), 'linestyle', 'none', 'FaceAlpha',0.3, 'HandleVisibility','off');
        plot(tt,  melt_mean, 'linewidth', 2, 'color', colmap(isc,:)); %SSP 5-8.5 includes 1850-2015
        plot([2100, 2100], [0, ylims(ish)], 'k--', 'LineWidth',1.3,'HandleVisibility','off' , 'Color', 0.5*[1,1,1]);


    end


    ax(ish).XLabel.String = 'time (years)';
    ax(ish).YLabel.String = 'mean melt rate (m/yr)';
    ax(ish).FontSize = 13;
    grid(ax(ish), 'on')
    if ish == 1; legend({'SSP1-1.9', 'SSP2-4.5', 'SSP5-8.5'}, 'location', 'northwest'); end
    xlim([2015, 2150])
    title(labels(ish), 'FontWeight','normal');
end
