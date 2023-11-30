% Plot the changes in melt flux for each scenario and for each ice shelf.

dx = 20e3;
dy = 20e3; %needs updating once we have the grid area data

nt = 300; %number of time points
shelves = ["Ross", "RonneFilchner", "Larsen", "Amery"];
melt_scenarios = nan(4, 3,10,nt); %4 shelves, 4 scenarios, 10 ensemble members
for ish = 1:4
    fname = strcat("Park2023_", shelves(ish), ".mat");
    indata = load(fname);
    for isc = 1:3
        for ie = 1:10
            for it = 1:nt
                %grid_area = 400e3*ones(280,280); %grid area, to be updated
                grid_area = park2023(isc,ie).darea(:,:,it);
                melt = park2023(isc,ie).melt(:,:,it);
                idx = (indata.IN) & (melt ~= 0) & (~isnan(melt)); %points to keep
                melt_keep = melt(idx);
                grid_keep = grid_area(idx);
                
                melt_total = sum(melt_keep .* grid_keep);
                area_total = sum(grid_keep);

                melt_scenarios(ish, isc,ie,it) = melt_total/area_total;


%                 melt(~indata.IN) = nan; %remove points outside shelf area specified by us
%                 melt(melt == 0) = nan;  %remove any non-floating points
%                 melt_scenarios(ish, isc,ie,it) = mean(mean(melt, 'omitnan'), 'omitnan');
            end
        end
    end
end

%%
figure(1); clf;
colmap = [0,1,1; %cyan: SSP1-1.9
    1, 0,1; %magenta SSP2-4.5
    1,0,0]; %red SSP5-8.5

melt_2100 = zeros(4, 3);
melt_2015 = zeros(4, 1);
for ish = 1:4
    subplot(2,2,ish); hold on; box on;
    for isc = 1:3
        melt_mean = squeeze(mean(melt_scenarios(ish, isc,:,:),3, 'omitnan')); %mean over 10 scenarios
        if isc < 3
            plot(2015:(2500-(486 - length(melt_mean))), melt_mean, 'linewidth', 2, 'color', colmap(isc,:));
            melt_2100(ish, isc) = melt_mean(86);
        else
            plot(1850:2149,  melt_mean, 'linewidth', 2, 'color', colmap(isc,:)); %SSP 5-8.5 includes 1850-2015
            melt_2100(ish, isc) = melt_mean(251);
        end
        
       
    end


    %get the 2015 mean
    melt_2015(ish) = melt_mean(164);


    xlabel('time (yrs)');
    ylabel('m/yr')
    if ish == 1; legend({'SSP1-1.9', 'SSP2-4.5', 'SSP5-8.5'}); end
    xlim([1850, 2150])
    title(shelves(ish));
end

% work out the enhancement
melt_2015 = repmat(melt_2015, [1,3]);
melt_enhancement = melt_2100./melt_2015;
perc_melt_enhancement = (melt_enhancement -1)*100;