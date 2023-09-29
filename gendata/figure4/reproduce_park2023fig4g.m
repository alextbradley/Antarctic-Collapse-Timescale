% reproduce figure 4g of Park et al. 2023
dx = 20e3;
dy = 20e3;

nt = 300; %number of time points
melt_scenarios = nan(3,10,nt);
for i = 1:3
    for j = 1:10
        for it = 1:nt
            melt = park2023(i,j).melt(:,:,it);
            melt_scenarios(i,j,it) = sum(sum(melt))*dx*dy;
        end
    end
end

%% make the plot
colmap = [0,1,1; %cyan: SSP1-1.9
          1, 0,1; %magenta SSP2-4.5
          1,0,0]; %red SSP5-8.5
figure(1); clf; hold on; box on
for is = 1:3
    melt_mean = squeeze(mean(melt_scenarios(is,:,:),2))/1e9; %mean over 10 scenarios
    if is < 3
        plot(2015:(2500-(486 - length(melt_mean))), melt_mean, 'linewidth', 2, 'color', colmap(is,:));
    else
        plot(1850:2149,  melt_mean, 'linewidth', 2, 'color', colmap(is,:)); %SSP 5-8.5 includes 1850-2015
    end

    
end
xlabel('time (yrs)');
ylabel('gt/yr2')
legend({'SSP1-1.9', 'SSP2-4.5', 'SSP5-8.5'})
xlim([1850, 2150])

