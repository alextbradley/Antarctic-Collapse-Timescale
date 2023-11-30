%Download the data from Park et al. 2023 (DOI: 10.1038/s41467-023-36051-9)

ssps = ["119", "245","585"]; %ssp scenarios
ensmembs = ["01","02","03","04","05","06","07","08","09","10"];

count  = 1;
park2023 = struct;
for i = 1:length(ssps)
    for j = 1:length(ensmembs)
        url = strcat("http://climatedata.ibs.re.kr:9090/dods/public-data/loveclip/SSP",ssps(i), "/loveclip-SSP",ssps(i),"_",ensmembs(j),"_ice_sh");
        park2023(i,j).lat = ncread(url, 'lat');
        park2023(i,j).lon = ncread(url, 'lon');
        park2023(i,j).melt = ncread(url, 'oceanmeltav');
        park2023(i,j).darea = ncread(url, 'darea');

        fprintf('Completed %g of %g \n', count, length(ssps)*length(ensmembs));

        count = count + 1;
    end
end

save('park2023.mat', 'park2023', '-v7.3');
