% Loop over different values of dm, generating the output each time.

shelf_name = "Ronne"; 
spatial_step = 1; 
%dm_timesonethousand = 1000*[0.0:0.1:0.9, 1:10];
dm_timesonethousand = 1000*[0.7:0.1:0.9, 1:10];
dm_timesonethousand = 1000*[0.1:0.1:0.4, 0.7, 0.8, 2:10];

num_cpu = 36; 
timestep = 1;
data_path = '../../data/ice_sheet_data.mat'; 
save_flag = 1; 
for i = 1:length(dm_timesonethousand)
    get_ave_shelf_collapse_time(shelf_name, spatial_step, dm_timesonethousand(i), timestep, num_cpu,  data_path,save_flag);
end
