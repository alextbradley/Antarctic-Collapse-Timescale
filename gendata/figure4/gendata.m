% Script to generate data for collapse time as a function of dM

data_path = "../../data/ice_sheet_data.mat";
spatial_step = 10;
timestep = 5; 
num_cpu = 2;
save_flag = 1;
dMs = 0.1:0.1:0.2;

shelf_names = ["Abbot","Amery","Borchgrevink","Brunt","Cook","Cosgrove","Crosson","Dotson","Filchner","FimbulJelbart",...
    "George6","Getz","KingBaudoin","Larsen",  "Nansen","PineIsland","PineIslandFast","PopeSmithKohler","RiiserLarsen", "Ronne","Ross",...
    "Shackleton", "Thwaites","TottenMoscow", "West","Wilkins"];

shelf_names = ["Abbot"];

for iM = 1:length(dMs)
    dM = dMs(iM);
    for is = 1:length(shelf_names)
        shelf_name = shelf_names(is);
  
        dm_timesonethousand = dM*1000;

        average_collapse_time = get_ave_shelf_collapse_time(shelf_name, spatial_step, dm_timesonethousand, timestep, num_cpu,  data_path, save_flag);

    end
end

