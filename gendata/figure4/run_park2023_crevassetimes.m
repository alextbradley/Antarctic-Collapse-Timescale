% "Ross", [256, 2301,3611]
% "Amery", [-278, 768, 3571]
% "Larsen", [-37, 2776, 4312]
% "Filchner", [223, 896, 2213]
% "Ronne", [223, 896, 2213]

shelf_name = "Ross";
for dm_timesonethousand = [256, 2301,3611]
average_collapse_time = get_ave_shelf_collapse_time(shelf_name, 10, dm_timesonethousand, 5, 2,  "../../data/ice_sheet_data.mat",1);
end

shelf_name = "Amery";
for dm_timesonethousand = [-278, 768, 3571]
average_collapse_time = get_ave_shelf_collapse_time(shelf_name, 10, dm_timesonethousand, 5, 2,  "../../data/ice_sheet_data.mat",1);
end

shelf_name = "Larsen";
for dm_timesonethousand = [-37, 2776, 4312]
average_collapse_time = get_ave_shelf_collapse_time(shelf_name, 10, dm_timesonethousand, 5, 2,  "../../data/ice_sheet_data.mat",1);
end

shelf_name = "Filchner";
for dm_timesonethousand = [223, 896, 2213]
average_collapse_time = get_ave_shelf_collapse_time(shelf_name, 10, dm_timesonethousand, 5, 2,  "../../data/ice_sheet_data.mat",1);
end


shelf_name = "Ronne";
for dm_timesonethousand = [223, 896, 2213]
average_collapse_time = get_ave_shelf_collapse_time(shelf_name, 10, dm_timesonethousand, 5, 2,  "../../data/ice_sheet_data.mat",1);
end

