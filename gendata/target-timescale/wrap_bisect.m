% Wrapper script for running loops over the bisection script and saving.
f  = load('../../data/ice_sheet_data.mat');
step = 10;
target_times = [127,227,327,377,427,477];

shelf_names = ["Abbot","Amery","Borchgrevink","Brunt","Cook","Cosgrove","Crosson","Dotson","Filchner","FimbulJelbart",...
    "George6","Getz","KingBaudoin","Larsen",  "Nansen","PineIsland","PineIslandFast","PopeSmithKohler","RiiserLarsen", "Ronne","Ross",...
    "Shackleton", "Thwaites_","TottenMoscow", "West","Wilkins"];

for it = 1:length(target_times)
    target_time = target_times(it);
    for is = 1:length(shelf_names)
        shelf = shelf_names(is);
        % run the bisection for this shelf and this value
        dM = bisect_shelf_collapse_time(shelf, step, target_time,f);


        % directory:
        folder = strcat("./step_", num2str(step), '/target_time_', num2str(target_time), '/');
        if ~exist(folder,'dir'); mkdir(folder); end
        fname = strcat(folder, shelf, ".mat");
        save(fname, 'dM');
    end
end


