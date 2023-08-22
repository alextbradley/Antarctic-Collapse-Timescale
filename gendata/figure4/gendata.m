% Script to generate data for collapse time as a function of dM

f  = load('../../data/ice_sheet_data.mat');
step = 10;
dMs = 0.1:0.1:0.9;

shelf_names = ["Abbot","Amery","Borchgrevink","Brunt","Cook","Cosgrove","Crosson","Dotson","Filchner","FimbulJelbart",...
    "George6","Getz","KingBaudoin","Larsen",  "Nansen","PineIsland","PineIslandFast","PopeSmithKohler","RiiserLarsen", "Ronne","Ross",...
    "Shackleton", "Thwaites_","TottenMoscow", "West","Wilkins"];

for iM = 1:length(dMs)
    dM = dMs(iM);
    for is = 1:length(shelf_names)
        shelf = shelf_names(is);
     
        % run the bisection for this shelf and this value d_m = guess;
%         s_epsxx = 1.7127e-04; 
%         s_dhdt = -0.0836;     %numbers based on all points around Antarctic 
        s_epsxx =  6.8109e-04; 
        s_dhdt = -0.4724;      %numbers based on regression of shelves (see supplementary figure D)

        d_m = dM;
        d_epsxx = d_m*s_epsxx;
        d_dhdt = d_m*s_dhdt;

        collapse_time = get_ave_shelf_collapse_time(shelf, step, d_m, d_epsxx,d_dhdt, f);


        % directory:
        folder = strcat("./step_", num2str(step), '/', shelf, '/');
        if ~exist(folder,'dir'); mkdir(folder); end
        fname = strcat(folder, strcat('dM_',num2str(dM)), ".mat");
        save(fname, 'collapse_time', 'dM');
    end
end

