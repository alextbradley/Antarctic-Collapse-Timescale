%% Plot the collpase time as a function of dM for different effects

%%
addpath('../../../functions');
f = load('../../../data/ice_sheet_data.mat', 'H', 'm'); %get the ice sheet thickness data (use this for shelf mask)
step = 10;        %grid resolution in km

shelf_names = ["Abbot","Amery","Borchgrevink","Brunt","Cook","Cosgrove","Crosson","Dotson","Filchner","FimbulJelbart",...
    "George6","Getz","KingBaudoin","Larsen",  "Nansen","PineIsland","PineIslandFast","PopeSmithKohler","RiiserLarsen", "Ronne","Ross",...
    "Shackleton", "Thwaites_","TottenMoscow", "West","Wilkins"];

shelf_info_folder = "../../../data/ice-shelves/all-shelves/"; %where to find co-ordinate files

effects = ["isolate_strainrate", "isolate_meltrate", "isolate_dhdt"];
%% Loop over the shelves. For each, make a an array with the dM and the mean collapse time
ct_shelves = cell(1,length(shelf_names));
dM_shelves = cell(1,length(shelf_names));
figure(1); clf; hold on; box on



for is = 1:length(shelf_names)
    clf; hold on; title(shelf_names(is))

    for ieff = 1:3

        data_folder = strcat(effects(ieff), '/step_', num2str(step));

        fpath = strcat(data_folder, '/', shelf_names(is));
        jdir = dir(strcat(fpath, '/*.mat')); %each of these is a different value of dM

        dMs = nan(1,length(jdir));
        cts = nan(1,length(jdir));
        for im  = 1:length(jdir)
            ff = load(strcat(fpath, '/', jdir(im).name));

            dMs(im) = ff.dM;
            cts(im) = ff.collapse_time;
        end
        %sort into increaseing dm order
        [dMs,I] = sort(dMs);
        cts = cts(I);

        ct_shelves{is} = cts;
        dM_shelves{is} = dMs;

        plot(dMs, cts)
    end
    legend('strain rate only', 'melt rate only', 'dhdt only')
    drawnow
    pause


end
