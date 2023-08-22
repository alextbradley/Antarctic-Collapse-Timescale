%Generate a text file (array-config.txt) containing informaution about
%array job. text file has six columns containing (1) array job number, (2)
%shelf name, (3) dm value (multiplied by 1000, so its an integer), (4)
%spatial step, (5) timestep, and (6) number of cpus (typically 2 on my
%macbook and 24/32 on the hpc)

shelf_names = ["Abbot","Amery","Borchgrevink","Brunt","Cook","Cosgrove","Crosson","Dotson","Filchner","FimbulJelbart",...
    "George6","Getz","KingBaudoin","Larsen",  "Nansen","PineIsland","PineIslandFast","PopeSmithKohler","RiiserLarsen", "Ronne","Ross",...
    "Shackleton", "Thwaites","TottenMoscow", "West","Wilkins"];

dm_timesonethousand = 100:100:900; %dm *1e-3

spacestep = 1; 
timestep  = 1; 
numcpu    = 32; 

num_entr = length(shelf_names)*length(dm_timesonethousand);

D = cell(num_entr, 5);
count = 1;
for i = 1:length(shelf_names)
    for j = 1:length(dm_timesonethousand)
        D{count,1} = count;
        D{count,2} = shelf_names(i);
        D{count,3} = dm_timesonethousand(j);
        D{count,4} = double(spacestep);
        D{count,5} = timestep;
        D{count,6} = numcpu;
        count = count + 1;
    end
end
%%

fido = fopen('array-config.txt', 'w');
for k1 = 1:size(D,1)
    fprintf(fido, '%d %s %d %d %d %d\n', D{k1,:});
end