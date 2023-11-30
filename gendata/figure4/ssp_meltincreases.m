% Compute the melt rate increase for the selected shelves shown in figure
% 4c. Shelves are:
% Ronne-Filcher (considered together)
% Ross
% Amery
% Larsen

% load in all data
%f = load('../../data/ice_sheet_data.mat');

%% Ross ice shelf
shelf = load('../../data/ice-shelves/all-shelves/Ross.mat');
melt = f.m;
melt = melt(shelf.IN);
melt(melt < 0) = 0;

fprintf('mean melt rate of Ross ice shelf is %.3f m/yr \n', mean(melt(~isnan(melt))))
fprintf('mean melt rate of Ross ice shelf under SSP1-1.9 (93 percent increase) is %.3f m/yr, i.e. dM = %.3f \n', mean(melt(~isnan(melt)))*1.93, mean(melt(~isnan(melt)))*0.93 )
fprintf('mean melt rate of Ross ice shelf under SSP2-4.5 (835 percent increase) is %.3f m/yr, i.e. dM = %.3f \n', mean(melt(~isnan(melt)))*9.35, mean(melt(~isnan(melt)))*8.35 )
fprintf('mean melt rate of Ross ice shelf under SSP5-8.5 (1311 percent increase) is %.3f m/yr, i.e. dM = %.3f \n', mean(melt(~isnan(melt)))*14.1, mean(melt(~isnan(melt)))*13.1 )

%% Amery ice shelf
shelf = load('../../data/ice-shelves/all-shelves/Amery.mat');
melt = f.m;
melt = melt(shelf.IN);
melt(melt < 0) = 0;

fprintf('mean melt rate of Amery ice shelf is %.3f m/yr \n', mean(melt(~isnan(melt))))
fprintf('mean melt rate of Amery ice shelf under SSP1-1.9  (-25 percent increase) is %.3f m/yr, i.e. dM = %.3f \n', mean(melt(~isnan(melt)))*0.75, mean(melt(~isnan(melt)))*(-0.25) )
fprintf('mean melt rate of Amery ice shelf under SSP2-4.5  (69 percent increase) is %.3f m/yr, i.e. dM = %.3f \n', mean(melt(~isnan(melt)))*1.69, mean(melt(~isnan(melt)))*0.69)
fprintf('mean melt rate of Amery ice shelf under SSP5-8.5  (321 percent increase) is %.3f m/yr, i.e. dM = %.3f \n', mean(melt(~isnan(melt)))*4.21, mean(melt(~isnan(melt)))*3.21 )

%% Larsen ice shelf
shelf = load('../../data/ice-shelves/all-shelves/Larsen.mat');
melt = f.m;
melt = melt(shelf.IN);
melt(melt < 0) = 0;

fprintf('mean melt rate of Larsen ice shelf is %.3f m/yr \n', mean(melt(~isnan(melt))))
fprintf('mean melt rate of Larsen ice shelf under SSP1-1.9  (-2 percent increase) is %.3f m/yr, i.e. dM = %.3f \n', mean(melt(~isnan(melt)))*0.98, mean(melt(~isnan(melt)))*(-0.02) )
fprintf('mean melt rate of Larsen ice shelf under SSP2-4.5  (150 percent increase) is %.3f m/yr, i.e. dM = %.3f \n', mean(melt(~isnan(melt)))*2.50, mean(melt(~isnan(melt)))*1.50 )
fprintf('mean melt rate of Larsen ice shelf under SSP5-8.5  (233 percent increase) is %.3f m/yr, i.e. dM = %.3f \n', mean(melt(~isnan(melt)))*3.33, mean(melt(~isnan(melt)))*2.33 )

%% FR ice shelf
Fshelf = load('../../data/ice-shelves/all-shelves/Filchner.mat');
Rshelf = load('../../data/ice-shelves/all-shelves/Ronne.mat');
melt = f.m;
melt = melt(Rshelf.IN | Fshelf.IN);
melt(melt < 0) = 0;

fprintf('mean melt rate of Filchner-Ronne ice shelf is %.3f m/yr \n', mean(melt(~isnan(melt))))
fprintf('mean melt rate of Filchner-Ronne ice shelf under SSP1-1.9  (45 percent increase) is %.3f m/yr, i.e. dM = %.3f \n', mean(melt(~isnan(melt)))*1.45, mean(melt(~isnan(melt)))*0.45 )
fprintf('mean melt rate of Filchner-Ronne ice shelf under SSP2-4.5  (181 percent increase) is %.3f m/yr, i.e. dM = %.3f \n', mean(melt(~isnan(melt)))*2.81, mean(melt(~isnan(melt)))*1.81 )
fprintf('mean melt rate of Filchner-Ronne ice shelf under SSP5-8.5  (447 percent increase) is %.3f m/yr, i.e. dM = %.3f \n', mean(melt(~isnan(melt)))*5.47, mean(melt(~isnan(melt)))*4.47 )
