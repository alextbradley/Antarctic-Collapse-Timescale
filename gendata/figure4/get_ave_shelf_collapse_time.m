function average_collapse_time = get_ave_shelf_collapse_time(shelf_name, spatial_step, dm_timesonethousand, timestep, num_cpu,  data_path,save_flag)
%return the average collapse time for the input variables:
% shelf_name: name of the ice shelf,
% spatial_step: the number of multiples of the spatial grid to use (e.g.
%            spatial_step = 2 corresponds to a resolution of 2*1km = 2km)
% dm_timesonesthousand: melt increase times 1000 (multiplication so that we
%            can have integer input)
% timestep: timestep in the solver
% num_cpu: number of cpus
% data_path: path to input data set (ice_shelf_data.mat)

% return the collapse time for the shelf inputted
poolobj = gcp('nocreate');
if ~isempty(poolobj);  delete(poolobj); end
poolobj = parpool('local',num_cpu);

shelf_name = string(shelf_name);
addpath('../../functions')

%whos data_path
f  = load(data_path);

%% Load the ice shelf data
%f  = load('../../data/ice_sheet_data.mat');
shelf = shelf_name;
fname = strcat('../../data/ice-shelves/all-shelves/' ,shelf, '.mat');
g = load(fname);

%% create restricted co-ordinates
[rId, cId] = find(g.IN) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns

%xs = f.xx(xminidx:spatial_step:xmaxidx);
%ys = f.yy(yminidx:spatial_step:ymaxidx); %restricted co-ordinates
hs = f.H(xminidx:spatial_step:xmaxidx, yminidx:spatial_step:ymaxidx);
ms = f.m(xminidx:spatial_step:xmaxidx, yminidx:spatial_step:ymaxidx);
dhdts = f.dMdtadj(xminidx:spatial_step:xmaxidx, yminidx:spatial_step:ymaxidx);
strains = f.eflow(xminidx:spatial_step:xmaxidx, yminidx:spatial_step:ymaxidx);
Bs = f.B(xminidx:spatial_step:xmaxidx, yminidx:spatial_step:ymaxidx);
sz = size(hs);

%% adhoc adjustments
%dhdts = -abs(dhdts);  
dhdts(:) = g.mean_dhdt;
ms(ms < 1e-1) = 1e-1; %set a minimum melt value

%% scale variables
d_m = dm_timesonethousand/1000; 
ms = ms + d_m; 
s_epsxx =  6.8109e-04; 
s_dhdt = -0.4724;      %numbers based on regression of shelves (see figure 2)
d_epsxx = d_m*s_epsxx;
d_dhdt = d_m*s_dhdt;

% power law relationship
d_epsxx = 0.0033*ms.^(0.4) - 0.0033*(ms - d_m).^(0.4);
d_dhdt  = -(0.66*ms.^(0.7) - 0.66*(ms - d_m).^(0.7));

dhdts = dhdts + d_dhdt; mean(mean(dhdts(~isnan(dhdts))));
strains = strains + d_epsxx; mean(mean(strains(~isnan(strains))));
fprintf('Mean melt rate is %.3f m/yr \n', mean(mean(ms(~isnan(ms)))))
fprintf('Current mean dhdt is %.3f m/yr \n', mean(mean(dhdts(~isnan(dhdts)))))
fprintf('Current mean strain rate is %.5f 1/yr \n', mean(mean(strains(~isnan(strains)))))

%% determine tags
% tags:
%    1: outside domain (no thickness data)
%    2: other missing data (no melt, strain rate, dhdt)
%    3: negative strain rate
%    4: positive thinning rate (currently stable)
%    5: positive thinning rate (currently unstable)
%    6: other (good data)

tags = 6*ones(size(ms));

tags(dhdts > 0) = 4;
tags(strains < 0) = 3;
tags(isnan(ms) | isnan(dhdts) | isnan(strains)) = 2;
tags(isnan(hs)) = 1;

%% Constant quantities
%constant parameters
Tb    = -2 + 273.15;     %basal temperature (kelvin)
Ts    = -22 + 273.15;    %surface temp
B0    = 1.928;  %viscosity constant
rhoi  = 918.0;  %ice density
g     = 9.81;   %gravitational acceleration
kappa = 36;     %ice diffusivity
n     = 3;      %glen flow coeff
frac_tough = 150*1e3;
ghf = 48; %geothermal heat flux

%store the params in a struct
pp = struct;
pp.Tb = Tb;
pp.Ts = Ts;
pp.B0 = B0;
pp.rhoi = rhoi;
pp.g = g;
pp.kappa = kappa;
pp.n = n;
pp.frac_tough = frac_tough;
pp.ghf = ghf;


collapse_time_square = nan(size(ms));
%tic
tmax = 1e4; %max time
dt = timestep; %rename timestep
parfor ix =  1:sz(1)
%for ix = 1:length(xs)

    %variable parameters
    row_H = hs(ix,:);
    row_dhdt = dhdts(ix,:);
    row_epsxx = strains(ix,:);
    row_tags = tags(ix,:);
    row_mdot = ms(ix,:);
    
    %timestepping parameters
    row_dt  = max(dt,abs(1./row_dhdt)); %timestep (larger timesteps for smaller dhdt?)
    row_tmax = tmax*ones(size(row_mdot));   %max time                  %maximum time
    
    %whos row_*

    %get the collapse time on this row
    collapse_time_row = get_collapse_time_row(row_H, row_dhdt, row_epsxx, row_mdot, row_tags,row_dt, row_tmax, ...
                                                Tb, Ts, B0, rhoi, g, kappa, n, frac_tough, ghf);

    collapse_time_square(ix,:) = collapse_time_row;
    %ix/sz(1)

end

% put this back into the array
collapse_time = nan(size(f.H));
collapse_time(xminidx:spatial_step:xmaxidx, yminidx:spatial_step:ymaxidx) = collapse_time_square;

idx = ~isnan(collapse_time_square) & (collapse_time_square < 0.95*tmax); %remove very large output points
collapse_time_count = collapse_time_square(idx); %the pts we count
kde = fitdist(collapse_time_count,'Kernel');
average_collapse_time = mean(kde);

if save_flag
        dM = d_m;
        % directory:
        folder = strcat("./step_", num2str(spatial_step), '/', shelf, '/');
        if ~exist(folder,'dir'); mkdir(folder); end
        fname = strcat(folder, strcat('dM_',num2str(dM)), ".mat");
        save(fname, 'collapse_time', 'dM', 'timestep', 'average_collapse_time');
end