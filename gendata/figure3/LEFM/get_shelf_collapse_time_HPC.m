function [collapse_time, collapse_time_square, tags] = get_shelf_collapse_time(shelf_name, step)
%return the collapse time (as part of the larger array -- collapse time
%square is the sub-array) for the shelf specified by shelf_name. Step is
%the step in the grid resolution (integer, set to 1 for whole grid)

% return the collapse time for the shelf inputted
poolobj = gcp('nocreate');
if ~isempty(poolobj);  delete(poolobj); end
num_cpu=56;
poolobj = parpool('local',num_cpu);

shelf_name = string(shelf_name);
addpath('../../../functions')

%% Load the ice shelf data
f  = load('../../../data/ice_sheet_data.mat');
shelf = shelf_name;
fname = strcat('../../../data/ice-shelves/all-shelves/' ,shelf, '.mat');
g = load(fname);

%% create restricted co-ordinates
[rId, cId] = find(g.IN) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns

%xs = f.xx(xminidx:step:xmaxidx);
%ys = f.yy(yminidx:step:ymaxidx); %restricted co-ordinates
hs = f.H(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
ms = f.m(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
dhdts = f.dMdtadj(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
strains = f.eflow(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
Bs = f.B(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
sz = size(hs);

%% adhoc adjustments
%dhdts = -abs(dhdts);  
dhdts = g.mean_dhdt;
ms(ms < 1e-1) = 1e-1; %set a minimum melt value

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
%tags(isnan(ms) | isnan(dhdts) | isnan(strains)) = 2;
tags(isnan(ms) | isnan(strains)) = 2;
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
tmax = 5*1e4; %max time
dt = 1; %timestep
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
    ix/sz(1)

end

% put this back into the array
collapse_time = nan(size(f.H));
collapse_time(xminidx:step:xmaxidx, yminidx:step:ymaxidx) = collapse_time_square;


%toc
savefolder = strcat('./step_', num2str(step));
if ~exist(savefolder, 'dir'); mkdir(savefolder); end
save(strcat(savefolder, '/collapse_time_', shelf,'_step', num2str(step), '.mat'), ...
    'collapse_time', 'collapse_time_square', 'tags', 'xminidx', 'xmaxidx', 'yminidx', 'ymaxidx', 'pp', 'ms', 'dhdts', 'strains', 'shelf_name', 'step');
end
