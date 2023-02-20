% Plot the LEFM predicted crevasse depth for (a) a linear profile and (b)
% boundary layer profile

% Preliminaries
addpath('../../functions');

% load data
f  = load('../../data/ice_sheet_data.mat');
h = f.H';
m = f.m';
eflow = f.eflow';
B = f.B';


% restrict co-ordinates
xminidx = 1; xmaxidx = length(f.xx);
yminidx = 1; ymaxidx = length(f.yy);
step = 4; %step in grid cell. Resolution becomes 500m*step
xs = f.xx(xminidx:step:xmaxidx);
ys = f.yy(yminidx:step:ymaxidx); %restricted co-ordinates
hs = h(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
ms = m(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
strains = eflow(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
Bs = B(xminidx:step:xmaxidx, yminidx:step:ymaxidx);

%% Boundary layer profile
dimless_crev_depth = nan(size(hs));
tic
parfor ix = 1:length(xs)
    row_H = hs(ix,:);
    row_epsxx = strains(ix,:);
    row_mdot = ms(ix,:);

    %constant parameters
    Tb    = -2 + 273.15;     %basal temperature (kelvin)
    Ts    = -18 + 273.15;    %surface temp
    B0    = 1.928;  %viscosity constant
    rhoi  = 918.0;  %ice density
    g     = 9.81;   %gravitational acceleration
    kappa = 36;     %ice diffusivity
    n     = 3;      %glen flow coeff
    frac_tough = 150*1e3;
    ghf = 48; %geothermal heat flux

   
    dim_crev_depth_row = get_dim_crev_depth_row(row_H, row_epsxx, row_mdot, ...
                                                Tb, Ts, B0, rhoi, g, kappa, n, frac_tough, ghf);


    dimless_crev_depth(ix,:) = dim_crev_depth_row;

    
    ix
    %toc
    
end

%% Repeat for the linear profile
dimless_crev_depth_linear = nan(size(hs));
tic
parfor ix = 1:length(xs)
    row_H = hs(ix,:);
    row_epsxx = strains(ix,:);
    row_mdot = ms(ix,:);

    %constant parameters
    Tb    = -2 + 273.15;     %basal temperature (kelvin)
    Ts    = -16 + 273.15;    %surface temp
    B0    = 1.928;  %viscosity constant
    rhoi  = 918.0;  %ice density
    g     = 9.81;   %gravitational acceleration
    kappa = 36;     %ice diffusivity
    n     = 3;      %glen flow coeff
    frac_tough = 150*1e3;
    ghf = 48; %geothermal heat flux

   
    dim_crev_depth_row = get_dim_crev_depth_row_linear(row_H, row_epsxx, row_mdot, ...
                                                Tb, Ts, B0, rhoi, g, kappa, n, frac_tough, ghf);


    dimless_crev_depth_linear(ix,:) = dim_crev_depth_row;

    
    ix
    %toc
    
end
