% Generate the data for figure 1 (see description therein for figure info).
%
% ATB (aleey@bas.ac.uk) 16/02/23, MIT licence
%
%% Preliminaries
clear
addpath('../../functions');
data = load('../../data/ice_sheet_data.mat'); %load in the data (large!)

savedat = 1; %flag to specify save (1) or not (0)

%% Constant parameters throughout
ghf = 48;
epsxx = 0.005; %strain rate
kappa = 36; %ice diffusivity
grav = 9.81;
Tb = -2; %basal temp
Ts = -22; %surface temp
B0 = 1.928; %viscosity constant
glen_n = 3;
frac_tough = 150e3;
rhoi = 918.0; 
dhdt = -1; %ice thickness rat eof change

%% Generate data for (a), (b)
ds = 1e3; %flowline spacing
fnames = ["Ross","PineIsland"];
nf = [24, 50]; %flowline number

for i = 1:2
    fname = strcat('../../data/ice-shelves/all-shelves/',fnames(i), '.mat');
    f = load(fname);

    %get flowline data (all flowlines)
    ss = get_flowline_data(fname, data, ds);
    s = ss(nf(i));

    %flowline co-ordinates
    flc = s.flowline;
    flc_all{i} = flc; %store for plotting later
    x = flc- flc(1,:);
    x = sqrt(x(:,1).^2 + x(:,2).^2);
    dzeta = 1e-3; %vertical dimensionless spacing
    zeta = dzeta:dzeta:(1-dzeta);
    S = 34.6 *ones(size(x)); %salinity along flowline
   

    %assemble profile
    TgfF = get_grounding_line_temp(ghf, Ts, s.h(1));
    [T,z,xx] = get_flowline_temp(x, zeta,s.h,s.speed,S, abs(s.melt), TgfF); %thsi gets the 'actual' Sergienko profile, but we use it for the co-ordinates and melt rate later
    T_exp = zeros(size(T));
    
    for ix = 1:length(x)

        mm = max(s.melt(ix), .1); 
        l  = kappa / mm / s.h(ix);
        anonT = @(z) ((Tb -  TgfF(z))*exp(-z/l)) ;
        TgfF = get_grounding_line_temp(ghf, Ts, s.h(ix));
        T_exp(:, ix) = TgfF(zeta) + (Tb - TgfF(zeta)).*exp(- zeta/l);
        ix/length(x) %for monitoring progress
    end
    Tflowline{i} = T_exp;
    Tflowline_act{i} = T;
    xflowline{i} = xx;
    zflowline{i} = z;
    mflowline{i} = s.melt;
    hflowline{i} = s.h;


end %end loop over shelves

%save
if savedat
save('figure1-data.mat', 'Tflowline','Tflowline_act',  'xflowline', 'zflowline', 'fnames', 'flc_all');
end
%% Generate inset data
figure1ab_inset = struct;
for i = 1:2
    fname = strcat('../../data/ice-shelves/all-shelves/',fnames(i), '.mat');
    f = load(fname);

    %restrict to this shelf
    in = f.IN;
    [rId, cId] = find(in) ;
    xminidx = min(rId); xmaxidx = max(rId);
    yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
    insetss(i).xs = data.xx(xminidx:xmaxidx);
    insetss(i).ys = data.yy(yminidx:ymaxidx); %restricted co-ordinates
    insetss(i).us = data.u(xminidx:xmaxidx, yminidx:ymaxidx); %velocity on this 'square'

end
if savedat
save('figure1-data.mat', "insetss", '-append');
end


%% Generate the data for (c), (d)
% (c) shows the crev depth as a function of time and (d) the temperature
% profile through time

% find where ice thickness passes thru a point
colthick = 600; %column thickness
for i = 1:2
    hh = cell2mat(hflowline(i));
    mm = cell2mat(mflowline(i));
    [~,idx] = min(abs(hh - colthick));
    m(i) = mm(idx); %melt rate we'll use
end
    
%% generate parameters for the collapse problem
pp = struct;
for i = 1:2
    pp(i).Tb    = Tb + 273.15;     %basal temperature (kelvin)
    pp(i).Ts    = Ts + 273.15;    %surface temp
    pp(i).B0    = B0;  %viscosity constant
    pp(i).rhoi  = rhoi;  %ice density
    pp(i).g     = grav;   %gravitational acceleration
    pp(i).epsxx = epsxx;  %strain rate
    pp(i).kappa = kappa;     %ice diffusivity
    pp(i).n     = glen_n;      %glen flow coeff
    pp(i).frac_tough = frac_tough; % ice fracture
    %pp(i).mdot  = m(i);
    pp(1).mdot = 0.1;
    pp(2).mdot = 48;
    pp(i).H0    = colthick;  %initial ice thickness
    pp(i).dhdt  = -1; %wlog (we'll normalize by dhdt)
    pp(i).dt    = 1; %timestep
    pp(i).ntmax = 1e6;  %max time
    pp(i).t     = 0; %seed the time
    pp(i).H     = pp(i).H0; %seed the thickness
   

    %compute dimensionless quantities
    pp(i).l      = pp(i).kappa/pp(i).H0/pp(i).mdot; %boundary layer lengthscale
    pp(i).lambda = 2*sign(pp(i).epsxx)*abs(pp(i).epsxx^(1/pp(i).n)) * pp(i).B0 / pp(i).rhoi / pp(i).g /  pp(i).H0;
    pp(i).F      = pp(i).frac_tough / (pp(i).H0)^(3/2) / pp(i).rhoi / pp(i).g;
    
end

%do the timestepping
ss = struct;
% pp(1).mdot = 0.1;
% pp(2).mdot = 90;
for ip = 1:2 %for each parameter set

    i = 1;
    dimless_crev_depth = 0; %seed the dimensionless crevasse depth
    TgfF = get_grounding_line_temp(ghf, pp(ip).Ts, pp(ip).H0);

    while (i < pp(ip).ntmax) && (dimless_crev_depth < 0.99)
        pp(ip).lambda;
        anonT = @(z) TgfF(z) + (pp(ip).Tb - TgfF(z)).*exp(-z/pp(ip).l); %with advected grounding line contribution
        [dimless_crev_depth, stress_intensity] = get_dimless_crev_depth(pp(ip), anonT);    
        dimless_crev_depth
        % store the solution
        ss(ip,i).t = pp(ip).t;
        ss(ip,i).H = pp(ip).H;
        ss(ip,i).dimless_crev_depth = dimless_crev_depth;
        ss(ip,i).stress_intensity = stress_intensity;
        zz = logspace(-4,0,1e2);
        ss(ip,i).zz = zz;
        ss(ip,i).temp_profile = anonT(zz);

        %update the thickness and dimensionless quantities
        pp(ip).H     = pp(ip).H + pp(ip).dt*pp(ip).dhdt;
        pp(ip).t      = pp(ip).t + pp(ip).dt;
        pp(ip).l      = pp(ip).kappa/pp(ip).H/pp(ip).mdot; %boundary layer lengthscale
        pp(ip).lambda = 2*sign(pp(ip).epsxx)*abs(pp(ip).epsxx^(1/pp(ip).n)) * pp(ip).B0 / pp(ip).rhoi / pp(ip).g /  pp(ip).H;
        pp(ip).F      = pp(ip).frac_tough / (pp(ip).H)^(3/2) / pp(ip).rhoi / pp(ip).g;
        %^^^ comment line out to expose only temperature profile changes
        %pause
        i = i+1;

    end

end

figc_data = ss;
figc_params = pp;
if savedat
save('figure1-data.mat', "figc_data", '-append');
save('figure1-data.mat', "figc_params", '-append');
end

%% Generate the data for (c)
HH = 40:10:1000;
mm = logspace(-2,2,7);

ct = zeros(length(HH), length(mm));
lambdac0 = zeros(length(HH), length(mm));
count = 1;
for im = 1:length(mm)
    for ih = 1:length(HH)
        % params
        tic
        pp = struct;
        pp.H0    = HH(ih);   %initial ice thickness
        pp.Tb    = Tb + 273.15;     %basal temperature (kelvin)
        pp.Ts    = Ts + 273.15;    %surface temp at grounding line
        pp.dhdt  = dhdt;      %rate of change of thickness
        pp.B0    = B0;  %viscosity constant
        pp.rhoi  = rhoi;  %ice density
        pp.g     = grav;   %gravitational acceleration
        pp.epsxx = epsxx;  %strain rate
        pp.kappa = kappa;     %diffusivity
        pp.mdot  = mm(im);     %melt rate
        pp.n     = glen_n;      %glen flow coeff
        pp.l     = pp.kappa/pp.H0/pp.mdot; %initial boundary layer lengthscale
        pp.frac_tough = frac_tough;
        pp.F = pp.frac_tough / (pp.H0)^(3/2) / pp.rhoi / pp.g;
        pp.ghf = ghf; %geothermal heat flux 

        tmax = 1000;
        dt = 1;

        ct(ih,im) = get_collapse_time_advect(pp, dt, tmax); 
        
        %compute the critical lambda
        lambdac0(ih,im) = get_critical_lambda(pp);

        fprintf('completed %.3f percent \n', 100* count/(length(HH)*length(mm)));
        toc
        count = count +1;
    end
end

figd_data = struct;
figd_data.ct = ct;
figd_data.h = HH;
figd_data.m = mm;
if savedat
save('figure1-data.mat', "figd_data", '-append');
end
