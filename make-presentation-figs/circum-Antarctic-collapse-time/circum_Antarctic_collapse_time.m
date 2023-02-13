function collapse_time_row = circum_Antarctic_collapse_time(nr,nrow)
% return the collapse time for the nth row of the partition, partitioned
% into nr pieces.
num_cpu=32;
poolobj = parpool('local',num_cpu);

addpath('../../functions')
f = load('../../data/ice_sheet_data.mat');
iidx = find(~isnan(f.H) & ~isnan(f.dhdtadj) & ~isnan(f.eflow) & ~isnan(f.m));

%nrow = 2; 
%nr = 20; %number of splits

np = ceil(length(iidx)/nr); %number of points per split

idx = zeros(np,nr);
idx(1:length(iidx)) = iidx;


% make array showing which run number
runsinfo = nan(size(f.H));
for ir = 1:nr
    idxrow = idx(:,ir);
    idxrow = idxrow(idxrow > 0);
    runsinfo(idxrow) = ir;
end

%%
h = f.H;
dhdt = f.dhdtadj;
epsxx = f.eflow;
mdot = f.m;

%% Constant quantities
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

%% loop over points in this 'row'
fidx = idx(:, nrow);
fidx = fidx(fidx ~=0);
if nrow == nr
    np = length(fidx);
end
collapse_time_row = nan(np,1);
H = h(idx(1:np, nrow));
DHDT = dhdt(idx(1:np, nrow));
EPSXX = epsxx(idx(1:np, nrow));
MDOT = mdot(idx(1:np, nrow));
parfor ip = 1:np

    %set up the parameters
    pp = struct;

    %variable quantities
    pp.H0    = H(ip);                   %initial ice thickness
    pp.dhdt  = -abs(DHDT(ip));          %rate of change of thickness (negative for thinning)
    pp.epsxx = EPSXX(ip);               %strain rate
    pp.mdot  = max([MDOT(ip), 1e-1]);   %melt rate

    %constant quantities
    pp.Tb    = Tb;     %basal temperature (kelvin)
    pp.Ts    = Ts;    %surface temp
    pp.B0    = B0;  %viscosity constant
    pp.rhoi  = rhoi;  %ice density
    pp.g     = g;   %gravitational acceleration
    pp.kappa = kappa;     %ice diffusivity
    pp.n     = n;      %glen flow coeff
    pp.frac_tough = frac_tough;
    pp.ghf   = ghf; %geothermal heat flux

    %dimensionless quantities
    pp.l      = pp.kappa/pp.H0/pp.mdot; %initial boundary layer lengthscale
    pp.F      = pp.frac_tough / (pp.H0)^(3/2) / pp.rhoi / pp.g;
    pp.lambda = 2*sign(pp.epsxx)*abs(pp.epsxx^(1/pp.n)) * pp.B0 / pp.rhoi / pp.g /  pp.H0;

    % timestepping quantities 
    tmax = 10000;
    dt = max([2,abs(1/pp.dhdt)]); %larger timesteps for smaller dhdt

    %get the relevant quantities
    if (~isnan(pp.dhdt) && ~isnan(pp.mdot) && ~isnan(pp.epsxx) && pp.epsxx > 0)
        collapse_time_row(ip) = get_collapse_time_advect(pp, dt, tmax);
    end
  % ip
end
%collapse_time_row
save(strcat('collapse_time_', num2str(nr),'_', num2str(nrow) , '.mat'), 'collapse_time_row', 'idx', 'nr', 'nrow');
end
