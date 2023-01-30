function collapse_time = get_collapse_time(pp, dt, tmax)
% Compute the collapse time of a column. pp are input parameters.

tnow = 0;
nsteps = floor(tmax/dt); %maximum number of timesteps
i = 0;

pp.l      = pp.kappa/pp.H0/pp.mdot; %boundary layer lengthscale
pp.lambda = 2*sign(pp.epsxx)*abs(pp.epsxx^(1/pp.n)) * pp.B0 / pp.rhoi / pp.g /  pp.H0;
pp.F      = pp.frac_tough / pp.H0 / pp.rhoi / pp.g;

anonT = @(z) (pp.Ts + (pp.Tb - pp.Ts)*exp(-z/pp.l));
dimless_crev_depth = get_dimless_crev_depth(pp, anonT); %time zero crevasse depth

while (i < nsteps) &&  (dimless_crev_depth < 1)
    an = @(z) (pp.Ts + (pp.Tb - pp.Ts)*exp(-z/pp.l));
    dimless_crev_depth = get_dimless_crev_depth(pp, anonT);

    %update thickness and time 
    pp.H0 = pp.H0 + dt*pp.dhdt;
    %pp.H0
    if dimless_crev_depth < 1
    tnow = tnow + dt;
    end

    %compute dimensionless quantities
    pp.l      = pp.kappa/pp.H0/pp.mdot; %boundary layer lengthscale
    pp.lambda = 2*sign(pp.epsxx)*abs(pp.epsxx^(1/pp.n)) * pp.B0 / pp.rhoi / pp.g /  pp.H0;
    pp.F = pp.frac_tough / pp.H0 / pp.rhoi / pp.g;

    i = i + 1;
end

collapse_time = tnow;
