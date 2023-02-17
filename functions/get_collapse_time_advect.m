function collapse_time = get_collapse_time_advect(pp, dt, tmax)
% Compute the collapse time of a column for the advected profile

tnow = 0;
nsteps = floor(tmax/dt); %maximum number of timesteps
i = 0;

pp.l      = pp.kappa/pp.H0/pp.mdot; %boundary layer lengthscale
pp.lambda = 2*sign(pp.epsxx)*abs(pp.epsxx^(1/pp.n)) * pp.B0 / pp.rhoi / pp.g /  pp.H0;
pp.F      = pp.frac_tough / (pp.H0)^(3/2) / pp.rhoi / pp.g;


TgfF = get_grounding_line_temp(pp.ghf, pp.Ts, pp.H0);
anonT = @(z) TgfF(z) + (pp.Tb- TgfF(z)).*exp(-z/pp.l); %with advected grounding line contribution
%anonT = @(z) min(anonT(z), 273.15);
dimless_crev_depth = get_dimless_crev_depth(pp, anonT); %time zero crevasse depth

while (i < nsteps) &&  (dimless_crev_depth < 1)
    anonT = @(z) TgfF(z) + (pp.Tb- TgfF(z)).*exp(-z/pp.l); %with advected grounding line contribution
    %anonT = @(z) min(anonT(z), 273.15);  %don't let it go above freezing
    dimless_crev_depth = get_dimless_crev_depth(pp, anonT);


    %update thickness and time
    pp.H0 = pp.H0 + dt*pp.dhdt;
    %pp.H0
    %fprintf('%.2f \n',pp.H0)
    if dimless_crev_depth < 0.99
        tnow = tnow + dt;
    end

    %compute dimensionless quantities
    pp.l      = pp.kappa/pp.H0/pp.mdot; %boundary layer lengthscale
    pp.lambda = 2*sign(pp.epsxx)*abs(pp.epsxx^(1/pp.n)) * pp.B0 / pp.rhoi / pp.g /  pp.H0;
    %pp.lambda
    pp.F = pp.frac_tough / (pp.H0)^(3/2) / pp.rhoi / pp.g;

    i = i + 1;
end

collapse_time = tnow;
