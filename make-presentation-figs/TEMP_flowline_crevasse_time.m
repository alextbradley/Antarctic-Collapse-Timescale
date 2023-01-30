
hh = ss(nf).h;
mm = ss(nf).melt;
eps_xx = ss(nf).eps_xx;
dhdt  = ss(nf).dhdt;

pp = struct;
pp.kappa = 36; 
pp.frac_tough = 150*1e3;
pp.rhoi = 918;
pp.g    = 9.81;
pp.n    = 3; 
pp.B0   = 1.928;
dt   = 1; 
tmax = 1000;

for ip = 14 %loop over flowline plts
    pp.H0 = hh(ip);
    pp.mdot = mm(ip);
    pp.epsxx = eps_xx(ip);
    pp.dhdt = dhdt(ip);
   
    
    %create anonymouse function of temp profile
    TT = T(:,ip);
    f = fit(zeta', TT, 'spline');
    anonT = @(z) f(z);

    %compute collapse time
    collapse_time = get_collapse_time_specTprof(pp, dt, tmax, anonT)


end
