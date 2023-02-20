function  collapse_time_row = get_collapse_time_row_Nye(row_H, row_dhdt, row_epsxx, row_mdot, row_tags, row_dt, row_tmax,...
                                Tb, Ts, B0, rhoi, g, kappa, n, frac_tough, ghf)
ny = length(row_H);
collapse_time_row = nan(1,ny);
for iy = 1:ny
    if row_tags(iy) == 6 
        %set up the parameters
        pp = struct;

        %variable quantities
        pp.H0    = row_H(iy);   %initial ice thickness
        pp.dhdt  = row_dhdt(iy);  %rate of change of thickness (negative for thinning)
        pp.epsxx = row_epsxx(iy);  %strain rate
        pp.mdot  = row_mdot(iy);     %melt rate

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
        pp.l     = pp.kappa/pp.H0/pp.mdot; %initial boundary layer lengthscale
        pp.F = pp.frac_tough / (pp.H0)^(3/2) / pp.rhoi / pp.g;
        pp.lambda = 2*sign(pp.epsxx)*abs(pp.epsxx^(1/pp.n)) * pp.B0 / pp.rhoi / pp.g /  pp.H0;


        %compute the collapse time
        %dt = abs(1/dhdts(ix,iy)); %larger timesteps for smaller dhdt
 
        dt = row_dt(iy); %larger timesteps for smaller dhdt
        tmax = row_tmax(iy);

        %iy
        nt = 100; %number of trial time output points
        collapse_time_row(iy) = get_collapse_time_advect_Nye(pp, nt, tmax);


        %fprintf('Completed %.3f percent of collase time points in square \n', count*100/ sum(sum(tags==6)))
        %toc
    end
end