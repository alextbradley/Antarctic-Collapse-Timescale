function dimless_crev_depth_row = get_dim_crev_depth_row_linear(row_H, row_epsxx, ...
    Tb, Ts, B0, rhoi, g, kappa, n, frac_tough, ghf)


dimless_crev_depth_row = nan(size(row_H));
ny = length(row_H);
for iy = 1:ny
    if (~isnan(row_H(iy)) && ~isnan(row_epsxx(iy)))
        %set up parameters


        pp = struct;
        %variable quantities
        pp.H0    = row_H(iy);   %initial ice thickness
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


        % dimensionless quantities
        pp.l      = pp.kappa/pp.H0/pp.mdot; %boundary layer lengthscale
        pp.lambda = 2*sign(pp.epsxx)*abs(pp.epsxx^(1/pp.n)) * pp.B0 / pp.rhoi / pp.g /  pp.H0;
        pp.F      = pp.frac_tough / pp.H0 / pp.rhoi / pp.g;


        %set up temp profile
        anonT = @(z) (Ts + (Tb - Ts)*z);

        dimless_crev_depth_row(iy) = get_dimless_crev_depth(pp, anonT);


    end
end
