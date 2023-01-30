function [T,z,xx] = get_flowline_temp(x, zeta,H,u,S, mdot, Tgf)
% Return the temperature T(x,z) along the flowline.
% Inputs: 
% x:    (nx x 1) horizontal co-ordinates (defines nx)
% zeta: (1 x ny) vertical dimensionless co-ordinates
% H:    (nx x 1) array of ice thickness along flowline
% u:    (nx x 1) array of ice speed along flowline        
% S:    (nx x 1) array of salinities along flowline
% mdot: (nx x 1) array of melt rates along flowline. NB: positive melting
% correspons to ice removal here!
% Tgf:  anonymous function of grounding line temp vertical profile 
      

%% Compute auxillary quantities
b = -918/1028*H; %ice shelf base
q = H.*u;  %ice flux
qg = q(1); %grounding line flux
kappai = 36;        %ice thermal diffusivity

%compute freezing temp along flowline
T0     = 8.32*1e-2; %freezing point offset
Lambda = 7.61*1e-4; %freezing point depth constant
Gamma  = 5.73*1e-2; %freezing point salinity constant
Tstar = T0 + Lambda*b - Gamma*S; %freezing temp

% remove negative melt rates for which theory breaks down (refreezing)
mdot(mdot<0) = nan;

%% Loop over flowline co-ordinates (x and zeta), constructing temp at each
%setup storage
z = nan(length(zeta), length(x));
T = nan(length(zeta), length(x));
xx = nan(length(zeta), length(x));

for izeta = 1:length(zeta)
    for ix = 1:length(x)

        % compute xi = 1 - (1-zeta)/qg * q
        xi = 1 - (1-zeta(izeta))/qg * q(ix);
        xi0 = 1 - (1-zeta(1))/qg * q(ix);    %xi at the grounding line

        %restrict xi, xi0 values to interval [0,1]
        if xi0 > 1
            xi0 = 1;
        elseif xi0 < 0
            xi0 = 0;
        end
        if xi > 1
            xi = 1;
        elseif xi < 0
            xi = 0;
        end
         
        z(izeta,ix) = zeta(izeta)*H(ix) + b(ix); 
        xx(izeta,ix) = x(ix);
        T(izeta, ix) = Tgf(xi) + (Tstar(ix) - Tgf(xi0))*exp(-mdot(ix) * H(ix)/kappai * zeta(izeta));
        if T(izeta, ix) >0 
            Tstar(ix)
            Tgf(xi0)
            mdot(ix)
             H(ix)
             kappai
             zeta(izeta)
           

            pause
        end

    end
end