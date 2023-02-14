% Make videos of the temperature profile, dimensionless temperature
% profile, crevasse depth and the stress intensity factor through time for
% (1) low melt rate, (2) high melt rate

%% Preliminaries
%clear
addpath('../..')
addpath('../../functions')
rmpath('functions-old/')

grounding_line_advection_profile = 1; %set to one to have the advected profile superimposed on the BL profile
%% Parameters
pp = struct;
for i = 1
    pp(i).Tb    = -5 + 273.15;     %basal temperature (kelvin)
    pp(i).Ts    = -20 + 273.15;    %surface temp
    pp(i).B0    = 1.928;  %viscosity constant
    pp(i).rhoi  = 918.0;  %ice den sity
    pp(i).g     = 9.81;   %gravitational acceleration
    pp(i).epsxx = 0.001;  %strain rate
    pp(i).kappa = 36;     %diffusivity
    pp(i).n     = 3;      %glen flow coeff
    pp(i).frac_tough = 150*1e3;

    % case specific quantities
    if i == 1 %warm case
        pp(i).mdot = 5;
        pp(i).fname = "high_melt_rate.mp4"; %save file name
        if grounding_line_advection_profile
            pp(i).fname = "high_melt_rate_GroundingLineAdvection.mp4";
        end
        pp(i).H0    = 400;    %initial ice thickness
        pp(i).dt = 0.2;
    else
        pp(i).mdot = 0.5;
        pp(i).fname = "low_melt_rate.mp4"; %save file name
        if grounding_line_advection_profile
            pp(i).fname = "low_melt_rate_GroundingLineAdvection.mp4";
        end
        
        pp(i).H0  = 330;  %initial ice thickness
        pp(i).dt = 10;
    end
    pp(i).dhdt = -pp(i).mdot;
    pp(i).H = pp(i).H0; %set thickness to initial thickness

    %compute dimensionless quantities
    pp(i).l      = pp(i).kappa/pp(i).H0/pp(i).mdot; %boundary layer lengthscale
    pp(i).lambda = 2*sign(pp(i).epsxx)*abs(pp(i).epsxx^(1/pp(i).n)) * pp(i).B0 / pp(i).rhoi / pp(i).g /  pp(i).H0;
    pp(i).F      = pp(i).frac_tough / (pp(i).H0)^(3/2) / pp(i).rhoi / pp(i).g;

    % timestepping quantities
    pp(i).t  = 0; %timenow
    %pp(i).ntmax = floor(abs(pp(i).H0/pp(i).dhdt));
    pp(i).ntmax = 1e6;
end

%% Generate data
ss = struct;
for ip = 1 %for each parameter set

    i = 1;
    dimless_crev_depth = 0; %seed the dimensionless crevasse depth
    ghf = 48; 
    TgfF = get_grounding_line_temp(ghf, pp(ip).Ts, pp(ip).H0);

    while (i < pp(ip).ntmax) && (dimless_crev_depth < 1)
        pp(ip).lambda;
        
        %anonT = @(z) (pp(ip).Ts + (pp(ip).Tb -pp(ip).Ts)*exp(-z/pp(ip).l)); %purely exponential 
        if grounding_line_advection_profile
        anonT = @(z) TgfF(z) + (pp(ip).Tb - TgfF(z)).*exp(-z/pp(ip).l); %with advected grounding line contribution
        end
        [dimless_crev_depth, stress_intensity] = get_dimless_crev_depth(pp(ip), anonT);
        DD(i) = dimless_crev_depth;
        
%         clf;plot_config(pp(ip), dimless_crev_depth, stress_intensity, anonT);
%         sgtitle(['t = ' num2str(pp(ip).t) ', H = ' num2str(pp(ip).H)], 'FontName', 'GillSans', 'FontSize', 20);
%         drawnow; %pause

        % store the solution
        ss(ip,i).t = pp(ip).t;
        ss(ip,i).H = pp(ip).H;
        ss(ip,i).dimless_crev_depth = dimless_crev_depth;
        ss(ip,i).stress_intensity = stress_intensity;

        %update the thickness and dimensionless quantities
        pp(ip).H     = pp(ip).H + pp(ip).dt*pp(ip).dhdt;
        pp(ip).t      = pp(ip).t + pp(ip).dt;
        pp(ip).l      = pp(ip).kappa/pp(ip).H/pp(ip).mdot; %boundary layer lengthscale
        pp(ip).lambda = 2*sign(pp(ip).epsxx)*abs(pp(ip).epsxx^(1/pp(ip).n)) * pp(ip).B0 / pp(ip).rhoi / pp(ip).g /  pp(ip).H;
        pp(ip).F      = pp(ip).frac_tough / (pp(ip).H)^(3/2) / pp(ip).rhoi / pp(ip).g;
        %^^^ comment line out to expose only temperature profile changes
        %pause
        i = i+1

    end

end

%profiling: default - 20s
%           change getG to explicit evaluation - 10s
%           remove row by row evaluation of T in get_stress_intensity - 3.2

% 15 s new, 55s old (already 2x speed up)
% --> down to 11s with vectorizing G in stress intensity
% --> down to 7s in fb
% --> down to 4s with vectorize G in gb




