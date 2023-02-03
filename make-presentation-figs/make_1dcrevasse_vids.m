% Make videos of the temperature profile, dimensionless temperature
% profile, crevasse depth and the stress intensity factor through time for
% (1) low melt rate, (2) high melt rate

%% Preliminaries
clear
addpath('..')
addpath('../functions')
make_vid = 1;
grounding_line_advection_profile = 1; %set to one to have the advected profile superimposed on the BL profile
%% Parameters
pp = struct;
for i = 1:2
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
        
        pp(i).H0    = 330;  %initial ice thickness
        pp(i).dt = 2;
    end
    pp(i).dhdt = -pp(i).mdot;
    pp(i).H = pp(i).H0; %set thickness to initial thickness

    %compute dimensionless quantities
    pp(i).l      = pp(i).kappa/pp(i).H0/pp(i).mdot; %boundary layer lengthscale
    pp(i).lambda = 2*sign(pp(i).epsxx)*abs(pp(i).epsxx^(1/pp(i).n)) * pp(i).B0 / pp(i).rhoi / pp(i).g /  pp(i).H0;
    pp(i).F      = pp(i).frac_tough / pp(i).H0 / pp(i).rhoi / pp(i).g;

    % timestepping quantities
    pp(i).t  = 0; %timenow
    %pp(i).ntmax = floor(abs(pp(i).H0/pp(i).dhdt));
    pp(i).ntmax = 1e6;
end

%% Generate data
ss = struct;
for ip = 1:2 %for each parameter set

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
        pp(ip).F      = pp(ip).frac_tough / pp(ip).H / pp(ip).rhoi / pp(ip).g;
        %^^^ comment line out to expose only temperature profile changes
        %pause
        i = i+1


    end

end

%% Make the video
t1 = [ss(1,:).t];
t2 = [ss(2,:).t];
tt = [length(t1),length(t2)];
nt = max([length(t1),length(t2)]);
colmap = zeros(2,nt,3);
cc = cmocean('matter', nt+20);
cc = cc(11:end-10,:);
colmap(1,:,:) = cc;

cc = cmocean('dense', nt+20);
cc = cc(11:end-10,:);
colmap(2,:,:) = cc;

for ip = 1:2
    if make_vid
        v = VideoWriter(pp(ip).fname, 'MPEG-4');
        v.FrameRate= 20;
        v.Quality =  100;
        open(v)
    end

    for it = 1:tt(ip) %loop over timesteps
        clf;

        l =  pp(ip).kappa /ss(ip,it).H / pp(ip).mdot;
        anonT = @(z) (pp(ip).Ts + (pp(ip).Tb - pp(ip).Ts)*exp(-z/l));
        if grounding_line_advection_profile
        anonT = @(z) TgfF(z) + (pp(ip).Tb- TgfF(z)).*exp(-z/l); %with advected grounding line contribution
        end

        %setup plot positions
        w = 0.14;
        gap = 0.08;
        startx = 0.1;
        shift2 = 0.06; %shift after second image
        starty = 0.12;
        h = 0.78; %height image

        positions = [0.1, starty, w, h;
            (w + startx + gap - shift2), starty, w, h;
            (2*w + startx + 2*gap - shift2), starty, w, h;
            (3*w + startx + 3*gap - shift2), starty, 2*w, h];

        fs = 16; %axis fontsize
        xl = pp(ip).Ts-1-273.15;
        xm = pp(ip).Tb + 1 - 273.15;
        zz = linspace(0,1,1e2);


        %plot the configuration in dimensional space
        ax(1) = subplot('Position', positions(1,:)); hold on; box on
        fill([-1,1,1,-1], ss(ip,it).H*[-0.9, -.9, -0.9 + ss(ip,it).dimless_crev_depth, -0.9 + ss(ip,it).dimless_crev_depth], 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2);
        fill([-1,1,1,-1], ss(ip,it).H*[-0.9 +  ss(ip,it).dimless_crev_depth, -0.9 +  ss(ip,it).dimless_crev_depth, 0.1, 0.1], 'b', 'LineStyle', 'none', 'FaceAlpha', 0.2);
        plot([-1,1], ss(ip,it).H*[.1,.1], 'k', 'linewidth', 2);
        plot([-1,1], ss(ip,it).H*[-0.9,-0.9], 'k', 'linewidth', 2);
        plot([-1,1], [0,0], 'k--', 'linewidth', 1.5);
        plot([0, 0], ss(ip,it).H*(-0.9 + [0, ss(ip,it).dimless_crev_depth]), 'k', 'linewidth', 2);

        ylabel('depth (m)', 'FontName', 'GillSans')
        ylim(pp(ip).H0*[-1,0.2]);
        xticks([]);
        yticks(-600:100:100)
        xlim([-1,1])
        ax(1).FontSize = fs;
        ax(1).FontName = 'GillSans';

        %plot the current dimensional temp profile
        ax(2) = subplot('Position', positions(2,:)); hold on; box on
        plot(100*[-1,1], ss(ip,it).H*[.1,.1], 'k', 'linewidth', 2);
        plot(100*[-1,1],  ss(ip,it).H*[-0.9,-0.9], 'k', 'linewidth', 2);
        plot(100*[-1,1], [0,0], 'k--', 'linewidth', 1.5);
        plot(anonT(zz)-273.15, zz* ss(ip,it).H + -0.9* ss(ip,it).H, 'linewidth', 1.75, 'color', colmap(ip,it,:));
        %ylabel('depth', 'FontName', 'GillSans')
        ylim(pp(ip).H0*[-1,0.2]);
        %xticks([]);
        yticks([]);
        xlim([xl*1.1,xm*0.9]);
        ax(2).FontSize = fs;
        ax(2).FontName = 'GillSans';
        xlabel('temp (C)', 'FontName', 'GillSans'); %ylabel('dimensionless depth', 'FontName', 'GillSans')


        % plot the temp profile in nondimensional depth for all times up to
        % now
        
        ax(3) = subplot('Position', positions(3,:)); box on

        hold on
        plot([xl,xm], [0,0], 'k', 'linewidth', 2);
        plot([xl,xm], [1,1], 'k', 'linewidth', 2);
        if it == 1 %save the time zero profile
            TT0 = anonT(zz)-273.15;         
        end
      
        plot(TT0, zz, 'k--','linewidth', 1.75);
        
        plot(anonT(zz)-273.15, zz, 'linewidth', 1.75,'Color',  colmap(ip,it, :));
        ylabel('dimensionless depth', 'FontName', 'GillSans')

        xlabel('temp (C)', 'FontName', 'GillSans'); %ylabel('dimensionless depth', 'FontName', 'GillSans')
        xlim([xl,xm])
        ylim([-.1,1.1]);
        ax(3).FontSize = fs;
        ax(3).FontName = 'GillSans';


        ax(4) = subplot('Position', positions(4,:)); hold on; box on
        xx = linspace(0,1, length(ss(ip,it).stress_intensity));
        plot(xx, ss(ip,1).stress_intensity, 'k--', 'linewidth', 2);
        plot(xx, ss(ip,it).stress_intensity, 'linewidth', 2,'color', colmap(ip,it, :));
        plot([0,1], [pp(ip).F, pp(ip).F], 'k--', 'linewidth', 2)
        xlabel('dimensionless basal crevasse depth', 'FontName', 'GillSans');
        ylabel('dimensionless stress intensity', 'FontName', 'GillSans')
        ylim([-1,1])
        %yticks([0,1])
        ax(4).FontSize = fs;
        ax(4).FontName = 'GillSans';

        fig = gcf;
        fig.Position(3:4) = [800, 420];
        fig.Color = [1,1,1];
        sgtitle(['t = ' num2str(ss(ip,it).t) 'yrs, H = ' num2str(ss(ip,it).H) 'm'] , 'FontName', 'GillSans', 'FontSize', 20);


%         pause
%         drawnow
        if make_vid
            frame = getframe(gcf);
            writeVideo(v,frame);
        end

    end %end loop over timesteps


    if make_vid
        close(v);
    end
end




%% functions
% plot the temperature profile and crack width
function plot_config(pp, dimless_crev_depth, stress_intensity, anonT)
%anonT = @(z) (pp.Ts + (pp.Tb - pp.Ts)*exp(-z/pp.l));

%setup plot positions
w = 0.14;
gap = 0.08;
startx = 0.1;
shift2 = 0.06; %shift after second image
starty = 0.12;
h = 0.78; %height image

positions = [0.1, starty, w, h;
    (w + startx + gap - shift2), starty, w, h;
    (2*w + startx + 2*gap - shift2), starty, w, h;
    (3*w + startx + 3*gap - shift2), starty, 2*w, h];

fs = 16; %axis fontsize
xl = pp.Ts-1-273.15;
xm = pp.Tb + 1 - 273.15;
zz = linspace(0,1,1e2);


%plot the configuration in dimensional space
ax(1) = subplot('Position', positions(1,:)); hold on; box on
fill([-1,1,1,-1], pp.H*[-0.9, -.9, -0.9 + dimless_crev_depth, -0.9 + dimless_crev_depth], 'r', 'LineStyle', 'none', 'FaceAlpha', 0.2);
fill([-1,1,1,-1], pp.H*[-0.9 + dimless_crev_depth, -0.9 + dimless_crev_depth, 0.1, 0.1], 'b', 'LineStyle', 'none', 'FaceAlpha', 0.2);
plot([-1,1], pp.H*[.1,.1], 'k', 'linewidth', 2);
plot([-1,1], pp.H*[-0.9,-0.9], 'k', 'linewidth', 2);
plot([-1,1], [0,0], 'k--', 'linewidth', 1.5);
plot([0, 0], pp.H*(-0.9 + [0, dimless_crev_depth]), 'k', 'linewidth', 2);

ylabel('depth (m)', 'FontName', 'GillSans')
ylim(pp.H0*[-1,0.2]);
xticks([]);
yticks(-600:100:100)
xlim([-1,1])
ax(1).FontSize = fs;
ax(1).FontName = 'GillSans';

%plot the dimensional temp profile
ax(2) = subplot('Position', positions(2,:)); hold on; box on
plot(100*[-1,1], pp.H*[.1,.1], 'k', 'linewidth', 2);
plot(100*[-1,1], pp.H*[-0.9,-0.9], 'k', 'linewidth', 2);
plot(100*[-1,1], [0,0], 'k--', 'linewidth', 1.5);
plot(anonT(zz)-273.15, zz*pp.H + -0.9*pp.H, 'linewidth', 1.75, 'color', [238, 75, 43]/255);
%ylabel('depth', 'FontName', 'GillSans')
ylim(pp.H0*[-1,0.2]);
%xticks([]);
yticks([]);
xlim([xl*1.1,xm*0.9]);
ax(2).FontSize = fs;
ax(2).FontName = 'GillSans';
xlabel('T', 'FontName', 'GillSans'); %ylabel('dimensionless depth', 'FontName', 'GillSans')



% plot the temp profile in nondimensional depth
ax(3) = subplot('Position', positions(3,:)); box on
hold on
plot([xl,xm], [0,0], 'k', 'linewidth', 2);
plot([xl,xm], [1,1], 'k', 'linewidth', 2);
plot(anonT(zz)-273.15, zz, 'linewidth', 1.75, 'color', [238, 75, 43]/255);
ylabel('dimensionless depth', 'FontName', 'GillSans')


xlabel('T', 'FontName', 'GillSans'); %ylabel('dimensionless depth', 'FontName', 'GillSans')
xlim([xl,xm])
ylim([-.1,1.1]);
ax(3).FontSize = fs;
ax(3).FontName = 'GillSans';


ax(4) = subplot('Position', positions(4,:)); hold on; box on
xx = linspace(0,1, length(stress_intensity));
plot(xx, stress_intensity, 'k', 'linewidth', 2);
plot([0,1], [pp.F, pp.F], 'k--', 'linewidth', 2)
xlabel('dimensionless basal crevasse depth', 'FontName', 'GillSans');
ylabel('dimensionless stress intensity', 'FontName', 'GillSans')
ylim([-1,1])
%yticks([0,1])
ax(4).FontSize = fs;
ax(4).FontName = 'GillSans';

fig = gcf;
fig.Position(3:4) = [800, 420];
fig.Color = [1,1,1];
end
