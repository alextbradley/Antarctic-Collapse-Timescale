%make a contour plot of crevasse timescale and dimensionless timescale as a
%function of initial thickness and melt rate and a plot of the dimensionless collapse
%timescale as a ice thickness
addpath('../functions')
HH = 40:40:1000;
%HH = 160:20:480;
mm = logspace(-1,2,13);
%mm = logspace(-1,2,10);

ct = zeros(length(HH), length(mm));
lambdac0 = zeros(length(HH), length(mm));
count = 1;
for im = 1:length(mm)
    for ih = 1:length(HH)
        % params
        tic
        pp = struct;
        pp.H0    = HH(ih);   %initial ice thickness
        pp.Tb    = -5 + 273.15;     %basal temperature (kelvin)
        pp.Ts    = -20 + 273.15;    %surface temp at grounding line
        pp.dhdt  = -1;      %rate of change of thickness
        pp.B0    = 1.928;  %viscosity constant
        pp.rhoi  = 918.0;  %ice density
        pp.g     = 9.81;   %gravitational acceleration
        pp.epsxx = 0.001;  %strain rate
        pp.kappa = 36;     %diffusivity
        pp.mdot  = mm(im);     %melt rate
        pp.n     = 3;      %glen flow coeff
        pp.l     = pp.kappa/pp.H0/pp.mdot; %initial boundary layer lengthscale
        pp.frac_tough = 150*1e3;
        pp.F = pp.frac_tough / pp.H0 / pp.rhoi / pp.g;
        pp.ghf = 48; %geothermal heat flux 

        tmax = 1000;
        dt = 1;

        ct(ih,im) = get_collapse_time_advect(pp, dt, tmax); 
        
        %
        % ct(ih,im) = get_collapse_time(pp, dt, tmax); %for only exponential profile

        %compute the critical lambda
        lambdac0(ih,im) = get_critical_lambda(pp);

        fprintf('completed %.3f percent \n', 100* count/(length(HH)*length(mm)));
        toc
        count = count +1;
    end
end

%% make plots
figure(1);clf;
subplot(1,3,1);
contourf(mm,HH, lambdac0, 50, 'linestyle', 'none');
c = colorbar; c.Label.String = '\lambda_{c,0}';
xlabel('melt rate');
set(gca, 'XScale', 'log');
ylabel('H_0');
ax = gca;
ax.FontSize = 12;

subplot(1,3,2);
ct(ct <= 0) = nan;
contourf(mm, HH, ct./abs(pp.dhdt), 50, 'linestyle', 'none');
c = colorbar; c.Label.String = 'collapse time/ dhdt';
xlabel('melt rate');
set(gca, 'XScale', 'log');
ylabel('H_0');
ax = gca;
ax.FontSize = 12;

subplot(1,3,3);
beta = 2*pp.epsxx^(1/pp.n) * pp.B0 / pp.rhoi / pp.g; %=lambda * H
Hr = repmat(HH', [1,length(mm)]);
xc = -(Hr - beta ./ lambdac0) / pp.dhdt; %x co-ordinates of collapsed data
xc(xc<0) = nan;
contourf(mm , HH, log10((ct./xc - 1) * 100), 50, 'linestyle', 'none');

c = colorbar; c.Label.String = '% enhancement over no melt change time';
xlabel('melt rate');
set(gca, 'XScale', 'log');
ylabel('H_0');
ax = gca;
ax.FontSize = 12;
clim([1,2]);
%c.Ticks = 0:0.5:1.5;
%c.TickLabels = []

%%
figure(2); clf; hold on; box on
colmap = flipud(cmocean('matter',length(mm)+2));
colmapcts = cmocean('matter');
colmapcts = colmapcts(40:end,:);
for im = 1:length(mm)
    yy = ct(:,im)./xc(:,im);
    %yy = smooth(yy,8);
    plot(HH, yy, 'linewidth', 1.5, 'color', colmap(im+2,:));
    %drawnow; pause
end
ax = gca;
plot([min(HH), max(HH)], [1,1], 'k--', 'linewidth', 1.5)
ax.FontSize = 16;
%ax.FontName = 'GillSans';
ylim([0.9, 2])
xlim([min(HH), max(HH)])
xlabel('initial ice thickness (m)', 'FontSize',16)%, 'FontName', 'GillSans')
cc = colorbar;
colormap(colmapcts)
%clabels = 10.^(-1:1:2);
cc.Ticks = 1/3*[0,1,2,3];
cc.TickLabels = {'10^{-1}','10^{0}','10^{1}','10^{2}'};
cc.Label.String = 'melt rate (m/yr)';
cc.Label.FontSize = 16;
%cc.Limits = [min(log10(clabels)), max(log10(clabels))];
%set(cc,'Ticks',log10(clabels),'TickLabels',clabels);

ax.XTick = (0:200:1000)
