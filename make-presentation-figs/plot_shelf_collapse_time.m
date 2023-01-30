% Plot the collapse time for a specific ice shelf using the scattered
% interpolant (advective approximation)
addpath('../functions/')
%% Load the ice shelf data
f  = load('../data/ice_sheet_data.mat');
shelf_name= 'PineIsland';
shelf = shelf_name;
fname = strcat('../data/ice-shelves/' ,shelf, '.mat');
g = load(fname);

%% create restricted co-ordinates
[rId, cId] = find(g.IN) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
xs = f.xx(xminidx:xmaxidx);
ys = f.yy(yminidx:ymaxidx); %restricted co-ordinates
hs = f.H(xminidx:xmaxidx, yminidx:ymaxidx);
ms = f.m(xminidx:xmaxidx, yminidx:ymaxidx);
dhdts = f.dhdtadj(xminidx:xmaxidx, yminidx:ymaxidx);
strains = f.eflow(xminidx:xmaxidx, yminidx:ymaxidx);
Bs = f.B(xminidx:xmaxidx, yminidx:ymaxidx);



ms(ms < 1e-1) = 1e-1; %set a minimum melt value

%% loop over grid points
% tags:
%    1: outside domain (no thickness data)
%    2: other missing data (no melt, strain rate, dhdt)
%    3: negative strain rate
%    4: positive thinning rate (currently stable)
%    5: positive thinning rate (currently unstable)
%    6: other (good data)

collapse_time_square = nan(size(ms));
tags = nan(size(ms));

% figure(1); clf; hold on
% ylim([min(xs), max(xs)])
% xlim([min(ys), max(ys)])

count = 1;
nruns = length(xs)*length(ys);
for ix =  1:length(xs)
    for iy = 1:length(ys)

        %set up the parameters
        pp = struct;
        pp.H0    = hs(ix,iy);   %initial ice thickness
        pp.Tb    = -2 + 273.15;     %basal temperature (kelvin)
        pp.Ts    = -18 + 273.15;    %surface temp
        pp.dhdt  = dhdts(ix,iy);      %rate of change of thickness (negativr for thinning)
        pp.B0    = 1.928;  %viscosity constant
        pp.rhoi  = 918.0;  %ice density
        pp.g     = 9.81;   %gravitational acceleration
        pp.epsxx = strains(ix,iy);  %strain rate
        pp.kappa = 36;     %ice diffusivity
        pp.mdot  = ms(ix,iy);     %melt rate
        pp.n     = 3;      %glen flow coeff
        pp.l     = pp.kappa/pp.H0/pp.mdot; %initial boundary layer lengthscale
        pp.frac_tough = 150*1e3;
        pp.F = pp.frac_tough / pp.H0 / pp.rhoi / pp.g;
        pp.ghf = 48; %geothermal heat flux

        tmax = 1000;
        dt = 1;

        %temp profile
        pp.lambda = 2*sign(pp.epsxx)*abs(pp.epsxx^(1/pp.n)) * pp.B0 / pp.rhoi / pp.g /  pp.H0;
        TgfF = get_grounding_line_temp(pp.ghf, pp.Ts, pp.H0);
        anonT = @(z) TgfF(z) + (pp.Tb- TgfF(z)).*exp(-z/pp.l); %with advected grounding line contribution

        %determine the tag
        if isnan(hs(ix,iy))
            tags(ix,iy) = 1;
        elseif isnan(ms(ix,iy)) || isnan(dhdts(ix,iy)) || isnan(strains(ix,iy))
            tags(ix,iy) = 2;
        elseif strains(ix,iy) < 0
            tags(ix,iy) = 3;
        elseif dhdts(ix,iy) > 0 %thickening
            %compute current crevasse depth
            dimless_crev_depth = get_dimless_crev_depth(pp, anonT);

            if dimless_crev_depth < 1
                tags(ix,iy) = 4;
            else
                tags(ix,iy) = 5;
            end

        else
            tags(ix,iy) = 6;

        end


    end
     fprintf('Completed x row. Progress: %.3f percentage of grid points in square \n', count*100/ length(xs))
     count = count + 1;

end

%%
figure(10); clf; contourf(ys,xs,flipud(tags), 20, 'linestyle', 'none')
c = colorbar;
clim([1,6])
colormap(lines(6))
c.Ticks = linspace(2-(7/12),5+(7/12),6);
c.TickLabels = {'1','2','3','4','5','6'};

%% compute the collapse time for those applicable
count = 1;
tic
for ix =  1:length(xs)
    for iy = 1:length(ys)
        if tags(ix,iy) == 6
        %set up the parameters
        %
        pp = struct;
        pp.H0    = hs(ix,iy);   %initial ice thickness
        pp.Tb    = -2 + 273.15;     %basal temperature (kelvin)
        pp.Ts    = -18 + 273.15;    %surface temp
        pp.dhdt  = dhdts(ix,iy);      %rate of change of thickness (negative for thinning)
        pp.B0    = 1.928;  %viscosity constant
        pp.rhoi  = 918.0;  %ice density
        pp.g     = 9.81;   %gravitational acceleration
        pp.epsxx = strains(ix,iy);  %strain rate
        pp.kappa = 36;     %ice diffusivity
        pp.mdot  = ms(ix,iy);     %melt rate
        pp.n     = 3;      %glen flow coeff
        pp.l     = pp.kappa/pp.H0/pp.mdot; %initial boundary layer lengthscale
        pp.frac_tough = 150*1e3;
        pp.F = pp.frac_tough / pp.H0 / pp.rhoi / pp.g;
        pp.ghf = 48; %geothermal heat flux

        tmax = 1000;
        dt = 1;

        %temp profile
        pp.lambda = 2*sign(pp.epsxx)*abs(pp.epsxx^(1/pp.n)) * pp.B0 / pp.rhoi / pp.g /  pp.H0;
        TgfF = get_grounding_line_temp(pp.ghf, pp.Ts, pp.H0);
        anonT = @(z) TgfF(z) + (pp.Tb- TgfF(z)).*exp(-z/pp.l); %with advected grounding line contribution

        %compute the collapse time
        %dt = abs(1/dhdts(ix,iy)); %larger timesteps for smaller dhdt
        dt = 1;
        dt = max([5,abs(1/dhdts(ix,iy))]); %larger timesteps for smaller dhdt

        %ix,iy%, pause

        [dimless_crev_depth,stress_intensity] = get_dimless_crev_depth(pp, anonT);
        collapse_time_square(ix,iy) = get_collapse_time_advect(pp, dt, tmax);
        
        fprintf('Completed %.3f percent of collase time points in square \n', count*100/ sum(sum(tags==6)))
        toc
        count = count + 1;
        end
    end
end

%% make the plot
figure(1); clf;  
ax(1) = gca; 
cmap1 = cmocean('matter',100);
%contourf(flipud(log10(collapse_time_square)), 50, 'linestyle', 'none')
pl = imagesc((log10(collapse_time_square)));
set(pl, 'AlphaData', ~isnan(collapse_time_square)); 
axis equal

c = colorbar;
colormap(ax(1), cmap1);
cmin = 1; cmax = 4;

clabels = 10.^(cmin:1:cmax); clim(log10(clabels([1 end]))); set(c,'Ticks',log10(clabels),'TickLabels',clabels); 
c.TickLabels = {'<10^1', '10^2','10^3', '>10^4'};
c.FontSize = 14;
c.FontName = 'GillSans';
c.Label.String = 'time to crevasse (yrs)';
c.Label.FontSize = 16;
ax(1).XTick = [];
ax(1).YTick = [];
ax(1).Visible = 'off';


% add the missing data
ax(2) = axes();
tags_missing = nan(size(tags));
tags_missing(tags == 2) = 1;
pl = imagesc(tags_missing);
set(pl, 'AlphaData', ~isnan(tags_missing)); 
colormap(ax(2), 0.6* [1,1,1])
axis equal


% add the negative strain rate as maximum
ax(3) = axes();
tags_negative_strain = nan(size(tags));
tags_negative_strain(tags == 3) = 1;
pl = imagesc( tags_negative_strain);
set(pl, 'AlphaData', ~isnan(tags_negative_strain)); 
colormap(ax(3), cmap1(end,:))
axis equal


% add the positive thinning rate as maximum
ax(4) = axes();
tags_pos_thin = nan(size(tags));
tags_pos_thin(tags == 4) = 1;
pl = imagesc(tags_pos_thin);
set(pl, 'AlphaData', ~isnan(tags_pos_thin)); 
colormap(ax(4), cmap1(end,:))
axis equal


% add positive thinning rate, currently unstable as minimum 
ax(5) = axes();
tags_pos_thin = nan(size(tags));
tags_pos_thin(tags == 5) = 1;
pl = imagesc( tags_pos_thin);
set(pl, 'AlphaData', ~isnan(tags_pos_thin)); 
colormap(ax(5), cmap1(1,:))
axis equal


% add zero collapse time as a minimum
ax(6) = axes();
zero_collapse = nan(size(collapse_time_square));
zero_collapse(collapse_time_square == 0) = 1;
pl = imagesc(zero_collapse);
set(pl, 'AlphaData', ~isnan(zero_collapse)); 
colormap(ax(6), cmap1(1,:))
axis equal



for i = 1:6
    ax(i).Position = ax(1).Position;
    ax(i).Visible = 'off';

end