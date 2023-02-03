% Plot the collapse time for a specific ice shelf using the scattered
% interpolant (advective approximation)
addpath('../functions/')

% parallel tasks
% delete(gcp('nocreate'))
% local_cluster = parcluster('local')
% parpool(local_cluster)
%% Load the ice shelf data
%f  = load('../data/ice_sheet_data.mat');
shelf_name= 'Ross';
shelf = shelf_name;
fname = strcat('../data/ice-shelves/' ,shelf, '.mat');
g = load(fname);

%% create restricted co-ordinates
[rId, cId] = find(g.IN) ; %rows and columns of non-zero entries
xminidx = min(rId); xmaxidx = max(rId);
yminidx = min(cId); ymaxidx = max(cId); %indices of min and max in rows and columns
step = 4; %step in grid cell. Resolution becomes 500m*step
xs = f.xx(xminidx:step:xmaxidx);
ys = f.yy(yminidx:step:ymaxidx); %restricted co-ordinates
hs = f.H(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
ms = f.m(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
dhdts = f.dhdtadj(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
strains = f.eflow(xminidx:step:xmaxidx, yminidx:step:ymaxidx);
Bs = f.B(xminidx:step:xmaxidx, yminidx:step:ymaxidx);


dhdts = -abs(dhdts);
ms(ms < 1e-1) = 1e-1; %set a minimum melt value

%% determine tags
% tags:
%    1: outside domain (no thickness data)
%    2: other missing data (no melt, strain rate, dhdt)
%    3: negative strain rate
%    4: positive thinning rate (currently stable)
%    5: positive thinning rate (currently unstable)
%    6: other (good data)

tags = nan(size(ms));
tags = 6*ones(size(ms));

tags(dhdts > 0) = 4;
tags(strains < 0) = 3;
tags(isnan(ms) | isnan(dhdts) | isnan(strains)) = 2;
tags(isnan(hs)) = 1;


% figure(1); clf; hold on
% ylim([min(xs), max(xs)])
% xlim([min(ys), max(ys)])

% count = 1;
% nruns = length(xs)*length(ys);
% for ix =  1:length(xs)
%     for iy = 1:length(ys)
% 
%         %set up the parameters
%         pp = struct;
%         pp.H0    = hs(ix,iy);   %initial ice thickness
%         pp.Tb    = -2 + 273.15;     %basal temperature (kelvin)
%         pp.Ts    = -18 + 273.15;    %surface temp
%         pp.dhdt  = dhdts(ix,iy);      %rate of change of thickness (negativr for thinning)
%         pp.B0    = 1.928;  %viscosity constant
%         pp.rhoi  = 918.0;  %ice density
%         pp.g     = 9.81;   %gravitational acceleration
%         pp.epsxx = strains(ix,iy);  %strain rate
%         pp.kappa = 36;     %ice diffusivity
%         pp.mdot  = ms(ix,iy);     %melt rate
%         pp.n     = 3;      %glen flow coeff
%         pp.l     = pp.kappa/pp.H0/pp.mdot; %initial boundary layer lengthscale
%         pp.frac_tough = 150*1e3;
%         pp.F = pp.frac_tough / pp.H0 / pp.rhoi / pp.g;
%         pp.ghf = 48; %geothermal heat flux
% 
%         tmax = 1000;
%         dt = 1;
% 
%         %temp profile
%         pp.lambda = 2*sign(pp.epsxx)*abs(pp.epsxx^(1/pp.n)) * pp.B0 / pp.rhoi / pp.g /  pp.H0;
%         TgfF = get_grounding_line_temp(pp.ghf, pp.Ts, pp.H0);
%         anonT = @(z) TgfF(z) + (pp.Tb- TgfF(z)).*exp(-z/pp.l); %with advected grounding line contribution
% 
%         %determine the tag
%         if isnan(hs(ix,iy))
%             tags(ix,iy) = 1;
%         elseif isnan(ms(ix,iy)) || isnan(dhdts(ix,iy)) || isnan(strains(ix,iy))
%             tags(ix,iy) = 2;
%         elseif strains(ix,iy) < 0
%             tags(ix,iy) = 3;
%         elseif dhdts(ix,iy) > 0 %thickening
%             %compute current crevasse depth
%             dimless_crev_depth = get_dimless_crev_depth(pp, anonT);
% 
%             if dimless_crev_depth < 1
%                 tags(ix,iy) = 4;
%             else
%                 tags(ix,iy) = 5;
%             end
% 
%         else
%             tags(ix,iy) = 6;
% 
%         end
% 
% 
%     end
%      fprintf('Completed x row. Progress: %.3f percentage of grid points in square \n', count*100/ length(xs))
%      count = count + 1;
% 
% end

%
figure(10); clf; contourf(ys,xs,flipud(tags), 20, 'linestyle', 'none')
c = colorbar;
clim([1,6])
colormap(lines(6))
c.Ticks = linspace(2-(7/12),5+(7/12),6);
c.TickLabels = {'1','2','3','4','5','6'};

%% compute the collapse time for those applicable
collapse_time_square = nan(size(ms));
%tic
parfor ix =  1:length(xs)
    %variable parameters
    row_H = hs(ix,:);
    row_dhdt = dhdts(ix,:);
    row_epsxx = strains(ix,:);
    row_tags = tags(ix,:);
    row_mdot = ms(ix,:);

    %constant parameters
    Tb    = -2 + 273.15;     %basal temperature (kelvin)
    Ts    = -18 + 273.15;    %surface temp
    B0    = 1.928;  %viscosity constant
    rhoi  = 918.0;  %ice density
    g     = 9.81;   %gravitational acceleration
    kappa = 36;     %ice diffusivity
    n     = 3;      %glen flow coeff
    frac_tough = 150*1e3;
    ghf = 48; %geothermal heat flux

    %get the collapse time on this row
    collapse_time_row = get_collapse_time_row(row_H, row_dhdt, row_epsxx, row_mdot, row_tags, ...
                                                Tb, Ts, B0, rhoi, g, kappa, n, frac_tough, ghf);

    collapse_time_square(ix,:) = collapse_time_row;
    ix

end

%% make the plot
figure(2); clf;  
ax(1) = gca; 
cmap1 = cmocean('matter',100);
%contourf(flipud(log10(collapse_time_square)), 50, 'linestyle', 'none')
pl = imagesc((log10(collapse_time_square)));
set(pl, 'AlphaData', ~isnan(collapse_time_square)); 
axis equal
%
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