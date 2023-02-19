% Make figure 2 of the ms:
% (a) contour plot of l = kappa/(mH), with ice shelves labelled
% (b) locations of ice shelves in (l, strain rate,  dhdt) space
%
% 17/02/23, ATB (aleey@bas.ac.uk), MIT licence
%% Preliminaries
%clear
addpath('../functions');
fig = figure(1);clf
positions = [0.05, 0.12, 0.39, 0.83;
             0.58, 0.12, 0.38, 0.83];
sz        = size(positions);
fig.Position(3:4) = [1040,430];
fs = 14; %fontsize
for i = 1:2
    ax(i) = subplot('Position', positions(i,:)); hold on
    ax(i).FontSize = fs;
end

%colourmap
cmap = cmocean('matter', 100);
cmap = flipud(cmap(1:end-60,:));
colormap(ax(1) ,cmap)

warm_col = [203,0,63]/255; %warm shelves colour
cold_col = [0, 63, 203]/255; %cold shelves colour
grey_col =  0.5*[1,1,1];
%% Make (a)

%
% Make the contour map
%
kappai = 36; %thermal diffusivity
H = linspace(100,800,100); %ice thicknesses
mdot = logspace(-1,log10(20),100);
[HH,mm] = meshgrid(H,mdot);
ell = kappai ./ HH ./ mm;

% contour plot
subplot('Position', positions(1,:))
%contourf(H, mdot, log10(ell), 100, 'linestyle', 'none');
imagesc(H,mdot, log10(ell));
set(ax(1), 'YDir', 'normal');
set(ax(1), 'YScale', 'log');

%labels
xlabel('ice thickness (m)', 'interpreter', 'none')
ylabel('melt rate (m/yr)', 'interpreter', 'none')

%colourbar
cc = colorbar;
clabels = 10.^(-3:1); 
clim(log10(clabels([1 end])));
set(cc,'Ticks',log10(clabels),'TickLabels',clabels);
cc.Label.Interpreter = 'latex';
cc.Label.String = '$\ell$';
cc.Label.FontSize = 15;
cc.Position(1) = 0.45;

% sort out limits
xl = [100,800];
yl = 10.^[min(log10(mdot)),max(log10(mdot))]; %x and y lims of plot
ax(1).XTick = (0:200:800);
ax(1).XLim = [xl(1)-1, xl(2)];
ax(1).YLim = [yl(1), yl(2)+0.1];
box(ax(1), 'on')

%
% add contours
%
hold on
levs = [-1.2, -1.5];
for i = 1:length(levs)
contour(H, mdot, log10(ell), levs(i)*[1,1], 'LineStyle','--', 'Color', grey_col, 'LineWidth', 2);
end


%
%%%%%% add shelf points %%%%%%%
%

% init storage
jdir = dir('../data/ice-shelves/major/*.mat'); %where shelf files are stored


m_ave = zeros(1,length(jdir));
h_ave = zeros(1,length(jdir));
mstd = zeros(1,length(jdir));
hstd = zeros(1,length(jdir));
percov = zeros(1,length(jdir)); %store the % coverage
shelf_type = zeros(1,length(jdir)); %flags for keeping the shelf or not (0), cold (1), warm (2)
epsxx_ave = zeros(1,length(jdir));
thinrate_ave = zeros(1,length(jdir));

% loop over shelves
shelf_names = strings;
for i = 1:length(jdir)
    shelf = jdir(i).name;
    shelf = strrep(shelf,'.mat',''); %strip the .mat at the end
    shelf_names(i) = shelf;
    
    fname = strcat('../data/ice-shelves/major/' ,shelf, '.mat');
    g = load(fname);

    % restrict data to those points in shelf with high enough melt rate
    m_shelf = f.m; m_shelf = m_shelf(g.IN); %only the points in this shelf
    h_shelf = f.H; h_shelf = h_shelf(g.IN);
    epsxx_shelf = f.eflow; epsxx_shelf = epsxx_shelf(g.IN);
    dhdt_shelf = f.dhdtadj; dhdt_shelf = dhdt_shelf(g.IN);
    idx =  (~isnan(h_shelf)) &  (m_shelf > 1e-6)  & (epsxx_shelf > 1e-6) &  (-dhdt_shelf > 1e-6) ; %points where we have point thickness,  melt rate > 0, strain > 0, thinning < 0
    idx_for_calc =  (~isnan(h_shelf)) &  (~isnan(m_shelf))  & (~isnan(epsxx_shelf)) &  (~isnan(-dhdt_shelf)); %points where we have data for 

    area_shelf(i) = length(h_shelf(~isnan(h_shelf))); %size in km^3 (grid size = 1e3 * 1e3)

    %idx =  (~isnan(h_iceshelf)) &  (~isnan(m_iceshelf)); %points where we have point thickness and melt rate > 0
    h_shelf = h_shelf(idx);
    m_shelf = m_shelf(idx); %arrays with points in particular shelf with both melt and thicknes
    epsxx_shelf = epsxx_shelf(idx);
    thinrate_shelf = -dhdt_shelf(idx);
    m_ave(i) = median((m_shelf));
    h_ave(i) = median((h_shelf));
    epsxx_ave(i) = median(epsxx_shelf);
    thinrate_ave(i) = median(thinrate_shelf);
    

    % plot point and add name
    shelf_total = area_shelf(i); %total number of points in shelf
    shelf_keep  = sum(idx_for_calc); %number of points with thickness data
    percov(i) = shelf_keep/shelf_total * 100;
    threshold_keep = 40;
    if (percov(i) > threshold_keep) || shelf == "Thwaites" %only plot if > threshold% coverage
        
        l = kappai / h_ave(i) / m_ave(i);
        if l > 10^(max(levs)) %cold shelves
            shelf_type(i) = 1;
            ptcol = cold_col;
            t(i) = text(h_ave(i)-81,m_ave(i),shelf, 'FontSize', fs, 'FontName', 'Arial');
        elseif l < 10.^(min(levs))
            shelf_type(i) = 2;
            ptcol = warm_col;
            
            t(i) = text(ax(1),h_ave(i)+11,m_ave(i),shelf, 'FontSize', fs, 'FontName', 'Arial');
        else 
             shelf_type(i) = 3;
             t(i) = text(ax(1),h_ave(i)+11,m_ave(i),shelf, 'FontSize', fs, 'FontName', 'Arial');
             ptcol =grey_col;

        end

        plot(ax(1), h_ave(i),m_ave(i), 'o', 'MarkerFaceColor', ptcol, 'MarkerEdgeColor', ptcol, 'MarkerSize', 6.5);


    end

    % add boxes around
%     mstd(i) = std(m_iceshelf);
%     hstd(i) = std(h_iceshelf);
%     xf = [hmean(i) - hstd(i), hmean(i) + hstd(i),  hmean(i) + hstd(i), hmean(i) - hstd(i)];
%     yf = [mmean(i) - mstd(i), mmean(i) - mstd(i),  mmean(i) + mstd(i), mmean(i) + mstd(i)];
%     yf = max(yf,0.1);
%     fill(xf, yf, 0.5*[1,1,1], 'FaceAlpha',0.3, 'EdgeAlpha',0.7);
     
end %end loop over shelves


%% Make (b)

ell_ave = kappai./ h_ave ./m_ave; %mean of the lengthscales
ell_ave(shelf_type == 0)= nan; %so we only plot those points in (a)

ms = 30; %markersize
% scatter3(ax(2), ell_ave(shelf_type == 1), epsxx_ave(shelf_type == 1), thinrate_ave(shelf_type == 1),ms, repmat(cold_col, sum(shelf_type == 1), 1), 'filled');
% scatter3(ax(2), ell_ave(shelf_type == 1), epsxx_ave(shelf_type == 1), thinrate_ave(shelf_type == 1),ms, repmat([0,0,0], sum(shelf_type == 1), 1)); %add the outside in black

% scatter3(ax(2),ell_ave(shelf_type == 2), epsxx_ave(shelf_type == 2), thinrate_ave(shelf_type == 2), ms,repmat(warm_col, sum(shelf_type == 2), 1), 'filled');
% scatter3(ax(2), ell_ave(shelf_type == 2), epsxx_ave(shelf_type == 2), thinrate_ave(shelf_type == 2),ms, repmat([0,0,0], sum(shelf_type == 2), 1)); %add the outside in black


ax(2).XLabel.String = '$\ell$'; ax(2).XLabel.Interpreter = 'latex';
ax(2).YLabel.String= 'strain rate (1/yr)';
ax(2).ZLabel.String = 'inverse thinning rate (yr/m)';
set(ax(2), 'XScale', 'log')
set(ax(2), 'YScale', 'log')
set(ax(2), 'ZScale', 'log')
ax(2).View = [45,36];
grid(ax(2), 'on')
ax(2).XLim = [0.0058    1.0];
ax(2).YLim = [8e-5    0.0140];
ax(2).ZLim = 1./flip([ 0.0361   10.0000]);
ax(2).ZLim = [0.15,40];
 hold(ax(2), 'on');
 box(ax(2), 'on')
%% project onto axes
% axinst = axes();
% axinst.Position = [0.5, 0.2, 0.2, 0.2];
% hold on
% box on% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% axinst.XLabel.String = 'strain rate (1/yr)';
% axinst.YLabel.String = 'thinning rate (m/yr)';


%project onto the l axis
scatter3(ax(2),min(ax(2).XLim)*ones(1,sum(shelf_type==1)), epsxx_ave(shelf_type == 1),1./thinrate_ave(shelf_type == 1),ms, repmat(cold_col, sum(shelf_type == 1), 1), 'Filled');
scatter3(ax(2),min(ax(2).XLim)*ones(1,sum(shelf_type==2)), epsxx_ave(shelf_type == 2),1./ thinrate_ave(shelf_type == 2), ms,repmat(warm_col, sum(shelf_type == 2), 1), 'Filled');
scatter3(ax(2),min(ax(2).XLim)*ones(1,sum(shelf_type==3)), epsxx_ave(shelf_type == 3),1./ thinrate_ave(shelf_type == 3), ms,repmat(grey_col, sum(shelf_type == 3), 1), 'Filled');

% project onto the epsxx axis
scatter3(ax(2), ell_ave(shelf_type == 1), max(ax(2).YLim)*ones(1,sum(shelf_type==1)),1./thinrate_ave(shelf_type == 1),ms, repmat(cold_col, sum(shelf_type == 1), 1), 'Filled');
scatter3(ax(2), ell_ave(shelf_type == 2), max(ax(2).YLim)*ones(1,sum(shelf_type==2)),1./ thinrate_ave(shelf_type == 2), ms,repmat(warm_col, sum(shelf_type == 2), 1), 'Filled');
scatter3(ax(2), ell_ave(shelf_type == 3), max(ax(2).YLim)*ones(1,sum(shelf_type==3)),1./ thinrate_ave(shelf_type == 3), ms,repmat(grey_col, sum(shelf_type == 3), 1), 'Filled');

%project onto the thinning rate axis
scatter3(ax(2), ell_ave(shelf_type == 1),  epsxx_ave(shelf_type == 1),min(ax(2).ZLim)*ones(1,sum(shelf_type==1)),ms, repmat(cold_col, sum(shelf_type == 1), 1), 'Filled');
scatter3(ax(2), ell_ave(shelf_type == 2),  epsxx_ave(shelf_type == 2), min(ax(2).ZLim)*ones(1,sum(shelf_type==2)), ms,repmat(warm_col, sum(shelf_type == 2), 1), 'Filled');
scatter3(ax(2), ell_ave(shelf_type == 3),  epsxx_ave(shelf_type == 3), min(ax(2).ZLim)*ones(1,sum(shelf_type==3)), ms,repmat(grey_col, sum(shelf_type == 3), 1), 'Filled');



%% final tidying
for i = 1:2
    ax(i).FontSize = fs;
    ax(i).XLabel.FontSize = fs+2;
    ax(i).YLabel.FontSize = fs+2;
    ax(i).FontName = 'Arial';
end

    ax(2).ZLabel.FontSize = fs+2;