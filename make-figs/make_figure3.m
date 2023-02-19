% Make figure 3 of the ms, showing the collapse time in (a) using the raw
% data and (b) using the absolute value of the thinning rate.
%
% ATB (aleey@bas.ac.uk), 17/02/2023. MIT Licence.
%
% (figure made as two individual figures)
%% Preliminaries
%clear
addpath('../functions');
fs = 14; %fontsize
cmap = cmocean('matter', 100);

%% Part (a)
figure(1);clf
ff = load('figure3-data.mat');
ct = ff.collapse_time; %collapse time data
tags = ff.tags;

% 1: collapse time
ax = gca;
pl = imagesc(ax, log10((ct))); 
set(pl, 'AlphaData', ~isnan(ct));
set(ax, "YDir", 'reverse')
axis equal
ax.Visible = 'off';
%pause

% 2: missing data
axlayer2 = axes();
tags_missing = nan(size(tags));
tags_missing(tags == 2) = 1;
pl = imagesc(tags_missing);
set(pl, 'AlphaData', ~isnan(tags_missing)); 
colormap(axlayer2, 0.6* [1,1,1])
axlayer2.Visible = 'off';
set(axlayer2, "YDir", 'reverse')

axis equal
axlayer2.Position = ax.Position;
%pause

% 3: negative strain rate
axlayer3 = axes();
tags_negative_strain =  nan(size(tags));
tags_negative_strain(tags == 3) = 1;
pl = imagesc( tags_negative_strain);
set(pl, 'AlphaData', ~isnan(tags_negative_strain)); 
colormap(axlayer3, cmap(end,:))
%colormap(ax(3),  0.3* [1,1,1])
axlayer3.Visible = 'off';
set(axlayer3, "YDir", 'reverse')
axis equal
axlayer3.Position = ax.Position;

% 4: negative dhdt
axlayer4 = axes();
tags_negative_strain =  nan(size(tags));
tags_negative_strain(tags == 4) = 1;
pl = imagesc( tags_negative_strain);
set(pl, 'AlphaData', ~isnan(tags_negative_strain)); 
colormap(axlayer4, cmap(end,:))
axlayer4.Visible = 'off';
set(axlayer4, "YDir", 'reverse')
axis equal
axlayer4.Position = ax.Position;

colormap(ax, cmap)
linkaxes([ax, axlayer2, axlayer3, axlayer4])
fig = gcf;
fig.Position(3:4) = [655,507];

%% Part (b) -- with approximation to thinning rate
figure(2);clf
clear ax axlayer2 axlayer3 axlayer4
ff = load('figure3-data.mat');
ct = ff.collapse_time; %collapse time data
tags = ff.tags;

% 1: collapse time
ax = gca;
pl = imagesc(ax, log10((ct))); 
c = colorbar(ax);
c.Limits = [1,4];
set(pl, 'AlphaData', ~isnan(ct));
set(ax, "YDir", 'reverse')
colormap(ax, cmap)
axis equal
ax.Visible = 'off';
%pause

% 2: missing data
axlayer2 = axes();
tags_missing = nan(size(tags));
tags_missing(tags == 2) = 1;
pl = imagesc(tags_missing);
set(pl, 'AlphaData', ~isnan(tags_missing)); 
colormap(axlayer2, 0.6* [1,1,1])
axlayer2.Visible = 'off';
set(axlayer2, "YDir", 'reverse')
ct(tags_missing == 1)=  nan;

axis equal
axlayer2.Position = ax.Position;
%pause

% 3: negative strain rate
axlayer3 = axes();
tags_negative_strain = nan(size(tags));
tags_negative_strain(tags == 3) = 1;
pl = imagesc( tags_negative_strain);
set(pl, 'AlphaData', ~isnan(tags_negative_strain)); 
colormap(axlayer3, cmap(end,:))
%colormap(ax(3),  0.3* [1,1,1])
axlayer3.Visible = 'off';
set(axlayer3, "YDir", 'reverse')
axis equal
axlayer3.Position = ax.Position;
ct(tags_negative_strain == 1)=  nan;

linkaxes([ax, axlayer2, axlayer3])


fig = gcf;
fig.Position(3:4) = [655,507];
c.Ticks = 1:4;
c.FontSize = fs;
c.TickLabels = {'10^1', '10^2', '10^3', '10^4'};

%% Make box plots of the shelves
% import the data from figure 2!

% jdir = dir('../data/ice-shelves/major/*.mat'); %where shelf files are stored
% warm_col = [203,0,63]/255; %warm shelves colour
% cold_col = [0, 63, 203]/255; %cold shelves colour
% bar_coldata = zeros(length(jdir), 3);
% 
% iskeep = zeros(1,length(jdir));
% 
% shelf_names = strings;
% 
% dhdtave = mean(mean(f.dMdtadj(~isnan(f.dMdtadj))));
% 
% for i = 1:length(jdir)
%     shelf = jdir(i).name;
%     shelf = strrep(shelf,'.mat',''); %strip the .mat at the end
%     shelf_names(i) = shelf;
%     
%     fname = strcat('../data/ice-shelves/major/' ,shelf, '.mat');
%     g = load(fname);
% 
%     % restrict to this shelf
%     ct_shelf = ct(g.IN); %only the points in this shelf
%     h = f.H;
%     hs = h(g.IN);
%     coverage(i) = sum(sum(~isnan(ct_shelf)))/sum(sum(~isnan(hs))) * 100;
% 
%     aa = (ct_shelf(~isnan(ct_shelf))); %non nan points
%     ct_ave(i) = median(aa(aa < 1e4));
% 
%     % compute the stuff to work out l
%     m_shelf = f.m; m_shelf = m_shelf(g.IN); %only the points in this shelf
%     h_shelf = f.H; h_shelf = h_shelf(g.IN);
%     epsxx_shelf = f.eflow; epsxx_shelf = epsxx_shelf(g.IN);
%     dhdt_shelf = f.dhdtadj; dhdt_shelf = dhdt_shelf(g.IN);
%     idx =  (~isnan(h_shelf)) &  (m_shelf > 1e-6)  & (epsxx_shelf > 1e-6) &  (-dhdt_shelf > 1e-6) ; %points where we have point thickness,  melt rate > 0, strain > 0, thinning < 0
%     idx_for_calc =  (~isnan(h_shelf)) &  (~isnan(m_shelf))  & (~isnan(epsxx_shelf)) &  (~isnan(-dhdt_shelf)); %points where we have data for 
%     area_shelf(i) = length(h_shelf(~isnan(h_shelf))); %size in km^3 (grid size = 1e3 * 1e3)
% 
% 
%     % plot point and add name
%     shelf_total = area_shelf(i); %total number of points in shelf
%     shelf_keep  = sum(idx_for_calc); %number of points with thickness data
%     percov(i) = shelf_keep/shelf_total * 100;
% 
%     h_shelf = h_shelf(idx);
%     m_shelf = m_shelf(idx); %arrays with points in particular shelf with both melt and thicknes
%     m_ave = median((m_shelf));
%     h_ave = median((h_shelf));
%     l(i) = kappai / h_ave / m_ave;
%     if (percov(i) > 50) || (shelf == "Thwaites")
%         iskeep(i) = 1;
%         levs = [-1.2, -1.5]; %%NB: check agrees w/ fig2!
%         if l(i) > 10^(max(levs)) %cold shelves
%             bar_coldata(i,:) = cold_col;
%         else
%             bar_coldata(i,:) = warm_col;
% 
%         end
%         %count = count +1;
%     end
% 
%     % compute the 'disappearance' timescale
%     dtc(i) = abs(h_ave / dhdtave);
% 
% end
% 
% %% make bar map
% figure(3); clf; 
% lk = l(iskeep==1);
% 
% 
% ct_avek = ct_ave(iskeep==1);
% [~,I] = sort(ct_avek);
% ct_avek = ct_avek(I); %put in order of increasing timesclae
% ct_avek(ct_avek<10) = 11; %just so it appears
% shelf_namesk = shelf_names(iskeep==1);
% shelf_namesk = shelf_namesk(I);
% 
% 
% bar_coldatak = bar_coldata(iskeep==1,:);
% bar_coldatak = bar_coldatak(I,:);
% 
% b = bar(ct_avek);
% b.FaceColor = 'flat';
% xticks(1:length(shelf_namesk))
% set(gca,'xticklabel',shelf_namesk)
% shg
% 
% b.CData = bar_coldatak;
% set(gca, 'YScale','linear')
% 
% %% repeat for the dimensionless timescale
% figure(4); clf; 
% lk = l(iskeep==1);
% 
% 
% 
% ct_avek = ct_ave(iskeep==1);
% [~,I] = sort(ct_avek);
% ct_avek = ct_avek(I); %put in order of increasing l
% ct_avek(ct_avek<10) = 11; %just so it appears
% 
% dtck = dtc(iskeep == 1);
% shelf_namesk = shelf_names(iskeep==1);
% shelf_namesk = shelf_namesk(I);
% 
% 
% bar_coldatak = bar_coldata(iskeep==1,:);
% bar_coldatak = bar_coldatak(I,:);
% 
% b = bar(ct_avek./dtck);
% b.FaceColor = 'flat';
% xticks(1:length(shelf_namesk))
% set(gca,'xticklabel',shelf_namesk)
% shg
% ylim([1e-3, 1e0])
% b.CData = bar_coldatak;
% set(gca, 'YScale','log')