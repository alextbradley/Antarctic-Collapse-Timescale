% Make supplementary figure A of the manuscript, showing violin plots of
% the density of thinning rate points.
%
% ATB (aleey@bas.ac.uk), 21/04/23, MIT licence.
%% Preliminaries
%clear
addpath('../functions');
%f = load('../data/ice_sheet_data.mat');
jdir = dir('../data/ice-shelves/all-shelves/*.mat'); %where shelf files are stored

%% Initialize plot
figure(1); clf; hold on;
x = linspace(-2,10,2e3); %dhdt points to sample in the distributions
w = 0.3;                 % relative width of the distributions
fillcol = [250,208,56]/255;  %filling colour


%% Loop over shelves
shelf_names = strings;
shelf_names_keep = strings;
count = 1;
for i = 1:length(jdir)
    shelf = jdir(i).name;
    shelf = strrep(shelf,'.mat',''); %strip the .mat at the end
    shelf_names(i) = shelf;

    fname = strcat('../data/ice-shelves/all-shelves/' ,shelf, '.mat');
    g = load(fname);

    % get the thickness, melt and strain rate so that we can consider only
    % points we would keep
    m_shelf = f.m; m_shelf = m_shelf(g.IN); %only the points in this shelf
    h_shelf = f.H; h_shelf = h_shelf(g.IN);
    epsxx_shelf = f.eflow; epsxx_shelf = epsxx_shelf(g.IN);
    dhdt_shelf = f.dMdtadj;
    dhdt_shelf = dhdt_shelf(g.IN);
    idx =  (~isnan(h_shelf))  &  (~isnan(-dhdt_shelf)) & (~isnan(m_shelf)) & (~isnan(epsxx_shelf))  ; %points where we have point thickness,  melt rate > 0, strain > 0, thinning < 0
    thinrate_shelf = -dhdt_shelf(idx);
    thinrate_ave(i) = median(thinrate_shelf);

    if any(idx) %if we have any points
        %fit a kde to it
        kde = fitdist(thinrate_shelf,'kernel');

        %evaluate it
        y = pdf(kde,x);

        %compute mean
        %meanval = median(kde);%+std(kde)/2;
        meanval = mean(kde)+std(kde); %1 std above mean
        mean_dhdt = -meanval;

        %scale the pdf
        y = y * w / max(y);

        % create fill array
        cline = count; %centreline of distribution
        xf = [cline + y, flip(cline - y)];
        yf = [x,flip(x)];

        %fill the data
        fill(xf, yf, fillcol, 'linewidth', 1, 'EdgeColor', 0.2*[1,1,1], 'FaceAlpha', 0.8)

        % add mean as a point
        plot(cline,meanval,'ko', 'markerfacecolor', 0*[1,1,1], 'markersize', 3,'LineWidth',1.1 )
        shelf_names_keep(count) = shelf;

        count = count + 1;
    else
        mean_dhdt = nan; %
    end

   % save(fname, 'mean_dhdt', '-append') %uncomment to save the mean point

end %end loop over shelves

%% tidy
ax = gca;
ax.XTick = 1:length(shelf_names_keep);
ax.XTickLabels = shelf_names;
ax.XLim = [0, length(shelf_names_keep)+1];

ax.YLim = [-2,10];
%ylim([-1,1])
box(ax, 'on');
ax.FontSize = 15;
ax.XTickLabels = shelf_names_keep;
ax.YLabel.String = 'thinning rate (m/yr)';