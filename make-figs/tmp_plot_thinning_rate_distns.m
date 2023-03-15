
%% Preliminaries
%clear
addpath('../functions');
%f = load('../data/ice_sheet_data.mat');
figure(3); clf; hold on;


%
%%%%%% add shelf points %%%%%%%
%

% init storage
jdir = dir('../data/ice-shelves/all-shelves/*.mat'); %where shelf files are stored

%we'll also work out shelf by shelf collapse data here
ff = load('figure3-data.mat');
ct = ff.collapse_time; %collapse time data
ct_Nye = ff.collapse_time_Nye;
tags = ff.tags;

%initialize storage
m_ave = zeros(1,length(jdir));
h_ave = zeros(1,length(jdir));
mstd = zeros(1,length(jdir));
hstd = zeros(1,length(jdir));
percov = zeros(1,length(jdir)); %store the % coverage of input data
collapse_coverage  =  zeros(1,length(jdir)); %per of points with non-nan output
shelf_type = zeros(1,length(jdir)); %flags for keeping the shelf or not (0), cold (1), warm (2)
epsxx_ave = zeros(1,length(jdir));
thinrate_ave = zeros(1,length(jdir));
shelf_counts = cell(1,length(jdir)); %for storing individual shelf counts
ct_ave = zeros(1,length(jdir)); %mean values of the collapse time for LEFM
ct_ave_Nye = zeros(1,length(jdir)); %mean values of the collapse time for Nye
ll = zeros(1,length(jdir)); %store the lengthscales l
shelf_cols = zeros(length(jdir), 3);%store plot colours
area_shelf = zeros(1,length(jdir)); %store shelf areas

x = linspace(-2,10,2e3); 
w = 0.3;


% loop over shelves
shelf_names = strings;
for i = 1:length(jdir)
    shelf = jdir(i).name;
    shelf = strrep(shelf,'.mat',''); %strip the .mat at the end
    shelf_names(i) = shelf;
    
    fname = strcat('../data/ice-shelves/all-shelves/' ,shelf, '.mat');
    g = load(fname);

    % restrict data to those points in shelf with high enough melt rate
    m_shelf = f.m; m_shelf = m_shelf(g.IN); %only the points in this shelf
    h_shelf = f.H; h_shelf = h_shelf(g.IN);
    epsxx_shelf = f.eflow; epsxx_shelf = epsxx_shelf(g.IN);
    err =  f.dMdterr;
    err(err > 10) = 10;
    
    
    dhdt_shelf = f.dMdtadj; 

    err = err(g.IN);
    dhdt_shelf = dhdt_shelf(g.IN);
    

    %dhdt_shelf = dhdt_shelf - err;
    relerr = err./dhdt_shelf;
  % dhdt_shelf(relerr > 1) = nan;


    idx =  (~isnan(h_shelf)) &  (m_shelf > 1e-6)  & (epsxx_shelf > 1e-6) &  (-dhdt_shelf > 1e-6) ; %points where we have point thickness,  melt rate > 0, strain > 0, thinning < 0
    %idx =  (~isnan(h_iceshelf)) &  (~isnan(m_iceshelf)); %points where we have point thickness and melt rate > 0
    idx_for_calc =  (~isnan(h_shelf)) &  (~isnan(m_shelf))  & (~isnan(epsxx_shelf)) &  (~isnan(-dhdt_shelf)); %points where we have data for 

    idx =  (~isnan(h_shelf))  &  (~isnan(-dhdt_shelf)) & (~isnan(m_shelf)) & (~isnan(epsxx_shelf))  ; %points where we have point thickness,  melt rate > 0, strain > 0, thinning < 0


    %compute the shelf area as a sanity check?
    area_shelf(i) = length(h_shelf(~isnan(h_shelf))); %size in km^3 (grid size = 1e3 * 1e3)

    %get data points which satisfy the idx criterion
    h_shelf = h_shelf(idx);
    m_shelf = m_shelf(idx); %arrays with points in particular shelf with both melt and thicknes
    epsxx_shelf = epsxx_shelf(idx);
    thinrate_shelf = -dhdt_shelf(idx);


    %compute averages of strain, thickness, melt, thining rate
%     
%     pd = fitdist(m_shelf,'kernel');  m_ave(i) = mean(pd);
%     pd = fitdist(h_shelf,'kernel');  h_ave(i) = mean(pd);
%     pd = fitdist(epsxx_shelf,'kernel');  epsxx_ave(i) = mean(pd);
%     pd = fitdist(thinrate_shelf,'kernel');   thinrate_ave(i) = mean(pd);
    
    m_ave(i) = median((m_shelf));
    h_ave(i) = median((h_shelf));
    epsxx_ave(i) = median(epsxx_shelf);
    thinrate_ave(i) = median(thinrate_shelf);


    if any(idx)
    %fit a kde to it
    kde = fitdist(thinrate_shelf,'kernel');

    %evaluate it
    y = pdf(kde,x);

    %compute mean
    %meanval = median(kde);%+std(kde)/2;
    meanval = mean(kde)+std(kde); %75th percentile
    mean_dhdt = -meanval;
    %meanval = mean(thinrate_shelf);

    %scale the pdf
    y = y * w / max(y); 

    % create fill array
    cline = i; %centreline of distribution
    xf = [cline + y, flip(cline - y)];
    yf = [x,flip(x)];

    %fill the data
    if meanval >0 
        fill(xf, yf, 'r', 'linewidth', 1, 'EdgeColor', 0.2*[1,1,1], 'FaceAlpha', 0.3)
    else
        fill(xf, yf, 'b', 'linewidth', 1, 'EdgeColor', 0.2*[1,1,1], 'FaceAlpha', 0.3)
    end


    % add mean as a point
    plot(cline,meanval,'ko', 'markerfacecolor', 0*[1,1,1], 'markersize', 6,'LineWidth',1.1 )
    else
        mean_dhdt = nan; %
    end
    save(fname, 'mean_dhdt', '-append')
    
      
end %end loop over shelves

ax = gca;
ax.XTick = 1:length(shelf_names);
ax.XTickLabels = shelf_names;

%ylim([-1,1])