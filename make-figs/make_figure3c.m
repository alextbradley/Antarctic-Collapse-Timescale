% have:
% lk : lengthscales of keep shelves
% shelf_namesk : names keen
% shelf_countsk: arrays of points for these shelves
% bar_coldatak : colours of points to keep
% shelf_typek : shelf types (1: cold, 2: warm);
% ct_avek : mean of distributions


shelf_namesk = shelf_names(shelf_type~=0);

% sort the data in increasing l
[~,I] = sort(ct_avek);
lks = lk(I);
shelf_namesks = shelf_namesk(I);
shelf_countsks = shelf_countsk(I);
bar_coldataks = bar_coldatak(I);
shelf_typeks = shelf_typek(I);
ct_aveks = ct_avek(I); 


%

clf; hold on; box on



w = 0.4; %wifth = the dist
for i = 1:length(shelf_typeks)

    %get the data as array
    counts = cell2mat(shelf_countsks(i));
    counts(counts < 10) = 10; %make the plotting a bit nicer

    %remove very long timescale pts
    %counts = counts(counts < 5*1e3);

    %fit a kde to it
    kde = fitdist(counts,'kernel');

    %evaluate it
    x = logspace(1,4);
    y = pdf(kde,x);

    %scale the pdf
    y = y * w / max(y); 

    % create fill array
    cline = i; %centreline of distribution
    xf = [cline + y, flip(cline - y)];
    yf = [x,flip(x)];

    %fill the data
    fill(xf, yf, bar_coldatak(i,:), 'linewidth', 1)

    % add mean as a point
    plot(cline,mean(kde),'ko', 'markerfacecolor', 'k', 'markersize', 5)
    %or maybe the point of 90 mass or so?
%     [countsort,~] = sort(counts);
%     frac = 0.75;
%     idx = round(length(countsort)*frac); %index frac way along counts
%     plot(cline, countsort(idx),'ko', 'markerfacecolor', 'k', 'markersize', 5)
    


end

xlim([0,length(shelf_namesk)+1]);
xticks(1:length(shelf_countsks));
xticklabels(shelf_namesks);
shg
set(gca, 'YScale', 'log')

fig = gcf;
fig.Position(3:4) = [1200, 400];