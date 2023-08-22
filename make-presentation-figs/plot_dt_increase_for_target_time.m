% Make figure 4 of the ms, showing timescale as a function of dM

% ATB (aleey@bas.ac.uk), 21/02/23. MIT licence.

%%
addpath('../functions');
f = load('../data/ice_sheet_data.mat', 'H', 'm'); %get the ice sheet thickness data (use this for shelf mask)
load('fig2_out_.mat');
step = 10;        %grid resolution in km

shelf_names = ["Abbot","Amery","Borchgrevink","Brunt","Cook","Cosgrove","Crosson","Dotson","Filchner","FimbulJelbart",...
    "George6","Getz","KingBaudoin","Larsen",  "Nansen","PineIsland","PineIslandFast","PopeSmithKohler","RiiserLarsen", "Ronne","Ross",...
    "Shackleton", "Thwaites_","TottenMoscow", "West","Wilkins"];

data_folder = strcat('../gendata/figure4_1/step_', num2str(step));
shelf_info_folder = "../data/ice-shelves/all-shelves/"; %where to find co-ordinate files

%% Loop over the shelves. For each, make a an array with the dM and the mean collapse time
% also plot the raw data
ct_shelves = cell(1,length(shelf_names));
dM_shelves = cell(1,length(shelf_names));
Pvals      = nan(2,length(shelf_names));
figure(1); clf; hold on; box on


ax(1) =    subplot(1,2,1) ;hold on; box on
ax(2) =    subplot(1,2,2) ;hold on; box on

fig2data = load('fig2_out_.mat'); 
for is = 1:length(shelf_names)
    fpath = strcat(data_folder, '/', shelf_names(is));
    jdir = dir(strcat(fpath, '/*.mat')); %each of these is a different value of dM

    dMs = nan(1,length(jdir));
    cts = nan(1,length(jdir));
    for im  = 1:length(jdir)
        ff = load(strcat(fpath, '/', jdir(im).name));
        
        dMs(im) = ff.dM;
        cts(im) = ff.collapse_time;
    end
    %sort into increaseing dm order
    [dMs,I] = sort(dMs);
    cts = cts(I);
    
    ct_shelves{is} = cts;
    dM_shelves{is} = dMs;

    %get the colour from figure 2 data
    idx = find(shelf_names(is) == fig2data.shelf_names);
    cc = fig2data.shelf_cols(idx,:);

    if shelf_type(idx) > 0
        xf = dMs(dMs<10);
        yf = cts(dMs<10);
        plot(ax(shelf_type(idx)),xf , yf,'o-','linewidth', 2, 'color', cc)

    end
end

for i =1:2 
    set(ax(i), 'YScale', 'log'); 
    ax(i).XLim = [0,2];
    ax(i).YLabel.String = 'Collapse time (yrs)';
    ax(i).XLabel.String = '$\Delta \dot{m}$';
    ax(i).XLabel.Interpreter = 'latex';
    ax(i).FontSize = 15;
    
end


%% Repeat for scaled data
%% Loop over the shelves. For each, make a an array with the dM and the mean collapse time
ct_shelves = cell(1,length(shelf_names));
dM_shelves = cell(1,length(shelf_names));
Pvals      = nan(2,length(shelf_names));
figure(2); clf; hold on; box on


ax(1) =    subplot(1,2,1) ;hold on; box on
ax(2) =    subplot(1,2,2) ;hold on; box on
  

fig2data = load('fig2_out_.mat'); 
for is = 1:length(shelf_names)
    fpath = strcat(data_folder, '/', shelf_names(is));
    jdir = dir(strcat(fpath, '/*.mat')); %each of these is a different value of dM

    dMs = nan(1,length(jdir));
    cts = nan(1,length(jdir));
    for im  = 1:length(jdir)
        ff = load(strcat(fpath, '/', jdir(im).name));
        
        dMs(im) = ff.dM;
        cts(im) = ff.collapse_time;
    end
    %sort into increaseing dm order
    [dMs,I] = sort(dMs);
    cts = cts(I);
    
    ct_shelves{is} = cts;
    dM_shelves{is} = dMs;

    %get the colour from figure 2 data
    idx = find(shelf_names(is) == fig2data.shelf_names);
    cc = fig2data.shelf_cols(idx,:);

    if shelf_type(idx) > 0
        xf = dMs(dMs<10);
        yf = cts(dMs<10);
        yy = (abs(log(yf) - log(yf(1))));  %if y ~ y0 exp (-a x^b), then yy = -(log(y) - log(y0)) ~ a x^b
        %and therefore log(yy)~log(a) + blog(x), i.e. it is linear in log
        %x.

        % if y ~ y0 x^b, then log (y) - log(y0) ~ b log(x), and a plot of
        % yy again log(x) would go thru the origin

        ii = yy>0;
        xxf = log(xf(ii));
        yyy = log(yy(ii));

        plot(ax(shelf_type(idx)),xxf , yyy,'o-','linewidth', 2, 'color', cc)

        % linear fit to the data
        P = polyfit(xxf,yyy,1);
        Pvals( :,is) = P;
         
        plot(ax(shelf_type(idx)),xxf , P(1)*xxf + P(2),'--','linewidth', 1.5, 'color', cc)
% 
  %  t = text(ax(shelf_type(idx)), xxf(1),yyy(1), shelf_names(is), 'FontSize', 14);
    end


end

for i = 1:2
%set(ax(i), 'YScale', 'log')
%set(ax(i), 'XScale', 'log')
%ax(i).XLim = [0,10];
%ax(i).XTick = 10.^[0, 0.5, 1];
%set(gca, 'YAxisLocation', 'right')

end
%ax(1).YLim = [0.5*1e2, 2*1e3];


%set(gca, 'XScale', 'log')

for i = 1:2
    ax(i).XLabel.String = 'log(melt increase)';
    ax(i).YLabel.String = 'log(\tau / \tau_0)';
    ax(i).FontSize = 14;
end

%% Plot the powers
figure(3); clf;
histogram(Pvals(1,shelf_type> 0), 10);
xlim([0,2])
xlabel('slope ');
ylabel('count')
ax = gca; ax.FontSize = 14;