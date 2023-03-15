% Make figure 4 of the ms, showing timescale as a function of dM

% ATB (aleey@bas.ac.uk), 21/02/23. MIT licence.

%%
addpath('../functions');
f = load('../data/ice_sheet_data.mat', 'H', 'm', 'dMdtadj'); %get the ice sheet thickness data (use this for shelf mask)
load('fig2_out_.mat');
step = 10;        %grid resolution in km

shelf_names = ["Abbot","Amery","Borchgrevink","Brunt","Cook","Cosgrove","Crosson","Dotson","Filchner","FimbulJelbart",...
    "George6","Getz","KingBaudoin","Larsen",  "Nansen","PineIsland","PineIslandFast","PopeSmithKohler","RiiserLarsen", "Ronne","Ross",...
    "Shackleton", "Thwaites_","TottenMoscow", "West","Wilkins"];

data_folder = strcat('../gendata/figure4_1/step_', num2str(step));
shelf_info_folder = "../data/ice-shelves/all-shelves/"; %where to find co-ordinate files


%% Plot the collapse time and prediction from upper bound
ct_shelves = cell(1,length(shelf_names));
dM_shelves = cell(1,length(shelf_names));
Pvals      = nan(2,length(shelf_names));


fig2data = load('fig2_out_.mat');
for is = 1:length(shelf_names)
    idx = find(shelf_names(is) == fig2data.shelf_names);
    if shelf_type(idx) > 0

        figure(2); clf; hold on; box on
        title(shelf_names(is))

        % load the data

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
        cc = fig2data.shelf_cols(idx,:);

        % restrict
        dM = dMs(dMs<=1);
        cts = cts(dMs<=1);

        plot(dM , cts,'o-','linewidth', 2.5, 'color', cc)

        % linear fit to the data
        P = polyfit(dM,log(cts),1);
        Pvals( :,is) = P;

        plot(dM , exp(P(1)*dM  + P(2)),'--','linewidth', 1.5, 'color',cc)
        %

        % get the upper bound fit
        alpha_hdot = -0.0836;
        fname = strcat(shelf_info_folder, shelf_names(is), '.mat');
        shelf_data = load(fname);
        idx = shelf_data.IN & ~isnan(f.H) & ~isnan(f.dMdtadj) ;
        dhdtbar = mean(mean(-abs(f.dMdtadj(idx))));

        plot(dM, cts(1)*exp(-dM*alpha_hdot/dhdtbar), 'k--', 'linewidth', 1.5);


        ax = gca; ax.FontSize = 14;
        ax.FontName = 'GillSans';
        ax.XLabel.String = '$\Delta \dot{m}$'; ax.XLabel.Interpreter = 'latex';
        ax.YLabel.String = 'collapse time';
        legend('raw data', 'linear fit', 'upper bound');

        set(gcf, 'color', 'w')
        xlim([0,1]);
        fig = gcf; fig.Position(3:4) = [450, 330];

        set(gca, 'YScale', 'log');
        drawnow; shg
        pause
    end


end


%% Make a plot showing (a) the collapse time for cold water shelves as a function of M and that data rescaled according to exponential fit

%% Plot the collapse time and prediction from upper bound
ct_shelves = cell(1,length(shelf_names));
dM_shelves = cell(1,length(shelf_names));
Pvals      = nan(2,length(shelf_names));

figure(1); clf; 

fig2data = load('fig2_out_.mat');
for is = 1:length(shelf_names)
    idx = find(shelf_names(is) == fig2data.shelf_names);
    if shelf_type(idx) ==  1

  

        % load the data

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
        cc = fig2data.shelf_cols(idx,:);

        % restrict
        dMsc = dMs(dMs<=1);
        ctsc = cts(dMs<=1);
        
        %raw data
        ax(1) = subplot(1,2,1); hold on; box on
        plot(dMs , cts,'linewidth', 2, 'color', cc)



        % linear fit to the data
        P = polyfit(dMc,log(ctsc),1);
        Pvals( :,is) = P;

        ax(2) = subplot(1,2,2); hold on; box on
        plot(dMs , smooth(cts),'linewidth', 2, 'color', cc)
        plot(dMs , exp(P(1)*dMs  + P(2)),'--','linewidth', 1.5, 'color',cc)
        %

    end


end


for i = 1:2
ax(i).FontSize = 16;
ax(i).FontName = 'GillSans';
ax(i).XLabel.String = '$\Delta \dot{m}$ (m/yr)'; ax(i).XLabel.Interpreter = 'latex';


end

ax(1).YLabel.String = 'collapse time (years)';
ax(2).YLabel.String = 'collapse time (years)';
%ax(2).YLim = [0,2];
ax(2).XLim = [0,1];
ax(1).XLim = [0,10];
set(ax(2), 'YScale', 'log');
set(gcf, 'color', 'w')

