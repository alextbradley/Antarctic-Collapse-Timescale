% Make figure 4 of the ms, showing timescale as a function of dM

% ATB (aleey@bas.ac.uk), 21/02/23. MIT licence.

%%
addpath('../functions');
load('fig2_out_.mat');
step = 10;        %grid resolution in km

shelf_names = ["Abbot","Amery","Borchgrevink","Brunt","Cook","Cosgrove","Crosson","Dotson","Filchner","FimbulJelbart",...
    "George6","Getz","KingBaudoin","Larsen",  "Nansen","PineIsland","PineIslandFast","PopeSmithKohler","RiiserLarsen", "Ronne","Ross",...
    "Shackleton", "Thwaites_","TottenMoscow", "West","Wilkins"];

data_folder = strcat('../gendata/figure4_1/step_', num2str(step));
shelf_info_folder = "../data/ice-shelves/all-shelves/"; %where to find co-ordinate files

% Loop over the shelves. For each, make a an array with the dM and the mean collapse time
% also plot the raw data
ct_shelves = cell(1,length(shelf_names));
dM_shelves = cell(1,length(shelf_names));
Pvals      = nan(2,length(shelf_names));


clf; hold on;
fig2data = load('fig2_out_.mat');
for is = 1:length(shelf_names)
     idx = find(shelf_names(is) == fig2data.shelf_names);

    if shelf_type(idx) == 1 %say == 1for cold shelves only , otherwise > 0 for both cold and warm shelves
        clf; hold on
        title(shelf_names(is))

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


        xf = dMs(dMs<10);
        yf = cts(dMs<10);
        plot(xf , log(yf),'o-','linewidth', 2, 'color', cc)

        % linear regression
        p = polyfit(xf(xf <1),log(yf(xf < 1)) ,1);

        plot(xf, p(2) + p(1)*xf, '--','linewidth', 2, 'color', cc)


        ax = gca;
        % set(ax, 'YScale', 'log');
        ax.XLim = [0,1];
        ax.YLabel.String = 'log(Collapse time) (log(yrs))';
        ax.XLabel.String = '$\Delta \dot{m}$';
        ax.XLabel.Interpreter = 'latex';
        ax.FontSize = 15;

        drawnow
        pause
    end
end


