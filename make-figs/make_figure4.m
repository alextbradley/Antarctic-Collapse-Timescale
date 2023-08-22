% Make figure 4 of the ms, showing collapse timescale as a function of dM
% on linear (a) and semilog (b) axes for cold ice shelves

% ATB (aleey@bas.ac.uk), 21/02/23. MIT licence.

%% Preliminaries
addpath('../functions');
%f = load('../data/ice_sheet_data.mat', 'H', 'm'); %get the ice sheet thickness data (use this for shelf mask)
fig2data = load('fig2_out_.mat');
step = 10;        %grid resolution in km

shelf_names = ["Abbot","Amery","Borchgrevink","Brunt","Cook","Cosgrove","Crosson","Dotson","Filchner","FimbulJelbart",...
    "George6","Getz","KingBaudoin","Larsen",  "Nansen","PineIsland","PineIslandFast","PopeSmithKohler","RiiserLarsen", "Ronne","Ross",...
    "Shackleton", "Thwaites","TottenMoscow", "West","Wilkins"];

data_folder = strcat('../gendata/figure4_1/step_', num2str(step));
shelf_info_folder = "../data/ice-shelves/all-shelves/"; %where to find co-ordinate files

%% Loop over the shelves and store the data
ct_shelves = cell(1,length(shelf_names));
dM_shelves = cell(1,length(shelf_names));
Pvals      = nan(2,length(shelf_names));


for is = 1:length(shelf_names)

    %get the path to data
    fpath = strcat(data_folder, '/', shelf_names(is));
    jdir = dir(strcat(fpath, '/*.mat')); %each of these is a different value of dM

    %loop over files in directory, each one is a different dM value
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

    %store the results
    ct_shelves{is} = cts;
    dM_shelves{is} = dMs;

end

%% Setup the plot
figure(1); clf; hold on; box on
fig= gcf; fig.Position(3:4) = [1300,500];

w = 0.35;
gapx = 0.1;
gapy = 0.05;
h = 0.4;
startx = (1 - (2*w + gapx))/2;
starty = 0.97 - (2*h + gapy);
positions = [startx,            starty,             w, h;
    startx,            starty + h + gapy,  w, h;
    startx + gapx + w, starty,             w, 2*h + gapy];
for i = 1:3
    ax(i) =    subplot('Position', positions(i,:)) ;
    hold(ax(i), 'on');
    box(ax(i), 'on');
end


%% Make plot
%cmap = lines(sum(fig2data.shelf_type == 1)); count = 1;
splittol = 1250; %where to chop the plots
for is =  1:length(shelf_names)

    %get the colour from figure 2 data
    idx = find(shelf_names(is) == fig2data.shelf_names);
    cc = fig2data.shelf_cols(idx,:);


    %get the data from above
    dMs = cell2mat(dM_shelves(is));
    cts = cell2mat(ct_shelves(is));

    if shelf_type(idx) == 1 %cold shelf
        if (~strcmp(shelf_names(is), "Wilkins"))  && (~strcmp(shelf_names(is), "Nansen")) %remove there
            xf = dMs(dMs<=10);
            yf = cts(dMs<=10);
           
            if yf(1) > splittol
                plot(ax(2),xf , yf,'-','linewidth', 2, 'color', cc)
                text(ax(2), xf(1), yf(1), shelf_names(is), 'FontSize', 13);
            else
                plot(ax(1),xf , yf,'-','linewidth', 2, 'color', cc)
                text(ax(1), xf(1), yf(1), shelf_names(is) , 'FontSize', 13);
            end


            xxf = xf(xf <= 1);
            yyf =  yf(xf <= 1)./yf(1);
            xxf = smooth(xxf);
            yyf = smooth(yyf);
            plot(ax(3),xxf, yyf,'-','linewidth', 2, 'color', cc)
            tt = text(ax(3), xxf(end), yyf(end), shelf_names(is),'FontSize', 13);
        end
    end

end


%%
for i =1:3
    %set(ax(i), 'YScale', 'log');
    %ax(i).XLim = [0,2];

    ax(i).XLabel.String = '$\Delta \dot{m}$ (m/yr)';
    ax(i).XLabel.Interpreter = 'latex';
    ax(i).YLabel.Interpreter = 'latex';
    ax(i).FontSize = 15;
    grid(ax(i), 'on')

end
ax(1).YLim = [0, splittol];
ax(2).YLim = [0, 3500];
 ax(2).XLabel.String = '';
 ax(2).XTick = ax(1).XTick;;
 ax(2).XTickLabel = {};
ax(3).YLabel.String = '$\tau(\Delta \dot{m})/ \tau(\Delta \dot{m} = 0)$';
ax(1).YLabel.String = 'collapse timescale, $\tau(\Delta \dot{m})$';
ax(2).YLabel.String = 'collapse timescale, $\tau(\Delta \dot{m})$';
set(ax(3), 'YScale', 'log');
ax(3).XLim = [0,1];
ax(3).YLim = [0.25,1];
%ax(1).YLim = [0,1e3];