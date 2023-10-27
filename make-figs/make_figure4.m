% Make figure 4 of the ms, showing collapse timescale as a function of dM
% on linear (a) and semilog (b) axes for cold ice shelves

% ATB (aleey@bas.ac.uk), 21/02/23. MIT licence.

%% Preliminaries
addpath('../functions');
%f = load('../data/ice_sheet_data.mat', 'H', 'm'); %get the ice sheet thickness data (use this for shelf mask)

step = 10;        %grid resolution in km

shelf_names = ["Abbot","Amery","Borchgrevink","Filchner","FimbulJelbart",...
    "KingBaudoin","Larsen","RiiserLarsen", "Ronne","Ross",...
    "Shackleton", "West"];

shelf_names_plot = ["ABB", "AME", "BOR", "FIL", "FIM", "BAU", "LAR", "RII",...
    "RON", "ROSS", "SHA", "WES"];

data_folder = strcat('../gendata/figure4/step_', num2str(step), '_CircumAntCoeffs');
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
fig= gcf; fig.Position(3:4) = [1150,800];

% setup panels (a)--(c)
w = 0.35;
gapx = 0.06;
gapy = 0.05;
h = 0.2;
startx = (1 - (2*w + gapx))/2;
starty = 0.97 - (2*h + gapy);
starty = 0.5;
positions = [startx,            starty,             w, h;
    startx,            starty + h + gapy,  w, h;
    startx + gapx + w, starty,             w, 2*h + gapy];
for i = 1:3
    ax(i) =    subplot('Position', positions(i,:)) ;
    hold(ax(i), 'on');
    box(ax(i), 'on');
end

% setup panels (d)--(g)
g2 = 0.06; %gap between these panels
w2 = (1 - 2*startx - 4*g2)/5; %width of these panels
h2 = 0.3; %height of these panels
starty2 = 0.1;
for i = 1:5
    positions(i+3,:) = [startx + (i-1)*g2 + (i-1)*w2, starty2, w2, h2];
    ax(i+3) =    subplot('Position', positions(i+3,:)) ;
    hold(ax(i+3), 'on');
    box(ax(i+3), 'on');
end




%% Panels (a)--(c) showing collapse times on linear scale (a,b) and log scale (c)
cmap = lines(length(shelf_names)); count = 1;
splittol = 1250; %where to chop the plots

for is =  1:length(shelf_names)


    %get the data from above
    dMs = cell2mat(dM_shelves(is));
    cts = cell2mat(ct_shelves(is));
    xf = dMs(dMs<=10);
    yf = cts(dMs<=10);

    if yf(1) > splittol
        plot(ax(2),xf , yf,'-','linewidth', 2, 'color', cmap(is,:))
        text(ax(2), xf(1), yf(1), shelf_names_plot(is), 'FontSize', 14);
    else
        plot(ax(1),xf , yf,'-','linewidth', 2, 'color', cmap(is,:))
        text(ax(1), xf(1), yf(1), shelf_names_plot(is) , 'FontSize', 14);
    end

    xxf = xf(xf <= 1);
    yyf =  yf(xf <= 1)./yf(1);
    xxf = smooth(xxf);
    yyf = smooth(yyf);
    plot(ax(3),xxf, yyf,'-','linewidth', 2, 'color', cmap(is,:))
    tt = text(ax(3), xxf(end), yyf(end), shelf_names_plot(is),'FontSize', 14);

end


for i = 1:3
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
ax(2).XTick = ax(1).XTick;
ax(2).XTickLabel = {};
ax(3).YLabel.String = '$\tau_C(\Delta \dot{m})/ \tau_C(\Delta \dot{m} = 0)$';
ax(1).YLabel.String = 'collapse timescale, $\tau_C(\Delta \dot{m}) \mathrm{(years)}$';
ax(1).YLabel.Position(2) = 1400;
%ax(2).YLabel.String = 'collapse timescale, $\tau(\Delta \dot{m})$';
set(ax(3), 'YScale', 'log');
ax(3).XLim = [0,1];
ax(3).YLim = [0.25,1];
%ax(1).YLim = [0,1e3];

%% Panels (d)--(g) showing histograms of the collapse time as a function of scenario by 2100

shelf_names = ["Ross"; "Ronne"; "Filchner"; "Amery"; "Larsen"];
dM          = [0.0001, 0.256, 2.301, 3.611;
    0.0001, 0.223, 0.896, 2.213;
    0.0001, 0.223, 0.896, 2.213;
    0.0001, -0.278, 0.768,3.571;
    0.0001, -0.037, 2.776, 4.312];

ct_ssp = nan(5,4);
% CData  = [70,150, 50; %green: SSP1-2.6
%           50, 120, 190; %blue: SSP2-4.5
%           120, 30, 120]/255; %purple: SSP5-8.5 (colours from SROCC https://www.ipcc.ch/site/assets/uploads/sites/3/2019/11/SROCC_FinalDraft_Chapter1-SM.pdf)

CData = linspecer(4);
for i = 1:5
    for is = 1:4
        if is > 1
            fpath = strcat('../gendata/figure4/step_', num2str(step), '/', shelf_names(i), '/dM_', num2str(dM(i, is)), '.mat');
            ct_ssp(i,is) = load(fpath, 'average_collapse_time').average_collapse_time;
        else
            fpath = strcat('../gendata/figure4/step_', num2str(step), '_CircumAntCoeffs/', shelf_names(i), '/dM_', num2str(dM(i, is)), '.mat');
            ct_ssp(i,is) = load(fpath, 'collapse_time').collapse_time;
        end


        b = bar(ax(i+3), 1:4,ct_ssp(i,:), 0.65);
        b.FaceColor = 'flat';
        b.CData = CData;

        box(ax(i+3), 'on');
        ax(i+3).XLim = [0.5, 4.5];
        ax(i+3).XTick = 1:4;
        ax(i+3).XTickLabel = {'Present day','SSP1-1.9', 'SSP2-4.5', 'SSP5-8.5'};
        ax(i+3).XTickLabelRotation = 45;
        ax(i+3).FontSize = 15;
        ax(i+3).YLabel.String = '$\tau_C \mathrm{(years)}$';
        ax(i+3).YLabel.Interpreter = 'latex';


    end
end
ax(5).YLim = [0, 3500];
shg


