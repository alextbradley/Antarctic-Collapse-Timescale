% Make figure 4 of the ms, showing (a) collapse timescale as a function of
% dM and (b)--(e) bar charts of collapse timescale under different SSP
% scenarios for selected ice shelves

% ATB (aleey@bas.ac.uk), 21/02/23. MIT licence.

%% Preliminaries
addpath('../functions');
%f = load('../data/ice_sheet_data.mat', 'H', 'm'); %get the ice sheet thickness data (use this for shelf mask)

step = 1;        %grid resolution in km

shelf_names = ["Abbot","Amery","Borchgrevink","Filchner","FimbulJelbart",...
    "KingBaudoin","RiiserLarsen", "Ronne","Ross",...
    "Shackleton", "West"];

shelf_names_plot = ["ABB", "AME", "BOR", "FIL", "FIM", "BAU", "RII",...
    "RON", "ROSS", "SHA", "WES"];

data_folder = strcat('../gendata/figure4/step_', num2str(step));
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
    is
    for im  = 1:length(jdir)
        ff = load(strcat(fpath, '/', jdir(im).name));

        dMs(im) = ff.dM;
        cts(im) = ff.average_collapse_time;
    end

    %sort into increaseing dm order
    [dMs,I] = sort(dMs);
    cts = cts(I);

    %store the results
    ct_shelves{is} = cts;
    dM_shelves{is} = dMs;


end

%% Setup figure
fig =figure(1);clf
fig.Position(3:4) = [1150, 450];

h = 0.4; %small plot size
gx = 0.05; %small plot size gap
gy = 0.05; 
w = 0.13; %small plot width
wbig = 2*w + 0.2;
pady = 0.1;
extragap = 0.1; %extra padding int he middle
startx = (1 - 2*w - wbig - gx - extragap)/2;
starty = (1 - 2*h - gy + pady)/2;
positions = [startx, starty, wbig,  2*h + gy;       
             startx + wbig + extragap, starty, w, h;
             startx + wbig + extragap + w + gy, starty, w, h;
             startx + wbig + extragap, starty + h + gy, w, h;
             startx + wbig + extragap + w + gy, starty + h + gy, w, h];

for i = 1:5
    ax(i) = subplot('Position', positions(i,:));
    hold on; box on;
    ax(i).FontSize = 14;
end

%% add (a)

% Panels (a)--(c) showing collapse times on linear scale (a,b) and log scale (c)
cmap = linspecer(length(shelf_names)); count = 1;

for is =  1:length(shelf_names)


    %get the data from above
    dMs = cell2mat(dM_shelves(is));
    cts = cell2mat(ct_shelves(is));
    xf = dMs(dMs<=10);
    yf = cts(dMs<=10);

    plot(ax(1),xf , yf/yf(1),'-','linewidth', 2, 'color', cmap(is,:))
    tt = text(ax(1), xf(end), yf(end)/yf(1), shelf_names_plot(is),'FontSize', 14);

end
set(ax(1), 'YScale', 'log');
ax(1).XLim = [0,10];
ax(1).YLim = [0.003,1];
grid(ax(1), 'on');
ax(1).XLabel.String = '$\Delta \dot{m}$ (m/yr)';
ax(1).YLabel.String = '$\tau_c(\Delta \dot{m})/\tau_c(\Delta \dot{m} = 0)$';
ax(1).XLabel.Interpreter = 'latex';
ax(1).YLabel.Interpreter = 'latex';

%% add (b)--(d)
shelf_names_ssp = ["Ross"; "Ronne"; "Filchner"; "Amery"];
dM          = [0.0001, 0.256, 2.301, 3.611;
    0.0001, 0.223, 0.896, 2.213;
    0.0001, 0.223, 0.896, 2.213;
    0.0001, -0.278, 0.768,3.571];

ct_ssp = nan(5,4);
% CData  = [70,150, 50; %green: SSP1-2.6
%           50, 120, 190; %blue: SSP2-4.5
%           120, 30, 120]/255; %purple: SSP5-8.5 (colours from SROCC https://www.ipcc.ch/site/assets/uploads/sites/3/2019/11/SROCC_FinalDraft_Chapter1-SM.pdf)

step = 10; 
CData = linspecer(4);
CData = linspecer(3);
CData = CData([3,1,2],:); %put green as ssp1, blue as ssp2 and red as ssp 3
CData = [254/255, 229/255, 134/255; CData]; %add the present day

for i = 1:4
    for is = 1:4
        if is > 1
            fpath = strcat('../gendata/figure4/step_', num2str(step), '/', shelf_names_ssp(i), '/dM_', num2str(dM(i, is)), '.mat');
            ct_ssp(i,is) = load(fpath, 'average_collapse_time').average_collapse_time;
        else
            fpath = strcat('../gendata/figure4/step_', num2str(step), '_CircumAntCoeffs/', shelf_names_ssp(i), '/dM_', num2str(dM(i, is)), '.mat');
            ct_ssp(i,is) = load(fpath, 'collapse_time').collapse_time;
        end


        b = bar(ax(i+1), 1:4,ct_ssp(i,:), 0.65);
        b.FaceColor = 'flat';
        b.CData = CData;

        box(ax(i+1), 'on');
        ax(i+1).XLim = [0.5, 4.5];
        ax(i+1).XTick = 1:4;
        ax(i+1).XTickLabel = {'Present day','SSP1-1.9', 'SSP2-4.5', 'SSP5-8.5'};
        ax(i+1).XTickLabelRotation = 45;
        ax(i+1).FontSize = 15;
        ax(i+1).YLabel.String = '$\tau_C \mathrm{(years)}$';
        ax(i+1).YLabel.Interpreter = 'latex';


    end
end
ax(3).YLim = [0, 3500];
shg

%% add the inset in (a)
% axnew = axes();
% hold(axnew, 'on');
%           
% axnew.Position = [0.31, 0.74, 0.2,0.2];
% box(axnew, 'on');
% for is =  1:length(shelf_names)
% 
% 
%     %get the data from above
%     dMs = cell2mat(dM_shelves(is));
%     cts = cell2mat(ct_shelves(is));
%     xf = dMs(dMs<=10);
%     yf = cts(dMs<=10);
% 
%     plot(axnew,xf , yf/yf(1),'-','linewidth', 2, 'color', cmap(is,:))
% 
% end    
%              
% axnew.XLim = [0,1];
% set(axnew, 'YScale', 'log');
% axnew.YLim = [0.2, 1];