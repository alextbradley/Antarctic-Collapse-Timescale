% Make supplementary figures J, showing the distributions of collapse
% timescale for different ice shelves and for different theories. 

%
% Load data
%
fig2data = load('fig2_out_.mat');
crevasse_times        = fig2data.crevasse_times;
crevasse_times_Nye    = fig2data.crevasse_times_Nye;
crevasse_times_ModNye = fig2data.crevasse_times_ModNye;
names = fig2data.shelf_names_adj;
avect = fig2data.ave_crevasse_time';
avect_Nye = fig2data.ave_crevasse_time_Nye';
avect_ModNye = fig2data.ave_crevasse_time_ModNye';

% 
% figure setup
%
figure(1);clf;
colors = [37,28,197; 
          197,37,28;
          28, 197, 38]/255;

%
% loop over shelves
%
count = 1;
%means = nan(length(names), 3);
for i = 1:length(names)
    if i ~= 25 %skip west?
    ax(count) = subplot(5,5,count); hold on; box on
    title(names(i))

    %
    % extract data
    %
    ct = cell2mat(crevasse_times(i));
    ct_Nye = cell2mat(crevasse_times_Nye(i));
    ct_ModNye = cell2mat(crevasse_times_ModNye(i));

    ct(ct == 0) = 1; %make the plotting a bit nicer on w/ log scale
    ct_Nye(ct_Nye == 0) = 1;
    ct_ModNye(ct_ModNye == 0) = 1; 


    maxval = max([ct; ct_Nye; ct_ModNye]);
    
    x = logspace(0,log10(maxval), 1e2); %where to evaluate the kde

    %
    % LEFM
    %
    kde = fitdist(ct,'kernel');
    y = pdf(kde,x);
    plot(x, y,'linewidth', 1.5, 'color', colors(1,:));
   % mean_val = mean(kde);
   % [~,idx] = min(abs(y-mean_val));
   % plot(x(idx), y(idx), 'o', 'markerfacecolor', colors(1,:), 'markeredgecolor', 'k', 'HandleVisibility','off');
   % means(i,1) = y(idx);

    %
    % Nye
    %
    kde = fitdist(ct_Nye,'kernel');
    y = pdf(kde,x);
    plot(x, y,'linewidth', 1.5, 'color', colors(2,:));
%     mean_val = mean(kde);
%     [~,idx] = min(abs(y-mean_val));
%     plot(x(idx), y(idx), 'o', 'markerfacecolor', colors(2,:), 'markeredgecolor', 'k','HandleVisibility','off');
%     means(i,2) = y(idx);

    
    %
    % Mod Nye
    %
    kde = fitdist(ct_ModNye,'kernel');
    y = pdf(kde,x);
    plot(x, y,'linewidth', 1.5, 'color', colors(3,:));
%     mean_val = mean(kde);
%     [~,idx] = min(abs(y-mean_val));
%     plot(x(idx), y(idx), 'o', 'markerfacecolor', colors(3,:), 'markeredgecolor', 'k','HandleVisibility','off');
%     means(i,3) = y(idx);


    ax(count).XLim = [1,maxval];
    ax(count).XLabel.String = 'crevasse timescale';
    ax(count).YLabel.String = 'density';
    set(ax(count), 'XScale', 'log');
    yl = ax(count).YLim;

    plot(avect(i)*[1,1], ax(count).YLim, '--', 'color', colors(1,:), 'linewidth', 1.5);
    plot(avect_Nye(i)*[1,1], ax(count).YLim, '--', 'color', colors(2,:), 'linewidth', 1.5);
    plot(avect_ModNye(i)*[1,1], ax(count).YLim, '--', 'color', colors(3,:), 'linewidth', 1.5);
    

    if count == 1
        legend('LEFM', 'Nye', 'Mod Nye');
    end
        

    count = count + 1;
    end


end
   

    