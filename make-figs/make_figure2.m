% Make figure 2 of the ms, showing:
% (a) Scatter plot of melt rate against strain rate, with regression 
% (b) Scatter plot of melt rate against thinning rate, with regression
%
% % 17/02/23, ATB (aleey@bas.ac.uk), MIT licence
%
% Note that this script also produces the violin plot data shown in figure
% 3b, and should be ran prior to the figure 3 script.
%% Preliminaries
clear
addpath('../functions');
saveout = 1; %flag to specify saving output. (Need output from this script to run figure3 and 4 scripts)
f = load('../data/ice_sheet_data.mat');
jdir = dir('../data/ice-shelves/major/*.mat'); %where shelf files are stored

%we'll also work out shelf by shelf collapse data here
fig3data = load('figure3-data.mat');
ct = fig3data.collapse_time; %collapse time data
ct_Nye = fig3data.collapse_time_Nye;
tags = fig3data.tags;

%% Fetch the data
% Loop over shelves, and for each, get the melt rate, thinning rate, strain
% rate, and collapse time. 

% init storage
shelf_names = strings;
melt_rates  = cell(1,length(jdir));
strain_rates  = cell(1,length(jdir));
thinning_rates = cell(1,length(jdir));
crevasse_times = cell(1,length(jdir));

ave_melt_rate = nan(1,length(jdir));
ave_strain_rate = nan(1,length(jdir));
ave_thinning_rate = nan(1,length(jdir));
ave_crevasse_time = nan(1,length(jdir));
ave_crevasse_time_Nye = nan(1,length(jdir));

std_melt_rate = nan(1,length(jdir));
std_strain_rate = nan(1,length(jdir));
std_thinning_rate = nan(1,length(jdir));

for i  = 1:length(jdir)

    %load data for this shelf
    shelf = jdir(i).name;
    shelf = strrep(shelf,'.mat',''); %strip the .mat at the end
    shelf_names(i) = shelf;
    fname = strcat('../data/ice-shelves/major/' ,shelf, '.mat');
    g = load(fname);

    % restrict data to those points in shelf
    m_shelf = f.m; 
    m_shelf = m_shelf(g.IN); %only the points in this shelf
    %m_shelf(m_shelf < 0) = 0; %remove any negative melt rate points?
    h_shelf = f.H; 
    h_shelf = h_shelf(g.IN);
    epsxx_shelf = f.eflow; 
    epsxx_shelf = epsxx_shelf(g.IN);
    thinrate_shelf = g.mean_dhdt*ones(size(epsxx_shelf));  %set the thinning rate to pre-determined value

    % restrict to where we have data
    idx =  (~isnan(h_shelf)) &  (~isnan(m_shelf))  & (~isnan(epsxx_shelf)); %points where we have data for melt, thichkness and strain rate
    h_shelf = h_shelf(idx);
    m_shelf = m_shelf(idx); %arrays with points in particular shelf with both melt and thicknes
    epsxx_shelf = epsxx_shelf(idx);
    epsxx_shelf = epsxx_shelf(epsxx_shelf > 0); % keep only positive strain rates

    % store
    melt_rates{i}     = m_shelf;
    %figure(2); clf; histogram(m_shelf); title(shelf); drawnow; xlim([-10,10]); shg; pause
    strain_rates{i}   = epsxx_shelf;
    thinning_rates{i} = thinrate_shelf;

    ave_melt_rate(i) = mean(m_shelf);
    ave_strain_rate(i) = mean(epsxx_shelf);
    ave_thinning_rate(i) = mean(thinrate_shelf);

    std_melt_rate(i) = std(m_shelf);
    std_strain_rate(i) = std(epsxx_shelf);
    std_thinning_rate(i) = std(thinrate_shelf);
 
    % compute the collapse timescale    
    ct_shelf = ct(g.IN); %only the points in this shelf
    ct_keep = (ct_shelf(~isnan(ct_shelf))); %all non nan and non infpoints in shelf
    ct_keep = ct_keep(~isinf(ct_keep));
    ct_keep = ct_keep(ct_keep < 5e4); %remove very long points from the distribution
    crevasse_times{i} = ct_keep;

    if ~isempty(ct_keep)
        pd = fitdist(ct_keep,'kernel');
        ave_crevasse_time(i) = mean(pd);
    end

    %compute mean for Nye
    ct_Nye_shelf = ct_Nye(g.IN);
    ct_Nye_keep = (ct_Nye_shelf(~isnan(ct_Nye_shelf))); %all non nan points in shelf
    ct_Nye_keep = ct_Nye_keep(~isinf(ct_Nye_keep));
    ct_Nye_keep  = ct_Nye_keep(ct_Nye_keep < 5e4);

    if ~isempty(ct_Nye_keep)
        pd = fitdist(ct_Nye_keep,'kernel');
        ave_crevasse_time_Nye(i) = mean(pd);
    end
   

end

% Adjust the names
shelf_names_adj = shelf_names;
shelf_names_adj(shelf_names_adj ==  "Thwaites") = "THW";
shelf_names_adj(shelf_names_adj ==  "Crosson") = "CRO";
shelf_names_adj(shelf_names_adj ==  "PopeSmithKohler") = "PSK";
shelf_names_adj(shelf_names_adj ==  "Dotson") = "DOT";
shelf_names_adj(shelf_names_adj ==  "PineIsland") = "PIG";
shelf_names_adj(shelf_names_adj ==  "Larsen") = "LAR";
shelf_names_adj(shelf_names_adj ==  "Wilkins") = "WIL";
shelf_names_adj(shelf_names_adj ==  "Brunt") = "BRU";
shelf_names_adj(shelf_names_adj ==  "George6") = "GVI";
shelf_names_adj(shelf_names_adj ==  "Getz") = "GET";
shelf_names_adj(shelf_names_adj ==  "Shackleton") = "SHA";
shelf_names_adj(shelf_names_adj ==  "West") = "WES";
shelf_names_adj(shelf_names_adj ==  "Cook") = "COO";
shelf_names_adj(shelf_names_adj ==  "Nansen") = "NAN";
shelf_names_adj(shelf_names_adj ==  "Cosgrove") = "COS";
shelf_names_adj(shelf_names_adj ==  "Borchgrevink") = "BOR";
shelf_names_adj(shelf_names_adj ==  "Abbot") = "ABB";
shelf_names_adj(shelf_names_adj ==  "SwinburneSulzbergerNickerson") = "SSN";
shelf_names_adj(shelf_names_adj ==  "FimbulJelbart") = "FIM";
shelf_names_adj(shelf_names_adj ==  "KingBaudoin") = "BAU";
shelf_names_adj(shelf_names_adj ==  "TottenMoscow") = "TOT";
shelf_names_adj(shelf_names_adj ==  "RiiserLarsen") = "RII";
shelf_names_adj(shelf_names_adj ==  "Amery") = "AME";
shelf_names_adj(shelf_names_adj ==  "Ronne") = "RON";
shelf_names_adj(shelf_names_adj ==  "Ross") = "ROSS ";
shelf_names_adj(shelf_names_adj ==  "Filchner") = "FIL";
%% Do the regression
% Power law
lm_thinning_rate_power = fitlm(log10(ave_melt_rate),log10(ave_thinning_rate));
lm_strain_rate_power = fitlm(log10(ave_melt_rate),log10(ave_strain_rate));

%get coeffs
p_thinning_rate_power = table2array(lm_thinning_rate_power.Coefficients(:,1));
p_strain_rate_power = table2array(lm_strain_rate_power.Coefficients(:,1));

fprintf('Regression coefficient between log10(melt rate) and log10(thinning rate) is %.4g \n', p_thinning_rate_power(2))
fprintf('Regression coefficient between log10(melt rate) and log10(strain rate) is %.4g \n', p_strain_rate_power(2))

% determine the correlation coefficienets
[R_dhdt,P_dhdt, RL_dhdt,RU_dhdt] = corrcoef(ave_melt_rate, ave_thinning_rate);
fprintf('Correlation coefficient between log10(melt rate) and log10(thinning rate) is %.3g \n', R_dhdt(1,2))
fprintf('p-value of correlation coefficient between log10(melt rate) and log10(thinning rate) is %.3g \n', P_dhdt(1,2))
fprintf('\n')

[R_epsxx,P_epsxx, RL_epsxx, RU_epsxx] = corrcoef(ave_melt_rate, ave_strain_rate);
fprintf('Correlation coefficient between log10(melt rate) and log10(strain rate) is %.3g \n', R_epsxx(1,2))
fprintf('p-value of correlation coefficient between log10(melt rate) and log10(strain rate) is %.3g \n', P_epsxx(1,2))


% Linear fit
lm_thinning_rate = fitlm((ave_melt_rate),(ave_thinning_rate));
lm_strain_rate = fitlm((ave_melt_rate),(ave_strain_rate));

%get coeffs
p_thinning_rate = table2array(lm_thinning_rate.Coefficients(:,1));
p_strain_rate = table2array(lm_strain_rate.Coefficients(:,1));
%% Setup the figure
fig = figure(1);clf

positions = [0.15, 0.08, 0.8, 0.41;
             0.15, 0.58, 0.8, 0.41];
sz        = size(positions);
fig.Position(3:4) = [430,700];
fs = 14; %fontsize
ms = 8; %markersize

for i = 1:2
    ax(i) = subplot('Position', positions(i,:)); hold on; box on; grid on
    ax(i).FontSize = fs;
end

%colourmap for melting
cmap = (cmocean('balance', 100)); cmap = cmap(20:end-20,:);

% Make the plot
% work out the map between log melt rates and colours: returns 1 for the lowest
% melt rate and length(cmap) for the highest melt rate, with linear interp
% in between
cindex = @(melt, melt_rates) (1 + (length(cmap)-1)*(log10(melt)-log10(min(melt_rates)))/(log10(max(melt_rates))-log10(min(melt_rates))));


% add the regression lines
m = logspace(-2,2);
plot(ax(2), m, 10^(p_strain_rate_power(1))*m.^(p_strain_rate_power(2)), 'k', 'linewidth', 1.5 );
plot(ax(1), m, 10^(p_thinning_rate_power(1))*m.^(p_thinning_rate_power(2)), 'k', 'linewidth', 1.5 );
% 
% linear regression
% plot(ax(2), m, p_strain_rate(1) + m*(p_strain_rate(2)), 'k--', 'linewidth', 1.5 );
% plot(ax(1), m, p_thinning_rate_power(1)+m* (p_thinning_rate_power(2)), 'k--', 'linewidth', 1.5 );

plot(ax(2), m, 0.0033*m.^(0.4), 'r--')
plot(ax(1), m, 0.66*m.^(0.4), 'r--')

shelf_colours = nan(length(ave_melt_rate), 3);
for i = 1:length(jdir)
    cidx = round(cindex(ave_melt_rate(i), ave_melt_rate));
    
    % plot the strain rate with errors
%     plot(ax(2), ave_melt_rate(i)*[1,1], [max(ave_strain_rate(i)- std_strain_rate(i), 0.0001), ave_strain_rate(i) + std_strain_rate(i)], 'k', 'linewidth', 1);
%     plot(ax(2), [ave_melt_rate(i)-std_melt_rate(i), ave_melt_rate(i)+std_melt_rate(i)], ave_strain_rate(i)*[1,1], 'k', 'linewidth', 1);
    plot(ax(2), ave_melt_rate(i), ave_strain_rate(i), 'ko', 'markersize', ms, 'markerfacecolor', cmap(cidx, :))

    % plot the thinning rate
    plot(ax(1), ave_melt_rate(i), ave_thinning_rate(i), 'ko', 'markersize', ms, 'markerfacecolor', cmap(cidx, :))

    % add the shelf name
    t(i) = text(ax(2),ave_melt_rate(i)*1.01, ave_strain_rate(i)*1.01,shelf_names_adj(i), 'FontSize', fs, 'FontName', 'Arial');
    t(i) = text(ax(1),ave_melt_rate(i)*1.01, ave_thinning_rate(i)*1.01,shelf_names_adj(i), 'FontSize', fs, 'FontName', 'Arial');

    % store the shelf colour
    shelf_colours(i,:) = cmap(cidx, :);

end

% Tidy the plot

set(ax(1), 'XScale', 'log')
set(ax(1), 'YScale', 'log')
ax(1).YLim = [0.03, 10];

set(ax(2), 'XScale', 'log')
set(ax(2), 'YScale', 'log')
ax(2).YLim = [0.0008, 0.04];

for i = 1:2
    ax(i).XLim = [0.05, 30];
end

ax(1).XLabel.String = 'melt rate (m/yr)';
ax(2).XLabel.String = 'melt rate (m/yr)';

ax(1).YLabel.String = 'thinning rate (m/yr)';
ax(2).YLabel.String = 'strain rate';

%% Save data for use in figures 3 and 4
if saveout
    save('fig2_out_.mat', 'shelf_names_adj', 'shelf_colours','ave_melt_rate', 'ave_crevasse_time', "ave_crevasse_time_Nye", "crevasse_times")
end




