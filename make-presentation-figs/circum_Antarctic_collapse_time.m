% Generate data for the circum Antartic collapse data
% output the solution row by row in the matrix, saved to folder
% 'circum-Antarctic-collapse-time'

%% Load the ice shelf data
addpath('../functions/')
f  = load('../data/ice_sheet_data.mat');

%% Generate tags mask
% tags:
%    1: outside domain (no thickness data)
%    2: other missing data (no melt, strain rate, dhdt)
%    3: negative strain rate
%    4: positive thinning rate (currently stable)
%    5: positive thinning rate (currently unstable)
%    6: other (good data)

xs = f.xx;
ys = f.yy;
ms = f.m;
hs = f.H;
dhdts = f.dhdtadj;
strains = f.eflow;
Bs = f.B;
ms(ms < 1e-1) = 1e-1; %set a minimum melt value
dhdts = -abs(dhdts);

tags = 6*ones(size(ms));

tags(dhdts > 0) = 4;
tags(strains < 0) = 3;
tags(isnan(ms) | isnan(dhdts) | isnan(strains)) = 2;
tags(isnan(hs)) = 1;

%
figure(1); clf; imagesc(ys,xs,tags)
c = colorbar;
clim([1,6])
colormap(lines(6))
c.Ticks = linspace(2-(7/12),5+(7/12),6);
c.TickLabels = {'1','2','3','4','5','6'};

