% Get the temperature profiles from Wang et al. figure 6.
%
% ATB (aleey@bas.ac.uk), 04/10/23, MIT licence


%% Panel a of Wang: profiles AM01, AM04, AM05
% show the image
I = imread("Wang2022Fig6a.png");
figure(1); clf; imshow(I); shg

%AM01
title('click along AM01')
[x_AM01,y_AM01] = ginput();
hold on; plot(x_AM01, y_AM01,  'ko-', 'markerfacecolor', 'k')

%AM04
title('click along AM04')
[x_AM04,y_AM04] = ginput();
hold on; plot(x_AM04, y_AM04,  'ko-', 'markerfacecolor', 'k')

%AM05
title('click along AM05')
[x_AM05,y_AM05] = ginput();
hold on; plot(x_AM05, y_AM05,  'ko-', 'markerfacecolor', 'k')



% calibrate:
title('click on [-20, -200]');
[xneg20, yneg200] = ginput(1);
plot(xneg20, yneg200, 'ro', 'markerfacecolor', 'r')

title('click on [-15, -300]');
[xneg15, yneg300] = ginput(1);
plot(xneg15, yneg300, 'ro', 'markerfacecolor', 'r')

% work out pixels to depth and temp
temp = @(pix) (-15 + (-20 + 15)/(xneg20 - xneg15)*(pix - xneg15));
depth = @(pix) (-300 + (-300 + 200)/(yneg300 - yneg200)*(pix - yneg300));

% get the profiles
temp_AM01 = temp(x_AM01);
depth_AM01 = depth(y_AM01);
temp_AM04 = temp(x_AM04);
depth_AM04 = depth(y_AM04);
temp_AM05 = temp(x_AM05);
depth_AM05 = depth(y_AM05);


figure(2); clf; hold on
plot(temp_AM01, depth_AM01, 'Color',[143, 182, 92]/255, 'linewidth', 1.5);
plot(temp_AM04, depth_AM04, 'Color',[118, 85, 147]/255, 'linewidth', 1.5);
plot(temp_AM05, depth_AM05, 'Color',[79, 162, 187]/255, 'linewidth', 1.5);

grid on
legend("AM01", "AM04", "AM05");

% save the data
save("Wang2022Fig6_profiles.mat",'temp_AM01', 'depth_AM01', 'temp_AM04', 'depth_AM04','temp_AM05', 'depth_AM05' )

%% Panel b of Wang: profiles AM02, AM03, AM06
I = imread("Wang2022Fig6b.png");
figure(1); clf; imshow(I); shg

%AM02
title('click along AM02')
[x_AM02,y_AM02] = ginput();
hold on; plot(x_AM02, y_AM02,  'ko-', 'markerfacecolor', 'k')

%AM03
title('click along AM03')
[x_AM03,y_AM03] = ginput();
hold on; plot(x_AM03, y_AM03,  'ko-', 'markerfacecolor', 'k')

%AM05
title('click along AM06')
[x_AM06,y_AM06] = ginput();
hold on; plot(x_AM06, y_AM06,  'ko-', 'markerfacecolor', 'k')

% calibrate:
title('click on [-20, -200]');
[xneg20, yneg200] = ginput(1);
plot(xneg20, yneg200, 'ro', 'markerfacecolor', 'r')

title('click on [-15, -300]');
[xneg15, yneg300] = ginput(1);
plot(xneg15, yneg300, 'ro', 'markerfacecolor', 'r')

% work out pixels to depth and temp
temp = @(pix) (-15 + (-20 + 15)/(xneg20 - xneg15)*(pix - xneg15));
depth = @(pix) (-300 + (-300 + 200)/(yneg300 - yneg200)*(pix - yneg300));

% get the profiles
temp_AM02 = temp(x_AM02);
depth_AM02 = depth(y_AM02);
temp_AM03 = temp(x_AM03);
depth_AM03 = depth(y_AM03);
temp_AM06 = temp(x_AM06);
depth_AM06 = depth(y_AM06);

%
figure(3);clf; hold on
plot(temp_AM02, depth_AM02, 'Color',[79, 115, 175]/255, 'linewidth', 1.5);
plot(temp_AM03, depth_AM03, 'Color',[179, 68, 71]/255, 'linewidth', 1.5);
plot(temp_AM06, depth_AM06, 'Color',[228, 170, 68]/255, 'linewidth', 1.5);

grid on
legend("AM02", "AM03", "AM06");
save("Wang2022Fig6_profiles.mat",'temp_AM02', 'depth_AM02', 'temp_AM03', 'depth_AM03','temp_AM06', 'depth_AM06','-append' )

% put into a struct
Wang2022data = struct;
Wang2022data(1).T  = temp_AM01;
Wang2022data(1).depth  = depth_AM01;
Wang2022data(1).label  = "AM01";

Wang2022data(2).T  = temp_AM02;
Wang2022data(2).depth  = depth_AM02;
Wang2022data(2).label  = "AM02";

Wang2022data(3).T  = temp_AM03;
Wang2022data(3).depth  = depth_AM03;
Wang2022data(3).label  = "AM03";

Wang2022data(4).T  = temp_AM04;
Wang2022data(4).depth  = depth_AM04;
Wang2022data(4).label  = "AM04";

Wang2022data(5).T  = temp_AM05;
Wang2022data(5).depth  = depth_AM05;
Wang2022data(5).label  = "AM05";

Wang2022data(6).T  = temp_AM06;
Wang2022data(6).depth  = depth_AM06;
Wang2022data(6).label  = "AM06";
save("Wang2022Fig6_profiles.mat",'Wang2022data');

%% Get the data for the other solutions
count = 1;
model_cmap = [233, 166, 60;
        41, 99, 174;
        205, 71, 40; 
        107, 166, 72]/255;

plotnos = [3,6,4,2,1,5]; %plot number from Wang 2022
for AM = 1:6
    figure(2); clf; I = imread(strcat("Wang2022Fig7_AM0", num2str(AM), ".png"));
    imshow(I);

    title('Click along BMB_ISMIP6 (yellow)', 'interpreter', 'none')
    [x_BMB_ISMIP6, y_BMB_ISMIP6] = ginput();

    title('Click along BMB_ROMS (blue)','interpreter', 'none')
    [x_BMB_ROMS, y_BMB_ROMS] = ginput();

    title('Click along BMB_CAL (red)','interpreter', 'none')
    [x_BMB_CAL, y_BMB_CAL] = ginput();

    title('Click along BMB_CAL2 (green)','interpreter', 'none')
    [x_BMB_CAL2, y_BMB_CAL2] = ginput();

    % calibrate:
    title('click on [-18, 0]');
    [xneg18, y0] = ginput(1);
    hold on; plot(xneg18, y0, 'ro', 'markerfacecolor', 'r')

    title('click on [-2, 1]');
    [xneg2, y1] = ginput(1);
    plot(xneg2, y1, 'ro', 'markerfacecolor', 'r')


%     % work out pixels to depth and temp
     temp = @(pix) (-18 + (-18 + 2)/(xneg18 - xneg2)*(pix - xneg18));
     depth = @(pix) (0 + (0 - 1)/(y0 - y1)*(pix - y0));

    % plot these
    figure(3); 
    ax(count) = subplot(2,3,plotnos(count)); hold on; box on; set(ax(count), 'YDir', 'reverse');
    plot(temp(x_BMB_ISMIP6), depth(y_BMB_ISMIP6), 'color', model_cmap(1,:))
    plot(temp(x_BMB_ROMS), depth(y_BMB_ROMS), 'color', model_cmap(2,:))
    plot(temp(x_BMB_CAL), depth(y_BMB_CAL), 'color', model_cmap(3,:))
    plot(temp(x_BMB_CAL2), depth(y_BMB_CAL2), 'color', model_cmap(4,:))
    plot(Wang2022data(count).T, -Wang2022data(count).depth /max(abs(Wang2022data(count).depth )), 'k')
    ax(count).YLim = [0,1];
    title(['AM_0', num2str(AM)], 'interpreter', 'none');

    % store 
    Wang2022data(count).T_BMB_ISMIP6 = temp(x_BMB_ISMIP6);
    Wang2022data(count).depth_BMB_ISMIP6 = -depth(y_BMB_ISMIP6)*max(abs(Wang2022data(count).depth));

    Wang2022data(count).T_BMB_ROMS = temp(x_BMB_ROMS);
    Wang2022data(count).depth_BMB_ROMS = -depth(y_BMB_ROMS)*max(abs(Wang2022data(count).depth));
    
    Wang2022data(count).T_BMB_CAL = temp(x_BMB_CAL);
    Wang2022data(count).depth_BMB_CAL = -depth(y_BMB_CAL)*max(abs(Wang2022data(count).depth));
    
    Wang2022data(count).T_BMB_CAL2 = temp(x_BMB_CAL2);
    Wang2022data(count).depth_BMB_CAL2 = -depth(y_BMB_CAL2)*max(abs(Wang2022data(count).depth));
    

count = count +1;
end



