% 3 theories, each with 3 surface temperatures assuming linear
% temperature profiles in the vertical

clear all;
close all;
tic;

% 3 temperatures: -2, -17, -32
Ts_c = [-2, -17, -32]; % celsius
Ts_k = Ts_c + 273.15;
% start by setting up common stresses, B(T), sampled ztilde

rhoi = 917; % ice density kg/m^3
g = 9.8; % m/s^2
rhow = 1028; % seawater density
Tb_choice = -2; % ice basal temperature celsius
Tb_k = 273.15 + Tb_choice; % kelvin

smallest = 1E-6;
ztilde = smallest:smallest:1; % z/H, =0 at base, =1 at surface

beta = 0:0.02:2; % buttressing reduction factor, large for Nye case.

% choosing strain rates
% 2 \langle B(T) \rangle \dot{\epsilon}_{xx}^{1/n} \approx \langle R_{xx} \rangle = \beta 0.5(1-rhoi/rhow) rhoi g H
% \dot{\epsilon}_{xx}^{1/n} = \beta (0.25/ \langle B(T) \rangle ) (1-rhoi/rhow) rhoi g H
Rxx = nan(length(Ts_c), length(ztilde), length(beta));
B = nan(length(Ts_c), length(ztilde), 1);
T = nan(length(Ts_c), length(ztilde));

intB_vec = nan(length(Ts_c), 1);

exx = nan(length(Ts_c), 1, length(beta));


B0 = 2.207;
T0 = 3155; % kelvin
c=0.16612;
Tr = 273.39;
k=1.17;

n=3;
H=300; % meters, arbitrary thickness value but representative

for i=1:length(Ts_c) 
    T(i,:) = 273.15 + Tb_choice + (Ts_c(i) - Tb_choice) .* ztilde; % kelvin
    B(i,:,1) = B0*exp(T0./(T(i,:))-c./((Tr-T(i,:)).^k)); % numbers from Hooke 1981

    Bfz = @(zt) B0*exp( T0 ./ (zt * (Ts_k(i) - Tb_k) + Tb_k) - c ./ (Tr - (zt * (Ts_k(i) - Tb_k) + Tb_k)).^k ); % linear
    intB_vec(i) = integral(Bfz, 0, 1);
        
    
    exx(i,:,:) = (beta .* 0.25*(1-rhoi/rhow)*rhoi*g*H/intB_vec(i)).^n; % 1/yr

    Rxx(i,:,:) = 2 .* B(i,:,1) .* exx(i,1,:) .^ (1/n);
end


overburden = - rhoi*g*H*(1-ztilde);
water_pressure = rhow*g*H*(rhoi/rhow - ztilde);

water_pressure(ztilde>=(rhoi/rhow)) = 0;

stress_from_vert_traction = Rxx + overburden + water_pressure; % dimensional
ND_stress_from_vert_traction = stress_from_vert_traction / (0.5*(1-rhoi/rhow)*rhoi*g*H);
toc;

for i=1:length(Ts_c)
    figure; hold on; title(append('$T_s$=', int2str(Ts_c(i))), Interpreter="latex");
    xline(0)
    xlabel('Nondimensional Stress Balance [1]', Interpreter='latex'); 
    ylabel('Nondimensional Height $z/H$ [1]', Interpreter='latex');
    for k=1:length(beta)
    plot(ND_stress_from_vert_traction(i,:,k), ztilde);
    end
end
toc;


% added 8/2/2023 for surface crevasses
surface_stress_from_vert_traction = Rxx + overburden; % dimensional
surface_ND_stress_from_vert_traction = surface_stress_from_vert_traction / (0.5*(1-rhoi/rhow)*rhoi*g*H);
toc;

for i=1:length(Ts_c)
    figure; hold on; title(append('$T_s$=', int2str(Ts_c(i))), Interpreter="latex");
    xline(0)
    xlabel('Nondimensional Ice Stress Balance [1]', Interpreter='latex'); 
    ylabel('Nondimensional Height $z/H$ [1]', Interpreter='latex');
    for k=1:length(beta)
    plot(surface_ND_stress_from_vert_traction(i,:,k), ztilde);
    end
end
toc;

%% Nye

% script it out
% sample 21 <Rxx>/RxxIT, 0:0.1:2
% 

stress_min = nan(length(Ts_c), 1, length(beta));
min_loc = nan(length(Ts_c), 1, length(beta));

db_nye_temp = nan(length(Ts_c), 1, length(beta));
% update 8/2/2023
surface_stress_min = nan(length(Ts_c), 1, length(beta)); 
surface_min_loc = nan(length(Ts_c), 1, length(beta));
ds_nye_temp = nan(length(Ts_c), 1, length(beta)); 
% update 8/2/2023

[~, zsl_ind] = min(abs( ztilde - rhoi/rhow )); % sea level index

dBdzt = nan(length(Ts_c), length(ztilde), length(beta));
% dBdzt_posonly = nan(length(Ts_c), length(ztilde), length(beta));
% eqn: exclude dBdzt > (rhow-rhoi)gH/ (2 \dot{\epsilon}_{xx}^{1/n})
expr = nan(length(Ts_c), 1, length(beta));

for i=1:length(Ts_c)
    for k=1:length(beta)
        % 8/2/2023 Surface Crevasse Update
        % You don't have to get fancy because there can only be one
        % crossing.
        [surface_stress_min(i,1,k),surface_min_loc(i,1,k)] = min(abs(surface_stress_from_vert_traction(i,:,k)));
        % 8/2/2023 Update

        [stress_min(i,1,k),min_loc(i,1,k)] = min(abs(stress_from_vert_traction(i,1:zsl_ind,k)));
        [stress_min_ns, min_loc_ns] = sort(abs(stress_from_vert_traction(i,1:zsl_ind,k))); % no save series

%         a little theory necessary for curvature condition : d sigma_xx^I
%         / dz > 0 needs to be excluded
%         num_deriv = diff(stress_from_vert_traction(i,1:zsl_ind,k)) ./ diff(ztilde(1:zsl_ind)); % numerical derivative for ease
%         [max_stress_diff, stress_diff_loc] = min(num_deriv(num_deriv>0)); % because you know the curvature you can do this

        % issue with this: write it out analytically as from whiteboard
        % then make sure you can get the derivative to sample below the
        % point of zero slope, not the points above.

% %         from mathematica:
%         coeff = B0*(c*k*(Tb_k-Ts_k(i))*(-Tb_k+Tr-(-Tb_k+Ts_k(i))*ztilde).^(-1-k) - (T0*(Tb_k+Ts_k(i)) ./ (Tb_k + (-Tb_k + Ts_k(i)) .*ztilde ).^2 ));
%         dBdzt(i,:,k) = coeff.*exp(-c*(-Tb_k+Tr-(-Tb_k+Ts_k(i)).*ztilde).^(-k) + (T0 ./ (Tb_k + (-Tb_k + Ts_k(i)).*ztilde) ) );
%         expr(i,1,k) = (rhow-rhoi)*g*H./(2*exx(i,1,k).^(1/n));
% 
%         dBdzt_posonly(i,:,k) = dBdzt(i,:,k);
%         dBdzt_posonly(dBdzt_posonly(i,:,k)-expr(i,1,k) < 0) = 10^9;
% 
%         [max_stress_diff, stress_diff_loc] = min(dBdzt_posonly(i,:,k)); % because you know the curvature you can do this


%         if  max_stress_diff > 0
        if round(ztilde(min_loc_ns(1)) - ztilde(min_loc_ns(2)), 3) > 0 % if the second lowest point is way off, take the smaller val
            % find a minimum for curvature negative, aka smallest basal
            % crevasse depth
%             [stress_min(i,1,k),min_loc(i,1,k)] = min(abs(stress_from_vert_traction(i,1:stress_diff_loc,k)));
            stress_min(i,1,k) = stress_min_ns(2);
            min_loc(i,1,k) = min_loc_ns(2);
        end

        if stress_min(i,1,k) < 1
            if k==1
                db_nye_temp(i,1,k) = 0; % 0 for beta=0.
                ds_nye_temp(i,1,k) = 1; % 0 for beta=0. 8/2/2023
            else
                db_nye_temp(i,1,k) = ztilde(min_loc(i,1,k)); % dimensionless
                ds_nye_temp(i,1,k) = ztilde(surface_min_loc(i,1,k)); % dimensionless. 8/2/2023
            end
        else
%             disp("No Stable Basal Crevasse Depth for:")
%             disp(Ts_c(i))
%             disp(beta(k))
        end
    end
end

db_nye_temp_2D = reshape(db_nye_temp, length(Ts_c), length(beta));
ds_nye_temp_2D = reshape(ds_nye_temp, length(Ts_c), length(beta)); % 8/2/2023 Update

figure; hold on;
lge = [];
colors = ["k", "r", "b", "c", "g", "m"];
for i=1:length(Ts_c)
    plot(beta, db_nye_temp_2D(i,:), colors(i), 'linewidth', 4);    
    lg = append("$T_s = $", int2str(Ts_c(i)), "$^{\circ}$C"); lge = [lge, lg];
end
yline(rhoi/rhow, "c:", 'linewidth', 2);
for i=1:length(Ts_c)
    [maxdbval, betaloc] = max(db_nye_temp_2D(i,:));
    plot(beta(betaloc), maxdbval, append(colors(i), "o"), 'markersize', 12);
    plot([beta(betaloc), beta(betaloc)], [maxdbval, ds_nye_temp_2D(i, betaloc)], append(colors(i), "--"), 'linewidth', 2);
end
lg = "Sea Level"; lge = [lge, lg];

% Surface Crevasse Update 8/2/2023
for i=1:length(Ts_c)
    [mindsval, betaloc] = min(ds_nye_temp_2D(i,:));
    plot(beta, ds_nye_temp_2D(i,:), append(colors(i), "-."), 'linewidth', 4);
    plot(beta(betaloc), mindsval, append(colors(i), "o"), 'markersize', 12);
%     lg = append("$T_s = $", int2str(Ts_c(i)), "$^{\circ}$C"); lge = [lge, lg];
end
% Update 8/2/2023

xlabel("$\bar{R}_{xx}/\bar{R}_{xx}^{IT}$ [1]", Interpreter="latex");
% ylabel('$d_b/H$ [1]', Interpreter='latex');
ylabel('$z/H$ [1]', Interpreter='latex'); % 8/2/2023 Update
ylim([0 1])
xlim([0 2.1])
title("Nye's Zero-Yield Stress", Interpreter='latex')
legend(lge, location='southeast', interpreter='latex')
set(gca, 'FontSize', 20)

toc;

%% save?
% saveas(gcf, append(pwd, "/fig2a.fig"))
% saveas(gcf, append(pwd, "/fig2a.jpg"))
saveas(gcf, append(pwd, "/fig2a_08042023.fig"))
saveas(gcf, append(pwd, "/fig2a_08042023.jpg"))

%% Modified Nye
% insertion from force_balance_with_Nye0Stress.m
clear all
close all
tic;
LinearRobinMarine = 0; % Linear 0 Robin 1 Marine 2. new as of 7/10/2023

% g=9.81; % taken out 5/8/2023 making sure same g across the board.
g=9.8;
rhoi=917;
rhow=1028;

% temperature
T=flip(-32:1:-2); % this will be what you vary
Tb=-2;
% Tofz = Tb + (T-Tb)*(y/H);

% intBTtil_vec = zeros(length(T), 1); % dimensionless

% gotta try
iind = 1:1:length(T);
% iind2 = 1:10:length(T); 
% iind2 = [1 21 29];
iind2 = [1 16 31]; % updated 4/5/2023 to reflect same temps as fig 2a
tds = [0:0.000001:(1-rhoi/rhow)+0.05, 1-rhoi/rhow+0.05]; % 8/3/2023 Update
tdb = [0:0.01:rhoi/rhow, rhoi/rhow]; % start with easy case for algorithm development
tdb = unique(tdb);
% tdb = rhoi/rhow;
tds_find = nan(length(iind), length(tdb));
beta = nan(length(iind), length(tdb));
sigma0 = nan(length(iind), length(tdb));

betamax = nan(size(iind));

betaterm_1 = nan(length(iind), length(tdb));
betaterm_2 = nan(length(iind), length(tdb));
betaterm_1_2 = nan(length(iind), length(tdb));
betaterm_3 = nan(length(iind), length(tdb));

tol=1E-5; % error tolerance
for i=iind
    breakflag = 0;
%     tildeB = @(zoH) exp(3155./((zoH.*(T(i)-Tb)+Tb)+273.15)-0.16612./((0.24-(zoH.*(T(i)-Tb)+Tb)).^1.17)); % fxn of z/H
    tildeB = @(zoH) exp(3155./((zoH.*(Tb-T(i))+T(i))+273.15)-0.16612./((0.24-(zoH.*(Tb-T(i))+T(i))).^1.17)); % fxn of z/H
    % new test 9.29.23
    % tildeB = @(zoH) exp(3155./((zoH.*(Tb-T(i))+T(i))+273.15));%-0.16612./((0.24-(zoH.*(Tb-T(i))+T(i))).^1.17)); % fxn of z/H
    % now let's try the taylor expansion
    % T0 = 3155; Tbk = 273.15+Tb; Tsk = T+273.15;
    % z0t = (T0*(1-Tsk(i)/Tbk)/Tbk)^(-1);
    % tildeB = @(zoH) exp(T0/Tbk)*exp((1-zoH)/z0t);
    % end new test

% new 7/10/2023
%     if LinearRobinMarine == 0
%         Bt = @(zt) exp( T0 ./ (zt * (Ts - Tb_k) + Tb_k) - c ./ (Tr - (zt * (Ts - Tb_k) + Tb_k)).^k ); % linear
%         tildeB = @(zoh) exp( T0 ./ (zt * (Ts - Tb_k) + Tb_k) - c ./ (Tr - (zt * (Ts - Tb_k) + Tb_k)).^k ); % differences: Kelvin, zt(Ts-Tb)+Tb. Flip both
%         tildeB = @(zoH) exp(3155./((zoH.*(Tb-T(i))+T(i))+273.15)-0.16612./((0.24-(zoH.*(Tb-T(i))+T(i))).^1.17)); % fxn of z/H
%         iBt_vec(i) = integral(Bt, 0, 1);
    if LinearRobinMarine == 1 % if you want robin solution
        Ts = 273.15 + T(i); % kelvin
        c=0.16612;
        k=1.17;
        B0=2.207;
        T0=3155; % kelvin
        Tr=273.39; % kelvin
        Geo=60*10^-3; % Watts / m^2. Geothermal heat flux at ice divide
        Hdivide=2400;%2400; % meters. Ice divide thickness
        kappa = 10^-6; % m^2/s, thermal diffusivity
        secinyear = 60*60*24*365.25;
        adot = 0.1 / secinyear; % m/s, snowfall rate at ice divide
        kcon = 2.1; % Watts /m /K, thermal conductivity
        erf1 = erf(sqrt(adot*Hdivide/(2*kappa)));
        % erf2 = erf(sqrt(adot/(2*Hdivide*kappa))*Hdivide.*zt);
        pre = (Geo/kcon)*sqrt(pi*Hdivide*kappa/(2*adot));
        
        Tb_k = -2+273.15; % Kelvin
        Gpre = (Tb_k - Ts) / (pre * erf1); % cheesing basal boundary condition
%         zt = 0:0.1:1;

%         check numerically
%         T_question = (Ts + Gpre * pre * (erf1 - erf(sqrt(adot/(2*Hdivide*kappa))*Hdivide.*zt)));
%         T_question2 = (Ts + Gpre * pre * (erf1 - (erf(sqrt(adot/(2*Hdivide*kappa))*Hdivide.*flip(zt))) )); % same thing just flip zt.

%         (Ts + Gpre * pre * (erf1 -
%         (erf(sqrt(adot/(2*Hdivide*kappa))*Hdivide.*(1 - zoh))) )) %final
%         T_question3 = Tb_k - Gpre * pre * (erf(sqrt(adot/(2*Hdivide*kappa))*Hdivide.*(zt)));

        Bt = @(zt) exp( T0 ./ (Ts + Gpre * pre * (erf1 - erf(sqrt(adot/(2*Hdivide*kappa))*Hdivide.*zt))) - c ./ (Tr - (Ts + Gpre * pre * (erf1 - erf(sqrt(adot/(2*Hdivide*kappa))*Hdivide.*zt)))).^k );
        tildeB = @(zoh) exp( T0./(Ts + Gpre * pre * (erf1 - (erf(sqrt(adot/(2*Hdivide*kappa))*Hdivide.*(1 - zoh))) )) - c./((Tr-(Ts + Gpre * pre * (erf1 - (erf(sqrt(adot/(2*Hdivide*kappa))*Hdivide.*(1 - zoh))) ))).^k) );
%         iBt_vec(i) = integral(Bt, 0, 1);
    end
% end new

%     intBT = integral(tildeB, 0, 1);
%     intBTtil_vec(i) = intBT;

%     tilBds = tildeB(tds); % please be cautious about z axis positive down
%     tilBdb = tildeB(1-tdb);
    for j=1:length(tdb)
%         expr = (rhoi/(rhow-rhoi))*tds - tdb(j) + 0.5*(rhoi/rhow)*(tilBdb(j) - tilBds)/intBTtil_vec(i);
        expr = (rhoi/(rhow-rhoi))*tds .* tildeB(1-tdb(j)) ./ tildeB(tds) - tdb(j);
        [minval, minloc] = min(abs(expr));
        if minval <= tol && tds(minloc) + tdb(j) <= 1
%             figure; plot(tds, expr, '-', 'linewidth', 4); xlabel("d_s/H [1]"); ylabel("Expr"); title('Roots of Expression?');
%             hold on; yline(0)
            tds_find(i,j) = tds(minloc);

%             plot(tds_find(j), 0, "o", 'markersize', 16);

            sigma0(i,j) = rhoi*g*tds_find(i,j) / tildeB(tds_find(i,j));

            % solving for beta
            betaterm_1(i,j) = (rhow/(rhow-rhoi))*tds_find(i,j).^2;
            betaterm_2(i,j) = (rhow/rhoi)*tdb(j).^2;
            betaterm_1_2(i,j) = betaterm_1(i,j) + betaterm_2(i,j);%(rhow/rhoi)*tdb(j).^2 + (rhow/(rhow-rhoi))*tds_find(i,j).^2;
            betaterm_3(i,j) = (2*sigma0(i,j)/(rhoi*g*(1-rhoi/rhow))) * integral(tildeB, tds_find(i,j), 1-tdb(j));
%             betaterm_3(i,j) = integral(tildeB, tds_find(i,j), 1 - tdb(j))/intBTtil_vec(i);
%             betaterm_4(i,j) = (1 - tdb(j) - tds_find(i,j)) * (2*(rhow/rhoi)*tdb(j) - tildeB(1-tdb(j))/intBTtil_vec(i));
%             beta(i,j) =  betaterm_1_2(i,j) + betaterm_3(i,j) + betaterm_4(i,j);
            beta(i,j) =  betaterm_1_2(i,j) + betaterm_3(i,j);
        elseif tds(minloc) + tdb(j) > 1
            disp(tds(minloc)) % surface crack depth
            disp(tdb(j)) % basal crack depth
            disp("Impossible Solution: ds+db > 1.")
            disp(i)
            disp(j)
            minval
            breakflag=1;
            break;
        else
            disp(i)
            disp(j)
            minval
            breakflag=1;
            break;
        end
    end
%     figure; plot(tds, expr, '-', 'linewidth', 4); xlabel("d_s/H [1]"); ylabel("Expr"); title('Roots of Expression?');
%     hold on; yline(0)
%     if breakflag == 0
        betamax(i) = max(beta(i,:));
%     end
end
toc;

colors = ["k", "r", "b", "c", "g", "m"];
figure; lge = [];
for ii=1; %:length(iind2)
    i=iind2(ii);
    hold on; 
%     plot(beta, tdb, '*', 'markersize', 14); % for yao
    plot(beta(i,:), tdb, colors(ii), 'linewidth', 4); % for you
%     betaboth = [beta(i,:), flip(beta(i,:))];
%     dbds = [tdb, flip(1-tds_find(i,:))];
%     plot(betaboth, dbds, '-', 'linewidth', 4); % for you    
    lg = append("$T_s$=", int2str(T(i)), "$^{\circ}$C"); lge = [lge, lg];
end
yline(rhoi/rhow, "c:", 'linewidth', 2);
for i=1:length(iind2)
    [maxbetaval, betaloc] = max(beta(iind2(i),:));
    plot(maxbetaval, tdb(betaloc), append(colors(i), "o"), 'markersize', 12);
    plot(maxbetaval, 1-tds_find(iind2(i),betaloc), append(colors(i), "o"), 'markersize', 12); % update 8/3/2023
    plot([maxbetaval, maxbetaval], [tdb(betaloc), 1-tds_find(iind2(i),betaloc)], append(colors(i), "--"), 'linewidth', 2); % update 8/3/2023
end
lg = "Sea Level"; lge = [lge, lg];
% Surface Crevasse Update 8/2/2023
for ii=1:length(iind2)
    i=iind2(ii);
    hold on;
    plot(beta(i,:), 1-tds_find(i,:), append(colors(ii), "-."), 'linewidth', 4); % 8/2/2023 Update
end
% 8/2/2023

xlabel("$\bar{R}_{xx}/\bar{R}_{xx}^{IT}$ [1]", Interpreter="latex");
% ylabel('$d_b/H$ [1]', Interpreter='latex');
ylabel('$z/H$ [1]', Interpreter='latex');
ylim([0 1])
xlim([0 2.1])
% title("Modified Nye's Zero-Yield Stress", Interpreter='latex')
% title("Modified Nye's", Interpreter='latex')
title("Nye's with $\Sigma F_x=0$", Interpreter='latex')
legend(lge, location='southeast', interpreter='latex')
set(gca, 'FontSize', 20)
toc;
figure; plot(T, betamax); ylim([0.9 1.1]); xlabel('Surface Temperature in Celsius'); ylabel("$\bar{R}_{xx}/\bar{R}_{xx}^{IT}$ [1]", Interpreter="latex");
% starts to drop down at -25 C, which you have seen before.
if LinearRobinMarine == 0
    title('Linear Temperature Profile');
elseif LinearRobinMarine == 1
    title('Robin Temperature Profile');
end
%% internal verification of what is going on
figure; hold on;
lge = [];
for ii=1:length(iind2)-1
    i=iind2(ii);
    hold on; 
%     plot(beta, tdb, '*', 'markersize', 14); % for yao
    colcount = (3*(ii-1));
    plot(betaterm_1(i,:), tdb, colors(colcount+1), 'linewidth', 4); % for you
    lg = append('$\beta_s$, $T_s=$', int2str(T(i)), "$^{\circ}$C"); lge = [lge, lg];
    plot(betaterm_2(i,:), tdb, colors(colcount+2), 'linewidth', 4);
    lg = append('$\beta_b$, $T_s=$', int2str(T(i)), "$^{\circ}$C"); lge = [lge, lg];
    plot(betaterm_3(i,:), tdb, colors(colcount+3), 'linewidth', 4);
    lg = append('$\beta_c$, $T_s=$', int2str(T(i)), "$^{\circ}$C"); lge = [lge, lg];
    plot(beta(i,:), tdb, append(colors(colcount+1), "--"), 'linewidth', 4);
    lg = append('$\beta$, $T_s=$', int2str(T(i)), "$^{\circ}$C"); lge = [lge, lg];
%     betaboth = [beta(i,:), flip(beta(i,:))];
%     dbds = [tdb, flip(1-tds_find(i,:))];
%     plot(betaboth, dbds, '-', 'linewidth', 4); % for you    
%     lg = append("$T_s$=", int2str(T(i)), "$^{\circ}$C"); lge = [lge, lg];
end
yline(rhoi/rhow, "k:", 'linewidth', 2);
lg = "Sea Level"; lge = [lge, lg];
ylabel('$d_b/H$ [1]', Interpreter='latex');
xlabel("$\beta$ and its terms [1]", Interpreter="latex");
title("Modified Nye's Zero-Yield Stress $\beta$ terms", Interpreter='latex')
leggy = legend(lge, location='best', interpreter='latex');
ylim([0 1])
set(gca, 'fontsize', 18, 'fontname', 'times');
% xlim([0 2.1])

% saveas(gcf, append(pwd, "/fig2b_components.fig"))
% saveas(gcf, append(pwd, "/fig2b_components.jpg"))
close;
%% save? Latest version 8/3/2023 Update
if LinearRobinMarine == 0
    saveas(gcf, append(pwd, "/fig2b_08032023.fig"))
    saveas(gcf, append(pwd, "/fig2b_08032023.jpg"))
    save("C:/Users/niall/OneDrive/Documents/YaoLab/Yao_Nature_Fractures/Data/Modified_Nye_rift_initiation_08032023.mat","T","betamax");
elseif LinearRobinMarine == 1
    saveas(gcf, append(pwd, "/Robinfig2b_08032023.fig"))
    saveas(gcf, append(pwd, "/Robinfig2b_08032023.jpg"))
    save("C:/Users/niall/OneDrive/Documents/YaoLab/Yao_Nature_Fractures/Data/Modified_Nye_Robin_rift_initiation_08032023.mat","T","betamax");
end
%% zoom in
set(leggy, 'visible', 'off');
xlim([0.98 1.02])
ylim([0.6 0.92])
if LinearRobinMarine == 0
    saveas(gcf, append(pwd, "/fig2b_08032023_zoom.fig"))
    saveas(gcf, append(pwd, "/fig2b_08032023_zoom.jpg"))
elseif LinearRobinMarine == 1
    saveas(gcf, append(pwd, "/Robinfig2b_08032023_zoom.fig"))
    saveas(gcf, append(pwd, "/Robinfig2b_08032023_zoom.jpg"))
end
%% save 8/2/2023 Update
if LinearRobinMarine == 0
    saveas(gcf, append(pwd, "/fig2b_08022023.fig"))
    saveas(gcf, append(pwd, "/fig2b_08022023.jpg"))
    % save("C:/Users/niall/OneDrive/Documents/YaoLab/Yao_Nature_Fractures/Data/Modified_Nye_rift_initiation.mat","T","betamax");
elseif LinearRobinMarine == 1
    saveas(gcf, append(pwd, "/Robinfig2b_08022023.fig"))
    saveas(gcf, append(pwd, "/Robinfig2b_08022023.jpg"))
    save("C:/Users/niall/OneDrive/Documents/YaoLab/Yao_Nature_Fractures/Data/Modified_Nye_Robin_rift_initiation_08022023.mat","T","betamax");
end
%% LEFM
% going to do this on della no longer! just going for it here. First cell
% takes 5 minutes

% please run cell 1 of this script (Nye's) to get the ordinary variables back

clearvars Rxx stress_from_vert_traction ND_stress_from_vert_traction exx % added exx 7/3/2023

tic;
% exx = reshape(exx, [length(Ts_c), length(beta)]);
K_Ic = 150000;

Ts = flip(Ts_c);

vert_grad_T = 1;
horiz_grad_T = 1;
depth_avg_T = 0;

dbt_takeaway = nan(length(Ts_c), length(beta(1:ceil(length(beta)/2))));

% calculate weight functions
A = 3.1416^0.5; B = 2/3.1416^0.5; rhoi=917; rhow=1028; g=9.8; 
%   was using this one  x=[0,logspace(-5,log10(0.9),100)]; x=[x,logspace(log10(0.9),log10(0.99),10)]; x=[x,logspace(log10(0.99),log10(0.996),10)]; x=[x,logspace(log10(0.996),log10(1),10)];% 
x=[0,logspace(-5,log10(0.9),100)]; x=[x,logspace(log10(0.9096),log10(0.99),9)]; x=[x,logspace(log10(0.9907),log10(0.996),9)]; x=[x,logspace(log10(0.9964),log10(1),9)];% 
x = [x, linspace(0.01, 0.99, 99)]; x = sort(x); % added every percent fracture depth

% 4/7/2023
x = 0:1E-3:1;
%     edited to make sure there are no repeat values so you can interpolate
syms y   % x=ds/H; y=z/ds;
FF= (2./x/pi.*tan(pi.*x/2)).^0.5.*(0.752 + 2.02.*x + 0.37*(1 - sin(pi.*x/2)).^3)./cos(pi.*x/2);
G1= 0.46 + 3.06.*x + 0.84.*(1 - x).^5 + 0.66.*x.^2.*(1 - x).^2;
G2= -3.52.*x.^2;
G3= 6.17 - 28.22.*x + 34.54.*x.^2 - 14.39.*x.^3 - (1 - x).^1.5 - 5.88.*(1 - x).^5 - 2.64.*x.^2.*(1 - x).^2;
G4= -6.63 + 25.16.*x - 31.04.*x.^2 + 14.41.*x.^3 + 2.*(1 - x).^1.5 +5.04.*(1 - x).^5 + 1.98.*x.^2.*(1 - x).^2;
G= G1 + G2.*y + G3.*y.^2 + G4.*y.^3; 
    % surface crack  % x=ds/H; y=z/ds;
    %ff= int((y).*G./((1 - x).^1.5.*(1 - y.^2).^0.5),y, 0, 1); %weight from ice
    %ff1= int((y-0.1).*G./((1 - x).^1.5.*(1 - y.^2).^0.5),y, 0.1, 1); %when water is 90% full%
    
    % bottom crack   % x=db/H; y=z/db;
fbi= int((1./x-y).*G./((1 - x).^1.5.*(1 - y.^2).^0.5),y,0,1);
toc;
%     fbw= int((rhoi/rhow./x-y).*G./((1 - x).^1.5.*(1 - y.^2).^0.5),y,0,1);
    %%%%%%% YOU ARE NOW CHANGING THE CODE TO YOUR OWN FORM BELOW
fbw = int((rhoi/rhow./x-y).*G./((1 - x).^1.5.*(1 - y.^2).^0.5),y,0,1);
toc;
for i=1:length(fbw) % x = d_b/H, so this should be 0 above sea level
    sealevel_ind = 0;
    if x(i) > rhoi/rhow
    
            % recalculate the G's
        G1p= 0.46 + 3.06*x(i) + 0.84.*(1 - x(i))^5 + 0.66*x(i)^2*(1 - x(i))^2;
        G2p= -3.52*x(i)^2;
        G3p= 6.17 - 28.22*x(i) + 34.54*x(i)^2 - 14.39*x(i)^3 - (1 - x(i))^1.5 - 5.88*(1 - x(i)).^5 - 2.64*x(i)^2*(1 - x(i))^2;
        G4p= -6.63 + 25.16*x(i) - 31.04*x(i)^2 + 14.41*x(i)^3 + 2*(1 - x(i))^1.5 +5.04*(1 - x(i))^5 + 1.98*x(i)^2*(1 - x(i))^2;
        Gp= G1p + G2p.*y + G3p.*y.^2 + G4p.*y.^3; 
    
        sealevel_ind = sealevel_ind + 1;
        upper_bound = rhoi/rhow/x(i); % again x = db/H
        fbw(i) = int((rhoi/rhow/x(i)-y)*Gp/((1 - x(i))^1.5.*(1 - y.^2).^0.5),y,0,upper_bound);
        i
    end
end
toc;

%% got most of the weight functions

% very high numerical resolution means this takes a while (30 min?)

Ks=[]; bpk=[]; xx=[];
if vert_grad_T == 0
    unstable = []; FbT_list = []; Ks_list = []; optimal_x = []; FbT_vec_list = []; gb_vec_list = []; fb_vec_list = [];
else
    T_min = floor(min(min(Ts))); T_max = ceil(max(max(Ts)));
    T_count = abs(T_max - T_min) + 1;
    unstable = [];%ones(T_count, 1); don't need the 1 there
    FbT_list = ones(T_count, 1); 
    Ks_list = ones(T_count, 1); 
    optimal_x = ones(T_count, 1); 
    FbT_vec_list = [];%ones(T_count, 1);
    gb_vec_list = [];%ones(T_count, 1); 
    fb_vec_list = [];%ones(T_count, 1); 
        
end    
% choose stress range
n=30;

T_min = round(min(min(Ts))); T_max = round(max(max(Ts)));
T_count = abs(T_max - T_min) + 1;
counter = 0;
for l=1:T_count
    Ks=[]; bpk=[]; xx=[];
    Ts_iter = T_min + l - 1;

    breakflagtop = 1;
    for j=1:length(Ts_c)
        if Ts_iter == Ts_c(j)
            breakflagtop = 0;
            counter = counter + 1; % to have an actual index for exx
        end
    end
    if breakflagtop == 1
        disp(int2str(Ts_iter))
        disp('brk');
        continue;
    end
    disp(Ts_iter)

    % temperature gradient
    if horiz_grad_T == 1 && vert_grad_T == 1
        Tb = Tb_choice; % this is where you're passing basal temperature from func defn, test -2 C
    elseif horiz_grad_T == 1 && vert_grad_T == 0
        if depth_avg_T == 0
            Tb = Ts_iter;
        else % if you are looking at depth-averaged tension case!!
            Tb = Tb_choice;
        end
    end

    if horiz_grad_T == 1 && vert_grad_T == 0
        xx = logspace(-3,-2,n); xx = [xx,linspace(0.0111,0.0416,n-1)];   %R_bar basal crack%
        %xx = logspace(-3,-2,n); xx = [xx,linspace(0.02, 0.1,n)];   %R_bar basal crack + DT%

    elseif horiz_grad_T == 1 && vert_grad_T == 1 % these are known from playing around
        xx = logspace(-3, -2, n); % for making the actual phase space plots with data
        if Ts_iter >= -10+Tb % works because og for Tb=0, now lin temp grad should allow us to add
            xx = [xx, linspace(0.015, 0.035, 5)];
            xx = [xx, linspace(0.04, 0.05, n)]; % this is confirmation that you are right - xx sample range increasing w/ incr |T|.

        elseif Ts_iter >= -14+Tb && Ts_iter < -10+Tb
            xx = [xx, linspace(0.015, 0.04, 6)];
            xx = [xx, linspace(0.045, 0.055, n)];

        elseif Ts_iter >= -17+Tb && Ts_iter < -14+Tb
            xx = [xx, linspace(0.015, 0.045, 7)];
            xx = [xx, linspace(0.05, 0.06, n)];

        elseif Ts_iter >= -20+Tb && Ts_iter < -17+Tb
            xx = [xx, linspace(0.015, 0.05, 8)];
            xx = [xx, linspace(0.055, 0.065, n)];

        elseif Ts_iter >= -22+Tb && Ts_iter < -20+Tb
            xx = [xx, linspace(0.015, 0.055, 9)];
            xx = [xx, linspace(0.06, 0.07, n)];

        elseif Ts_iter >= -24+Tb && Ts_iter < -22+Tb
            xx = [xx, linspace(0.015, 0.06, 10)];
            xx = [xx, linspace(0.065, 0.075, n)];
                
        elseif Ts_iter >= -26+Tb && Ts_iter < -24+Tb
            xx = [xx, linspace(0.015, 0.065, 11)];
            xx = [xx, linspace(0.07, 0.08, n)];

        elseif Ts_iter >= -28+Tb && Ts_iter < -26+Tb
            xx = [xx, linspace(0.015, 0.07, 12)];
            xx = [xx, linspace(0.075, 0.085, n)];

        elseif Ts_iter >= -31+Tb && Ts_iter < -28+Tb
            xx = [xx, linspace(0.015, 0.075, 13)];
            xx = [xx, linspace(0.08, 0.09, n)];

        else
            xx = [xx, linspace(0.015, 0.08, 14)];
            xx = [xx, linspace(0.085, 0.095, n)];
        end
    end

%     overwrite the whole thing haha
% 4/5/2023
%     xx = (2*intB_vec(end+1-counter)*exx(length(Ts_c)+1-counter, 1:ceil(length(beta)/2)).^(1/3)) ./ (rhoi*g*H); % depth-avg Rxx. going to update Bs below to make it work
%     xx = 0.5*(1-rhoi/rhow)*beta(1:ceil(length(beta)/2));
    xx = 0.5*(1-rhoi/rhow)*[0:0.01:2]; % this is just numerical sampling, don't read into this.

    %figure('visible','off'); hold on;
    figure; hold on;
    % time to see if it is larger than K_Ic tilde
    K_Ic_tilde = K_Ic / (rhoi*g*H^(3/2));
    plot(x, K_Ic_tilde*ones(size(x)), 'r-', 'linewidth', 2);
%     xlabel('$d_b/H$ [1]', Interpreter='latex'); 
    xlabel('$z/H$ [1]', Interpreter='latex'); 
    ylabel('$K_{Ic}/\rho_i g H^{3/2}$ [1]', Interpreter='latex');
    for j=1:length(xx) % for each ND extensional stress
        j
        a = xx(j);  %R_bar%
            
            % %%%%%%%%%%%%% Considering Vertical Variation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if depth_avg_T == 0
            ffirn=[0]; xfirn=[0]; FBfirn=[0]; FsT_=[0]; FbT_=[0]; FB0firn=[0]; 
            for xi=1e-5:0.001:1
                G1= 0.46 + 3.06.*xi + 0.84.*(1 - xi).^5 + 0.66.*xi.^2.*(1 - xi).^2;
                G2= -3.52.*xi.^2;
                G3= 6.17 - 28.22.*xi + 34.54.*xi.^2 - 14.39.*xi.^3 - (1 - xi).^1.5 - 5.88.*(1 - xi).^5 - 2.64.*xi.^2.*(1 - xi).^2;
                G4= -6.63 + 25.16.*xi - 31.04.*xi.^2 + 14.41.*xi.^3 + 2.*(1 - xi).^1.5 +5.04.*(1 - xi).^5 + 1.98.*xi.^2.*(1 - xi).^2;

                xfirn=[xfirn,xi];
                
                
                Bs = 2.207*exp(3155./(Ts_iter+273.15)-0.16612./((0.24-Ts_iter).^1.17)); % When Ts=?
%                 Bs = intB_vec(end+1-counter);% 2.207*exp(3155./(Ts_iter+273.15)-0.16612./((0.24-Ts_iter).^1.17)); % depth-average 4/6/2023

                fun4 = @(y) ((2.207*exp(3155./((xi.*y.*(Ts_iter-Tb)+Tb)+273.15)-0.16612./((0.24-(xi.*y.*(Ts_iter-Tb)+Tb)).^1.17)))./Bs.*(G1 + G2.*y + G3.*y.^2 + G4.*y.^3)./((1 - xi).^1.5.*(1 - y.^2).^0.5)); %temperature from Hooke
                q4 = integral(fun4,0,1);
                FbT_=[FbT_,q4];

            end
                % firn = interp1(xfirn,ffirn,x);
                % Ffirn = interp1(xfirn,FBfirn,x);
                % FsT = interp1(xfirn,FsT_,x);
            FbT = interp1(xfirn,FbT_,x);

            Ks  = B*a*x.^0.5.*FbT - B*x.^1.5.*fbi + rhow/rhoi*B.*x.^1.5.*fbw;   %basal crack + Temperature gradient (unstable at xx = 0.104)

        elseif depth_avg_T == 1
            Ks  = A*a*x.^0.5.*FF - B*x.^1.5.*fbi + rhow/rhoi*B.*x.^1.5.*fbw;   %basal crack
        end

        
        Kdiff = Ks - K_Ic_tilde;
        [maxKval, maxKloc] = max(Kdiff);
        Kdiffpos = Kdiff;
        Kdiffpos(Kdiffpos < 0) = 1E10;

        if maxKval < 0
            dbt_takeaway(counter, j) = 0; % no fracture
        else

%             dK_I_dx = diff(Ks(3:end-3)) ./ diff(x(3:end-3)); % find where this is 0, sample values >= than this
        % while derivative is 0 or negative, find the intersection between SIF and fracture toughness.
%             [mindK, mindKloc] = min(abs(dK_I_dx)); % want slope to be 0

            [mindK, mindKloc] = min(Kdiffpos);
            
            if mindKloc == 1
                min_db_tilde = x(mindKloc);
            else
                min_db_tilde = x(mindKloc-1); % placeholder for now, ~d before maxima, eventually ~d "at" maxima.
            end

            for ll=mindKloc:length(x)
                if Ks(ll) >= K_Ic_tilde % if we're above fracture toughness
                    min_db_tilde = x(ll); % stepping forwards
                else
                    break;
                end
            end
            if min_db_tilde > 0.99
                dbt_takeaway(counter,j) = 1; % full thickness!
            else
                dbt_takeaway(counter, j) = min_db_tilde;
                Ks_circle = Ks(ll-1); % minus 1 because current ll broke
            end
            
        end

        plot(x, Ks, 'linewidth', 4); 
        set(gca, 'xscale', 'log');
        set(gca, 'yscale', 'log');
        if dbt_takeaway(counter,j) ~= 0 && dbt_takeaway(counter,j) ~= 1
            plot(dbt_takeaway(counter,j), Ks_circle, "o", 'markersize', 12);
        end

        if dbt_takeaway(counter, j) == 0
            disp('No Fracture');
        elseif dbt_takeaway(counter, j) == 1
            disp('Unstable Fracture');
%             disp(exx(length(Ts_c)+1-counter, j));
            disp(j);
        else
            disp(append('Stable Fracture, d_b/H = ', num2str(dbt_takeaway(counter,j))));
        end
        % time to save 

    end
    title(append('$T_s = $', int2str(Ts_iter)), Interpreter='latex');

end
% okay let's get it! 4/7/2023 Shoutout Forester for Undercover
Ts_k = Ts+273.15;
counter = 0;
for l=1:T_count
    Ts_iter = T_min + l - 1;

    breakflagtop = 1;
    for j=1:length(Ts_c)
        if Ts_iter == Ts_c(j)
            breakflagtop = 0;
            counter = counter + 1; % to have an actual index for exx
        end
    end
    if breakflagtop == 1
        disp(int2str(Ts_iter))
        disp('brk');
        continue;
    end
    disp(Ts_iter)
    Bfz = @(zt) B0*exp( T0 ./ (zt * (Ts_k(counter) - Tb_k) + Tb_k) - c ./ (Tr - (zt * (Ts_k(counter) - Tb_k) + Tb_k)).^k ); % linear
    intB_vec(counter) = integral(Bfz, 0, 1);
    beta_for_LEFM(counter,:) = intB_vec(counter) * xx ./ (Bfz(1) *0.5*(1-rhoi/rhow) ); % B at surface
end
% beta_for_LEFM = nan(size(xx)); % in the end you'll make it iBt*xx/Bt(surface), or Bbar * Rxx(surf) / B(surf).

% find when these first equal 1, and then index location before that draw
% the circle and dashed line
ubs = nan(size(Ts)); % upper bounds
for i=1:length(Ts)
    iind3 = length(Ts)+1-i;
    for j=1:length(xx)
        if dbt_takeaway(iind3,j) == 1
            ubs(i) = j; % include the first 1, don't include the 2nd 
            break;
        end
    end
end

colors = ["k", "r", "b", "c", "g", "m"];
figure; lge = [];
hold on;
for i=1:length(Ts)
    iind3 = length(Ts)+1-i;
    plot(beta_for_LEFM(iind3,1:ubs(i)), dbt_takeaway(iind3, 1:ubs(i)), colors(i), 'linewidth', 4);
    lg = append("$T_s$=", int2str(Ts(iind3)), "$^{\circ}$C"); lge = [lge, lg];
end
yline(rhoi/rhow, "c:", 'linewidth', 2);
lg = "Sea Level"; lge = [lge, lg];

% for i=1:length(Ts) % would use this if you want circles and dashed lines
%     asdf
% end

xlabel("$\bar{R}_{xx}/\bar{R}_{xx}^{IT}$ [1]", Interpreter="latex");
ylabel('$d_b/H$ [1]', Interpreter='latex');
ylim([0 1])
xlim([0 2.1])
title("LEFM", Interpreter='latex')
legend(lge, location='southeast', interpreter='latex')
set(gca, 'FontSize', 20)
%% when updating old figure
set(gca, 'FontSize', 20)
h = get(gca, 'Children');
set(h(3), 'Color', 'r');
set(h(2), 'Color', 'b');
%% save?
saveas(gcf, append(pwd, "/fig2c.fig"))
saveas(gcf, append(pwd, "/fig2c.jpg"))
close

% Calculation of Rxx
function Rxx = recalc_stress(T, exx, H) % this T should be in Kelvin
    %exx was furst's strain rate
    %now, exx is whichever strain rate you are using. Wearing or Furst.

    rhoi=917; g=9.8; %H=300;
    Bhooke_surf=2.207*exp(3155./(T)-0.16612./((273.39-T).^1.17)); %Kelvin, Hooke temperature, T=RACMO
    Rxx = 2.*Bhooke_surf.*abs(exx).^(1/3)/rhoi/g./H;   Rxx(exx<=0) = -Rxx(exx<=0);
    % this is already non-dimensionalized
end