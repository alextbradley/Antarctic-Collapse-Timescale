function [Rxx_found, ds_found, db] = spacesearch_ModifiedNye(anonT)
% Return the dimensionless dimensionless resistive stress Rxx_found and
% corresponding surface and basal crevasse depths ds and bd for the
% Modified Nye theory (Coffey et al. 2023)

%constants
T0 = 3155;
Tr = 273.39;
C = 0.166;
k = 1.17;
rhoi = 918.0;
rhow = 1028.0;
g = 9.81;


% compute depth averaged viscosity 
T = anonT;
dz = 1e-2; 
zz = dz:dz:(1-dz);
B = exp(T0./T(zz) - C./(Tr - T(zz)).^k);
Bf = @(z) exp(T0./T(z) - C./(Tr - T(z)).^k); %anonymous functional form of B
%Bbar= 2*sum(B)*dz;


%depth averaged resitive stress
%Rxx = 2*Bbar.* sign(pp.epsxx).*(abs(pp.epsxx)).^(1/n);
%Rxx_IT = 0.5*(1-rhoi/rhow)*rhoi*g*pp.H; 

tol = 1e-5;
tds = [0:0.000001:(1-rhoi/rhow)+0.05, 1-rhoi/rhow+0.05]; 
tdb = [0:0.01:rhoi/rhow, rhoi/rhow];
tds_find = nan(1, length(tdb));
beta = nan(1, length(tdb));
Bf = @(z) exp(T0./(anonT(1-z))-C./((Tr-(anonT(1-z))).^k)); % fxn of z/H

%Bf = @(zoH) exp(3155./((zoH.*(Tb-Ts)+Ts)+273.15)-0.16612./((0.24-(zoH.*(Tb-Ts)+Ts)).^1.17)); % fxn of z/H

for j=1:length(tdb)
    %         expr = (rhoi/(rhow-rhoi))*tds - tdb(j) + 0.5*(rhoi/rhow)*(tilBdb(j) - tilBds)/intBTtil_vec(i);
    expr = (rhoi/(rhow-rhoi))*tds .* Bf(1-tdb(j)) ./ Bf(tds) - tdb(j);
    [minval, minloc] = min(abs(expr));
    if minval <= tol && tds(minloc) + tdb(j) <= 1

        tds_find(j) = tds(minloc);

        %             plot(tds_find(j), 0, "o", 'markersize', 16);

        sigma0 = rhoi*g*tds_find(j) / Bf(tds_find(j));

        % solving for beta
        betaterm_1 = (rhow/(rhow-rhoi))*tds_find(j).^2;
        betaterm_2= (rhow/rhoi)*tdb(j).^2;
        betaterm_1_2 = betaterm_1 + betaterm_2;%(rhow/rhoi)*tdb(j).^2 + (rhow/(rhow-rhoi))*tds_find(i,j).^2;
        betaterm_3 = (2*sigma0/(rhoi*g*(1-rhoi/rhow))) * integral(Bf, tds_find(j), 1-tdb(j));
    beta(j) =  betaterm_1_2 + betaterm_3;
    elseif tds(minloc) + tdb(j) > 1

        break;
    else

        break;
    end
end

db = tdb;
ds_found = tds_find;
Rxx_found = beta;