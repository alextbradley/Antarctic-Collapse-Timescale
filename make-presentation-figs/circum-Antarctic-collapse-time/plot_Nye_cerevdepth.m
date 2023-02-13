% compute the current crevasse depth over the whole of antarctica.
clear
f = load('../../data/ice_sheet_data.mat');

%% compute the viscosity
dz = 1e-2; 
zz = dz:dz:(1-dz);
idx = find(~isnan(f.H)); %non-zero thickness points
zzM = repmat(zz', [1,size(idx)]); % array with rows corresponding to depths and columns to non-zero points
hidx = f.H(idx)';
hM  = repmat(hidx, [length(zz), 1]); %icethicknesses 

Ts = -18 + 273.15; %surface temp
Tb = -2 + 273.15;

ghf = 48; %geothermal heat flux


k = 1.18*1e-6;
Adot = 3.17*1e-9;
K = 5.3e-3*4184 * 10^2;
dthetadhbot = ghf/K;

L = sqrt(2*hM*k/Adot); %lengthscale

%grounding line profile 
Tgf = Ts + dthetadhbot*L.*(erf(hM./L) - erf(zzM.*hM./L))*sqrt(pi)/2;

Tgf(Tgf > Tb) = Tb; %effectively pressure melting
%% boundary layer profile
kappa = 36;
midx = f.m(idx)';

mM  = repmat(midx, [length(zz), 1]); %ice melt rate
mM(mM < 1e-3) = 1e-3;
l = kappa./hM./mM;
T = Tgf +  (Tb - Tgf).*exp(-zzM./l);

n = 3; %glen flow coef
T0 = 3155;
Tr = 273.39;
C = 0.166;
k = 1.17;
B = exp(T0./T - C./(Tr - T).^k); %ice viscosity with depth 

% depth average the viscoity
Bbar = 2*sum(B,1)*dz;

%%
strainidx = f.eflow(idx)';
Rxx = 2*Bbar.* sign(strainidx).*(abs(strainidx)).^(1/n);
rhoi = 918.0;
rhow = 1028.0;
g = 9.18; 
db = Rxx / g / (rhow - rhoi);

%% repopulate the matrix
db(db < 0) = 0; %negative strain
crev_dep = nan(size(f.H));
crev_dep(idx) = db;
%crev_dep(crev_dep) >


dcv = crev_dep ./ (f.H);
dcv(dcv > 1) = 1;
dcv(dcv < 0) = 1e-3; 

figure(1); clf; p = imagesc((dcv)); 
set(p, 'AlphaData', ~isnan(dcv));
colorbar
ax = gca;

%% timescale
tcrev = (f.H - crev_dep)./abs(f.dhdtadj);
tcrev(tcrev < 0) = 0;
figure(2); p = imagesc(log10((tcrev)));
set(p, 'AlphaData', ~isnan(tcrev));
c = colorbar;
ax2 = gca;
%ax2.Colormap = ax.Colormap;
ax2.Colormap = cmocean('Matter');
clim([1,4])
axis equal
c.Ticks = 1:4;
c.FontName = 'GillSans';
c.FontSize = 16;
ax2.Visible = 'off';
