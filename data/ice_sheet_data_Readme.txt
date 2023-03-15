===================================
The attached file includes

T: unit [K] temperature from RACMO, same as "T_mask"
eflow: unit [1/year] along-flow strain rate from Wearing
strain: unit [1/year] along-flow strain rate from Furst; strain = (along_flow_stress./2./B).^3;
Rxxs: unit [NaN] surface resistive stress, Furst's strain rate + B(T), calculated details below*
Rxx: unit [NaN] resistive stress, Furst's along-flow stress from ice-flow inversion
B: unit [Pa year^1/3] ice viscosity from Furst's inversion, same as "visc" from Furst
H: unit [m] ice thickness, same as "thick_mask"
S: unit [m] ice surface slope
u: unit [m/year] ice velocity along the x direction
v: unit [m/year] ice velocity along the y direction
m: unit [m/year] basal melt rate from Adusumilli et al. 2020 (doi: 10.1038/s41561-020-0616-z)
tf_max: unit [C] maximum thermal forcing between 200 and 800m in the water column Adusumilli et al. 2020 (doi: 10.1038/s41561-020-0616-z)
bed: unit [m] bed depth from Bedmachine V3
x: unit [m] x position
y: unit [m] y position
xx: unit [m] x position (in 1 x nx array)
yy: unit [m] y position (in 1 x ny array)
dhdt: unit [m/yr] thinning rate dh/dt from Smith et al. 2020 (doi: 10.1126/science.aaz5845), downscaled from 5km grid of Smith et al. to 1km grid
dhdtadj: unit [m/yr] thinning rate 'dhdt' from Smith et al. 2020 (doi: 10.1126/science.aaz5845), downscaled from 5km grid of Smith et al. to 1km grid. 
                     WAIS section replaced by more complete data provided by B. Smith (pers comm).
dmdt: unit [m/yr] thinning rate 'dmdt' from Smith et al. 2020 including firn-air correction. 
dmdtadj: unit [m/yr] thinning rate 'dmdt' from Smith et al. 2020 including firn-air correction with ad-hoc WAIS correction. 
dmdterr: unit [m/yr] RMSE dmdt error from Smith et al. 2020.

MOA2009_1km: unit [NaN] satellite image
bb0448974g_2_1.h5

=== * Calculation of Rxxs ========
exx is strain rate from Furst
rhoi=917; rhow=1028; g=9.8; 
Bhooke_surf=2.207*exp(3155./(T)-0.16612./((273.39-T).^1.17)); %Hooke temperature, T=RACMO (in Kelvin)
Rxx = 2.*Bhooke_surf.*abs(exx).^(1/3)/rhoi/g/H;   Rxx(exx<=0) = -Rxx(exx<=0);


===================================
To plot in matlab, use:

imagesc(variable_name); caxis([minvalue maxvalue])
