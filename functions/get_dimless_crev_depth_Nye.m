function dimless_crev_depth = get_dimless_crev_depth_Nye(pp, anonT)
%return the dimensionless crevasse depth for the Nye crack.
dz = 1e-2; 

%constants
n = 3; %glen flow coef
T0 = 3155;
Tr = 273.39;
C = 0.166;
k = 1.17;

% compute depth averaged viscosity 
T = anonT;
zz = dz:dz:(1-dz);
B = exp(T0./T(zz) - C./(Tr - T(zz)).^k);
Bbar= 2*sum(B)*dz;


Rxx = 2*Bbar.* sign(pp.epsxx).*(abs(pp.epsxx)).^(1/n);
rhoi = 918.0;
rhow = 1028.0;
g = 9.18; 
dimless_crev_depth = Rxx / g / (rhow - rhoi) / pp.H0;

end
