function dimless_crev_depth = get_dimless_crev_depth_ModifiedNye(pp, anonT)
%return the dimensionless crevasse depth for the modified Nye case
dz = 1e-2;
%search the space
[Rxx_nondim, ds, db] = spacesearch_ModifiedNye(anonT);

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
g = 9.81;
RxxIT = 0.5 * (1-rhoi/rhow)*rhoi*g*pp.H0;

if Rxx > RxxIT
    dimless_crev_depth = 1;
else %find interxections 
    [~,idx] = min(abs(Rxx_nondim - Rxx/RxxIT)); 
    dimless_crev_depth = ds(idx) + db(idx);

end