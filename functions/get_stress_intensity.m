function [stress_intensity,term1,term2,term3] = get_stress_intensity(dz,db,lambda,T)
%Return the stress intensity as a function of db. NB: T must be input in Kelvin

%constants
n = 3; %glen flow coef
T0 = 3155;
Tr = 273.39;
C = 0.166;
k = 1.17;

zz = dz:dz:(1-dz);

% initialize terms in stress intensity
term1 = zeros(1,length(db));
term2 = zeros(1,length(db));
term3 = zeros(1,length(db));

    for idb = 1:length(db)
     %   G = zeros(size(zz));
     %   for iG = 1:length(G)
     %       G(iG) = get_G(zz(iG),db(idb));
     %   end
       G = get_G(zz,db(idb));

        integrand = G .* (1 - db(idb)).^(-3/2) .* (1 - zz.^2) .^(-1/2);
       % BB = zeros(size(zz));
        BB = exp(T0./T(zz) - C./(Tr - T(zz)).^k);
       %  for iz = 1:length(zz)
       %      BB(iz) = exp(T0/T(zz(iz)) - C/(Tr - T(zz(iz)))^k);
       % end
        intgrl = sum(BB .* integrand * dz);

        term1(idb) = 2/sqrt(pi) * lambda* db(idb)^(1/2) * intgrl; %with depth dependent viscosity
        term2(idb) = -2/sqrt(pi) * db(idb)^(3/2) *fb(db(idb));    
        term3(idb) = 2/sqrt(pi) * 1028/918 * db(idb)^(3/2) * gb(db(idb));
        
%                 %surface crevasses
%         term1(idb) = sqrt(pi)*Rxx/(rhoi*grav*H)*db(idb)^(1/2)* F(db(idb));
%         term2(idb) = -2/sqrt(pi) * db(idb)^(3/2) *fs(db(idb));
%         term3(idb) = 2/sqrt(pi) * 1028/918 * db(idb)^(3/2) * gs(db(idb));
        

    end
    stress_intensity = term1 + term2 + term3;
end