function [dimless_crev_depth, stress_intensity] = get_dimless_crev_depth(pp, anonT)
dz = 1e-2; 
db = 1e-2:1e-2:(1-1e-2);  %grid for the stress intensity algoritm

[stress_intensity,~,~,~] = get_stress_intensity(dz,db,pp.lambda,anonT);


dbc = intersections([0,1], [pp.F,pp.F], db, stress_intensity);

if isempty(dbc) %no root
    if all(stress_intensity > pp.F)
    dimless_crev_depth = 1;
    else 
    dimless_crev_depth = 0;
    end

elseif length(dbc) == 1 %must be increasing so only one root
    dimless_crev_depth = 1;
elseif length(dbc) == 2
    dimless_crev_depth = max(dbc);
end
end
