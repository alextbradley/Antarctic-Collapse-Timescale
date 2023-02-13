function g = gb(d)
    dz = 1e-2;
    zz = dz:dz:(1-dz);
%     G = zeros(size(zz));
%     for iG = 1:length(G)
%         G(iG) = get_G(zz(iG),d);
%     end
    G =  get_G(zz,d);
    g = sum((918/1028 /d - zz) .*  G .* (1 - d)^(-3/2) .* (1- zz.^2).^(-1/2))*dz;
end
