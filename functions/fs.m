function f = fs(d)
    dz = 1e-3;
    zz = dz:dz:(1-dz);
    G = zeros(size(zz));
    for iG = 1:length(G)
        G(iG) = get_G(zz(iG),d);
    end
    f = sum(zz .*  G .* (1 - d)^(-3/2) .* (1- zz.^2).^(-1/2))*dz;
end
