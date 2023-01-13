function g = gs(d)
  
    dw = 0.1;
    zz = linspace((1-dw),(1-1e-3));
    dz = diff(zz); dz = dz(1);
    G = zeros(size(zz));
    for iG = 1:length(G)
        G(iG) = get_G(zz(iG),d);
    end
    g = sum((zz- (1-dw)) .*  G .* (1 - d)^(-3/2) .* (1- zz.^2).^(-1/2))*dz;
end