function Tgf = get_grounding_line_temp(ghf, Ts,H)
% return an anoymous function of grounding line temperature with geothermal
% heat flux ghf and surface temperature Ts

k = 1.18*1e-6;
Adot = 3.17*1e-9;
K = 5.3e-3*4184 * 10^2;
dthetadhbot = ghf/K;
%dthetadhbot = 0.1;


L = sqrt(2*H*k/Adot); %lengthscale
Tgf = @(z) Ts + dthetadhbot*L*(erf(H/L) - erf(z*H/L))*sqrt(pi)/2;
