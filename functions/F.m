function FF = F(d)
    FF = sqrt(2/pi/d * tan(pi/ 2 * d))* ((0.752 + 2.02*d + 0.37*(1-sin(pi/2 *d))^3)/cos(pi/2 *d));

end