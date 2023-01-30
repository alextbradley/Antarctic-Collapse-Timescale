function lambda = get_critical_lambda(pp)
% Apply a bisection method to determine the critical strain rate as a
% function of input parameters.
db = (1e-2):1e-2:(1-1e-2);
dz = 1e-2;
anonT = @(z) (pp.Ts + (pp.Tb - pp.Ts)*exp(-z/pp.l));
nitermax = 30; %maximum number of iterations in the bisection
errmax = 1e-11;


% Set initial guess and seed the bound booleans
lb = 1e-9;
is_lb = 0;
ub = 1e-6;
is_ub = 0; %

% Check that they are indeed lower and upper bounds and adjust if not
while ~is_lb
    [stress_intensity,~,~,~] = get_stress_intensity(dz,db,lb,anonT);

    dbc = intersections([0,1], [pp.F,pp.F], db, stress_intensity);
    if ((length(dbc) == 1) || all(stress_intensity > pp.F)) %single root or always bigger than stress intensity
        %= have upper bound
        lb = lb/10;
      
    else %
        is_lb = 1;
    end
    
end

while ~is_ub
    [stress_intensity,~,~,~] = get_stress_intensity(dz,db,ub,anonT);

    dbc = intersections([0,1], [pp.F,pp.F], db, stress_intensity);
    if ((length(dbc) == 1) || all(stress_intensity > pp.F)) %single root or always bigger than stress intensity
        %= have upper bound
        is_ub = 1;
      
    else %
        ub = ub*10;
    end
end

%ub, lb

% Now bisect
err = ub - lb;
niter = 0;
while ( err > errmax) && (niter < nitermax)
    guess = (ub + lb)/2;
     [stress_intensity,~,~,~] = get_stress_intensity(dz,db,guess,anonT);

    dbc = intersections([0,1], [pp.F,pp.F], db, stress_intensity);
    if ((length(dbc) == 1) || all(stress_intensity > pp.F)) %single root or always bigger than stress intensity
        %= have upper bound
        ub = guess;
      
    else %
        lb = guess;
    end
    niter = niter + 1;
    ub;
    lb;
    err = ub - lb;
end

lambda = guess;
end