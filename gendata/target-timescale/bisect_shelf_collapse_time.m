function dM = bisect_shelf_collapse_time(shelf, step, target_time,f)

%return the melt increase dM required to achieve a median shelf collapse
%time of target_time. Step is the grid step. f is the input ice sheet data

% do the bisection for the current conditions
fprintf(strcat("Getting present day collapse time for ", shelf, "\n"))
average_collapse_time = get_ave_shelf_collapse_time(shelf, step,0, 0, 0,f);
fprintf('Collapse time for shelf under present day is %.3f \n', average_collapse_time)

if average_collapse_time < target_time
    fprintf('present day collapse time is below target time, so I am done... \n')
    dM = 0;
    return
end

% do the regression to determine increase parameters
%[s_epsxx, s_dhdt] = get_regression_coefficients(shelf);

% or use the circum-Antarctic values
s_epsxx = 1.7127e-04;
s_dhdt = -0.0836;

% if you make it this far, 0 is a lower bound. Now find an upper bound.
init_lb = 0; %initial lower bound
inc_step = 10;  %how much to increase in this initial stage

have_init_ub = 0; %do we have an initial upper bound?
guess = init_lb + inc_step;
fprintf('Finding initial upper bound \n')
while ~have_init_ub
    d_m = guess;
    d_epsxx = guess*s_epsxx;
    d_dhdt = guess*s_dhdt;
    average_collapse_time = get_ave_shelf_collapse_time(shelf, step, d_m, d_epsxx,d_dhdt,f);
    
    if average_collapse_time < target_time
        fprintf('Guess of d_m = %.3f gave collapse time %.3f, so I have an upper bound. \n \n', guess, average_collapse_time);
        init_ub = guess;
        have_init_ub = 1;

    else
        init_lb = guess;
        guess = guess + inc_step;
        fprintf('Guess of d_m = %.3f gave collapse time %.3f, so not an upper bound. \n   New guess: %.3f \n \n', init_lb, average_collapse_time, guess);

    end
       
end

% do the bisection
lb = init_lb;
ub = init_ub;
guess = (ub + lb)/2;
err = ub - lb;
tol = 0.1;
fprintf('Performing the bisection... \n')

while err > tol
    d_m = guess;
    d_epsxx = guess*s_epsxx;
    d_dhdt = guess*s_dhdt;
    average_collapse_time = get_ave_shelf_collapse_time(shelf, step, d_m, d_epsxx,d_dhdt,f);
    fprintf('Guess of d_m = %.3f gave collapse time %.3f \n \n', guess, average_collapse_time);


    if average_collapse_time < target_time %have a new upper bound
        ub = guess;
        fprintf('Setting new upper bound (upper bound = %.3f, lower bound = %.3f) \n', ub, lb)
    else %have a new lower bound
        lb = guess;        
        fprintf('Setting new upper bound (upper bound = %.3f, lower bound = %.3f) \n', ub, lb)

    end
    guess = (ub + lb)/2;
    err = (ub -lb);
end

dM = guess;