function runsim_shelf_hpc(shelf_name, step)
    %wrapper script for running sequential shelf files
    [~, ~, ~] = get_shelf_collapse_time_HPC(shelf_name, step);
end