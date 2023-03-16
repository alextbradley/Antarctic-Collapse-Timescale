function runsim_shelf_hpc(shelf_num, step)
    %wrapper script for running sequential shelf files

    jdir = dir('../../../data/ice-shelves/all-shelves/*.mat');
    shelf_name = strrep(jdir(shelf_num).name,'.mat','');
    if ~(jdir(shelf_num).name == "NonClassifiedShelves.mat") %don't run on non-classified, for which code not yet set up right!
        [~, ~, ~] = get_shelf_collapse_time_HPC(shelf_name, step);
    end
end