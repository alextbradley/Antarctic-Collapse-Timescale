% launch simulations for circum Antarctic. This code produces individual
% collapse files, held within directory indexed by 'step'

jdir = dir('../../../data/ice-shelves/all-shelves/*.mat');

tic
for i = 1:length(jdir)
    shelf_name = strrep(jdir(i).name,'.mat','');
    step = 10;
    if ~(jdir(i).name == "NonClassifiedShelves.mat") %don't run on non-classified, for which code not yet set up right!
    [collapse_time, collapse_time_square, tags] = get_shelf_collapse_time(shelf_name, step);
    end
    tt(i) = toc;
end