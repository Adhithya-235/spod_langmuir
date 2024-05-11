function [x, y, z, X, Y, Z] = get_space_data(folder_name, probe_toggle, wrap)

% This function reads HDF5 data produced by dedalus and extracts the spatial grid. 
% Specify date-based folder name. Also, wrap=1 means the last gridpoint in 
% periodic grids is kept. Useful for integration purposes. Finally, if probe_toggle = 1, 
% probe data will be collected instead of more finely sampled field snapshots.


%% FILENAME 
if probe_toggle == 1
    fname = sprintf('../%s/field_probes/field_probes_s1.h5', folder_name);
else if probe_toggle == 0    
    fname = sprintf('../%s/field_snapshots/field_snapshots_s1.h5', folder_name);
end

%% READ x, y AND z DATA

x = h5read(fname,'/scales/x/1.0');
y = h5read(fname,'/scales/y/1.0');
z = h5read(fname,'/scales/z/1.0');

%% WRAP 

if wrap == 1
    x = cat(1, x, 2*x(end)-x(end-1));
    y = cat(1, y, 2*y(end)-y(end-1));
end

%% MESHGRID

[Y, X, Z] = meshgrid(y, x, z);

end