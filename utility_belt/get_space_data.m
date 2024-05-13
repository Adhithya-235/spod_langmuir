function [x, y, z, X, Y, Z] = get_space_data(folder_name, file_name, wrap)

% This function reads HDF5 data produced by dedalus and extracts the spatial grid. 
% Specify date-based folder name and data file name. Also, wrap=1 means the last 
% gridpoint in periodic grids is kept. Useful for integration purposes.

%% FILENAME 

fname = sprintf('../%s/%s/%s_s1.h5', folder_name, file_name, file_name);

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