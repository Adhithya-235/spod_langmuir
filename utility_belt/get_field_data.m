function [t, U, V, W, Psi, nf] = get_field_data(folder_name, file_name, svec, wrap)

% This function reads HDF5 data produced by dedalus and extracts temporal
% grid as well as 4-dimensional (x, y, z, t) primitive variable fields. 
% Specify date-based folder name, data file name as well as the desired 
% data series numbers to be plotted. Also, wrap=1 means the last gridpoint 
% in periodic grids is kept. Useful for integration purposes.


%% FILENAME

fname = string.empty;
maxs = length(svec);
for s = 1:maxs
    fname(s) = sprintf('../%s/%s/%s_s%d.h5', folder_name, file_name, file_name, svec(s));
end

%% GET DATA FROM FILE

U   = [];
V   = [];
W   = [];
Psi = [];
t   = [];

for s = 1:maxs
   t   = [t; h5read(fname(s),'/scales/sim_time')];
   U   = cat(4, U, h5read(fname(s), '/tasks/U'));
   V   = cat(4, V, h5read(fname(s), '/tasks/V'));
   W   = cat(4, W, h5read(fname(s), '/tasks/W'));
   Psi = cat(4, Psi, h5read(fname(s), '/tasks/Psi'));
 end


%% PERMUTATION

U   = permute(U,[3,2,1,4]);
V   = permute(V,[3,2,1,4]);
W   = permute(W,[3,2,1,4]);
Psi = permute(Psi,[3,2,1,4]);

%% WRAP

if wrap == 1
    U   = cat(1, U, U(1, :, :, :));
    U   = cat(2, U, U(:, 1, :, :));
    V   = cat(1, V, V(1, :, :, :));
    V   = cat(2, V, V(:, 1, :, :));
    W   = cat(1, W, W(1, :, :, :));
    W   = cat(2, W, W(:, 1, :, :));
    Psi = cat(1, Psi, Psi(1, :, :, :));
    Psi = cat(2, Psi, Psi(:, 1, :, :));
end

%% DETERMINE TIMESERIES LENGTH

nf = length(t);

end