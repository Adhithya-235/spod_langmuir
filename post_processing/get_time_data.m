function [t, nf, ke, cwke] = get_time_data(folder_name, maxs)

% This function reads HDF5 data produced by dedalus and extracts temporal
% grid as well as timeseries of the total and cross-wind kinetic energies. 
% Specify date-based folder name as well as the maximum number of data series produced. 
% Also, maxs is the maximum series number.
    
%% FILENAME
    
fname = string.empty;
    
for s = 1:maxs
    fname(s) = sprintf('../%s/energy_timeseries/energy_timeseries_s%d.h5', folder_name, s);
end

%% READ DATA
    
t = [];
ke = [];
cwke = [];

for s = 1:maxs
    t = [t; h5read(fname(s),'/scales/sim_time')];
    tempen = h5read(fname(s),'/tasks/KE');
    tempec = h5read(fname(s),'/tasks/CWKE');
    ke = [ke; tempen(:)];
    cwke = [cwke; tempec(:)];
end

nf = length(t);

end