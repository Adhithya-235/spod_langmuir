%=======================================================================%
%   Animate snapshots for visualization purposes. Best practice to run
%   from command line on a unix environment is to use the following
%   syntax:
%
%        matlab -batch "clear;close all;clc; ...
%               folder_name='$folder_name'; svec=[1:3]; wrap=0; animate_fields"
%
%=======================================================================%

%% READ DATA

[x, y, z, X, Y, Z]    = get_space_data(folder_name, file_name, wrap);
[t, ~, ~, ~, Psi, nf] = get_field_data(folder_name, file_name, svec, wrap);


%% DEFINE SLICES

xs = x(1);
ys = y(end);
zs = z(33);

%% MAKE MOVIE AND COLLECT FRAMES

MOV = help_animate_fields(X,Y,Z,Psi,xs,ys,zs,'$\Psi$',t,nf);

%% CREATE AVI FILE

filename = sprintf('../%s/plots/movies/streamfunction.avi', folder_name);
v = VideoWriter(filename); 
v.FrameRate = 16;
open(v)
for i = 1:nf
   writeVideo(v, MOV(i)) 
end
close(v)
