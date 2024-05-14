%=======================================================================%
%   Plot field snapshots for visualization purposes. Best practice to run
%   from command line on a unix environment is to use the following
%   syntax:
%
%        matlab -batch "clear;close all;clc; ...
%               folder_name='$folder_name'; svec=[1:3]; wrap=0; plot_fields"
%
%=======================================================================%

%% READ DATA

[x, y, z, X, Y, Z] = get_space_data(folder_name, file_name, wrap);
[t, U, V, W, Psi, nf] = get_field_data(folder_name, file_name, svec, wrap);


%% DEFINE SLICES

xs = x(1);
ys = y(1);
zs = z(end);

%% INITIALIZE FIGURE

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])

%% MAKE AND SAVE PLOT

for i = 1:nf

    %% PLOT U

    subplot(221)
    clim = [min(min(min(min(U)))), max(max(max(max(U))))];
    help_plot_fields(X,Y,Z,U(:,:,:,i),xs,ys,zs,clim,'U',t(i))
    title(sprintf('Time = %3.3f', t(i)))

    %% PLOT V

    subplot(222)
    clim = [min(min(min(min(V)))), max(max(max(max(V))))];
    help_plot_fields(X,Y,Z,V(:,:,:,i),xs,ys,zs,clim,'V',t(i))
    title(sprintf('Time = %3.3f', t(i)))

    %% PLOT W

    subplot(223)
    clim = [min(min(min(min(W)))), max(max(max(max(W))))];
    help_plot_fields(X,Y,Z,W(:,:,:,i),xs,ys,zs,clim,'W',t(i))
    title(sprintf('Time = %3.3f', t(i)))

    %% PLOT Psi

    subplot(224)
    clim = [min(min(min(min(Psi)))), max(max(max(max(Psi))))];
    help_plot_fields(X,Y,Z,Psi(:,:,:,i),xs,ys,zs,clim,'$\Psi$',t(i))
    title(sprintf('Time = %3.3f', t(i)))

    %% DRAW FRAME

    drawnow
    
    %% SAVE PLOT
     
    saveas(f, sprintf('../%s/plots/frames/frame_time_%3.3f.png', folder_name, t(i))) 

end