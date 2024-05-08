%=======================================================================%
%   Plot a timeseries of cross-wind kinetic energy. Best practice to run
%   from command line on a unix environment is to use the following
%   syntax:
%
%        matlab -batch "clear;close all;clc; ...
%               folder_name='$folder_name'; maxs=1; plot_cwke"
%
%=======================================================================%

%% READ DATA

[t, ~, ~, cwke] = get_timeseries_data(folder_name, maxs);

%% PLOT TIMESERIES

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
plot(t, cwke, '-', 'linewidth', 4)
xlabel('$T$', 'interpreter', 'latex')
ylabel('$CWKE$', 'interpreter', 'latex')
set(gca, 'fontsize', 30)
xlim([t(1), t(end)])
grid on
box on
set(gca, 'linewidth', 5)

%% SAVE VARIABLES AND PLOT

saveas(f, sprintf('../%s/plots/timeseries/cwke_timeseries.png', folder_name)) 

