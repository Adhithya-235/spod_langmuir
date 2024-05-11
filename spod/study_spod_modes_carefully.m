%========================================================================%
%   Study SPOD modes produced from CL velocity fields, use the state
%   vector w, and a weight-matrix based on the trapezoidal rule.
%========================================================================%

clear
close all
clc

%% WRITE La

La  = 0.001;
Lat = 0;

%% LOAD FILES

cl   = load('import/la1e-3_new/wrclv.mat', 'W', 'x', 'y', 'z', 't'); 
data = permute(cl.W, [4,1,2,3]);
dt = mean(diff(cl.t));

%% GET WEIGHT MATRIX AND DEFINE WINDOW LENGTH

weight_xyz = calc_3Dtrapzweights(cl.x, cl.y, cl.z);
Nf         = 50;

%% PERFORM SPOD 

[L,P,f,Lc] = spod(data,Nf,weight_xyz,[],dt); 
ang_f      = 2*pi*f;

%% PLOT MODE ENERGY

figure(1)
hold on
patch([ang_f(2:end) fliplr(ang_f(2:end))],...
    [Lc(2:end,1,1)' fliplr(Lc(2:end,1,2)')], 'g')
plot(ang_f, L, 'b-','LineWidth', 3)
xline(2*pi/25, 'k--', 'LineWidth', 3)
title(sprintf('SPOD SPECTRUM, La = %1.3f', La))
xlabel('$\omega$', 'interpreter', 'latex')
ylabel('$E_{mode}$', 'interpreter', 'latex')
xlim([1e-1 1e1])
ylim([1e-6 1e1])
box on
set(gca, 'boxstyle', 'full', 'linewidth', 3, 'fontsize', 20,...
        'yscale', 'log', 'xscale', 'log')

%% PLOT SELECT MODES

[Y, X, Z] = meshgrid(cl.y, cl.x, cl.z);
mode = real(squeeze(P(2,:,:,:,1)));
Cmax = max(max(max(mode)));
Cmin = min(min(min(mode)));
figure(4);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
plot_slices(X,Y,Z,mode,cl.x(end),cl.y(end),cl.z(33),[Cmin,Cmax],'spod mode',[]);

%% CREATE SOLUTION FIELD

time     = 0:dt:100;
dom_mode = squeeze(P(4,:,:,:,1));
dom_evol = zeros([size(dom_mode), length(time)]);
for ti = 1:length(time)
    dom_evol(:,:,:,ti) = dom_mode*exp(1i*ang_f(4)*time(ti));
end

%% SPACE-TIME PLOT

probe_w   = squeeze(cl.W(1,:,33,:));
probe_spw = squeeze(dom_evol(1,:,33,:));
figure(4)
[xg, tg] = meshgrid(cl.y,cl.t);
subplot(121)
pcolor(xg,tg,probe_w')
colormap parula
c = colorbar;
shading interp
xlabel('$y$', 'interpreter', 'latex')
ylabel('$t$', 'interpreter', 'latex')
title(sprintf('La = %0.3f, La_t = %0.3f', La, Lat))
c.Label.Interpreter = 'latex';
c.FontSize = 30;
axis tight
axis square
box on
set(gca, 'fontsize', 30, 'linewidth', 3, 'boxstyle', 'full')

subplot(122)
pcolor(xg,tg,real(probe_spw)')
colormap parula
c = colorbar;
shading interp
xlabel('$y$', 'interpreter', 'latex')
title(sprintf('La = %0.3f, La_t = %0.3f', La, Lat))
c.Label.Interpreter = 'latex';
c.Label.String = '$w_{spod}$';
c.FontSize = 30;
axis tight
axis square
box on
set(gca, 'fontsize', 30, 'linewidth', 3, 'boxstyle', 'full')