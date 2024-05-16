%========================================================================%
%   Study SPOD modes produced from CL velocity fields, use the state
%   vector w, and a weight-matrix based on the trapezoidal rule.
%========================================================================%

%% READ DATA

[x, y, z, ~, ~, ~]  = get_space_data(folder_name, file_name, wrap);
[t, ~, ~, ~, Psi, nf] = get_field_data(folder_name, file_name, svec, wrap);
data = permute(Psi, [4,1,2,3]);
dt = mean(diff(t));

%% GET WEIGHT MATRIX AND DEFINE WINDOW LENGTH

weight_xyz = calc_3Dtrapzweights(x, y, z);
Nf         = 250;

%% PERFORM SPOD 

[L,P,f,Lc]       = spod(data,Nf,weight_xyz,[],dt); 
ang_f            = 2*pi*f;
[maxrow, maxcol] = find(L == max(max(L)))

%% PLOT MODE ENERGY

f1 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
hold on
patch([ang_f(2:end) fliplr(ang_f(2:end))],...
    [Lc(2:end,1,1)' fliplr(Lc(2:end,1,2)')], 'g')
plot(ang_f, L, 'b-','LineWidth', 3)
xline(2*pi/25, 'k--', 'LineWidth', 3)
xline(2*pi, 'k--', 'LineWidth', 3)
title(sprintf('SPOD SPECTRUM, La = %1.3f', La))
xlabel('$\omega$', 'interpreter', 'latex')
ylabel('$E_{mode}$', 'interpreter', 'latex')
xlim([1e-1 1e2])
ylim([1e-15 1e1])
box on
set(gca, 'boxstyle', 'full', 'linewidth', 3, 'fontsize', 20,...
        'yscale', 'log', 'xscale', 'log')
saveas(f1, sprintf('../%s/plots/spod/spod_eigenvalues.png', folder_name))

%% PLOT SELECT MODES

[Y, X, Z] = meshgrid(y, x, z);
mode = real(squeeze(P(maxrow,:,:,:,maxcol)));
Cmax = max(max(max(mode)));
Cmin = min(min(min(mode)));
f2 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
help_plot_fields(X,Y,Z,mode,x(end),y(end),z(33),[Cmin,Cmax],'dominant spod mode',[]);
saveas(f2, sprintf('../%s/plots/spod/spod_dom_mode.png', folder_name))

%% CREATE SOLUTION FIELD

time     = t(1):dt:t(end);
time     = time - time(1);
dom_mode = squeeze(P(maxrow,:,:,:,maxcol));
dom_evol = zeros([size(dom_mode), length(time)]);
for ti = 1:length(time)
    dom_evol(:,:,:,ti) = dom_mode*exp(1i*ang_f(maxrow)*time(ti));
end

%% SPACE-TIME PLOT

probe_w   = squeeze(Psi(1,:,33,:));
probe_spw = squeeze(dom_evol(1,:,33,:));
f3 = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
[xg, tg] = meshgrid(y,t);
subplot(121)
pcolor(xg,tg,probe_w')
colormap parula
c = colorbar;
shading interp
xlabel('$y$', 'interpreter', 'latex')
ylabel('$t$', 'interpreter', 'latex')
title(sprintf('La = %0.3f', La))
c.Label.Interpreter = 'latex';
c.FontSize = 20;
axis tight
axis square
box on
set(gca, 'fontsize', 20, 'linewidth', 3, 'boxstyle', 'full')

subplot(122)
pcolor(xg,tg,real(probe_spw)')
colormap parula
c = colorbar;
shading interp
xlabel('$y$', 'interpreter', 'latex')
title(sprintf('La = %0.3f', La))
c.Label.Interpreter = 'latex';
c.Label.String = '$\Psi_{spod}$';
c.FontSize = 20;
axis tight
axis square
box on
set(gca, 'fontsize', 20, 'linewidth', 3, 'boxstyle', 'full')
drawnow
saveas(f3, sprintf('../%s/plots/spod/spod_dom_trace.png', folder_name))