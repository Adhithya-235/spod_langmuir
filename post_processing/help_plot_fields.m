function [] = help_plot_fields(X,Y,Z,PHI,xs,ys,zs,clim,clabel,time)

%   Plot slices of PHI. Inputs are self-explanatory (heheh), 
%   apart from clim, which is to be entered as a two-element 
%   increasing vector. 

hold off
slice(Y,X,Z,PHI,ys,xs,zs)  
caxis manual
caxis(clim)
colormap turbo
shading interp
c1 = colorbar;
xlabel('$y$', 'interpreter', 'latex')
ylabel('$x$', 'interpreter', 'latex')
zlabel('$z$', 'interpreter', 'latex')
title(sprintf('Time = %.2f', time))
c1.Label.Interpreter = 'latex';
c1.Label.String = clabel;
c1.Limits = clim;
c1.FontSize = 20;
c1.Location = 'eastoutside';
hold on
set(gca, 'fontsize', 20)
axis tight
box on
set(gca, 'boxstyle', 'full', 'linewidth', 2)
view([-15, 7])
%view([90, 90])
%view([0, 0])
daspect([1 1 1])

end