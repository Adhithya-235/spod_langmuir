function MOV = animate_slices(X,Y,Z,PHI,xs,ys,zs,clabel,time,nf)

%   Animates slices of the given data. 

%% INITIALIZE FIGURE

f = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])

%% PREALLOCATE FRAMES

MOV(nf) = struct('cdata',[],'colormap',[]);

%% COMPUTE COLOUR AXIS LIMITS

Cmax = max(max(max(max(PHI))));
Cmin = min(min(min(min(PHI))));

%% ANIMATE SLICES

for i=1:nf
   help_plot_fields(X, Y, Z, squeeze(PHI(:, :, :, i)), xs, ys, zs,...
       [Cmin,Cmax], clabel, time(i));
   drawnow
   MOV(i) = getframe(f); 
end

end