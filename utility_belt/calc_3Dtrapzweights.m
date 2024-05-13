function [weight_xyz] = calc_3Dtrapzweights(x,y,z)

% GET_3DTRAPZWEIGHTS Integration weight tensor for 3D cartesian coordinates 
% using trapezoidal rule
% AS, 2023

%% SET UP SOLUTION ARRAYS

weight_x = zeros(length(x),1);
weight_y = zeros(length(y),1);
weight_z = zeros(length(z),1);

%% X DIRECTION WEIGHTS

weight_x(1)   = (x(2)-x(1))/2;
weight_x(end) = (x(end)-x(end-1))/2;

for i = 2:length(x)-1
    weight_x(i) = (x(i)-x(i-1))/2 + (x(i+1)-x(i))/2; 
end

%% Y DIRECTION WEIGHTS

weight_y(1)   = (y(2)-y(1))/2;
weight_y(end) = (y(end)-y(end-1))/2;

for i = 2:length(y)-1
    weight_y(i) = (y(i)-y(i-1))/2 + (y(i+1)-y(i))/2; 
end

%% Z DIRECTION WEIGHTS

weight_z(1)   = (z(2)-z(1))/2;
weight_z(end) = (z(end)-z(end-1))/2;

for i = 2:length(z)-1
    weight_z(i) = (z(i)-z(i-1))/2 + (z(i+1)-z(i))/2; 
end

%% TENSOR PRODUCT

weight_xyz = squeeze(tensorprod(squeeze(tensorprod(weight_x,weight_y)),weight_z));

end