close all; 
clear all; 
clc; 

[xi_x, xi_y] = meshgrid(0:0.1:1, 0:0.1:1); 
xi_x = xi_x(:); 
xi_y = xi_y(:); 
xi = [xi_x, xi_y]; 
xi = xi((xi_x + xi_y) <= 1, :); 
nxi = size(xi, 1); 

N = zeros(nxi, 3); 
for i = 1:nxi
    N(i, :) = ShapeFunTri(xi(i, :)); 
end


pltN = @(i, s, c) plot3(xi(:, 1), xi(:, 2), N(:, i), 'marker', s, 'linestyle', 'none', ...
    'color', c, 'markersize', 10, 'markerfacecolor', c); 

figure; 
hold on;
pltN(1, 'o', 'r'); 
pltN(2, 'o', 'g'); 
pltN(3, 'o', 'b');
view(16, 15); 
axis equal;


