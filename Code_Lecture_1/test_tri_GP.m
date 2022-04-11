close all; 
clear all; 
clc; 

x  = [0 0; 
    1 0; 
    0 1]; 
ix = [1 2 3]; 

nint = 1:3; 
sty  = {'x', 'o', '^'}; 
clr  = {'r', 'b', 'k'}; 
figure; 
patch('vertices', x, 'faces', ix, 'facecolor', 'none', 'edgecolor', 'g'); 
hold on; 
for i = 1:length(nint)
    [g, w] = tri_GP(nint(i)); 
    plot(g(:, 1), g(:, 2), 'marker', sty{i}, 'color', clr{i}, 'linestyle', 'none', ...
        'markersize', 16); 
end