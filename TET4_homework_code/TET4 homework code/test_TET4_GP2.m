close all;
clear all;
clc;
% Draw tetrahedron and its gauss points

x = [0 0 0; 
     1 0 0; 
     0 1 0;
     0 0 1];
% Tetrahedron vertex coordinates in isoparametric space

ix=[1 2 3 4 1 3 4 2];
% Tetrahedral connection order

nint=1:3;
% Order
for k=1:length(nint)
    [g, w] = TET4_GP(k);
    % Gauss points and weighs
    figure;
    patch('vertices', x, 'faces', ix, 'facecolor', 'none', 'edgecolor', 'b');
    % Draw tetrahedron in isoparametric space
    hold on;
    plot3(g(:, 1), g(:, 2), g(:, 3), 'marker', 'x', 'color', 'r', 'linestyle', 'none', ...
        'markersize', 16); 
    % Draw gauss point of corresponding order (nint)
    view(43, 22);
    % Fixed viewing angle
end
% Contributed by Xiong
