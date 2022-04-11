% Draw tetrahedron and its gauss points
% Contributed by Xiong
close all;
clear all;
clc;

% Tetrahedron vertex coordinates in isoparametric space
x = [0 0 0; 
     1 0 0; 
     0 1 0;
     0 0 1];

% Tetrahedral connection order
ix=[1 2 3 4 1 3 4 2];

% Order
nint=1:3;

for k=1:length(nint)
    % Gauss points and weighs
    [g, w] = TET4_GP(k);
    % Draw tetrahedron in isoparametric space
    figure;
    patch('vertices', x, 'faces', ix, 'facecolor', 'none', 'edgecolor', 'b');
    % Draw gauss point of corresponding order (nint)
    hold on;
    plot3(g(:, 1), g(:, 2), g(:, 3), 'marker', 'x', 'color', 'r', 'linestyle', 'none', ...
        'markersize', 16); 
    % Fixed viewing angle
    view(43, 22);
end

