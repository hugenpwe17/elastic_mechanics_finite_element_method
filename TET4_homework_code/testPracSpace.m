% Draw the tetrahedron and Gauss point in real space
% Contributed by OuYang, Xiong
clear;
close all;
clc;

% Order
nint = 1:3;

% Draw a specific tetrahedron
PraCor = [0.2511,0.3517,0.5497,0.7572;
        0.6160,0.8308,0.9172,0.7537;
        0.4733,0.5853,0.2858,0.3804];

% Find the dimensionality and node's number in input data
[D,nnde] = size(PraCor);

% Tetrahedral connection order
ix   =[1 2 3 4 1 3 4 2];


 for k=1:length(nint)
    % Draw the tetrahedron in isoparametric space
    figure;
    patch('vertices', PraCor', 'faces', ix, 'facecolor', 'none', 'edgecolor', 'b')
    hold on
    
    % Calculate the coordinates (g) and number (ngp) of gauss points
    [g, w] = TET4_GP(k);
    ngp    = size(g,1);
    
    % Define gauss point matrix in real space
    GusCor = zeros(ngp,D);
    
    for i = 1:ngp
        % Calculate shape function value in the gauss point
        [N,~] = ShapeFun(g(i,:));
        % Interpolate nodes to find Gauss points in real space
        GusCor(i,:) = (sum(N.*PraCor,2))';
    end
    
   % Draw the position of the Gauss point in the tetrahedron in real space
    scatter3(GusCor(:,1),GusCor(:,2),GusCor(:,3),50)
    view(72,35);
    
 end
