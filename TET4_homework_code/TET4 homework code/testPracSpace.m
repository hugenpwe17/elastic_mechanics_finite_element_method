clear;
close all;
clc;
% Draw the tetrahedron and Gauss point in real space

nint = 1:3;
% Order

PraCor = [0.2511,0.3517,0.5497,0.7572;
        0.6160,0.8308,0.9172,0.7537;
        0.4733,0.5853,0.2858,0.3804];
% Draw a specific tetrahedron

[D,nnde] = size(PraCor);
% Find the dimensionality and node's number in input data

ix   =[1 2 3 4 1 3 4 2];
% Tetrahedral connection order

 for k=1:length(nint)
    figure;
    patch('vertices', PraCor', 'faces', ix, 'facecolor', 'none', 'edgecolor', 'b')
    % Draw the tetrahedron in isoparametric space
    hold on
    
    [g, w] = TET4_GP(k);
    ngp    = size(g,1);
    % Calculate the coordinates (g) and number (ngp) of gauss points
    GusCor = zeros(ngp,D);
    % Define gauss point matrix in real space
     
    for i = 1:ngp
        [N,~] = ShapeFun(g(i,:));
        % Calculate shape function value in the gauss point 
        GusCor(i,:) = (sum(N.*PraCor,2))';
        % Interpolate nodes to find Gauss points in real space
    end
   
    scatter3(GusCor(:,1),GusCor(:,2),GusCor(:,3),50)
    view(72,35);
    % Draw the position of the Gauss point in the tetrahedron in real space
 end
 % % Contributed by OuYang, Xiong