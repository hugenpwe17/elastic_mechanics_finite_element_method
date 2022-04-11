% show shape function in isoparametric space
% Contributed by OuYang
close all;
clear all;
clc;

% creat grid points
[Xi_x,Xi_y,Xi_z]...
        =   meshgrid(0:0.1:1,0:0.1:1,0:0.1:1);
Xi_x    =   Xi_x(:);
Xi_y    =   Xi_y(:);
Xi_z    =   Xi_z(:);
Xi      =   [Xi_x(:),Xi_y(:),Xi_z(:)];

% select the points inside and on the tetrahedron
Xi      =   Xi((Xi_x + Xi_y + Xi_z) <= 1,:);

% find the shape function of each node 
NumXi   =   size(Xi,1);
N       =   zeros(NumXi,4);
for i   =   1:NumXi
    N(i,:)...
        =   ShapeFun(Xi(i,:));
end

% Draw the image of the shape function in the parameter function space
figure(1);
scatter3(Xi(:,1),Xi(:,2),Xi(:,3),20,N(:, 1),'filled')
view(31,38);
colorbar;
figure(2);
scatter3(Xi(:,1),Xi(:,2),Xi(:,3),20,N(:, 2),'filled')
view(31,38);
colorbar;
figure(3);
scatter3(Xi(:,1),Xi(:,2),Xi(:,3),20,N(:, 3),'filled')
view(31,38);
colorbar;
figure(4);
scatter3(Xi(:,1),Xi(:,2),Xi(:,3),20,N(:, 4),'filled')
view(31,38);
colorbar;