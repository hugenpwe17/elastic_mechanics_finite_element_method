%% Draw tetrahedron and its gauss points
% Contributed by Xiong
close all;
clear;
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

% define the gauss points cors and weight cell g and w
g = cell(1,3);
w = cell(1,3);

% Calculate Gauss points of different orders
for k=1:length(nint)
    [g{k}, w{k}] = Tet4Gp(k);
end

% Draw tetrahedron in isoparametric space
for l=1:length(nint)
    subplot(1,3,l)
    patch('vertices', x, 'faces', ix, 'facecolor', 'none', 'edgecolor', 'b');
    str = ['Gauss point of order ' num2str(l)];
    MyPatch(x, ix, 'none', 'b', str, 'x', 'y', 'z');
    % Draw gauss point of corresponding order (nint)
    hold on;
    plot3(g{l}(:,1), g{l}(:,2), g{l}(:,3), 'marker', 'x', 'color', 'r', 'linestyle', 'none', ...
        'markersize', 16); 
    % Fixed viewing angle
    view(43, 22);
end
set(figure(1),'Position',[0,0,1800,500]);

%% show shape function in isoparametric space
% Contributed by OuYang
close all;
clear;
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
for i=1:4
subplot(2,2,i)
scatter3(Xi(:,1),Xi(:,2),Xi(:,3),20,N(:, i),'filled')
view(31,38);
colorbar()
end
% set(figure(1),'Position',[0,0,800,600]);


%% Draw the tetrahedron and Gauss point in real space
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
[D,~] = size(PraCor);

% Tetrahedral connection order
ix   =[1 2 3 4 1 3 4 2];


 for k=1:length(nint)
    % Draw the tetrahedron in isoparametric space
    figure;
    patch('vertices', PraCor', 'faces', ix, 'facecolor', 'none', 'edgecolor', 'b')
    hold on
    
    % Calculate the coordinates (g) and number (ngp) of gauss points
    [g, w] = Tet4Gp(k);
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

%% jacobi test 
clc

dN = [1 0 0;
      0 1 0;
      0 0 1;
     -1 -1 -1];
 
% real TET4
% x  = rand(4,3);
x = [0 0 0; 
     1 0 0; 
     0 1 0;
     0 0 1];

% calculate the volume of real TET4
V  = x-x(1,:);
V(1,:) = [];
Vol= abs(dot(V(1,:),cross(V(2,:),V(3,:))))/6;
str1 = ['Vol = ' num2str(Vol)];
disp(str1)

% from other 
[J, detJ] = ShapeFunJacob(dN, x);
str2 = ['detJ = ' num2str(detJ*(1/6))];
disp(str2)

if (abs(Vol)-abs(detJ)) <10^-9
    disp('right');
else
    disp('wrong');
end

%% test
% material and model properties
D    = 3;
nnde = 4;
nint = 3;
E    = 10;
nu   = 0.33;

%element info 
x = [0.6948 0.4387 0.1869
    0.3171 0.3816 0.4898
    0.9502 0.7655 0.4456
    0.0344 0.7952 0.6463];

ix = [1 2 3 4];
fx = vertoface(x,ix);

%% Analytical solution
% fx = VerToFace(x(ix(2,:),:),ix(2,:));

% n = 3;
% [~,fx] = VerToFace(x(ix(n,:),:),ix(n,:));
% patch('vertices', x(ix(n,:),:), 'faces', fx, ...
%'facecolor', 'none', 'edgecolor', 'b');

% Fix position 
x(:,1) = x(:,1) - 0.5;
x(:,2) = x(:,2) - 0.5;

% Draw grid
for i = 1: length(ix)
n = i;
[~,fx] = VerToFace(x(ix(n,:),:),ix(n,:));
patch('vertices', x(ix(n,:),:), 'faces', fx, ...
    'facecolor', 'none', 'edgecolor', 'k');
hold on
end
hold on 
axis equal

% h1 = mypatch(x, ix, 'b', 'k', []); 

% material parameters
rho = 0.33;
E   = 100;
nu  = 0.99;

% Geometric parameters
Lz  = max(x(:,3)) - min(x(:,3));
Ly  = max(x(:,2)) - min(x(:,2));
Lx  = max(x(:,1)) - min(x(:,1));
% physical parameters
g   = 9.8;

% Analytical solution value
sig_z = rho* g* x(:,3);
h = mypatch(x, ix, 'interp', 'none', sig_z, 0.5); 

eps_xx = -(nu* rho* g* x(:,3))/E;
eps_yy = -(nu* rho* g* x(:,3))/E;
eps_zz = (rho* g* x(:,3))/E;

h = mypatch(x, ix, 'interp', 'none', eps_xx, 0.5); 
h = mypatch(x, ix, 'interp', 'none', eps_yy, 0.5); 
h = mypatch(x, ix, 'interp', 'none', eps_zz, 0.5); 

u_x = -(nu* rho* g).* x(:,1).* x(:,3) /E;
u_y = -(nu* rho* g).* x(:,2).* x(:,3) /E;
u_z = (rho* g/(2*E)).* (x(:,3).^2- Lz^2+ nu*(x(:,1).^2+x(:,2).^2));

x1 = x+ [u_x,zeros(length(u_x),1),zeros(length(u_x),1)];
x2 = x+ [zeros(length(u_y),1),u_y,zeros(length(u_y),1)];
x3 = x+ [zeros(length(u_z),1),zeros(length(u_z),1),u_z];

xx = x+ [u_x,u_y,u_z];

h = mypatch(x1, ix, 'interp', 'none', u_x, 0.5); 
h = mypatch(x2, ix, 'interp', 'none', u_y, 0.5); 
h = mypatch(x3, ix, 'interp', 'none', u_z, 0.5); 


%% optimize 
for i = 1 : length(ix)
[fxs([(4*i-3):(4*i)],:),~] = VerToFace(x(ix(i,:),:),ix(i,:));
end

%% draw the 
h1 = patch('vertices', x, 'faces', fxs, ...
'facecolor', 'b', 'edgecolor', 'k');
hold on 
axis equal