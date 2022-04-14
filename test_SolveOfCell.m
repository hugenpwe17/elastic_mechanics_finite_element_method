%% Finite element solution for a single element
% Contributed by OuYang,Xiong

%% definte parameters
% material parameters
% young modulus
E    = 10; 
% poisson's ratio
nu   = 0.33;
% model parameters
% dimensionality
D    = 3; 
% number of nodes
nnde = 4; 
% order of numerical integration
nint = 3;
% all degrees of freedom
nf   = D* nnde;

%% define element info
% coordinates
x = [-0.1968   -0.1969         0
     -0.5000         0         0
     -0.5000   -0.5000         0
     -0.5000   -0.1934    0.3260]*100;
% node number
ix = [1 2 3 4]; 

%% define numerical integration parameters
% get parameters
SF = GenerateShapeFunction(D,nnde,nint);
% calculate elast tensor
CC = ElastTensor(E,nu);

%% take a numerical integration
[K, M]=IntKMLoc(SF, CC, x(ix, :));

%% Show K matrix and its eigenvalues
disp('K = : '); 
disp(K); 
disp('eig = : ')
disp(eig(K)); 

%% initialization displacement and force vector
u    = zeros(D*nnde,1);
fext = zeros(D*nnde,1);

%% define constraint point and others
pu2 = [1 2 3 4 8 12 11];
pu1 = setxor((1:nf), pu2); 

%% impose constraints
u(pu2)  = 0;
fext(:) = 0;
% extrude the 4th node a distance 
u(12)   = 0.1;

%% define Known displacement and force 
u2 = u(pu2);
f1 = fext(pu1);

%% calculate Unknown displacement and force
u(pu1)    = K(pu1, pu1) \ (f1 - K(pu1, pu2) * u2); 
fext(pu2) = K(pu2, :) * u; 

%% Draw displacement distribution and force distribution
% reshape force and displacement matrix  
fext0 = reshape(fext,[3,4])';
u0    = reshape(u,[3,4])';
% get the current location of node
xu    = x+u0;
% convert vertex matrix to face matrix
fx    = VerToFace(x,ix);

% difine title
ftl ={'Force in x', 'Force in y', 'Force in z'};
% draw all force distribution 
for i =1:3
figure
% draw the frame
patch('vertices', x, 'faces', fx, 'facecolor', 'none', 'edgecolor', 'b');
hold on
% draw force distribution
mypatch(x, fx, 'interp', 'none', fext0(:,i), 0.5, 0,ftl{i},'x axis','y axis','z axis',0,3*10^4);
% mark node number
gui_label(3, x, fx, 4, 1, ix);
end

%define title
dtl ={'Displacement in x', 'Displacement in y', 'Displacement in z'};
% draw all displacement distribution
for i =1:3
figure
% draw frame
patch('vertices', x, 'faces', fx, 'facecolor', 'none', 'edgecolor', 'b');
hold on
% draw displacement distribution
mypatch(xu, fx, 'interp', 'none', u0(:,i), 0.5, 1,dtl{i},'x axis','y axis','z axis',0,0.15);
% mark node number
gui_label(3, x, fx, 4, 1, ix);
end
