close all; 
clear ; 
clc; 

% mesh grid
Lx = 1; 
Ly = 1; 
Lz = 3.5; 

Nx = 19; 
Ny = 19; 
Nz = 38; 

mesher = MyRegular_Mesher; 
mesher.create(5, [Lx, Ly, Lz], [Nx, Ny, Nz], [-Lx/2 -Ly/2 0]); 
x   = mesher.x; 
nx  = mesher.nx; 
ix  = mesher.ix; 
nix = mesher.nix; 

% figure; 
% plot_mesh(3, x, ix, nnde * ones(nix, 1), 'none', 'k', 1); 
% axis equal
% xlabel('x'); 
% ylabel('y'); 
% zlabel('z'); 
% view(19, 23);

% material and model properties
nnde  = 4; 
D     = 3; 
nintx = 3;    % int. order for K, M 
nints = 3;    % int. order for recovery matarices (Mr, Reps, Rsig)
E     = 2000; 
nu    = 0.33; 
%　density
rho = 3;
% Geometric parameters
% length in x,y,z
Lz  = max(x(:,3)) - min(x(:,3));
Ly  = max(x(:,2)) - min(x(:,2));
Lx  = max(x(:,1)) - min(x(:,1));
% physical force
% gravity
g   = -9.8;

% shape function structures
SFx = GenerateShapeFunction(D,nnde,nintx); 
SFs1 = GenerateShapeFunction(D,nnde,nints); 
SFs2 = GenerateShapeFunction(D,nnde,1);
% elastic tensor 
CC = ElastTensor(E, nu);
Dc = size(CC, 1); 

% integrate K, M
 tic;
[K, M] = IntKM(x, ix, SFx, CC);
 disp(['integration of K, M finished for :', num2str(toc), 'sec']);

% integrate Mr, Reps, Rsig
 tic;
[Mr, Reps, Rsig] = intRecoverMat(x, ix, SFs1, SFs2, CC);
disp(['integration of Mr, Reps, Rsig finished for :', num2str(toc), 'sec']);

% Boundary conditions
topID   = find(abs(x(:, 3) - Lz) < 1e-5); 
topLeft = find(abs(x(:, 3) - Lz) < 1e-5 & abs(x(:, 1) - 0) < 1e-5); 
topBack = find(abs(x(:, 3) - Lz) < 1e-5 & abs(x(:, 2) - 0) < 1e-5); 

% constrain
% pu2    = [topID * D - 2, topID * D - 1,topID * D]'; 
pu2    = [topID * D; topLeft * D - 2; topBack * D - 1]; 
pu2    = pu2(:); 
u      = zeros(nx * D, 1); 
u(pu2)       = 0;

% external forces
puf  = topID * D; 
fext = zeros(nx, D); 
fext(:,3) = -rho* g* Lx* Ly* Lz / (nx);
fext = fext';
fext = fext(:);

% solve for displacement
[u, fext] = solveLin(K, u, fext, pu2);

% recover stress and strain
sigVec = Mr \ (Rsig * u); 
sig    = reshape(sigVec, Dc, nx)';
epsVec = Mr \ (Reps * u); 
eps    = reshape(epsVec, Dc, nx)'; 

% deformed configuration; 
u0 = reshape(u, D, nx)';
x1 = x + u0; 
nnd = nnde *  ones(size(ix, 1), 1);


h1 =figure(1); 
 mypatch2(x1,ix,nnd,'interp','none',u0(:,3),1,1,'u_z','x axis','y axis','z axis',-0.1,0);

h2 =figure(2);
 mypatch2(x1,ix,nnd,'interp','none',u0(:,2),1,1,'u_y','x axis','y axis','z axis',-8*10^(-3),8*10^(-3));

h3 =figure(3);
 mypatch2(x1,ix,nnd,'interp','none',u0(:,1),1,1,'u_x','x axis','y axis','z axis',-8*10^(-3),8*10^(-3));

h4 =figure(4);
 mypatch2(x1,ix,nnd,'interp','none',eps(:,1),1,1,'\epsilon_{xx}','x axis','y axis','z axis',-18*10^(-3),0);

h5 =figure(5);
 mypatch2(x1,ix,nnd,'interp','none',eps(:,2),1,1,'\epsilon_{yy}','x axis','y axis','z axis',-18*10^(-3),0);

h6 =figure(6);
 mypatch2(x1,ix,nnd,'interp','none',eps(:,3),1,1,'\epsilon_{zz}','x axis','y axis','z axis',0,0.05);

h7 =figure(7);
 mypatch2(x1,ix,nnd,'interp','none',eps(:,4),1,1,'\epsilon_{xy}','x axis','y axis','z axis',-15*10^(-4),15*10^(-4));

h8 =figure(8);
 mypatch2(x1,ix,nnd,'interp','none',eps(:,5),1,1,'\epsilon_{xz}','x axis','y axis','z axis',-3*10^(-3),5*10^(-3));

h9 =figure(9);
 mypatch2(x1,ix,nnd,'interp','none',eps(:,6),1,1,'\epsilon_{yz}','x axis','y axis','z axis',-3*10^(-3),5*10^(-3));

h10 =figure(10);
 mypatch2(x1,ix,nnd,'interp','none',sig(:,1),1,1,'\sigma_{xx}','x axis','y axis','z axis',-2,3);

h11=figure(11);
  mypatch2(x1,ix,nnd,'interp','none',sig(:,2),1,1,'\sigma_{yy}','x axis','y axis','z axis',-2,3);

h12=figure(12);
mypatch2(x1,ix,nnd,'interp','none',sig(:,3),1,1,'\sigma_{zz}','x axis','y axis','z axis',0,110);

h13=figure(13);
 mypatch2(x1,ix,nnd,'interp','none',sig(:,4),1,1,'\sigma_{xy}','x axis','y axis','z axis',-1.5,0.3);

h14=figure(14);
 mypatch2(x1,ix,nnd,'interp','none',sig(:,5),1,1,'\sigma_{xz}','x axis','y axis','z axis',-2,4);

h15 =figure(15);
 mypatch2(x1,ix,nnd,'interp','none',sig(:,6),1,1,'\sigma_{yz}','x axis','y axis','z axis',-2,4);

% disp('Saving...');
% print(h1,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/u_z','-dpng','-r600')
% print(h2,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/u_y','-dpng','-r600')
% print(h3,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/u_x','-dpng','-r600')
% print(h4,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/epsilon_{xx}','-dpng','-r600')
% print(h5,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/epsilon_{yy}','-dpng','-r600')
% print(h6,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/epsilon_{zz}','-dpng','-r600')
% print(h7,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/epsilon_{xy}','-dpng','-r600')
% print(h8,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/epsilon_{xz}','-dpng','-r600')
% print(h9,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/epsilon_{yz}','-dpng','-r600')
% print(h10,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/sigma_{xx}','-dpng','-r600')
% print(h11,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/sigma_{yy}','-dpng','-r600')
% print(h12,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/sigma_{zz}','-dpng','-r600')
% print(h13,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/sigma_{xy}','-dpng','-r600')
% print(h14,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/sigma_{xz}','-dpng','-r600')
% print(h15,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/FEM/sigma_{yz}','-dpng','-r600')

disp('done');

