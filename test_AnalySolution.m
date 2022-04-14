% Analytical solution
% Contributed by OuYang,Xiong
close all;
clear;
clc;
% define parameters
% material parameters
% young modulus
E   = 2000;
% poisson's ratio
nu  = 0.33;
%　density
rho = 3;
% load mesh model 
% x  : coordinates (3* total nnde)
% ix : node number matrix (nnde in cell* cell number) 
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
nnde = 4; 

% figure; 
% plot_mesh(3, x, ix, nnde * ones(nix, 1), 'none', 'k', 1); 
% axis equal
% xlabel('x'); 
% ylabel('y'); 
% zlabel('z'); 
% view(19, 23);


% physical force
%gravity
g   = 9.8;

% calculate analytical solution value
% stress (x y z)
sig_zz  = rho* g* x(:,3);
sig_xx  = zeros(size(sig_zz));
sig_yy  = zeros(size(sig_zz));
sig_xy  = zeros(size(sig_zz));
sig_yz  = zeros(size(sig_zz));
sig_xz  = zeros(size(sig_zz));
% strain (x y z)
eps_xx = -(nu* rho* g* x(:,3))/E;
eps_yy = -(nu* rho* g* x(:,3))/E;
eps_zz = (rho* g* x(:,3))/E;
eps_xy  = zeros(size(eps_zz));
eps_yz  = zeros(size(eps_zz));
eps_xz  = zeros(size(eps_zz));
% displacement (x y z) 
u_x = -(nu* rho* g).* x(:,1).* x(:,3) /E;
u_y = -(nu* rho* g).* x(:,2).* x(:,3) /E;
u_z = (rho* g/(2*E)).* (x(:,3).^2- Lz^2+ nu*(x(:,1).^2+x(:,2).^2));
%current configuration's coordinates
xx = x+ [u_x,u_y,u_z];
u = [u_x,u_y,u_z];
sig=[sig_xx,sig_yy,sig_zz,sig_xy,sig_xz,sig_yz];
eps = [eps_xx,eps_yy,eps_zz,eps_xy,eps_xz,eps_yz];

% draw the figures
nnd = nnde *  ones(size(ix, 1), 1);

h1 =figure(1); 
 mypatch2(x,ix,nnd,'interp','none',u(:,3),1,1,'u_z','x axis','y axis','z axis',-0.1,0);

h2 =figure(2);
 mypatch2(x,ix,nnd,'interp','none',u(:,2),1,1,'u_y','x axis','y axis','z axis',-8*10^(-3),8*10^(-3));

h3 =figure(3);
 mypatch2(x,ix,nnd,'interp','none',u(:,1),1,1,'u_x','x axis','y axis','z axis',-8*10^(-3),8*10^(-3));

h4 =figure(4);
 mypatch2(x,ix,nnd,'interp','none',eps(:,1),1,1,'\epsilon_{xx}','x axis','y axis','z axis',-18*10^(-3),0);

h5 =figure(5);
 mypatch2(x,ix,nnd,'interp','none',eps(:,2),1,1,'\epsilon_{yy}','x axis','y axis','z axis',-18*10^(-3),0);

h6 =figure(6);
 mypatch2(x,ix,nnd,'interp','none',eps(:,3),1,1,'\epsilon_{zz}','x axis','y axis','z axis',0,0.05);

h7 =figure(7);
 mypatch2(x,ix,nnd,'interp','none',eps(:,4),1,1,'\epsilon_{xy}','x axis','y axis','z axis',-15*10^(-4),15*10^(-4));

h8 =figure(8);
 mypatch2(x,ix,nnd,'interp','none',eps(:,5),1,1,'\epsilon_{xz}','x axis','y axis','z axis',-3*10^(-3),5*10^(-3));

h9 =figure(9);
 mypatch2(x,ix,nnd,'interp','none',eps(:,6),1,1,'\epsilon_{yz}','x axis','y axis','z axis',-3*10^(-3),5*10^(-3));

h10 =figure(10);
 mypatch2(x,ix,nnd,'interp','none',sig(:,1),1,1,'\sigma_{xx}','x axis','y axis','z axis',-2,3);

h11=figure(11);
  mypatch2(x,ix,nnd,'interp','none',sig(:,2),1,1,'\sigma_{yy}','x axis','y axis','z axis',-2,3);

h12=figure(12);
mypatch2(x,ix,nnd,'interp','none',sig(:,3),1,1,'\sigma_{zz}','x axis','y axis','z axis',0,110);

h13=figure(13);
 mypatch2(x,ix,nnd,'interp','none',sig(:,4),1,1,'\sigma_{xy}','x axis','y axis','z axis',-1.5,0.3);

h14=figure(14);
 mypatch2(x,ix,nnd,'interp','none',sig(:,5),1,1,'\sigma_{xz}','x axis','y axis','z axis',-2,4);

h15 =figure(15);
 mypatch2(x,ix,nnd,'interp','none',sig(:,6),1,1,'\sigma_{yz}','x axis','y axis','z axis',-2,4);
 
% disp('Saving...');
% print(h1,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/u_z','-dpng','-r600')
% print(h2,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/u_y','-dpng','-r600')
% print(h3,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/u_x','-dpng','-r600')
% print(h4,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/epsilon_{xx}','-dpng','-r600')
% print(h5,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/epsilon_{yy}','-dpng','-r600')
% print(h6,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/epsilon_{zz}','-dpng','-r600')
% print(h7,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/epsilon_{xy}','-dpng','-r600')
% print(h8,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/epsilon_{xz}','-dpng','-r600')
% print(h9,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/epsilon_{yz}','-dpng','-r600')
% print(h10,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/sigma_{xx}','-dpng','-r600')
% print(h11,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/sigma_{yy}','-dpng','-r600')
% print(h12,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/sigma_{zz}','-dpng','-r600')
% print(h13,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/sigma_{xy}','-dpng','-r600')
% print(h14,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/sigma_{xz}','-dpng','-r600')
% print(h15,'/Users/xiongzhihao/Documents/FiniteElementMethod-final/Result/Theory/sigma_{yz}','-dpng','-r600')

disp('done');