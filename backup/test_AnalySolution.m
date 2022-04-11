%% parameters
% material parameters
E   = 600;
nu  = 1.33;
rho = 1;
% Geometric parameters
Lz  = max(x(:,3)) - min(x(:,3));
Ly  = max(x(:,2)) - min(x(:,2));
Lx  = max(x(:,1)) - min(x(:,1));
% outer field
g   = 9.8;

%% Analytical solution value
% sigma (only z,x and y is zero)
sig_z  = rho* g* x(:,3);
sig_x  = zeros(size(sig_z));
sig_y  = zeros(size(sig_z));
%epsilon (x y z)
eps_xx = -(nu* rho* g* x(:,3))/E;
eps_yy = -(nu* rho* g* x(:,3))/E;
eps_zz = (rho* g* x(:,3))/E;
%replacement 
u_x = -(nu* rho* g).* x(:,1).* x(:,3) /E;
u_y = -(nu* rho* g).* x(:,2).* x(:,3) /E;
u_z = (rho* g/(2*E)).* (x(:,3).^2- Lz^2+ nu*(x(:,1).^2+x(:,2).^2));
%current configuration's coordinates
xx = x+ [u_x,u_y,u_z];

%% Convert vertex matrix to face matrix
for i = 1 : length(ix)
[fx([(4*i-3):(4*i)],:),~] = VerToFace(x(ix(i,:),:),ix(i,:));
end

%% draw the figures
% frame
% h1 = patch('vertices', x, 'faces', fx, ...
% 'facecolor', 'none', 'edgecolor', 'r');
% hold on 
% axis equal

% sigma figure 
figure(1)
sig=[sig_x,sig_y,sig_z];
for i =1:3
subplot(1,3,i)
patch('vertices', x, 'faces', fx, ...
'facecolor', 'none', 'edgecolor', 'r');
hold on 
axis equal
mypatch(x, ix, 'interp', 'none', sig(:,i), 0.5); 
view(50,25)
colorbar
end

% epsilon figure 
eps = [eps_xx,eps_yy,eps_zz];
figure(2)
for i =1:3
subplot(1,3,i)
patch('vertices', x, 'faces', fx, ...
'facecolor', 'none', 'edgecolor', 'r');
hold on 
axis equal
mypatch(x, ix, 'interp', 'none', eps(:,i), 0.5); 
view(50,25)
colorbar
end

% replacement figure
u = [u_x,u_y,u_z];
figure(3)
for i =1:3
subplot(1,3,i)
patch('vertices', x, 'faces', fx, ...
'facecolor', 'none', 'edgecolor', 'r');
hold on 
axis equal
mypatch(xx, ix, 'interp', 'none', u(:,i), 0.5); 
view(50,25)
colorbar
end 