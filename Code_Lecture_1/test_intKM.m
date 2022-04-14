close all; 
clear all; 
% material and model properties
nnde = 3; 
D    = 2; 
nint = 2;
E    = 10; 
nu   = 0.33; 

% mesh grid
L  = 100; 
nL = 10; 

xL = linspace(0, L, nL); 

[px, py] = meshgrid(xL, xL); 

x  = [px(:), py(:)]; 
nx = size(x, 1); 
ix = delaunay(x(:, 1), x(:, 2)); 
nix = size(ix, 1);

figure; 
h1 = mypatch(x, ix, 'none', 'k', []); 
set(h1, 'linewidth', 2);
% Gauss. Pt
[g, wg] = tri_GP(nint); 
ngp     = length(wg); 

% Shape function info.
SF.N  = zeros(ngp, nnde); 
SF.dN = zeros(nnde, D, ngp); 
for i = 1:ngp
    [SF.N(i, :), SF.dN(:, :, i)] = ShapeFunTri(g(i, :)); 
end
SF.wg = wg; 
SF.Vc = 0.5; 

% elastic tensor 
CC = elastTensor(D, E, nu);

%
tic;
[K, M] = intKM(x, ix, SF, CC);
toc;

% Boundary condition
% constrain;
lowerID = find(abs(x(:, 2) - 0) < 1e-5); 
topID   = find(abs(x(:, 2) - L) < 1e-5); 

hold on; 
plot(x(lowerID, 1), x(lowerID, 2), 'ob', 'markerfacecolor', 'b'); 
plot(x(topID, 1), x(topID, 2), 'ob', 'markerfacecolor', 'r'); 

pu2 = [lowerID * D - 1, lowerID * D]'; 
pu2 = pu2(:); 

puf = topID * D; 

% disp. fext
u = zeros(nx * D, 1); 
u(pu2) = 0; 

fext = zeros(nx * D, 1); 
fext(puf) = 100 / length(topID); 

% solve
[u, fext] = solveLin(K, u, fext, pu2);

% 
x1 = x + reshape(u, D, nx)'; 

h2 = mypatch(x1, ix, 'interp', 'none', u(1:2:end)); 
%set(h2, 'linewidth', 2);
% 
disp('done');