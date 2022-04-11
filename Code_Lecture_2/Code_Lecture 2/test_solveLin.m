close all; 
clear all; 
% material and model properties
nnde = 3; 
D    = 2; 
nint = 2;
E    = 10; 
nu   = 0.33; 

% % mesh grid
% L  = 100; 
% nL = 10; 
% 
% xL = linspace(0, L, nL); 
% 
% [px, py] = meshgrid(xL, xL); 
% 
% x  = [px(:), py(:)]; 
% nx = size(x, 1); 
% ix = delaunay(x(:, 1), x(:, 2)); 
% nix = size(ix, 1);
% 
% showforce = @(x, f, id) quiver(x(id, 1), x(id, 2), f(id, 1), f(id, 2)); 
% 
% figure; 
% h1 = mypatch(x, ix, 'none', 'k', []); 
% set(h1, 'linewidth', 2);
x=[0 0;
    1,0;
    0,1;];
ix=[1 2 3];
nx = size(x, 1);
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
[K, M] = intKMLoc(SF, CC, x);

% Boundary condition
% constrain;
% lowerID = find(abs(x(:, 2) - 0) < 1e-5); 
% topID   = find(abs(x(:, 2) - L) < 1e-5); 
% 
% hold on; 
% plot(x(lowerID, 1), x(lowerID, 2), 'ob', 'markerfacecolor', 'b'); 
% plot(x(topID, 1), x(topID, 2), 'ob', 'markerfacecolor', 'r'); 
% 
% pu2 = [lowerID * D - 1, lowerID * D]'; 
% pu2 = pu2(:); 
%
pu2 = [1 2 3 4];
% puf = topID * D; 
% 
% % disp. fext
u = zeros(nx * D, 1); 
u(pu2) = 0; 
% 
fext = zeros(nx * D, 1); 
%fext(puf) = 100 / length(topID); 




% solve
[u, fext] = solveLin(K, u, fext, pu2);

% 
% x1 = x + reshape(u, D, nx)'; 
% 
% h2 = mypatch(x1, ix, 'interp', 'none', u(1:2:end), 0.8); 
% set(h2, 'linewidth', 2);
% fextv = reshape(fext, D, nx)';
% showforce(x1, fextv, lowerID); 
% showforce(x1, fextv, topID); 
% 
% u_1    = zeros(nx * D, 1); 
% fext_1 = zeros(nx * D, 1); 
% 
% pu2_1b = [lowerID * D - 1, lowerID * D]'; 
% pu2_1t = topID * D; 
% u_1(pu2_1t) = 10; 
% pu2_1  = [pu2_1t(:); pu2_1b(:)]; 
% 
% [u_1, fext_1] = solveLin(K, u_1, fext_1, pu2_1);
% 
% x1_1    = x + reshape(u_1, D, nx)'; 
% fext_1v = reshape(fext_1, D, nx)'; 
% 
% 
% 
% figure; 
% hold on;
% mypatch(x, ix, 'none', 'k', []); 
% mypatch(x1_1, ix, 'interp', 'none', u_1(1:2:end), 0.8); 
% showforce(x1_1, fext_1v, lowerID); 
% showforce(x1_1, fext_1v, topID); 
% 
% 
disp('done');