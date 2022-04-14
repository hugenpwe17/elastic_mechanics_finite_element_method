close all; 
clear; 
clc; 

% material and model properties
nnde  = 3; 
D     = 2; 
nintx = 2;    % int. order for K, M 
nints = 3;    % int. order for recovery matarices (Mr, Reps, Rsig)
E     = 10; 
nu    = 0.33; 

% mesh grid
L  = 100; 
nL = 10; 
xL = linspace(0, L, nL); 
[px, py] = meshgrid(xL, xL); 
x  = [px(:), py(:)]; 
nx = size(x, 1); 
ix = delaunay(x(:, 1), x(:, 2)); 
nix = size(ix, 1);

% shape function structures
SFx = generate_SF_tri(nintx); 
SFs1 = generate_SF_tri(nints); 
SFs2 = generate_SF_tri(1); 
% elastic tensor 
CC = elastTensor(D, E, nu);
Dc = size(CC, 1); 

% integrate K, M
tic;
[K, M] = intKM(x, ix, SFx, CC);
disp(['integration of K, M finished for           :', num2str(toc), 'sec']);

% integrate Mr, Reps, Rsig
tic;
[Mr, Reps, Rsig] = intRecoverMat(x, ix, SFs1, SFs2, CC);
disp(['integration of Mr, Reps, Rsig finished for :', num2str(toc), 'sec']);

% Boundary conditions
lowerID = find(abs(x(:, 2) - 0) < 1e-5); 
topID   = find(abs(x(:, 2) - L) < 1e-5); 
leftID  = find(abs(x(:, 1) - 0) < 1e-5); 
rightID = find(abs(x(:, 1) - L) < 1e-5); 

o = min(x, [], 1); 
F = [1.02, 0.015; 
    0.075, 0.98]; 

bcids = unique([lowerID; topID; leftID; rightID]); 
xbc   = x(bcids, :); 
ubc   = (xbc - o) * F' + o - xbc;


% constrain
pu2    = [bcids * D - 1, bcids * D]'; 
pu2    = pu2(:); 
u      = zeros(nx * D, 1); 
u(pu2) = reshape(ubc', length(pu2), 1);

% external forces
fext = zeros(nx * D, 1); 

% solve for displacement
[u, fext] = solveLin(K, u, fext, pu2);

% recover stress and strain
sigVec = Mr \ (Rsig * u); 
sig    = reshape(sigVec, Dc, nx)';
epsVec = Mr \ (Reps * u); 
eps    = reshape(epsVec, Dc, nx)'; 

% deformed configuration; 
x1 = x + reshape(u, D, nx)'; 


figure; 
subplot(2, 4, 1); 
mypatch(x1, ix, 'interp', 'none', u(1:2:end), 0.8); 
colorbar; 
title('u_x'); 

subplot(2, 4, 5); 
mypatch(x1, ix, 'interp', 'none', u(2:2:end), 0.8); 
colorbar; 
title('u_y'); 

subplot(2, 4, 2); 
mypatch(x1, ix, 'interp', 'none', eps(:, 1), 0.8); 
colorbar; 
title('\epsilon_{xx}'); 

subplot(2, 4, 3); 
mypatch(x1, ix, 'interp', 'none', eps(:, 2), 0.8); 
colorbar; 
title('\epsilon_{yy}'); 

subplot(2, 4, 4); 
mypatch(x1, ix, 'interp', 'none', eps(:, 3) / 2.0, 0.8); 
colorbar; 
title('\epsilon_{xy}'); 

subplot(2, 4, 6); 
mypatch(x1, ix, 'interp', 'none', sig(:, 1), 0.8); 
colorbar; 
title('\sigma_{xx}'); 

subplot(2, 4, 7); 
mypatch(x1, ix, 'interp', 'none', sig(:, 2), 0.8); 
colorbar; 
title('\sigma_{yy}'); 

subplot(2, 4, 8); 
mypatch(x1, ix, 'interp', 'none', sig(:, 3), 0.8); 
colorbar; 
title('\sigma_{xy}'); 

disp('done');
