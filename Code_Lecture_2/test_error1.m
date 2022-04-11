close all; 
clear all; 
clc; 
% inline functions
showforce = @(x, f, id) quiver(x(id, 1), x(id, 2), f(id, 1), f(id, 2)); 

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

% constrain
pu2    = [lowerID * D - 1, lowerID * D]'; 
pu2    = pu2(:); 
u      = zeros(nx * D, 1); 
u(pu2)       = 0;


% external forces
puf  = topID * D; 
fext = zeros(nx * D, 1); 
fext(puf) = 100 / length(topID); 
% solve for displacement
[u, fext] = solveLin(K, u, fext, pu2);

% recover stress and strain
sigVec = Mr \ (Rsig * u); 
sig    = reshape(sigVec, Dc, nx)';
epsVec = Mr \ (Reps * u); 
eps    = reshape(epsVec, Dc, nx)'; 



[errAbs, errRel, eng] = recoverError(eps, u, SFx, CC, x, ix); 

errAbsTot = sum(errAbs); 
errRelTot = errAbsTot / sum(eng); 

disp(['Total absolute Error: ', num2str(errAbsTot)]); 
disp(['Total relative Error: ', num2str(errRelTot)]); 
disp(['Total Elastic Energy: ', num2str(0.5 * sum(eng))]); 

% deformed configuration; 
xu = x + reshape(u, D, nx)'; 

%
bisec = MyRecursiveBisect; 


figure; 
subplot(2, 2, 1); 
hold on; 
mypatch(x, ix, 'none', 'k', [], 1); 
mypatch(xu, ix, 'interp', 'interp', u(1:2:end), 0.8); 
colorbar; 
title('u_x'); 

subplot(2, 2, 2); 
hold on; 
mypatch(x, ix, 'none', 'k', [], 1); 
mypatch(xu, ix, 'interp', 'interp', u(2:2:end), 0.8); 
colorbar; 
title('u_y'); 

subplot(2, 2, 3); 
mypatch(x, ix, 'none', 'k', [], 1); 
mypatch(xu, ix, 'flat', 'none', eng, 0.8); 
title('2 x Energy'); 
colorbar;

subplot(2, 2, 4); 
errThred = 0.05; 
elidBad  = errRel > errThred; 
mypatch(x, ix, 'none', 'k', [], 1); 
mypatch(xu, ix, 'flat', 'none', errRel, 0.8); 
mypatch(xu, ix(elidBad, :), 'none', 'r', errRel(elidBad), 1); 
title('Relative Error'); 
colorbar;
caxis([0, 0.1]); 

[ix1, x1, u1, fext1] = bisec.refine(find(elidBad), ix, x, reshape(u, D, nx)', ...
    reshape(fext, D, nx)');  %#ok<*FNDSB>
xu1 = x1 + u1; 

figure;
subplot(1, 2, 1); 
hold on; 
mypatch(xu, ix, 'interp', 'k', u(1:2:end), 0.8); 
colorbar; 
title('u_x'); 

subplot(1, 2, 2); 
hold on; 

mypatch(xu1, ix1, 'interp', 'k', u1(:, 1), 0.8); 
mypatch(xu, ix(elidBad, :), 'none', 'r', errRel(elidBad), 1); 
colorbar; 
title('u_x'); 



nR  = 5; 
ixR  = cell(nR, 1); 
xR   = cell(nR, 1); 
uR   = cell(nR, 1); 
fR   = cell(nR, 1); 
errR = cell(nR, 1); 
epsR = cell(nR, 1); 
sigR = cell(nR, 1); 
engR = cell(nR, 1); 

disp('done');
