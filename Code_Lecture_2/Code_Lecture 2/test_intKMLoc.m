close all; 
clear all; 
clc; 
% material and model properties
nnde = 3; 
D    = 2; 
nint = 3;
E    = 10; 
nu   = 0.33; 
% element info. 

% x = [0.8147, 0.9134
%     0.9058, 0.6324
%     0.4270, 0.3975] * 100; 
% ix = [2 1 3]; 

x = [1 0 ; 0 1 ; 0 0];
ix = [2 1 3];

% plot_mesh(D, x, ix, nnde, 'b', 'b', 1);


% Gauss. Pt
[g, wg] = tri_GP(nint); 
ngp     = length(wg); 
SF = generate_SF_tri(nint);
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

% K, M
[K, M] = intKMLoc(SF, CC, x(ix, :)); 

% 
disp('K = : '); 
disp(K); 
disp(eig(K)); 
% 
disp('done'); 