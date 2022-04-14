close all; 
clear all; 
clc; 
% material and model properties
nnde = 4; 
D    = 3; 
nint = 3;
E    = 10; 
nu   = 0.33; 
% element info. 

x = [0.9157    0.0357    0.7577;
    0.7922    0.8491    0.7431;
    0.9595    0.9340    0.3922;
    0.6557    0.6787    0.6555]*100;
ix = [1 2 3 4]; 

% x =[0 0 0; 1 0 0; 0 1 0; 0 0 1];
% ix = [1 2 3 4]; 

% Gauss. Pt
[g, w] = Tet4Gp(nint);
ngp    = length(w); 

SF = GenerateShapeFunction(D,nnde,nint);

% elastic tensor 
CC = ElastTensor(E,nu);

% K, M
[K, M]=IntKMLoc(SF, CC, x(ix, :));

% 
disp('K = : '); 
disp(K); 
disp('eig = : ')
disp(eig(K)); 
% 
disp('done'); 