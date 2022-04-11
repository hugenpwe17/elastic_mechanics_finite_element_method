x =[0 0 0; 1 0 0; 0 1 0; 0 0 1];
ix = [1 2 3 4]; 

% x = [0.9157    0.0357    0.7577;
%     0.7922    0.8491    0.7431;
%     0.9595    0.9340    0.3922;
%     0.6557    0.6787    0.6555]*100;
% ix = [1 2 3 4]; 

%
nnde = 3; 
D    = 2; 
nint = 3;
E    = 10; 
nu   = 0.33; 

%
u = zeros(D*nnde,1);
fext = zeros(D*nnde,1);
%
nf  = length(u);
pu2 = [1 3 5];
pu1 = setxor((1:nf), pu2); 
%
u(pu2) = 0;
%
u2 = u(pu2); 
f1 = fext(pu1);
% replace inv(A)*b to A\b
u(pu1)    = K(pu1, pu1) \ (f1 - K(pu1, pu2) * u2); 
fext(pu2) = K(pu2, :) * u; 

disp('done')