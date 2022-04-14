% x =[0 0 0; 1 0 0; 0 1 0; 0 0 1];
% ix = [1 2 3 4]; 

x = [0.9157    0.0357    0.7577;
    0.7922    0.8491    0.7431;
    0.9595    0.9340    0.3922;
    0.6557    0.6787    0.6555]*100;
ix = [1 2 3 4]; 

%
nnde = 4; 
D    = 3; 
nint = 3;
E    = 10; 
nu   = 0.33; 

%
u = zeros(D*nnde,1);
fext = zeros(D*nnde,1);
%
nf  = length(u);
pu2 = [1 2 3 4 5 7];
pu1 = setxor((1:nf), pu2); 
%
u(pu2) = 0;
fext(8) = 0.5;
%
u2 = u(pu2); 
f1 = fext(pu1);
% replace inv(A)*b to A\b
u(pu1)    = K(pu1, pu1) \ (f1 - K(pu1, pu2) * u2); 
fext(pu2) = K(pu2, :) * u; 

fext1 = reshape(fext,[3,4])';
u1    = reshape(u,[3,4])';

xu   = x+u1;
fx   = VerToFace(x,ix);
patch('vertices', x, 'faces', fx, 'facecolor', 'none', 'edgecolor', 'b');
hold on
quiver3(x(:,1),x(:,2),x(:,3),fext1(:,1),fext1(:,2),fext1(:,3));

disp('done')