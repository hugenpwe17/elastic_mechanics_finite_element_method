function [u, fext] = solveLin(K, u, fext, pu2)
% [u, fext] = solveLin(K, u, fext, pu)
% solve constrained linear system
% u   : vector form of nodal displacement (nx * D by 1)
% fext: vector form of external forces on nodes (nx * D by 1)
% K   : stiffness matrix (nx * D by nx * D)
% pu2 : constrained dofs

nf = length(u); 
pu1 = setxor((1:nf), pu2);

u2 = u(pu2); 
f1 = fext(pu1); 


u(pu1)    = K(pu1, pu1) \ (f1 - K(pu1, pu2) * u2); 
fext(pu2) = K(pu2, :) * u; 
end