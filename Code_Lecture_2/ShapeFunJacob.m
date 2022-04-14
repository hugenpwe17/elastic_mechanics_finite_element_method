function [J, detJ] = ShapeFunJacob(dN, x)
% [J, detJ] = ShapeFunJacob(dN, x)
% J   : dx/dg
% detJ: determinant of J
% dN  : dN(i, j) = dN(i) / dg(j)
% x   : nodal coords. (number of nodes by D)
[nnde, D] = size(x); 
J         = zeros(D, D);
for i = 1:nnde
    J = J + x(i, :)' * dN(i, :); 
end
detJ = det(J); 
end