function [J, detJ] = ShapeFunJacob(dN, x)
% Jacobi matrix and det
% [J, detJ] = ShapeFunJacob(dN, x)
% Contributed by Xiong

% dN  : dN(i, j) = dN(i) / dg(j)
% x   : nodal coords. (number of nodes by D)

% J   : dx/dg
% detJ: determinant of J

% % Calculate jacobi matrix
% J = x'*dN; 
% 
% % Calculate det of jacobi matrix
% detJ  = det(J); 

[nnde, D] = size(x); 
J         = zeros(D, D);
for i = 1:nnde
    J = J + x(i, :)' * dN(i, :); 
end
detJ = det(J); 
end
