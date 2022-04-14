function [J, detJ] = ShapeFunJacob(dN, x)
% Jacobi matrix and det
% [J, detJ] = ShapeFunJacob(dN, x)

% dN  : dN(i, j) = dN(i) / dg(j)
% x   : nodal coords. (number of nodes by D)

% J   : dx/dg
% detJ: determinant of J

[D,nnde] = size(x); 
% Find the dimensionality and node's number in input data

J = zeros(D, D);
% Define jacobi matrix

for i = 1:nnde
    J = J + x(:,i)*dN(i, :); 
end
% Calculate jacobi matrix
detJ  = det(J); 
% Calculate det of jacobi matrix
end
% Contributed by Xiong