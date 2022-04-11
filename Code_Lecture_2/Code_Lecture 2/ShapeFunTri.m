function [N, dN] = ShapeFunTri(g)
% [N, dN] = ShapeFunTri(g)
% N: shape functions [N1, N2, N3]
% dN: dN / dg dN(i, j) = dN(i)/g(j)
N  = [g(1), g(2), 1 - g(1) - g(2)]; 
dN = [1 0;
    0 1; 
    -1 -1];
end