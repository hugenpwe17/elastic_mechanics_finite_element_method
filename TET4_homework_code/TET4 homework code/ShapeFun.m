function [N,dN] = ShapeFun(Xi)
% Shape function of tetrahedron
% Syntax: [N,dN] = ShapeFun(Xi)

%	N : Shape Function N=[N1,N2,N3,N4]
%   dN: dN(i,j)=dN(i)/dXi(j)
   
N  = [Xi(1),Xi(2),Xi(3),1-Xi(1)-Xi(2)-Xi(3)];
% Shape function    
dN = [1 0 0;
      0 1 0;
      0 0 1;
     -1 -1 -1];
% Divergence of all component of shape function in natural coordinates
end
% Contributed by OuYang