function [N,dN] = ShapeFun(Xi)
% Shape function of tetrahedron
% Syntax: [N,dN] = ShapeFun(Xi)
% Contributed by OuYang

%	N : Shape Function N=[N1,N2,N3,N4]
%   dN: dN(i,j)=dN(i)/dXi(j)
  
% Shape function    
N  = [Xi(1),Xi(2),Xi(3),1-Xi(1)-Xi(2)-Xi(3)];

% Divergence of all component of shape function in natural coordinates
dN = [1 0 0;
      0 1 0;
      0 0 1;
     -1 -1 -1];

end
