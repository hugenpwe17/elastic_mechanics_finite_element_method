function SF = GenerateShapeFunction(D,nnde,nint)
% Generate Shape Function TET4
% SF = GenerateShapeFunction(nint)
% Contributed by Xiong

% nint : Gauss piont int order
% SF   : 
%       .N  : shape function values matrix (ngp by nnde)
%       .dN : derivatives of N w.r.t. iso. coords (nnde by D by ngp)
%       .Vc : volume coeff.
%       .w : weighting coeff. of GPs 
SF.nnde  = nnde;
SF.D     = D;
[g,w] = Tet4Gp(nint);
ngp   = length(w);

SF.N  = zeros(ngp, SF.nnde); 
SF.dN = zeros(SF.nnde, SF.D, ngp); 

for i = 1:ngp
    [SF.N(i, :), SF.dN(:, :, i)] = ShapeFun(g(i, :));
end
SF.w = w; 
SF.Vc = 1/6; 
end