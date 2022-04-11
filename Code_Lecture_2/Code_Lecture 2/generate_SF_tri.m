function SF = generate_SF_tri(nint)
% SF = generate_SF_tri(nint)
% generate shape function data structure
% nint: Gauss. int. order
% SF : struct, 
%      .N  : shape function values matrix (ngp by nnde)
%      .dN : derivatives of N w.r.t. iso. coords (nnde by D by ngp)
%      .Vc : volume coeff.
%      .wg : weighting coeff. of GPs
nnde = 3; 
D    = 2; 
% Gauss. Pt
[g, wg] = tri_GP(nint); 
ngp     = length(wg); 

% Shape function info.
SF.N  = zeros(ngp, nnde); 
SF.dN = zeros(nnde, D, ngp); 
for i = 1:ngp
    [SF.N(i, :), SF.dN(:, :, i)] = ShapeFunTri(g(i, :)); 
end
SF.wg = wg; 
SF.Vc = 0.5; 
end