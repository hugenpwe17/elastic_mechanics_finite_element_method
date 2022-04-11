function [K, M]=IntKMLoc(SF, CC, x)
% integrate stiffness matrix and mass matrix for each element
% [K, M]=IntKMLoc(SF, CC, x)
% Contributed by Xiong

% x  : nodal coords. (nnde by D)
% SF   : 
%  .N  : shape function values matrix (ngp by nnde)
%  .dN : derivatives of N w.r.t. iso. coords (nnde by D by ngp)
%  .Vc : volume coeff.
%   .w : weighting coeff. of GPs
% CC : elastic tensor matrix (Ce by Ce)

% K  : stiffness matrix ( (nnde * D) by (nnde * D))
% M  : mass matrix (same size K)

[nnde, D] = size(x);
Dc        = size(CC, 1); 
ngp       = size(SF.N, 1);
nf        = nnde * D;

K = zeros(nf, nf);
M = zeros(nf, nf);

for i = 1:ngp        
    [J, detJ] = ShapeFunJacob(SF.dN(:, :, i), x);
    B         = UpdateB(J, SF.dN(:, :, i));
    N         = kron(SF.N(i, :), eye(D));
    dV = SF.w(i) * detJ * SF.Vc;
    K  = K + (B' * (CC * B)) * dV; 
%         K  = (B' * (CC * B)) * dV; 
    M  = M + (N' * N) * dV; 
%         M  = (N' * N) * dV; 
end

end