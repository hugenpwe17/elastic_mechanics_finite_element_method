function [K, M] = intKMLoc(SF, CC, x)
% [K, M] = intKMLoc(SF, CC, x)
% integrate stiffness matrix and mass matrix for each element
% x  : nodal coords. (nnde by D)
% SF : struct, 
%      .N  : shape function values matrix (ngp by nnde)
%      .dN : derivatives of N w.r.t. iso. coords (nnde by D by ngp)
%      .Vc : volume coeff.
%      .wg : weighting coeff. of GPs
% CC : elastic tensor matrix (Ce by Ce)

% K  : stiffness matrix ( (nnde * D) by (nnde * D))
% M  : mass matrix (same size K)

[nnde, D] = size(x);
ngp       = size(SF.N, 1); 
nf        = nnde * D; 

K = zeros(nf, nf); 
M = zeros(nf, nf); 

for i = 1:ngp        
    [J, detJ] = ShapeFunJacob(SF.dN(:, :, i), x);
    B         = updateB(J, SF.dN(:, :, i));
    N         = kron(SF.N(i, :), eye(D)); 
    dV = SF.wg(i) * detJ * SF.Vc; 
    K  = K + (B' * (CC * B)) * dV; 
    M  = M + (N' * N) * dV; 
end


end