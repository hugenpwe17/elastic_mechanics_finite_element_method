function [Mr, Reps, Rsig] = intRecoverMatLoc(SF1, SF2, CC, x)
% [Mr, Reps, Rsig] = intRecoverMatLoc(SF1, SF2, CC, x)
% integrate matrices required for stress/strain recovery
% x  : nodal coords. (nnde by D)
% SF : struct, (SF1 for Mr, SF2 for Reps and Rsig)
%      .N  : shape function values matrix (ngp by nnde)
%      .dN : derivatives of N w.r.t. iso. coords (nnde by D by ngp)
%      .Vc : volume coeff.
%      .wg : weighting coeff. of GPs
% CC : elastic tensor matrix (Ce by Ce)
% Mr    : int(N'N)
% Resp  : int(N'B)
% Rsig  : int(N'CB)

[nnde, D] = size(x);
ngp1      = size(SF1.N, 1);
ngp2      = size(SF2.N, 1); 
Dc        = size(CC, 1); 
nfx       = nnde * D; 
nfs       = nnde * Dc; 


Mr   = zeros(nfs, nfs);
Rsig = zeros(nfs, nfx); 
Reps = zeros(nfs, nfx); 

for i = 1:ngp1        
    [~, detJ] = ShapeFunJacob(SF1.dN(:, :, i), x);    
    N         = kron(SF1.N(i, :), eye(Dc)); 
    
    dV   = SF1.wg(i) * detJ * SF1.Vc; 
    Mr   = Mr + (N' * N) * dV;     
end
for i = 1:ngp2
    [J, detJ] = ShapeFunJacob(SF2.dN(:, :, i), x);
    B         = updateB(J, SF2.dN(:, :, i));
    N         = kron(SF2.N(i, :), eye(Dc));    
    dV   = SF2.wg(i) * detJ * SF2.Vc; 
    Rsig = (N' * (CC * B)) * dV; 
    Reps = (N' * B) * dV;     
end

end