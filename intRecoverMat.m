function [Mr, Reps, Rsig] = intRecoverMat(x, ix, SF1, SF2, CC)

% [Mr, Reps, Rsig] = intRecoverMat(x, ix, SF, CC)
% x : nodal coords. (nx by D)
% ix: mesh connectivity (nix by nnde)
% SF : struct, (SF1 for Mr, SF2 for Reps and Rsig)
%      .N  : shape function values matrix (ngp by nnde)
%      .dN : derivatives of N w.r.t. iso. coords (nnde by D by ngp)
%      .Vc : volume coeff.
%      .wg : weighting coeff. of GPs
% CC : elastic tensor matrix (Dc by Dc)
% Mr    : int(N'N)
% Resp  : int(N'B)
% Rsig  : int(N'CB)

[nx, D]     = size(x); 
[nix, nnde] = size(ix); 
Dc          = size(CC, 1); 

nfx    = nx * D; 
nfs    = nx * Dc;

nfLocx  = nnde * D; 
nfLocs  = nnde * Dc; 

nf2Mr   = nfLocs * nfLocs; 
nf2R    = nfLocs * nfLocx;  

% Mr
rowidsMr = zeros(nix * nf2Mr, 1); 
colidsMr = zeros(nix * nf2Mr, 1); 
valMr    = zeros(nix * nf2Mr, 1); 

% Rsig, Reps
rowidsR  = zeros(nix * nf2R, 1); 
colidsR  = zeros(nix * nf2R, 1); 
valRsig  = zeros(nix * nf2R, 1); 
valReps  = zeros(nix * nf2R, 1); 

idsMr = (1:nf2Mr)'; 
idsR  = (1:nf2R)'; 

pid2dof = @(pid, D) pid * D + (-(D-1):0)'; 

for i = 1:nix
    pid  = ix(i, :);
    xloc = x(pid, :);             
    dofsx = pid2dof(pid, D); 
    dofsx = dofsx(:); 
    dofss = pid2dof(pid, Dc); 
    dofss = dofss(:); 
    
            
    [Mrloc, Repsloc, Rsigloc] = intRecoverMatLoc(SF1, SF2, CC, xloc);
    
    % Mr
    tmpMr           = repmat(dofss, 1, nfLocs);     
    rowidsMr(idsMr) = reshape(tmpMr, nf2Mr, 1); 
    colidsMr(idsMr) = reshape(tmpMr', nf2Mr, 1);
    valMr(idsMr)    = Mrloc(:); 
    idsMr           = idsMr + nf2Mr; 
    
    % Reps, Rsig
    tmpRrow       = repmat(dofss, 1, nfLocx);     
    tmpRcol       = repmat(dofsx', nfLocs, 1); 
    rowidsR(idsR) = tmpRrow(:); 
    colidsR(idsR) = tmpRcol(:); 
    valReps(idsR) = Repsloc(:); 
    valRsig(idsR) = Rsigloc(:);
    idsR          = idsR + nf2R; 
    
%     if mod(i, 100) == 0
%         disp([num2str(i), '/', num2str(nix)]); 
%     end
end

Mr   = sparse(rowidsMr, colidsMr, valMr, nfs, nfs); 
Reps = sparse(rowidsR, colidsR, valReps, nfs, nfx); 
Rsig = sparse(rowidsR, colidsR, valRsig, nfs, nfx); 

end