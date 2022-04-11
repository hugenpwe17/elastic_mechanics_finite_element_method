function [K, M] = intKM(x, ix, SF, CC)
% [K, M] = intKM(x, ix, SF, CC)
% x : nodal coords. (nx by D)
% ix: mesh connectivity (nix by nnde)
% SF : struct, 
%      .N  : shape function values matrix (ngp by nnde)
%      .dN : derivatives of N w.r.t. iso. coords (nnde by D by ngp)
%      .Vc : volume coeff.
%      .wg : weighting coeff. of GPs
% CC : elastic tensor matrix (Ce by Ce)
% K  : stiffness matrix ( (nx * D) by (nx * D)), sparse
% M  : mass matrix (same size K), sparse
[nx, D]     = size(x);      % [number of nodes , dimensionality]
[nix, nnde] = size(ix);     % [number of cells , number of nodes in cell]

nf     = nx * D;            % degrees of freedom of all nodes
nfLoc  = nnde * D;          % degrees of freedom of all nodes where one cell
nfLoc2 = nfLoc * nfLoc; 

rowids = zeros(nix * (nnde * D)^2, 1); 
colids = rowids; 
Kval   = rowids; 
Mval   = rowids; 
ids    = (1:nfLoc2)'; 

for i = 1:nix
    pid  = ix(i, :);
    xloc = x(pid, :);             
    dofs = pid * D + (-(D-1):0)'; 
    dofs = dofs(:); 
    tmp  = repmat(dofs, 1, nfLoc); 
            
    [Kloc, Mloc] = intKMLoc(SF, CC, xloc); 
    
    rowids(ids) = reshape(tmp, nfLoc2, 1); 
    colids(ids) = reshape(tmp', nfLoc2, 1);
    Kval(ids)   = Kloc(:); 
    Mval(ids)   = Mloc(:);
    ids         = ids + nfLoc2; 
    if mod(i, 10000) == 0
        disp([num2str(i), '/', num2str(nix)]); 
    end
end

K = sparse(rowids, colids, Kval, nf, nf); 
M = sparse(rowids, colids, Mval, nf, nf); 

end