function [K, M] = IntKM(x, ix, SF, CC)

% parameters
% nx : nodes number 
% D  : dimensionality
[nx, D]     = size(x);
% nix  : cell number
% nnde : nodes number in a cell
[nix, nnde] = size(ix); 
% nf : all degrees of freedom
nf    = nx * D; 
% nfLoc : all degrees of freedom in a cell
nfLoc  = nnde * D;
nfLoc2 = nfLoc * nfLoc; 

% assmble global matrix
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
            
    [Kloc, Mloc] = IntKMLoc(SF, CC, xloc); 
    
    rowids(ids) = reshape(tmp, nfLoc2, 1); 
    colids(ids) = reshape(tmp', nfLoc2, 1);
    Kval(ids)   = Kloc(:); 
    Mval(ids)   = Mloc(:);
    ids         = ids + nfLoc2; 
    if mod(i, 10) == 0
        disp([num2str(i), '/', num2str(nix)]); 
    end
end

K = sparse(rowids, colids, Kval, nf, nf); 
M = sparse(rowids, colids, Mval, nf, nf); 

end
     
    