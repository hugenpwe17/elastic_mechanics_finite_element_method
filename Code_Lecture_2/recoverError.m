function [errAbs, errRel, eng] = recoverError(epsNode, u, SF, CC, x, ix)
% [errAbs, errRel] = recoverError(epsNode, u, SF, CC, x)
% errAbs : elementwise absolute error, nix by 1
% errRel : elementwise relative error, nix by 1
% eng    : elementwise elastic energy, nix by 1
% epsNode: recovered nodal strain, nx by Dc
% CC     : elastic tensor matrix, Dc by Dc
% x      : nodal coords., nx by D
% ix     : mesh connectivity, nix by nnde
% SF : struct, shape function info 
%      .N  : shape function values matrix (ngp by nnde)
%      .dN : derivatives of N w.r.t. iso. coords (nnde by D by ngp)
%      .Vc : volume coeff.
%      .wg : weighting coeff. of GPs
D   = size(x, 2); 
nix = size(ix, 1); 
ngp = size(SF.N, 1); 

errAbs = zeros(nix, 1); 
eng    = zeros(nix, 1); 

pid2dof = @(pid) pid(:) * D + (-(D-1):0); 

for i = 1:nix
    pid  = ix(i, :); 
    dof  = pid2dof(pid)'; 
    dof  = dof(:);
    xloc   = x(pid, :); 
    uloc   = u(dof); 
    epsLoc = epsNode(pid, :);  % nnde by Dc
    % int((eps - eps1)^T CC (eps-eps1))
    errAbsLoc = 0; 
    engLoc    = 0; 
    for j = 1:ngp
        [J, detJ] = ShapeFunJacob(SF.dN(:, :, j), xloc);
        B         = updateB(J, SF.dN(:, :, j));
        dV        = SF.wg(j) * SF.Vc * detJ; 
        eps1      = B * uloc;  % Dc by 1
        eps       = SF.N(j, :) * epsLoc; % 1 by Dc
        deps      = eps' - eps1; 
        engLoc    = engLoc + (eps1' * (CC * eps1)) * dV; 
        errAbsLoc = errAbsLoc + (deps' * (CC * deps)) * dV; 
    end
    eng(i)    = engLoc; 
    errAbs(i) = errAbsLoc; 
end

eng    = sqrt(eng); 
errAbs = sqrt(errAbs); 
errRel = errAbs ./ eng; 

end