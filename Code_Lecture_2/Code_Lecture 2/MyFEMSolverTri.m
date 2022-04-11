function dat = MyFEMSolverTri(dat)
% dat
% inputs:
%   .x   : nodal coords.; nx by D
%   .ix  : connectivity
%   .E   : Young's moudulus
%   .nu  : Poisson's ration
%   .pu2 : 
%   .u   : 
%   .fext:
%   .nintx: 
%   .nints: 
%
% output
%   .nx  : number of nodes
%   .D   : dimension

%   .SFx : shape function 
%   .SFs : shape function for recovery
%   .u   : displacement, nx * D by 1; 
%   .fext: nodal force, 

% 
[dat.nx, dat.D]     = size(dat.x); 
[dat.nix, dat.nnde] = size(dat.nix); 

% shape function structures
dat.SFx = generate_SF_tri(dat.nintx); 
dat.SFs1 = generate_SF_tri(dat.nints); 
dat.SFs2 = generate_SF_tri(1); 
% elastic tensor 
dat.CC = elastTensor(dat.D, dat.E, dat.nu);
dat.Dc = size(dat.CC, 1); 

% integrate K, M
[dat.K, dat.M] = intKM(dat.x, dat.ix, dat.SFx, dat.CC);

% integrate Mr, Reps, Rsig
[dat.Mr, dat.Reps, dat.Rsig] = intRecoverMat(dat.x, dat.ix, dat.SFs1, ...
    dat.SFs2, dat.CC);


% solve for displacement
[dat.u, dat.fext] = solveLin(dat.K, dat.u, dat.fext, dat.pu2);

% recover stress and strain
sigVec  = dat.Mr \ (dat.Rsig * dat.u); 
dat.sig = reshape(sigVec, dat.Dc, dat.nx)';
epsVec  = dat.Mr \ (dat.Reps * dat.u); 
dat.eps = reshape(epsVec, dat.Dc, dat.nx)'; 



[dat.errAbs, dat.errRel, dat.eng] = recoverError(dat.eps, dat.u, dat.SFx, ...
    dat.CC, dat.x, dat.ix); 

end