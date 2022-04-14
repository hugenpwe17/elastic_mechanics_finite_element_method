function [B, dNx] = UpdateB(J, dN)
% calculate B and dNx
% [B, dNx] = updateB(J, dN)
% Contributed by OuYang

%   B  : displacement-strain matrix (Ce x (nnde x D), Ce = 6 (3D)
%   dNx: dNx(i, j) = dN(i)/dx(j) (nnde x D)

% Calculate the node number and dimensions of a element
[nnde, D] = size(dN); 

% calculate dN
dNx       = dN / J;

%define B 
Ce = 6; 
B  = zeros(Ce, nnde * D); 

% calculater B
for i = 1:nnde
    B(:, ((i-1)*D + 1):(i*D)) ...
      =    [dNx(i, 1), 0, 0;
           0, dNx(i, 2), 0;
           0, 0, dNx(i, 3);
           dNx(i, 2), dNx(i, 1), 0;
           dNx(i, 3), 0, dNx(i, 1);
           0, dNx(i, 3), dNx(i, 2)]; 
end

end