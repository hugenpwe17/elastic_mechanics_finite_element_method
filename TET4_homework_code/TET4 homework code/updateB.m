function [B, dNx] = updateB(J, dN)
% [B, dNx] = updateB(J, dN)
% B:   displacement-strain matrix (Ce x (nnde x D), Ce = 3 (2D) or 6 (3D)
% dNx: dNx(i, j) = dN(i)/dx(j) (nnde x D)
[nnde, D] = size(dN); 
dNx       = dN / J;

    Ce = 6; 
    B  = zeros(Ce, nnde * D); 
    for i = 1:nnde
        B(:, ((i-1)*D + 1):(i*D)) = [dNx(i, 1), 0, 0;
            0, dNx(i, 2), 0;
            0, 0, dNx(i, 3);
            dNx(i, 2), dNx(i, 1), 0;
            dNx(i, 3), 0, dNx(i, 1);
            0, dNx(i, 3), dNx(i, 2)]; 
    end
end