function [g, w] = tri_GP(nint)
% [g, w] = tri_GP(nint)
% g: iso coords of Gauss. points
% w: weighting coeffs.
% nint: integration order
assert(nint > 0 || nint < 4,'integration order must be within 1 and 4');
switch nint 
    case 1 
        g = [1/3, 1/3]; 
        w = 1;
    case 2
        g = [1/2, 1/2; 
            0, 1/2; 
            1/2, 0]; 
        w = 1/3 * ones(1, 3); 
    case 3
        g = [1/3, 1/3;
            0.6, 0.2;
            0.2, 0.6;
            0.2, 0.2]; 
        w = [-27/48, 25/48, 25/48, 25/48];     
end
end