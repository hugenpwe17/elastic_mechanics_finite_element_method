function [g, w] = TET4_GP(nint)
% Gauss points and wights
% [g,w] = TET4_GP(nint)

% nint: order of gauss points

% g : gauss point coordinates
% w : corresponding weigh

switch nint
    case 1
        g = [1/4, 1/4, 1/4];
        w = 1;
    case 2
        a = 0.58541020;
        b = 0.13819660;
        g = [a, b, b;
             b, a, b; 
             b, b, a;
             b, b, b];
         w = 1/4 * ones(1, 4);
    case 3
        g = [1/4, 1/4, 1/4;
             1/2, 1/6, 1/6;
             1/6, 1/2, 1/6;
             1/6, 1/6, 1/2;
             1/6, 1/6, 1/6];
         w = [-4/5, 9/20, 9/20, 9/20, 9/20];
end
% Store gauss point data of order one to three
end
% Contributed by Xiong