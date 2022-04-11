function CC = elastTensor(D, E, nu)
% CC = elastTensor(D, E, nu)
% CC: elastic tensor (Ce by Ce, Ce = 3 (2D) or 6 (3D)
% D : dimensionality
% E : Young's modulus
% nu: Poisson's ratio
mu = E / 2 / (1 + nu);
lm = E * nu / (1 + nu) / (1 - 2 * nu); 
switch D
    case 2
        CC = [...
            2 * mu + lm, lm, 0;
            lm, 2 * mu + lm, 0; 
            0, 0, mu]; 
    case 3
        CC = [...
            2 * mu + lm, lm, lm, 0, 0, 0;
            lm, 2 * mu + lm, lm, 0, 0, 0; 
            lm, lm, 2 * mu + lm, 0, 0, 0;
            0, 0, 0, mu, 0, 0;
            0, 0, 0, 0, mu, 0;
            0, 0, 0, 0, 0, mu];
    otherwise
        error('illegal dimensionality');
end
end