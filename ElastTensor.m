function CC = ElastTensor(E,nu)
% Elastic Tensor of 3D
% Syntax: CC = ElastTensor(E,nu)
% Contributed by OuYang

%   E  : Young's modulus
%   nu : Poisson's ratio

%   CC : ElastTensor (Ce by Ce),Ce = 6 (3D)

% Lame constants
mu = E/(2*(1+nu));
lm = E*nu/((1+nu)*(1-2*nu)); 

% Calculate elastic tensor
CC = [...
     2*mu+lm, lm, lm, 0, 0, 0;
     lm, 2*mu+lm, lm, 0, 0, 0; 
     lm, lm, 2*mu+lm, 0, 0, 0;
     0, 0, 0, mu, 0, 0;
     0, 0, 0, 0, mu, 0;
     0, 0, 0, 0, 0, mu];
 
end
