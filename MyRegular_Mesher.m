classdef MyRegular_Mesher < hgsetget
    properties(SetAccess = private)
        x;
        nx;
        ix;
        nix;
    end
    methods
        function obj = MyRegular_Mesher
        end
        function create(obj, ID, L, N, o)
            [obj.x, obj.nx, obj.ix, obj.nix] = regular_mesh( ID, L, N, o);
        end        
    end
end

function [x, nx, ix, nix] = regular_mesh( el_type, L, N, o)
% [x, nx, ix, nix] = regular_mesh( el_type, L, N, o)
%-------------------------------------------------------------------------
% generating regular mesh:
%       1D: BAR2
%       2D: TRI3 or QUAD4
%       3D: TET4 or HEX8
% Copyright @ W. Shan, 08.2013
%-------------------------------------------------------------------------
% inputs:
%   el_type: element type (string) 'BAR2', 'TRI3', 'QUAD4', 'TET4' or
%               'HEX8'
%   L: dimension of the domain [Lx, Ly, Lz]
%   N: number of seeds in each direction [Nx, Ny, Nz]
%   o: lower-left corner of the domain [ox, oy, oz]
% outputs
%   x: nodal coordinates, row vectors
%   nx: number of nodes
%   ix: element connectivity, row vectors
%   nix: number of elements


D = length( L);
seed = cell(1, D); % create seeds for each direction
for i = 1:D
    seed{i} = linspace( 0, 1, N(i));
end

p = cell(1, D);
if D == 1
    p{1} = seed{1};
else
    if D == 2
        [p{[2 1]}] = meshgrid( seed{[2 1]}); % generating nodal coordinates
    else
        [p{[2 1 3]}] = meshgrid( seed{[2 1 3]});
    end
end
nx = numel( p{i});
x  = zeros( nx, D); % formatting coordinates into matrix
for i = 1:D
    x(:, i) = p{i}(:);
end
clear p;

Nx = N(1);
if D > 1
    Ny = N(2);
    if D > 2
        Nz = N(3);
    end
end

is_tri = false;
area_fun = @(i) i;
switch el_type
    case 2
        err_D(D, 1);
        nix = Nx - 1;
        ix  = [1:nx-1; 2:nx]';
    case 3
        err_D(D, 2);
        ix = delaunay( x(:, 1), x(:, 2));
        nix = size( ix, 1);
        area_fun = @(i) tri_area(ix(i, :), x);
        is_tri   = true;
    case 4
        % Connectivity
        %           2---------1
        %           |         |
        %           |         h
        %           |         |
        %         o 3----w----4
        err_D(D, 2);
        nix = (Nx-1)*(Ny - 1);
        ix  = zeros( nix, 4);
        for j = 1:Ny-1
            tmp1 = (j-1)*(Nx-1);
            tmp2 = (j-1)*Nx;
            for i = 1:Nx-1
                id = i + tmp1;
                id1 = i + tmp2;
                ix( id, :) = [id1, id1 + 1, id1 + Nx + 1, id1 + Nx];
            end
        end
        ix = ix(:, [3, 4, 1, 2]);
    case 5
        err_D(D, 3);
        ix = delaunay( x(:, 1), x(:, 2), x(:, 3));
        nix = size( ix, 1);
        area_fun = @(i) tet_vol(ix(i, :), x);
        is_tri = true;
    case 8
        %               4-------------3
        %              /:-1,-1,1     /|1,-1,1
        %             / :           / |
        %            /-1,1,1       /1,1,1                  ^ z
        %           8-------------7   |                    |
        %           |   1.........|...2                    |
        %           |  / -1,-1,-1 |  / 1,-1,-1             o----> x
        %           | /           | /                     /
        %           |/-1,1,-1     |/                    |_ y
        %           5-------------6 1,1,-1
        err_D(D, 3);
        nelx = Nx - 1; nely = Ny - 1; nelz = Nz - 1;
        nix  = nelx*nely*nelz;
        ix   = zeros( nix, 8);
        for hei = 1:nelz
            tmp2 = (hei-1)*Nx*Ny;
            tmp3 = hei*Nx*Ny;
            tmp6 = (hei-1)*nelx*nely;
            for row = 1:nely
                tmp1 = (row-1)*Nx;
                tmp5 = (row-1)*nelx;
                for col = 1:nelx
                    iel = col + tmp5 + tmp6;
                    % Connectivity
                    tmp4 = [ col col+1 col+nelx+2 col+nelx+1] + tmp1;
                    ix(iel,1:4) = tmp4 + tmp2;
                    ix(iel,5:8) = tmp4 + tmp3;
                end
            end
        end
        %         ix = ix(:, [1, 2, 6, 5, 4, 3, 7, 8]);
    otherwise
        error( 'unknown element!');
end

if is_tri
    if D == 2
        swap = [1 3 2];
    else
        swap = [1 2 4 3];
    end
    id = 1:nix;
    a  = arrayfun(area_fun, id);
    is_swap = a < 0;
    ix(is_swap, :) = ix(is_swap, swap);
end

% scaling and shifting
for i = 1:D
    x(:, i) = x(:, i)*L(i) + o(i);
end
end

function err_D(D, i)
if D ~= i
    error( 'dimensionality dosen''t match element type');
end
end

function a = tri_area(ix, x)
col0 = ones(3, 1);
col12 = x(ix, :);
a = det([col0, col12]);
end

function v = tet_vol(ix, x)
col0 = ones(4, 1);
col13 = x(ix, :);
v = det([col0, col13]);
end