classdef MyRecursiveBisect < hgsetget
    methods
        function [ix, x, u, r] = refine(~, id, ix, x, u, r)
            nnde = size(ix, 2); 
            nix  = size(ix, 1); 
            nx   = size(x, 1);             
            if isempty(id)
                return;
            end
            switch nnde
                case 3
                    [ix, ~, x, ~, u, r] = recr_bisc_tri3(id, ix, nix, x, nx, u, r); 
                case 4
                    [ix, ~, x, ~, u, r] = recr_bisc_tet4(id, ix, nix, x, nx, u, r); 
                otherwise
                    error('Recursive refinement is only valid for TRI3 or TET4 elements!'); 
            end
        end
    end
end

function [ix, nix, x, nx, u, r] = recr_bisc_tet4( id, ix, nix, x, nx, u, r)
% [ix, nix, x, nx, u, r] = recr_bisc_tet4( id, ix, nix, x, nx, u, r)
% ------------------------------------------------------------------
%   Recursively bisect the longest edge of targeted elements
%   WZ. Shan, 03/03/2009
%   ----------------------------------------------------------------
%       Input:
%           id: indices of targeted elements
%           ix: element connectivity, [p1, p2, p3, p4]
%           nix : number of elements
%           x   : nodal coordinates, [x, y, z]
%           nx  : number of nodes
%           u   : nodal displacements, [ux, uy, uz]
%           r   : nodal force, [rx, ry, rz]
%       Output:
%           updated input value, See above

nid = length( id);
disp( ['Refine ', num2str( nid), ' elements...']);
tic;
% find longest edges of targeted elements
A = zeros( length( id), 2);
edg_id = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
ix_tmp = ix( id, :);
for i = 1:length( id)
    l = x( ix_tmp(i, edg_id(:, 1)), :) - x( ix_tmp(i, edg_id(:, 2)), :);
    dl = l(:, 1).^2 + l(:, 2).^2 + l(:, 3).^2;
    [~, I] = max( dl);
    A( i, :) = ix_tmp( i, edg_id( I, :));
end
A = sort(A, 2);

set( 0, 'RecursionLimit', max(size(A, 1)*2, 2000));
% recursively refine point set A
[ix, nix, x, nx, u, r] = recursive( A, ix, nix, x, nx, u, r);

t = toc;
disp( ['Refinement finished for ', num2str( t), ' secs.']);
set( 0, 'RecursionLimit', 500);
end

function [ix, nix, x, nx, u, r] = recursive( A, ix, nix, x, nx, u, r)
while ~isempty( A)
    % sharing elements of A(1, :)
    ID = ismember( ix, A(1, :));
    ID = sum( ID, 2);
    B  = find( ID == 2);
    nB = length( B);
    % check for imcompatible elements and recursively refine them
    edg_id = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    lmax      = zeros( nB, 1);
    lmax_id   = zeros( nB, 1);
    A1        = zeros( nB, 2);
    
    for i = 1:nB
        l = x( ix( B(i), edg_id(:, 1)), :) - x( ix( B(i), edg_id(:, 2)), :);
        dl = l(:, 1).^2 + l(:, 2).^2 + l(:, 3).^2;
        [lmax(i), lmax_id(i)] = max( dl);
        A1(i, :)              = ix( B(i), edg_id( lmax_id(i), :));
    end
    
    id_add = abs(lmax - norm( x( A(1, 1), :) - x( A(1, 2), :))^2) > 1E-5;
    if any( id_add)
        % recursive call
        A      = [unique(sort(A1( id_add, :), 2), 'rows'); A]; %#ok<*AGROW>
%         [ix, nix, x, nx, u, r] = recursive( A, ix, nix, x, nx, u, r);
        continue
    else
        % bisection
        x_m = .5 * ( x( A(1, 1), :) + x( A(1, 2), :));
        if ~isempty( u)
            u_m = .5 * ( u( A(1, 1), :) + u( A(1, 2), :));
            u   = [u; u_m];
        end
        if ~isempty( r)
            r_m = .5 * ( r( A(1, 1), :) + r( A(1, 2), :));
            r   = [r; r_m];
        end
        % --bisect all neighbors
        elad = zeros( 2*nB, 4);
        for i = 1:nB
            nodes = ix( B( i), :);
            % left child
            elad( i*2-1, :) = nodes;
            elad( i*2-1, nodes == A(1, 2)) = nx + 1;

            % right child
            elad( i*2, :) = nodes;
            elad( i*2, nodes == A(1, 1))   = nx + 1;
        end
        % update
        x = [x; x_m];        
        nx = nx + 1;
        ix( B, :) = [];
        ix        = [ix; elad];
        nix       = nix + nB;
        % remove A(1, :)
        rm_id = (A(:, 1) == A( 1, 1)) & (A(:, 2) == A(1, 2));
        A(rm_id, :) = [];    
    end
end
end

function [ix, nix, x, nx, u, r] = recr_bisc_tri3( id, ix, nix, x, nx, u, r)

nid = length( id);
disp( ['Refine ', num2str( nid), ' elements...']);
tic;
% find longest edges of targeted elements
A = zeros( length( id), 2);
edg_id = [1 2; 1 3; 2 3];
ix_tmp = ix( id, :);
for i = 1:length( id)
    l = x( ix_tmp(i, edg_id(:, 1)), :) - x( ix_tmp(i, edg_id(:, 2)), :);
    dl = l(:, 1).^2 + l(:, 2).^2;
    [~, I] = max( dl);
    A( i, :) = ix_tmp( i, edg_id( I, :));
end
A = sort(A, 2);

set( 0, 'RecursionLimit', max(size(A, 1)*2, 2000));
% recursively refine point set A
[ix, nix, x, nx, u, r] = recursive_2D( A, ix, nix, x, nx, u, r);

t = toc;
disp( ['Refinement finished for ', num2str( t), ' secs.']);
set( 0, 'RecursionLimit', 500);
end

function [ix, nix, x, nx, u, r] = recursive_2D( A, ix, nix, x, nx, u, r)
while ~isempty( A)
    % sharing elements of A(1, :)
    ID = ismember( ix, A(1, :));
    ID = sum( ID, 2);
    B  = find( ID == 2);
    nB = length( B);
    % check for imcompatible elements and recursively refine them
    edg_id = [1 2; 1 3; 2 3];
    lmax      = zeros( nB, 1);
    lmax_id   = zeros( nB, 1);
    A1        = zeros( nB, 2);
    
    for i = 1:nB
        l = x( ix( B(i), edg_id(:, 1)), :) - x( ix( B(i), edg_id(:, 2)), :);
        dl = l(:, 1).^2 + l(:, 2).^2;
        [lmax(i), lmax_id(i)] = max( dl);
        A1(i, :)              = ix( B(i), edg_id( lmax_id(i), :));
    end
    
    id_add = abs(lmax - norm( x( A(1, 1), :) - x( A(1, 2), :))^2) > 1E-5;
    if any( id_add)
        % recursive call
        A      = [unique(sort(A1( id_add, :), 2), 'rows'); A]; %#ok<*AGROW>
%         [ix, nix, x, nx, u, r] = recursive( A, ix, nix, x, nx, u, r);
        continue
    else
        % bisection
        x_m = .5 * ( x( A(1, 1), :) + x( A(1, 2), :));
        if ~isempty( u)
            u_m = .5 * ( u( A(1, 1), :) + u( A(1, 2), :));
            u   = [u; u_m];
        end
        if ~isempty( r)
            r_m = .5 * ( r( A(1, 1), :) + r( A(1, 2), :));
            r   = [r; r_m];
        end
        % --bisect all neighbors
        elad = zeros( 2*nB, 3);
        for i = 1:nB
            nodes = ix( B( i), :);
            % left child
            elad( i*2-1, :) = nodes;
            elad( i*2-1, nodes == A(1, 2)) = nx + 1;

            % right child
            elad( i*2, :) = nodes;
            elad( i*2, nodes == A(1, 1))   = nx + 1;
        end
        % update
        x = [x; x_m];        
        nx = nx + 1;
        ix( B, :) = [];
        ix        = [ix; elad];
        nix       = nix + nB;
        % remove A(1, :)
        rm_id = (A(:, 1) == A( 1, 1)) & (A(:, 2) == A(1, 2));
        A(rm_id, :) = [];    
    end
end
end