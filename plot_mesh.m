function h = plot_mesh( D, x, ix, nnde, cf, ce, fa)
% h = plot_mesh( D, x, ix, nnde, cf, ce, fa)
%-------------------------------------------------
% plot mesh of following elements
%   D1N2: bar
%   D2N3: triangle
%   D2N4: Quad
%   D3N4: Tetrahedron
%   D3N8: Cube
% Copyright @ W. Shan, 07.2013
%-------------------------------------------------
% input:
%   D: dimensionality, 1,2,3
%   x: nodal coordinates row vectors
%   ix: connectivity
%   nnde: number of nodes
%   cf  : face color
%   ce  : edge color
%   fa  : face alpha
% output:
%   h: handle of the mesh

switch D
    case 1
        D1  = size(x, 2); 
        nix = size(ix, 1);
        h   = zeros( nix, 1);
        hold on; 
        switch D1
            case 1
                for i = 1:nix
                    h(i) = plot( x(ix(i, :)), zeros(1, 2), 'o-', ...
                        'markeredgecolor', ce, 'markerfacecolor', cf, 'color', ce);
                end
            case 2
                for i = 1:nix
                    h(i) = plot( x(ix(i, :), 1), x(ix(i, :), 2), 'o-', ...
                        'markeredgecolor', ce, 'markerfacecolor', cf, 'color', ce);
                end
            case 3
                for i = 1:nix
                    h(i) = plot3( x(ix(i, :), 1), x(ix(i, :), 2), x(ix(i, :), 3), 'o-', ...
                        'markeredgecolor', ce, 'markerfacecolor', cf, 'color', ce);
                end
        end
        if ischar(cf)
            if strcmpi( cf, 'none')
                return
            end
        end
    case {2, 3}
        nnde_val = unique( nnde);
        ntype    = length( nnde_val);        
        h         = zeros( ntype, 1);
        for i = 1:ntype
            tmp = ix( nnde == nnde_val(i), 1:nnde_val(i));
            switch D
                case 2
                    faces = tmp;
                case 3
                    switch nnde_val(i)
                        case 4
                            faces = [tmp(:, [1 2 3]);
                                tmp(:, [1 3 4]);
                                tmp(:, [1 4 2]);
                                tmp(:, [2 3 4])];
                        case 8
                            faces = [tmp(:, [1 4 3 2]);
                                tmp( :, [1 2 6 5]);
                                tmp( :, [2 6 7 3]);
                                tmp( :, [3 4 8 7]);
                                tmp( :, [4 1 5 8]);
                                tmp( :, [5 6 7 8])];
                    end
            end
            hold on;
            h(i) = patch('vertices', x, 'faces', faces, 'facecolor', cf, ...
                'edgecolor', ce, 'facealpha', fa);            
        end
        if D == 3
            view( -30, 20);
        end
end