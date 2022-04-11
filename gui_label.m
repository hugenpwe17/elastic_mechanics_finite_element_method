function h = gui_label( ndim, x, con, nnde, type, id)
% h = gui_label( ndim, x, ix, nnde, typ, id)
%--------------------------
% plot labels             %
% Wenzhe Shan, 14/01/2009 %
% add id, 21/01/2009      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input:
%   x   : coordinates, in row vectors
%   con : connectivity matrix
%   ndim: dimension
%   nnde: number of nodes
%   id  : indices
%   type: 1: node, 2: edge, 3:face: 4: element
% output:
%   h : graphic handles

nid            = length( id);
str            = num2str( id(:));
if type > 1
    xmitte = zeros(nid, 1);    
    ymitte = zeros(nid, 1);
    if ndim > 2
        zmitte = zeros(nid, 1);
    end
    
end
if type == 1
    valign = 'Bottom';
else
    valign = 'Middle';
end    
switch type
    case 1
        % nodes
        xmitte = x(id, 1);
        if ndim > 1
            ymitte = x( id, 2);
            if ndim > 2
                zmitte = x(id, 3);
            end
        end        
    case 2
        % edges
        h = zeros( nid, 1);
        p1 = con(id, 1);
        p2 = con(id, 2);
        xmitte = (x(p1, 1) + x( p2, 1))/2;
        if ndim > 1
            ymitte = (x(p1, 2) + x( p2, 2))/2;
            if ndim > 2
                zmitte = (x(p1, 3) + x( p2, 3))/2;
            end
        end        
    case {3, 4}
        % faces, elements
        h = zeros( nid, 1);
        array_id = (1:nid)';
        if ndim == 1
            x = x(:);
        end
        xmid = @(i) mean( x( con( id(i), 1:nnde(i)), 1));
        xmitte = cell2mat( arrayfun(xmid, array_id, 'uniformoutput', false));        
        if ndim > 1
            ymid = @(i) mean( x( con( id(i), 1:nnde(i)), 2));
            ymitte = cell2mat( arrayfun(ymid, array_id, 'uniformoutput', false));
            if ndim > 2
                zmid = @(i) mean( x( con( id(i), 1:nnde(i)), 3));
                zmitte = cell2mat( arrayfun(zmid, array_id, 'uniformoutput', false));
            end
        end
end
switch ndim
    case {1, 2}
        if ndim == 1
            ymitte = zeros(size(xmitte));
        end
        h = text( xmitte, ymitte, str);
    case 3
        h = text( xmitte, ymitte, zmitte, str);
end
set(h, 'Color', 'b', 'FontSize', 10, 'HorizontalAlignment', 'Center',...
                    'VerticalAlignment', valign);
                