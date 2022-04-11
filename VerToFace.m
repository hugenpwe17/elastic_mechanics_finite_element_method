function [fx,fxs] = VerToFace(x,ix)
% vertex to face matrix
% [fx,fxs] = VerToFace(x,ix)
% Contributed by OuYang

% x:  point's coordinates
% ix: point's name code

% fxs: area matrix in cell
% fx : area matrix in global

% ix number 
innde = length(ix);
% get all possible point set that can form a plane
tmp    = perms((1:1:innde));
% delete a column as three points form a plane 
tmp(:,1) = [];
% a counter 
count = 1;
% midpoint of a cell
ap = sum(x,1)/length(ix);
% find the node set that make up the positive surface
% method:
% subtract the selected node set from the midpoint to form a vector set,
% and perform the mixed product in order.if it is positive, it is the node 
% set we need.
for i = 1:length(tmp)
    vt = x(tmp(i,:),:)-ap;
    jug = dot(vt(1,:),cross(vt(2,:),vt(3,:)));
    if jug > 0
        tmp2(count,:)=tmp(i,:);
        count =count +1;
    end
end
% delete duplicates
[~,n] = unique(sort(tmp2,2),'row');
% area matrix in a single cell
fxs = tmp2(n,:);
% area matrix in global
fx  = ix(tmp2(n,:));
end