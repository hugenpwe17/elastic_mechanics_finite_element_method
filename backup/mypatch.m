function h = mypatch(x, ix, cf, ce, cd, fa)
if nargin < 6
    fa = 1;
end
h = patch('vertices', x, 'faces', ix, 'facecolor', cf, 'edgecolor', ce, 'cdata', cd, 'facealpha', fa); 
end