function h = mypatch(x, ix, cf, ce, cd)
h = patch('vertices', x, 'faces', ix, 'facecolor', cf, 'edgecolor', ce, 'cdata', cd); 
end