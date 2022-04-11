function h = mypatch(x, ix, cf, ce, cd, fa, g, tl, xl, yl, zl,a,b)
% h = mypatch(x, ix, cf, ce, cd, fa, g, tl, xl, yl, zl)
% Contributed by Xiong

% x : coord. of point
% ix: order of point
% cf: facecolor
% ce: edgecolor
% fa: facealpha

% g : =1:x,y,zgrid
%     =0:none

% tl:title
% xl:xlabel
% yl:ylabel
% zl:zlabel
% a :colorbar upper limit
% b :colorbar lower limit

[~,D]=size(x(1,:));
h = patch('vertices', x, 'faces', ix, 'facecolor', cf, 'edgecolor', ce, 'cdata', cd, 'facealpha', fa);
axis equal;


xlabel(xl,'FontSize',13);
ylabel(yl,'FontSize',13);
title(tl,'FontSize',17);

if D==3
    zlabel(zl,'FontSize',13);
end

if g==1
    ha=gca;
    set(ha,'xgrid','on')
    set(ha,'ygrid','on')
    if D==3
        set(ha,'zgrid','on')
    end
end

if D==3
    view(35,20)
end
colorbar
if nargin >11
    
    caxis([a,b])
end
end