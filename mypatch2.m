function h = mypatch2(x, ix,nnd, cf, ce,cd, fa, g, tl, xl, yl, zl,a,b)
% h = mypatch2(x, ix,nnd, cf, ce,cd, fa, g, tl, xl, yl, zl,a,b)
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
% nnd = nnde * ones(nix,1);
% h = patch('vertices', x, 'faces', ix, 'facecolor', cf, 'edgecolor', ce, 'cdata', cd, 'facealpha', fa);
h = plot_mesh1(3, x, ix,nnd, cf, ce,cd, fa);
axis([-0.5,0.5,-0.5,0.5,0,3.5])
set(gca,'FontSize',13);
axis equal;


xlabel(xl,'FontSize',15);
ylabel(yl,'FontSize',15);
title(tl,'FontSize',25);

if D==3
    zlabel(zl,'FontSize',15);

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
    view(36,16);
end
colorbar('FontSize',15);
if nargin >12    
    caxis([a,b]);
end
end