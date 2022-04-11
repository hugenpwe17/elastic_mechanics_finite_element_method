close all;
clear all;
clc;

%选择高斯点阶数
%nint = 1;
%绘制四面体
x  = [0 0 0; 
    1 0 0; 
    0 1 0;
    0 0 1];%四面体定点坐标

ix=[1 2 3 4 1 3 4 2];%四面体连线顺序

 figure; 
 patch('vertices', x, 'faces', ix, 'facecolor', 'none', 'edgecolor', 'b');%四面体绘制

 nint=1:3;
 %绘制高斯点
sty  = {'x', 'o', '^'}; 
clr  = {'r', 'b', 'k'}; 
hold on;%将三阶图的高斯点画一起
for i = 1:length(nint)
    [g, w] = TET4_GP(nint(i)); 
    plot3(g(:, 1), g(:, 2), g(:, 3), 'marker', sty{i}, 'color', clr{i}, 'linestyle', 'none', ...
        'markersize', 16); 
end
view(10, 45);%视图旋转