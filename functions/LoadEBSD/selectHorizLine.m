function [X,Y] = selectHorizLine(x,y,image,V)

Xvec = x(1,:);
Yvec = y(:,1);

H1 = figure;
hold on
imagesc(Xvec,Yvec,image)
set(gca,'YDir','normal')
axis equal
xlim([min(Xvec) max(Xvec)])
ylim([min(Yvec) max(Yvec)])
colorbar
str = ['Select Horizontal Line (',V.type,'_',num2str(V.i),'_',num2str(V.j),')'];
title(str)

xlabel(['x-position [',V.unit,']']);
ylabel(['y-position [',V.unit,']']);
colormap jet

[X,Y] = ginput(2);
plot(X,Y)
hold off

