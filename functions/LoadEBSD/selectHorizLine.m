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
set(H1,'Name',['Select Horizontal Line (',V.type,'_',num2str(V.i),'_',...
    num2str(V.j),')'],'NumberTitle','off');
title('Click the Crack start and tip')
set(gca,'CLim',[min(min(image)) max(max(image))*5/8]); 
xlabel(['x-position [',V.unit,']']);
ylabel(['y-position [',V.unit,']']);
colormap jet

pos = get(gcf,'position');          set(gcf,'position',[100 100 pos(3:4)*2]) 

[X,Y] = ginput(2);
plot(X,Y)
hold off
close all

