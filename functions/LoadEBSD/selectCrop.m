function [Xcrop,Ycrop,fig] = selectCrop(X,Y,Z,lineX,lineY,V)
% Plots Z (defined at co-ordinates [X],[Y]) and the associated horizontal
% line defined by lineX and lineY

% Asks for user input to crop a region

Xvec = X(1,:);
Yvec = Y(:,1);
fig = figure;
hold on
imagesc(Xvec,Yvec,Z)
set(gca,'YDir','normal')
axis equal
xlim([min(Xvec) max(Xvec)])
ylim([min(Yvec) max(Yvec)])
colorbar
str = ['Select crop boundary (',V.type,'_',num2str(V.i),'_',num2str(V.j),'^n^e^w)'];
title(str)
xlabel(['x-position [',V.unit,']']);
ylabel(['y-position [',V.unit,']']);
plot(lineX,lineY,'color','y')
hold off
colormap jet

set(fig,'Name','Crop required data','NumberTitle','off');
[Xcrop,Ycrop] = ginput(2);
Xcrop = [min(Xcrop);max(Xcrop)];
Ycrop = [min(Ycrop);max(Ycrop)];
hold on
plot([Xcrop(1) Xcrop(2) Xcrop(2) Xcrop(1) Xcrop(1)],...
    [Ycrop(1) Ycrop(1) Ycrop(2) Ycrop(2) Ycrop(1)],'color','r')
hold off
set(fig,'Name','Cropped region shown in red','NumberTitle','off');