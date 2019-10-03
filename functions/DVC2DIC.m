function [alldata] = DVC2DIC(inDir,inFile)
%infli is in .dat fromat
%% import DVC data
dataPath = [ inDir '\' inFile];
DVCresult = importdata(dataPath);
%% average non-zero data points 
% import x,y,z,ux and average ux data (non-zero)
C(:,1)=DVCresult.data(:,1);
C(:,2)=DVCresult.data(:,2);
C(:,3)=DVCresult.data(:,3);
C(:,4)=DVCresult.data(:,4);
indices = find(C(:,4)==0);
C(indices,4) = NaN;
[xy,~,j] = unique(C(:,[1 2]),'rows');
xyv = [xy accumarray(j,C(:,4),[size(xy,1),1],@nanmean)];
%import x,y,z,uy and average uy data (non-zero)
D(:,1)=DVCresult.data(:,1);
D(:,2)=DVCresult.data(:,2);
D(:,3)=DVCresult.data(:,3);
D(:,4)=DVCresult.data(:,5);
indices = find(D(:,4)==0);
D(indices,4) = NaN;
[ab,~,k] = unique(D(:,[1 2]),'rows');
abc = [ab accumarray(k,D(:,4),[size(ab,1),1],@nanmean)];
%% generate output data with x,y,z,average ux, average uy
outputdata(:,1)=xyv(:,1);
outputdata(:,2)=xyv(:,2);
outputdata(:,3)=xyv(:,3);
outputdata(:,4)=abc(:,3);
%% figure displacement map 
% 120 and 97 are unique size of average ux and uy, you need to convert x and y into matrix to obtain size of the matrix
%check the size of x and y by using size(unique(x); and size(unique(y));
%(in this case, the x=97 and y=120)
x=outputdata(:,1);
y=outputdata(:,2);
ux=outputdata(:,3);
uy=outputdata(:,4);
% xr=reshape(x,length(unique(y)),length(unique(x)));
xr=reshape(x,90,104);%(120,97 need to be changed accordingly based on your own data size)% 
% yr=reshape(y,length(unique(y)),length(unique(x)));
yr=reshape(y,90,104);
% uxr=reshape(ux,length(unique(y)),length(unique(x)));
uxr=reshape(ux,90,104);
% uyr=reshape(uyr,length(unique(y)),length(unique(x)));
uyr=reshape(uy,90,104);
% figure uyr (if imagesc doesn't work use function: scatter)
imagesc(x,y,uxr);
c = colorbar('FontSize',20);
xlabel('x (mm)','FontSize',20)
ylabel('y (mm)','FontSize',20)
c.Label.String = 'ux displacement (mm)';
title('U_x')
saveas(gcf,[inDir 'Av. Ux.png']); close all
imagesc(x,y,uyr);
c = colorbar('FontSize',20);
xlabel('x (mm)','FontSize',20)
ylabel('y (mm)','FontSize',20)
c.Label.String = 'uy displacement (mm)';
title('U_y')
saveas(gcf,[inDir 'Av. Uy.png']); close all
%% Added by Abdo to save DATA
alldata = [xr(:) yr(:) uxr(:) uyr(:)];