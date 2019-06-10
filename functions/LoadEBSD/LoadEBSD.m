function [ EBSDdata,input ] = LoadEBSD( input )
%LOADEBSD Summary of this function goes here
%   Detailed explanation goes here

% LoadEBSD.m
% Prepares EBSD data for import into JMAN

load(input.fullpath)

UserSelectHorizLine = true;
UserCropData = true;
DisplayFigure = true;

fprintf(1,'Original Data Size:%4.0f x%4.0f (%4.0f data points)\n',...
    size(Maps.X,1),size(Maps.X,2),size(Maps.X,1)*size(Maps.X,2));

if input.isXEBSD
    % Displacement Gradient Tensor where A(i,j) = du_i/dx_j
    % A = E (strain tensor) + W (rotation tensor)
    Maps.A{1,1} = Maps.E11_F;
    Maps.A{1,2} = Maps.W12_F1+Maps.E12_F;
    Maps.A{1,3} = Maps.W13_F1+Maps.E13_F;
    Maps.A{2,1} = Maps.W21_F1+Maps.E12_F;
    Maps.A{2,2} = Maps.E22_F;
    Maps.A{2,3} = Maps.W23_F1+Maps.E23_F;
    Maps.A{3,1} = Maps.W31_F1+Maps.E13_F;
    Maps.A{3,2} = Maps.W32_F1+Maps.E23_F;
    Maps.A{3,3} = Maps.E33_F;
    
    if strcmp(input.StressUnits,'MPa')
        Maps.S = {Maps.S11_F Maps.S12_F Maps.S13_F;...
            Maps.S12_F Maps.S22_F Maps.S23_F;...
            Maps.S13_F Maps.S23_F Maps.S33_F};
    elseif strcmp(input.StressUnits,'GPa') % If GPa need to convert to MPa
        sf = 1000;
        Maps.S = {Maps.S11_F.*sf Maps.S12_F.*sf Maps.S13_F.*sf;...
            Maps.S12_F.*sf Maps.S22_F.*sf Maps.S23_F.*sf;...
            Maps.S13_F.*sf Maps.S23_F.*sf Maps.S33_F.*sf};
    end
    
    Maps.E = {Maps.E11_F Maps.E12_F Maps.E13_F;...
        Maps.E12_F Maps.E22_F Maps.E23_F;...
        Maps.E13_F Maps.E23_F Maps.E33_F};
    
    % Maps.W = 1/2{SijEij}
    Maps.W = (1/2).*(Maps.S{1,1}.*Maps.E{1,1} + Maps.S{2,1}.*Maps.E{2,1} ...
        + Maps.S{3,1}.*Maps.E{3,1} +...
        Maps.S{2,1}.*Maps.E{2,1} + Maps.S{2,2}.*Maps.E{2,2} + ...
        Maps.S{2,3}.*Maps.E{2,3} +...
        Maps.S{3,1}.*Maps.E{3,1} + Maps.S{3,2}.*Maps.E{3,2} +...
        Maps.S{3,3}.*Maps.E{3,3});
    
end
%% Find Stepsize
% Scale the X and Y axes to get lengths in [um]
input.stepsize = Data.Stepsize; % um
stepsize = input.stepsize * 0.001; % convert um into mm

V.unit = 'mm'; %units for plotting
Maps.X = Maps.X .* stepsize;
Maps.Y = Maps.Y .* stepsize;

%% DEFINE HORIZONTAL LINE

% Which component to visualise for user input? A{Vi,Vj}
V.i = 1; V.j = 2; V.type = 'A';

if UserSelectHorizLine
    % Take user input to define horizontal line
    if V.type == 'A'
        [ptX,ptY]=selectHorizLine(Maps.X,Maps.Y,Maps.A{V.i,V.j},V);
    elseif V.type == 'S'
        [ptX,ptY]=selectHorizLine(Maps.X,Maps.Y,Maps.S{V.i,V.j},V);
    else
        error('Incorrect visualisation type');
    end
else
    % OR use default values
    ptX = [24.5378;88.0487];
    ptY = [29.4905;58.7832];
end

theta = atan((ptY(2)-ptY(1))/(ptX(2)-ptX(1)));
if input.rotate180
    theta = theta + pi; % rotates by 180 degrees
end
deg = theta*180/pi;

% X,Y are the x and y co-ordinates of the original data points, defined
% with respect to the original frame of reference
X = Maps.X(:);
Y = Maps.Y(:);

a = DirectionCosine(theta);

%% Rotate the locations of the data points

% Define a new set of points, with respect to the original axes, but whose
% locations are rotated by an angle theta about the origin.
Xnew = a(1,1).*X + a(1,2)*Y;
Ynew = a(2,1).*X + a(2,2)*Y;
Xnew = double(Xnew); Ynew = double(Ynew);
ptXnew = a(1,1).*ptX + a(1,2)*ptY;
ptYnew = a(2,1).*ptX + a(2,2)*ptY;


% create a coarse grid of points to visualise data
gridCoarseness = stepsize;
i = (min(Xnew):gridCoarseness:max(Xnew));
j = (min(Ynew):gridCoarseness:max(Ynew));
% xq and yq are the co-ordinates of a grid parellel to the crack,
% defined with respect to the original frame of reference.
[xq,yq] = meshgrid(i,j);
xq = double(xq); yq = double(yq);
% interpolate data onto this new grid of xq,yq to allow a contour plot to
% be made
[Maps] = interpolateData(Xnew,Ynew,Maps,xq,yq);

%% Perform tensor transformation to new, rotated, co-ordinate system

% Displacement Gradient Tensor:
[Maps.txf.A] = componentTensorTransform(Maps.rot.A,a);
% Stress Tensor:
[Maps.txf.S] = componentTensorTransform(Maps.rot.S,a);
% Strain Tensor:
[Maps.txf.E] = componentTensorTransform(Maps.rot.E,a);

%% CROP THE DATA

% Which component to visualise for user input? A{Vi,Vj}
V.i = 2; V.j = 2; V.type = 'S';

if UserCropData
    % Get User Input to crop the required data
    if V.type == 'A'
        [selXcrop,selYcrop,h2] = selectCrop(xq,yq,Maps.txf.A{V.i,V.j},...
            ptXnew,ptYnew,V);
    elseif V.type == 'S'
        [selXcrop,selYcrop,h2] = selectCrop(xq,yq,Maps.txf.S{V.i,V.j},...
            ptXnew,ptYnew,V);
    else
        error('Incorrect visualisation type');
    end
else
    % or use default crop values...
    selXcrop = [99,101];
    selYcrop = [4,5.5];
end

% Generate crop Mask
cropMask = xq;
cropMask(xq<min(selXcrop))=0;
cropMask(xq>max(selXcrop))=0;
cropMask(yq<min(selYcrop))=0;
cropMask(yq>max(selYcrop))=0;
cropMask(cropMask~=0)=1;

[Ycrop, Xcrop] = find(cropMask~=0);
Xcrop = [min(Xcrop) max(Xcrop)];
Ycrop = [min(Ycrop) max(Ycrop)];

% Crop the data
% Crop the x-coordinates
Maps.crop.X = xq(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop));
% Crop the y-coordinates
Maps.crop.Y = yq(min(Ycrop):max(Ycrop),min(Xcrop):max(Xcrop)); 
for i = 1:3
    for j = 1:3
        [Maps.crop.A{i,j}] = Crop(Maps.txf.A{i,j},Xcrop,Ycrop);
        [Maps.crop.S{i,j}] = Crop(Maps.txf.S{i,j},Xcrop,Ycrop);
        [Maps.crop.E{i,j}] = Crop(Maps.txf.E{i,j},Xcrop,Ycrop);
    end
end
Maps.crop.W = (1/2).*(Maps.crop.S{1,1}.*Maps.crop.E{1,1} +...
    Maps.crop.S{2,1}.*Maps.crop.E{2,1} + Maps.crop.S{3,1}.*Maps.crop.E{3,1} +...
    Maps.crop.S{2,1}.*Maps.crop.E{2,1} + Maps.crop.S{2,2}.*Maps.crop.E{2,2} +...
    Maps.crop.S{2,3}.*Maps.crop.E{2,3} +...
    Maps.crop.S{3,1}.*Maps.crop.E{3,1} + Maps.crop.S{3,2}.*Maps.crop.E{3,2} +...
    Maps.crop.S{3,3}.*Maps.crop.E{3,3});

% Shift the X and Y axes to give values starting at zero
Maps.crop.Xsc = Maps.crop.X - min(Maps.crop.X(:));
Maps.crop.Ysc = Maps.crop.Y - min(Maps.crop.Y(:));

fprintf(1,'Cropped Data Size:%4.0f x%4.0f (%4.0f data points)',...
    size(Maps.crop.X,1),size(Maps.crop.X,2),size(Maps.crop.X,1)*...
    size(Maps.crop.X,2));

%% PLOT FIGURE SHOWING CROPPED DATA
if DisplayFigure
    Xvec = Maps.crop.Xsc(1,:);
    Yvec = Maps.crop.Ysc(:,1);
    figure
    if V.type == 'A'
        imagesc(Xvec,Yvec,Maps.crop.A{V.i,V.j})
    elseif V.type == 'S'
        imagesc(Xvec,Yvec,Maps.crop.S{V.i,V.j})
    else
        error('Incorrect visualisation type');
    end
    % hold on
    % plot([min(Maps.crop.X(:));ptXnew(2)],ptYnew,'color','y') 
    % plots horizontal line
    % hold off
    set(gca,'YDir','normal')
    axis equal
    str = ['Cropped Data [Displaying the ',V.type,'_',num2str(V.i),...
        '_',num2str(V.j),'^n^e^w component, with respect to the new axes)'];
    title(str)
    xlim([min(Xvec) max(Xvec)])
    ylim([min(Yvec) max(Yvec)])
    colorbar
    xlabel(['x-position [',V.unit,']']);
    ylabel(['y-position [',V.unit,']']);
else
end

StressPlot(Maps)
input.Save = fullfile(input.results,[ num2str(input.fillname) '.' ...
    num2str(input.Trail_number) ' Strain Map.png']);
saveas(gcf,input.Save); close all
%% BUILD EBSDdata ARRAY
% Build the 'alldata' matrix containing all displacement gradient
% information.
% [x y A11 A12 A13 A21 A22 A23 A31 A32 A33 S11 S12 S13 S22 S23 S33]

EBSDdata = [Maps.crop.Xsc(:) Maps.crop.Ysc(:) ...
    Maps.crop.A{1,1}(:) Maps.crop.A{1,2}(:) Maps.crop.A{1,3}(:)...
    Maps.crop.A{2,1}(:) Maps.crop.A{2,2}(:) Maps.crop.A{2,3}(:)...
    Maps.crop.A{3,1}(:) Maps.crop.A{3,2}(:) Maps.crop.A{3,3}(:)...
    Maps.crop.S{1,1}(:) Maps.crop.S{1,2}(:) Maps.crop.S{1,3}(:)...
    Maps.crop.S{2,2}(:) Maps.crop.S{2,3}(:) Maps.crop.S{3,3}(:)];
EBSDdata(isnan(EBSDdata))=0; %Replace NaN values with 0

end

