function [mesh] =  loadUd(pixel_size,scan_ID, scan_dir, input_unit,Operation)
load(fullfile(scan_dir,[num2str(scan_ID) '.mat']));
mesh.Operation=Operation;
mesh.operation='Norm';
xx=unique(alldata(:,1),'first'); 
mesh.x.length=length(xx);
mesh.x.xi=(max(xx))-min(xx)/2;

if mesh.Operation=='Dis'
% Convert [x y] to metres
switch input_unit
    case 'mm'
        % convert mm --> m
        alldata = alldata.* (pixel_size * 10^-3); 
    case 'um'
        % convert um --> m
        alldata = alldata.* (pixel_size * 10^-6); % convert pixels --> m
    case 'm'
        alldata = alldata.* (pixel_size); 
    otherwise
        error('Invalid input_unit')
end

dispsU_1 = alldata(:,1:2);
dispsd_1 = alldata(:,3:4);
tmp(1:2,:)=dispsU_1';
tmp(3:4,:)=dispsd_1';
%Guess the dimensions of the image
mesh.winDIC = GetSize(tmp(2,:));
% sort the dataset
tmp = sortrows(tmp',[1 2])';
mesh.UDIC =tmp(1:2,:);
mesh.dDIC =tmp(3:4,:);
mesh.winFE  = [2 2; mesh.winDIC(1)-1 mesh.winDIC(2)-1];
% mesh.winFE  = [1 1; mesh.winDIC(1) mesh.winDIC(2)];

%% Prepare strain field data
elseif mesh.Operation=='Str'
    % Convert [x y] to metres
switch input_unit
    case 'mm'
        % convert mm --> m
        alldata(:,1:2) = alldata(:,1:2).* 10^-3; 
    case 'um'
        % convert um --> m
        alldata(:,1:2) = alldata(:,1:2).* 10^-6; % convert pixels --> m
    case 'm'
        %nothing
    otherwise
        error('Invalid input_unit')
end

alldata(isnan(alldata))=0;
dispsU_1 = alldata(:,1:2);

straindata = alldata(:,3:end);

%Interpolate strain data on a regular grid.
str_limits = [min(dispsU_1);max(dispsU_1)];
nbpts = [length(unique(dispsU_1(:,1))) length(unique(dispsU_1(:,2)))];
pas = (str_limits(2,:)-str_limits(1,:))./nbpts;

xrdXq=str_limits(1,1):pas(1):str_limits(2,1);
xrdYq=str_limits(1,2):pas(2):str_limits(2,2);
[Xq, Yq] = meshgrid(xrdXq,xrdYq);

Fxx = scatteredInterpolant(dispsU_1(:,1), dispsU_1(:,2), straindata(:,1));
Fyy = scatteredInterpolant(dispsU_1(:,1), dispsU_1(:,2), straindata(:,2));
Fxy = scatteredInterpolant(dispsU_1(:,1), dispsU_1(:,2), straindata(:,3));
intpExx = Fxx(Xq,Yq); intpEyy = Fyy(Xq,Yq); intpExy = Fxy(Xq,Yq);

straindata = [Xq(:) Yq(:) intpExx(:) intpEyy(:) -1.*intpExy(:)];

% Build the grid where strain points are gauss points
baseXq=[xrdXq(1)-(pas(1)/2) xrdXq+(pas(1)/2)];
baseYq=[xrdYq(1)-(pas(2)/2) xrdYq+(pas(2)/2)];
[Xbase, Ybase] = meshgrid(baseXq,baseYq);

dispsU_1 = [Xbase(:) Ybase(:)];
dispsd_1 = nan(size(dispsU_1));

%Guess the dimensions of the image
tmp(1:2,:)=dispsU_1'; tmp(3:4,:)=dispsd_1';

mesh.winDIC = GetSize(tmp(2,:));

%sort the dataset
tmp = sortrows(tmp',[1 2])';

mesh.UDIC = tmp(1:2,:);
mesh.dDIC = tmp(3:4,:);
mesh.straindata = straindata;

mesh.winFE  = [2 2; mesh.winDIC(1)-1 mesh.winDIC(2)-1];
end
