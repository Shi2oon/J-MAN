function [files,Operation,unit] = Westergaard(StressIntensityFactor)
% close all; clear; clc

state = [{'plane_stress'},{'plane_strain'},{'plane_stress'}];

minGrid = -4E-3; % [m]
gridStep = 0.1E-3; % [m]
maxGrid = 4E-3; % [m]

xvec = minGrid : gridStep : maxGrid; % [m]
yvec = minGrid : gridStep : maxGrid; % [m]

[x,y] = meshgrid(xvec,yvec); % [m]
% StressIntensityFactor = 30; % [MPa m^0.5]
fprintf('preparing synthetic Westergaard Solution Data .. ');
K_I = StressIntensityFactor * 1E6; % Stress intensity factor [Pa m^0.5]
E = 210E9; % Young's Modulus [Pa]
nu = 0.3; % poisson ratio
mu = E/(2.*(1+nu)); % Shear Modulus [Pa]
for iv =1:3
    switch state{iv}
        case 'plane_strain'
            kappa = 3 - (4 .* nu); % [/]
        case 'plane_stress'
            kappa = (3 - nu)./(1 + nu); % [/]
    end
[theta,r] = cart2pol(x,y);
ux = (K_I./(2.*mu)).*sqrt(r/(2*pi)).*cos(theta/2).*(kappa - 1 + 2.*(sin(theta/2)).^2); % Anderson p99 A2.44a
uy = (K_I./(2.*mu)).*sqrt(r/(2*pi)).*sin(theta/2).*(kappa + 1 - 2.*(cos(theta/2)).^2); % Anderson p99 A2.44b
[dux_dx,dux_dy] = gradient(ux,gridStep);
[duy_dx,duy_dy] = gradient(uy,gridStep);
exx = dux_dx;
eyy = duy_dy;
exy = 0.5*(dux_dy + duy_dx);
if iv == 1
    %% save disps data
    Operation{1} = 'Dis';
%     alldata = [x(:) y(:) ux(:) uy(:)]; % [m]
%     save([num2str(StressIntensityFactor) 'MPa_[m]_DISP.mat'],'alldata')
    alldata = [x(:) y(:) ux(:) uy(:)].*1000; % [mm]
    save([num2str(StressIntensityFactor) 'MPa_[mm]_DISP.mat'],'alldata')
    files{1} = [num2str(StressIntensityFactor) 'MPa_[mm]_DISP'];
    unit{1}  = 'mm';
elseif iv ==2
    %% save strain data
    Operation{2} = 'Str';
    % alldata = [x(:) y(:) exx(:) eyy(:) exy(:)]; % [m]
    % save([num2str(StressIntensityFactor) 'MPa_[m]_Strain.mat'],'alldata')
    alldata = [x(:).*1000 y(:).*1000 exx(:) eyy(:) exy(:)]; % [mm]
    save([num2str(StressIntensityFactor) 'MPa_[mm]_Strain.mat'],'alldata')
    files{2} = [num2str(StressIntensityFactor) 'MPa_[mm]_Strain'];
    unit{2} = 'mm';    
elseif iv == 3
%% save EBSD data 
    Operation{3} = 'xED';
    stiff(1,1) = E*(1-nu)/((1-2*nu)*(1+nu))*1e-9; 
    stiff(1,2) = nu*stiff(1,1)/(1-nu);
    stiff(1,3) = stiff(1,1)-stiff(1,2);
    stiff = [stiff(1,1), stiff(1,2), stiff(1,2), 0,          0,          0; 
             stiff(1,2), stiff(1,1), stiff(1,2), 0,          0,          0; 
             stiff(1,2), stiff(1,2), stiff(1,1), 0,          0,          0; 
             0,          0,          0,          stiff(1,3), 0,          0; 
             0,          0,          0,          0,          stiff(1,3), 0; 
             0,          0,          0,          0,          0,          stiff(1,3)];
    Maps.GND = (10^15.5-10^14).*rand(length(eyy),length(eyy),1) + 10^14;   % from xEBSD  
    % rotation
    Maps.W11 = zeros(size(exx));      Maps.W12 = exx;                  Maps.W13 = eyy;           
    Maps.W21 = -exx;                  Maps.W22 = zeros(size(exx));     Maps.W23 = exy;    
    Maps.W31 = -eyy;                  Maps.W32 = -exy;                 Maps.W33 = zeros(size(exx)); 
    %strain
    Maps.E11 = exx;                   Maps.E12 =exy;                   Maps.E13 = exx-exy;          	
    Maps.E22 = eyy;                   Maps.E23 = eyy-exy;              Maps.E33 = zeros(size(exx));
    %stress in GPa
    Maps.S11 = Maps.E11*stiff(1,1);	  Maps.S12 = Maps.E12*stiff(1,2);  Maps.S13 = Maps.E13*stiff(1,3);         	
    Maps.S22 = Maps.E22*stiff(2,2);   Maps.S23 = Maps.E23*stiff(2,3);  Maps.S33 = Maps.E33*stiff(3,3);
    %% stepsize
    Maps.Stiffness  = stiff;
    Maps.X   = x.*1e6;                 Maps.Y   = y.*1e6;
    Maps.stepsize =(abs(Maps.X(1,1)-Maps.X(1,2)));
    save([num2str(StressIntensityFactor) 'MPa_[um]_xEBSD.mat'],'Maps')
    files{3} = [num2str(StressIntensityFactor) 'MPa_[um]_xEBSD'];
    unit{3} = 'um';    
end
end
fprintf ('DONE\n\n');
%% Plot
% % figure
% % quiver(x,y,ux,uy)
% 
% figure
% % subplot(2,4,1)
% % imagesc(xvec,yvec,x)
% % colorbar
% % colormap(jet)
% % title('x')
% % xlabel('x')
% % ylabel('y')
% % axis('image')
% % set(gca,'YDir','normal')
% % 
% % subplot(2,4,2)
% % imagesc(xvec,yvec,y)
% % colorbar
% % title('y')
% % xlabel('x')
% % ylabel('y')
% % axis('image')
% % set(gca,'YDir','normal')
% 
% subplot(2,4,1)
% imagesc(xvec,yvec,r)
% colorbar
% title('r')
% xlabel('x')
% ylabel('y')
% axis('image')
% set(gca,'YDir','normal')
% 
% subplot(2,4,2)
% imagesc(xvec,yvec,theta)
% colorbar
% title('theta')
% xlabel('x')
% ylabel('y')
% axis('image')
% set(gca,'YDir','normal')
% 
% subplot(2,4,3)
% C = ux;
% imagesc(xvec,yvec,C)
% colorbar
% title('ux')
% xlabel('x')
% ylabel('y')
% axis('image')
% set(gca,'YDir','normal')
% 
% subplot(2,4,4)
% C = uy;
% imagesc(xvec,yvec,C)
% colorbar
% title('uy')
% xlabel('x')
% ylabel('y')
% axis('image')
% set(gca,'YDir','normal')
% 
% subplot(2,4,5)
% C = exx;
% clim = [0 0.01];
% imagesc(xvec,yvec,C,clim)
% colorbar
% title('exx')
% xlabel('x')
% ylabel('y')
% axis('image')
% set(gca,'YDir','normal')
% 
% subplot(2,4,6)
% C = eyy;
% clim = [0 0.01];
% imagesc(xvec,yvec,C,clim)
% colorbar
% title('eyy')
% xlabel('x')
% ylabel('y')
% axis('image')
% set(gca,'YDir','normal')
% 
% subplot(2,4,7)
% C = exy;
% clim = [-0.01 0.01];
% imagesc(xvec,yvec,C,clim)
% colorbar
% title('exy')
% xlabel('x')
% ylabel('y')
% axis('image')
% set(gca,'YDir','normal')
% 
% 
% % figure
% % mesh(x,y,uy)
% % zlabel('uy')
% % xlabel('x')
% % ylabel('y')