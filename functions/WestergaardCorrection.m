function [Results,error] =WestergaardCorrection (Results,mat,operation,maxGrid,elements)
Ktrue=1E-6.*((mean(Results(:,2))+max(Results(:,2)))/2);

%if operation == 'Dis'
xvec = linspace(-maxGrid,maxGrid,elements); % [m]
yvec = linspace(-maxGrid,maxGrid,elements); % [m]
gridStep = abs(yvec(1)-yvec(2)); % [m]

[x,y] = meshgrid(xvec,yvec); % [m]
K_I = Ktrue * 1E6; % Stress intensity factor [Pa m^0.5]
mu = mat.E/(2.*(1+mat.nu)); % Shear Modulus [Pa]
state = mat.stressstate;
switch state
    case 'plane_strain'
        kappa = 3 - (4 .* mat.nu); % [/]
    case 'plane_stress'
        kappa = (3 - mat.nu)./(1 + mat.nu); % [/]
end
[theta,r] = cart2pol(x,y);
ux = (K_I./(2.*mu)).*sqrt(r/(2*pi)).*cos(theta/2).*(kappa - 1 +...
    2.*(sin(theta/2)).^2); % Anderson p99 A2.44a
uy = (K_I./(2.*mu)).*sqrt(r/(2*pi)).*sin(theta/2).*(kappa + 1 - ...
    2.*(cos(theta/2)).^2); % Anderson p99 A2.44b
[dux_dx,dux_dy] = gradient(ux,gridStep);
[duy_dx,duy_dy] = gradient(uy,gridStep);
exx = dux_dx;
eyy = duy_dy;
exy = 0.5*(dux_dy + duy_dx);

% %save disps data
if operation == 'Dis'
    alldata = [x(:) y(:) ux(:) uy(:)]; % [m]   
elseif operation=='Str'
    alldata = [x(:) y(:) exx(:) eyy(:) exy(:)]; % [m]    
else
    fprintf('You need to specify the calculation mode\nDis (for Displacement) or Str (for Strain)')
end

[KJ_MAN] = Westegraard_JMAN(mat,alldata,maxGrid,operation);

clc;    error=Ktrue/KJ_MAN;
        Results(:,3:4)=Results*error;