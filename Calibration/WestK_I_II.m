StressIntensityFactor = 30;
state = [{'plane_stress'},{'plane_strain'}];

%%
minGrid = -4E-3;    % [m]
gridStep = 0.1E-3;  % [m]
maxGrid = 4E-3;     % [m]

xvec  = minGrid : gridStep : maxGrid;   % [m]
yvec  = minGrid : gridStep : maxGrid;   % [m]
[x,y] = meshgrid(xvec,yvec);            % [m]
K_I = StressIntensityFactor * 1E6; % Stress intensity factor [Pa m^0.5]
E = 210E9; % Young's Modulus [Pa]
nu = 0.3; % poisson ratio
mu = E/(2.*(1+nu)); % Shear Modulus [Pa]
for iv=1:2
    switch state{iv}
        case 'plane_strain'
            kappa = 3 - (4 .* nu); % [/]
        case 'plane_stress'
            kappa = (3 - nu)./(1 + nu); % [/]
    end

    [theta,r] = cart2pol(x,y);
    ux = (K_I./(2.*mu)).*sqrt(r/(2*pi)).*cos(theta/2).*(kappa - 1 + 2.*(sin(theta/2)).^2); % Anderson p99 A2.44a
    uy = (K_I./(2.*mu)).*sqrt(r/(2*pi)).*sin(theta/2).*(kappa + 1 - 2.*(cos(theta/2)).^2); % Anderson p99 A2.44b
    for K =1:2
        if K ==1 && iv==1
            alldata = [x(:) y(:) ux(:) uy(:)].*1000; % [mm]
            save('30MPa_[mm]_DISP_I_N.mat','alldata')
            % near zero ux
            UX = ux/1e8;
            alldata = [x(:) y(:) UX(:) uy(:)].*1000; % [mm]
            save('30MPa_[mm]_DISP_I_NZ.mat','alldata')
            % zero ux
            UX =ux*0;
            alldata = [x(:) y(:) UX(:) uy(:)].*1000; % [mm]
            save('30MPa_[mm]_DISP_I_Z.mat','alldata')
            
        elseif K ==1 && iv==2
            [exx,eyy,exy] = Calc_Str(ux,uy,gridStep);
            alldata = [x(:).*1000 y(:).*1000 exx(:) eyy(:) exy(:)]; % [mm]
            save('30MPa_[mm]_Strain_I_N.mat','alldata')
            % near zero ux
            UX = ux/1e8;
            [exx,eyy,exy] = Calc_Str(UX,uy,gridStep);
            alldata = [x(:).*1000 y(:).*1000 exx(:) eyy(:) exy(:)]; % [mm]
            save('30MPa_[mm]_Strain_I_NZ.mat','alldata')
            % zero ux
            UX =ux*0;
            [exx,eyy,exy] = Calc_Str(UX,uy,gridStep);
            alldata = [x(:).*1000 y(:).*1000 exx(:) eyy(:) exy(:)]; % [mm]
            save('30MPa_[mm]_Strain_I_Z.mat','alldata')
            
        elseif K==2  && iv==1
            alldata = [x(:) y(:) uy(:) ux(:)].*1000; % [mm]
            save('30MPa_[mm]_DISP_II_N.mat','alldata')
            % near zero ux
            UX = ux/1e8;
            alldata = [x(:) y(:) uy(:) UX(:)].*1000; % [mm]
            save('30MPa_[mm]_DISP_II_NZ.mat','alldata')
            % zero ux
            UX =ux*0;
            alldata = [x(:) y(:) uy(:) UX(:)].*1000; % [mm]
            save('30MPa_[mm]_DISP_II_Z.mat','alldata')
            
        elseif K==2  && iv==2
            [exx,eyy,exy] = Calc_Str(uy,ux,gridStep);
            alldata = [x(:).*1000 y(:).*1000 exx(:) eyy(:) exy(:)]; % [mm]
            save('30MPa_[mm]_Strain_II_N.mat','alldata')
            % near zero ux
            UX = ux/1e8;
            [exx,eyy,exy] = Calc_Str(uy,UX,gridStep);
            alldata = [x(:).*1000 y(:).*1000 exx(:) eyy(:) exy(:)]; % [mm]
            save('30MPa_[mm]_Strain_II_NZ.mat','alldata')
            % zero ux
            UX =ux*0;
            [exx,eyy,exy] = Calc_Str(uy,UX,gridStep);
            alldata = [x(:).*1000 y(:).*1000 exx(:) eyy(:) exy(:)]; % [mm]
            save('30MPa_[mm]_Strain_II_Z.mat','alldata')
        end
    end
end
    
    function [exx,eyy,exy] = Calc_Str(ux,uy,gridStep)
        [dux_dx,dux_dy] = gradient(ux,gridStep);
        [duy_dx,duy_dy] = gradient(uy,gridStep);
        exx = dux_dx;
        eyy = duy_dy;
        exy = 0.5*(dux_dy + duy_dx);
    end