function [Values] = loadingXEBSD(DirxEBSD)
load(DirxEBSD,'Maps','GND','iPut','Data','C_rotated');      	

if      exist('iPut') == 1 ||     exist('C_rotated') == 1       
    Values.GND     = GND.total;   % from xEBSD  
    % rotation
    Values.W11 = Maps.W11_F1;   Values.W12 = Maps.W12_F1;	Values.W13 = Maps.W13_F1;           
    Values.W21 = Maps.W21_F1;   Values.W22 = Maps.W22_F1;	Values.W23 = Maps.W23_F1;    
    Values.W31 = Maps.W31_F1;  	Values.W32 = Maps.W32_F1;   Values.W33 = Maps.W33_F1;  

    if     exist('C_rotated') == 0  % old version of xEBSD
        % stress
        Values.S11 = Maps.S11_F;	Values.S12 = Maps.S12_F;    Values.S13 = Maps.S13_F;         	
        Values.S22 = Maps.S22_F;    Values.S23 = Maps.S23_F;   	Values.S33 = Maps.S33_F;  
        %strain
        Values.E11 = Maps.S11_F; 	Values.E12 = Maps.S12_F;    Values.E13 = Maps.S13_F;          	
        Values.E22 = Maps.S22_F;    Values.E23 = Maps.S23_F;  	Values.E33 = Maps.S33_F;
        %
        Values.Stiffness  = iPut.stiffnessvalues;
        Values.X   = Data.X;        Values.Y   = Data.Y;
    
    elseif exist('C_rotated') == 1 %NEW xEBSD version
        load(DirxEBSD,'Map_stress_sample','Map_strain_sample');
        % stress
        Values.S11 = Map_stress_sample(:,:,1,1);    Values.S12 = Map_stress_sample(:,:,1,2);
        Values.S13 = Map_stress_sample(:,:,1,3);    Values.S22 = Map_stress_sample(:,:,2,2);
        Values.S23 = Map_stress_sample(:,:,2,3);    Values.S33 = Map_stress_sample(:,:,3,3);
        % strain
        Values.E11 = Map_strain_sample(:,:,1,1);    Values.E12 = Map_strain_sample(:,:,1,2);
        Values.E13 = Map_strain_sample(:,:,1,3);    Values.E22 = Map_strain_sample(:,:,2,2);
        Values.E23 = Map_strain_sample(:,:,2,3);    Values.E33 = Map_strain_sample(:,:,3,3); 
        %
        Values.Stiffness = C_rotated;   
        Values.X   = Data.XSample;               	Values.Y   = Data.YSample;
    end   
    % stepsize
    Values.stepsize =(abs(Values.X(1,1)-Values.X(1,2)));
    
elseif exist('iPut') ~= 1 && exist('C_rotated') ~= 1  && exist('Maps') ~= 1
    Values = Maps;
else %if using my code
    [filepath,~,~] = fileparts(DirxEBSD);
    [~,name,~]     = fileparts(filepath);
    load([filepath '\' name '.mat'],'Maps');
    Values = Maps;
end
