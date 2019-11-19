checks.SIF = [0.1, 0.5, 1, 5, 10, 50, 100];        
addpath(genpath('C:\Users\scro3511\OneDrive - Nexus365\Documents\GitHub\J-MAN')); 
for ii =3:length(checks.SIF)
    [checks.files,checks.Operation,checks.unit] = Westergaard(checks.SIF(ii));
    for i=1:3
        clc; close all; clearvars -except i checks ii
        Dir.fullpath        = pwd;
        Dir.fillname        = checks.files{i};
        Operation           = checks.Operation{i};
        Dir.input_unit      = checks.unit{i};
        
        main_JMAN
        
        contrs = length(Results);           contrs = contrs - round(contrs/5);
        checks.Ktrue(ii,i)  = 1E-6.*((mean(Results(contrs:end,2))+max(Results(contrs:end,2)))/2);
        checks.KCorr(ii,i)  = Ktrue;
        checks.error(ii,i)  = ER;
    end
end

close all
CF.Dis = 1-mean(checks.Ktrue(:,1)./checks.SIF); %DIC correction factor
CF.Str = 1-mean(checks.Ktrue(:,2)./checks.SIF); %DIC correction factor
CF.xED = 1-mean(checks.Ktrue(:,3)./checks.SIF); %DIC correction factor
% save('Holy_Correction_Factors.mat','CF'); %if any 
