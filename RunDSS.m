% addpath('P:\Abdo\GitHub\mtex-5.2.beta2');          startup_mtex;  
addpath('C:\Users\scro3511\OneDrive - Nexus365\Documents\GitHub\My-Mtex-Code\Data');
for xii = 1:9
    for xiX =1:3
%% Data
    clearvars -except xii xiX path Ok
    close all;          colordef white;              clc;        warning off;
    switch xii
        case 1;     B_34_191021_DSSB_100hr_v06;      % 19-10-21
        case 2;     B_34_191020_DSSB_100hr_v06;      % 19-10-20
        case 3;     B_34_191019_DSSB_100hr_v06;      % 19-10-19
        case 4;     B_34_191018_DSSB_100hr_v06;      % 19-10-18
        case 5;     B_34_191018_DSSB_100hr_v06;      % 19-10-14
        case 6;     B_34_191013_DSSB_100hr_v11;      % 19-10-13
        case 7;     B_34_191012_DSSB_100hr_v11;      % 19-10-12
        case 8;     B_34_191011_DSSB_100hr_v01;      % 19-10-11
        case 9;     AM190930DSSB_100hr;              % 19-09-30
        case 10;    AM190929DSSB_100hr;              % 19-09-29
        case 11;    AM190916DSSB_100hr;              % 19-09-16
        case 12;    AM190914DSSB_100hr;              % 19-09-14
        case 13;    AM190913DSSB_100hr;              % 19-09-13
        case 14;    AM190912DSSB_100hr;              % 19-09-12
        case 15;    AM190909DSSB_100hr;              % 19-09-09
        case 16;    AM190906DSSB_100hr;              % 19-09-06
        case 17;    AM190903DSSB_100hr;              % 19-09-03
        case 18;    AM190902DSSB_100hr;              % 19-09-02
        case 19;    AM190408DSSB_100hr;              % 19-04-08
        case 20;    AM190401DSSB_100hr;              % 19-04-01
            %AM190728DSSB_100hr
    end
    [~,fname,~]           = fileparts(fname);
    fprintf('%s with crack no. %d\n\n',fname,xiX);
    Ok.file{xiX,xii}      = fname;  

%% J-intergal
    main_JMAN     
    Ok.Results{xiX,xii}   = Results;
    Ok.Ktrue(xiX,xii)     = Ktrue;
    Ok.Kdiv(xiX,xii)      = Kdiv;
    Ok.nu(xiX,xii)        = Dir.nu;
    Ok.E(xiX,xii)         = Dir.E*1e-6;
    Ok.y(xiX,xii)         = Dir.yield*1e-9;
    Tg = table( Results(:,1),Results(:,2),Results(:,3),Results(:,4),...
                'VariableNames',{'J', 'K_Pa','J_corr', 'K_Pa_corr'}); % 
    Excel_path = 'P:\Abdo\EBSD Data\Alex xEBSD\D-SS.xlsx'; 
    sheet2 = [num2str(xii) '.' num2str(xiX)];
    writetable(Tg,Excel_path,'Sheet',sheet2);
    Tg = table(Dir.Stiffness.*1e-6,'VariableNames',{'Stiffness_tensors_GPa'}); 
    writetable(Tg,Excel_path,'Sheet',sheet2,'Range','F2:K9')
    Ok.Crack(xiX,xii)     = input('What was that, Twin (T), Crack (C) or Slip (S)?  ', 's');
    end
    %% Traces
    addpath('C:\Users\scro3511\OneDrive - Nexus365\Documents\GitHub\Trace Analysis');
    DataLoc = 'C:\Users\scro3511\OneDrive - Nexus365\Documents\GitHub\My-Mtex-Code\Data';
    [Traces{xii}] = SlipTrace(fullfile(DataLoc,Sd));
    writetable(Traces{xii},Excel_path,'Sheet',sheet2,'Range','F11:N18')
end
Tg = table(Ok.file(:),Ok.Crack(:),Ok.nu(:),Ok.E(:),Ok.y(:),Ok.Ktrue(:),Ok.Kdiv(:),...
    'VariableNames',{'file', 'Crack','Possion_Ratio','Young_Modulus_GPa',...
    'Yield_Stress_GPa' 'K_MPa', 'Div_MPa'}); %
writetable(Tg, Excel_path,'Sheet','All');