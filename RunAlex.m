addpath('P:\Abdo\GitHub\mtex-5.2.beta2');          startup_mtex;  
addpath('C:\Users\scro3511\OneDrive - Nexus365\Documents\GitHub\My-Mtex-Code\Data');

%%
for xii = 9:9
    for xiX =1:3
    clearvars -except xii xiX path Ok
    close all;          colordef white;         clc;            warning off;
 
    [pname,fname,ebsd,CS] = Alex_Data(xii);
    [~,fname,~]           = fileparts(fname);
    fprintf('%s with crack no. %d\n\n',fname,xiX);
    Ok.file{xiX,xii}      = fname;
    Ok.Crack(xiX,xii)     = xiX;

%% J-intergal
    main_JMAN 
    
    Ok.Results{xiX,xii}   = Results;
    Ok.Ktrue(xiX,xii)     = Ktrue;
    Ok.Kdiv(xiX,xii)      = Kdiv;
    Ok.nu(xiX,xii)        = Dir.nu;
    Ok.E(xiX,xii)         = Dir.E*1e-9;
    Ok.y(xiX,xii)         = Dir.yield*1e-9;

    Tg = table( Results(:,1),Results(:,2),Results(:,3),Results(:,4),...
                'VariableNames',{'J', 'K_Pa','J_corr', 'K_Pa_corr'}); % 
    Excel_path = 'P:\Abdo\EBSD Data\Alex xEBSD\AlexData_v5.xlsx'; 
    sheet2 = [num2str(xii) '.' num2str(xiX)];
    writetable(Tg,Excel_path,'Sheet',sheet2);
    
    Tg = table(Dir.Stiffness.*1e-6,'VariableNames',{'Stiffness_tensors_GPa'}); 
    writetable(Tg,Excel_path,'Sheet',sheet2,'Range','F2:K9')
    end
%% Traces   
    addpath('C:\Users\scro3511\OneDrive - Nexus365\Documents\GitHub\Trace Analysis');
    DataLoc = 'C:\Users\scro3511\OneDrive - Nexus365\Documents\GitHub\My-Mtex-Code\Data';
    if      xii==1;             Sd = 'Alex_6H_750_irr_1000nm';
    elseif  xii==2;             Sd = 'Alex_6H_Ne_1000nm';
    elseif  xii==3;             Sd = 'Alex_6H_Si_300_high_unirr_1000nm';
    elseif  xii==4;             Sd = 'Alex_Si_indent_21_EBSD';
    elseif  xii==5;             Sd = 'Alex_Ne_1000nm_indent_25';
    elseif  xii==6;             Sd = 'Alex_Si_1000nm_unirr_32';
    elseif  xii==7;             Sd = 'Alex_unirradiatedSi300high';
    elseif  xii==8;             Sd = 'Alex_unirr_6H_1000nm';
    elseif  xii==9;             Sd = 'Alex_indent_25';
    end
    [Traces{xii}] = SlipTrace(fullfile(DataLoc,Sd));
    writetable(Traces{xii},Excel_path,'Sheet',sheet2,'Range','F11:N18')
end

Tg = table(Ok.file(:),Ok.Crack(:),Ok.nu(:),Ok.E(:)*1000,Ok.y(:),Ok.Ktrue(:),Ok.Kdiv(:),...
    'VariableNames',{'file', 'Crack','Possion_Ratio','Young_Modulus_GPa',...
    'Yield_Stress_GPa' 'K_MPa', 'Div_MPa'}); %
writetable(Tg, Excel_path,'Sheet','All');