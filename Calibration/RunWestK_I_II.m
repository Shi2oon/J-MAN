files{1,1}  = '30MPa_[mm]_DISP_I_N';          files{1,2}  = 'Dis';
files{2,1}  = '30MPa_[mm]_DISP_I_NZ';         files{2,2}  = 'Dis';
files{3,1}  = '30MPa_[mm]_DISP_I_Z';          files{3,2}  = 'Dis';

files{4,1}  = '30MPa_[mm]_DISP_II_N';         files{4,2}  = 'Dis';
files{5,1}  = '30MPa_[mm]_DISP_II_NZ';        files{5,2}  = 'Dis';
files{6,1}  = '30MPa_[mm]_DISP_II_Z';         files{6,2}  = 'Dis';

files{7,1}  = '30MPa_[mm]_Strain_I_N';        files{7,2}  = 'Str';
files{8,1}  = '30MPa_[mm]_Strain_I_NZ';       files{8,2}  = 'Str';
files{9,1}  = '30MPa_[mm]_Strain_I_Z';        files{9,2}  = 'Str';

files{10,1} = '30MPa_[mm]_Strain_II_N';       files{10,2} = 'Str';
files{11,1} = '30MPa_[mm]_Strain_II_NZ';      files{11,2} = 'Str';
files{12,1} = '30MPa_[mm]_Strain_II_Z';       files{12,2} = 'Str';

for ig = 1:length(files)
    clearvars -except files ig KK_true
    main_JMAN       %fix values for crack locaiton (mask values)
    if ig < 7
        KK_true(1,ig)   = round(Ktrue,1);
    else
        KK_true(2,ig-6) = round(Ktrue,1);
    end
end
k = [1 1 1 2 2 2];
table(k(:),KK_true(1,1:6), KK_true(1,7:12), 'VariableNames',...
    {'Mode','Displacement','Strain'})
    