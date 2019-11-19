function [Ktrue2,Operation]=Westegraard_JMAN(mat,alldata,maxGrid,Operation)
close all; clc
[mesh] = loadUdCopy(alldata,Operation,mat);                %loads data
mesh.maxGrid=maxGrid;
[el2,mesh] = meshDIC(mesh);
%runs isoparematric FE anlysis to slove for strain and stress at Gauss points
[el2] = FEanalysisStrains(el2,mat,mesh);
%Interpolates FE results into a global domain at nodal co-ordinates
[gl2] = makeglobal(el2, mesh); 
[Jint2]=getMask(mesh,gl2,el2,mat); 
close all % masking
[el2,mesh,gl2,Jint2]=StrainInegration(mat,mesh,el2,Jint2,gl2);% STRAIN INTEGRATION
close all;  
[Results2,Jint2] =J_IntegralCacl(Jint2,mat,mesh,el2); % J-INTEGRAL CALCULATION 

contrs = length(Results2);           contrs = contrs - round(contrs/5);
Ktrue2 = 1E-6.*((mean(Results2(contrs:end,2))+max(Results2(contrs:end,2)))/2);
% 
% %% checks
% Logix (1) = max(max(gl.Ux-gl2.Ux))          == 0;
% Logix (2) = max(max(gl.Uy-gl2.Uy))          == 0;
% Logix (3) = max(max(gl.dx-gl2.dx))          == 0;
% Logix (4) = max(max(gl.dy-gl2.dy))          == 0;
% Logix (5) = max(max(gl.Uxdef-gl2.Uxdef))    == 0;
% Logix (6) = max(max(gl.Uydef-gl2.Uydef))    == 0;
% Logix (7) = max(max(gl.smises-gl2.smises))  == 0;
% clc;
% if max(Logix) == 1
%     fprintf ('\nnothing is wrong :)\n');
% else
%     fprintf('\nyou may need to check each function one bey one :(\n');
% end

