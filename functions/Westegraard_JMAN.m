function [Ktrue,Operation]=Westegraard_JMAN(mat,alldata,maxGrid,Operation)
close all; clc
[mesh] = loadUdCopy(alldata,Operation,mat);                %loads data
mesh.maxGrid=maxGrid;
[el,mesh] = meshDIC(mesh);
%runs isoparematric FE anlysis to slove for strain and stress at Gauss points
[el] = FEanalysisStrains(el,mat,mesh);
%Interpolates FE results into a global domain at nodal co-ordinates
[gl] = makeglobal(el, mesh); 
[Jint]=getMask(mesh,gl,el,mat); 
close all % masking
[el,mesh,gl,Jint]=StrainInegration(mat,mesh,el,Jint,gl);% STRAIN INTEGRATION
close all;  
[Results,Jint] =J_IntegralCacl(Jint,mat,mesh,el); % J-INTEGRAL CALCULATION 

Ktrue=1E-6.*(mean(Results(:,2))+max(Results(:,2)))/2; 