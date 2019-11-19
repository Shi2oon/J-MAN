%% FE code to evaluate J-integral from DIC, EDXD and HR-EBSD 
% DIC IMAGE:        <--------------mesh.winsize(1)---------->
%                  ...........................................
%             |    :  FE----------------------------------   :                   
%             |    :  |  Jo----------------------------  |   :
%             |    :  |  |  Ji----------------------  |  |   :
% mesh.winsize(2)  :========> crack in pos x dir.  |  |  |   :
%             |    :  |  |  ----------------------Ji  |  |   :
%             |    :  |  ----------------------------Jo  |   :
%             |    :  ----------------------------------FE   :                  
%                  :.........................................:
%VARIABLES:             % gl.   - Global interpolation of FE analysis from Gauss points
% mat.  - material and element properties       % Jint. - J-integral elements, variables and results
% el.   - FE elements, variables and results    % mesh. - FE mesh 
% tic; java.lang.Thread.sleep(Duration_in_Sec * 1000); toc; 
close all;       clc;clear;      colordef white;     restoredefaultpath; 
warning off;     	addpath(genpath(pwd));              printCredit;        
set(0,'defaultAxesFontSize',20);    set(0,'DefaultLineMarkerSize',12)   
 
%% Setting the scene
Dir.fullpath     = 'P:\Abdo\GitHub\xEBSD2ABAQUS\Calibration\12 WS'; % filr Directory
Dir.fillname     = '12MPa_[mm]_DISP';       % all in .mat w/ data in alldata var except DVC in .dat format
Dir.Trail_number = 1;           % file number %try to do at least 4
Dir.pixel_size   = 1;           % if DIC/DVC values are in pixel, 1 if in physical units;
Dir.input_unit   = 'm';        % meter (m) or milmeter (mm) or micrometer(um);
Operation        = 'Dis';       % calculation mode\n Dis = DIC, Str = Strain, xED = xEBSD, DVC
Dir.results      = fullfile(Dir.fullpath ,['J-MAN5_' Dir.fillname]);          mkdir(Dir.results);
% INPUT MATERIAL PROPERTIES AND DATA
% Poisson's ratio,          Young's Modulus [Pa],           Yield Stress [Pa] 
  Dir.nu    = 0.3;          Dir.E  = 210E9;                 Dir.yield = 4E9;
  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END of USER INTERFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%% Load, MESHING & STRAIN/STRESSES CALCULATIONS
% loadUD output is in meter
[mesh,Dir] = loadUd(Dir.pixel_size,Dir.fullpath, Dir.fillname, Dir.input_unit,Operation,Dir);
[mat]      = matprop(Dir.E,Dir.nu,Dir.yield,Dir.stressstat,Dir.input_unit); % material props.
[el,mesh]  = meshDIC(mesh);
%runs isoparematric FE anlysis to slove for strain and stress at Gauss points
[el]      = FEanalysisStrains(el,mat,mesh);
%Interpolates FE results into a global domain at nodal co-ordinates
[gl]      = makeglobal(el, mesh); 
% [mesh,mat,el,gl] = Crack_align(mesh,mat,el,gl);   %Crack Aligement

%% J_intergal Mask.
[Jint]    = getMask(mesh,gl,el,mat);                        % masking
saveas(gcf,[Dir.results '/' Dir.fillname '.' num2str(Dir.Trail_number) '_selection.png']); close all
[el,mesh,gl,Jint] = StrainInegration(mat,mesh,el,Jint,gl);	% STRAIN INTEGRATION, gl.dy & dx
if mesh.Operation == 'Str';         saveas(gcf,fullfile(Dir.results,[Dir.fillname '.' ...
   num2str(Dir.Trail_number) '_integrated_displacements.fig']));    close all;      end

%% J Calc
[Results,Jint]  = J_IntegralCacl(Jint,mat,mesh,el); % J-INTEGRAL CALCULATION
[Results,ER]    = WestergaardCorrection(Results,mat,mesh.Operation,mesh.x.xi,mesh.x.length);
[Jtrue,Jdiv,Ktrue,Kdiv] = PlotAll(Results,Dir,ER,mesh,el,gl,Jint); % Plotting 
save([Dir.results '\' Dir.fillname '.' num2str(Dir.Trail_number) 'JMAN.mat']);