%% FE code to evaluate J-integral from DIC images, by Thorsten Becker - b.thorsten@gmail.com
%% and Selim Barhli - selim.barhli@ensiacet.fr

%DIC IMAGE:         <--------------mesh.winsize(1)---------->
%                  ...........................................
%             |    :  FE----------------------------------   :                   
%             |    :  |  Jo----------------------------  |   :
%             |    :  |  |  Ji----------------------  |  |   :
% mesh.winsize(2)  :========> crack in pos x dir.  |  |  |   :
%             |    :  |  |  ----------------------Ji  |  |   :
%             |    :  |  ----------------------------Jo  |   :
%             |    :  ----------------------------------FE   :                  
%                  :.........................................:

% FE => mesh.winFE(x&y,point)-window size of image considered for the FE analysis
% Jo=>Jint.nout(x&y,point)-outer Jintegral path     % Ji=>Jint.nin(x&y,point)-inner Jintegral path
%VARIABLES:             % gl.   - Global interpolation of FE analysis from Gauss points
% mat.  - material and element properties    % Jint. - J-integral elements, variables and results
% el.   - FE elements, variables and results % mesh. - FE mesh 

close all;  clear;   colordef white;   clc;  
addpath([pwd '\functions']);       printCredit
 
%% Setting the scene
Dir.scan_dir  = 'X:\EE 13579 Diamond Light Source\Diamond EE 13579\Strain Analysis 6\66881 to 66887\J-MAN 66881_66887';
Dir.scan_ID   = 66881066887;   
Dir.Trail_number = 1;     %file number %try to do at least 4
pixel_size    = 3;      % 0.31; %0.54; %[um/px] from Davis :: Pixelsize=[value]um :: pixel_size_um=1/value;
input_unit    = 'mm';   % meter (m) or milmeter (mm) or micrometer(um) ;
Operation     = 'Str';    % calculation mode\n Dis (for Displacement) or Str (for Strain)'
Dir.results_dir=Dir.scan_dir;
addpath([pwd '\functions']);

%% INPUT MATERIAL PROPERTIES AND DATA
% Poisson's ratio,  Young's Modulus [Pa],   Yield Stress [Pa] 
mat.nu    = 0.3;       mat.E  = 210E9;         mat.yield = 820E6;
mat.stressstate = 'plane_strain';    %plane_strain (bulk) or plane_stress (surface)
[mat]     = matprop(mat);               %sets material peoperties: E,nu,etc...

%% Load, MESHING & STRAIN/STRESSES CALCULATIONS
[mesh]    = loadUd(pixel_size,Dir.scan_ID, Dir.scan_dir, input_unit,Operation) ;  %loads DIC data
[el,mesh] = meshDIC(mesh);
%runs isoparematric FE anlysis to slove for strain and stress at Gauss points
[el]      = FEanalysisStrains(el,mat,mesh);
%Interpolates FE results into a global domain at nodal co-ordinates
[gl]      = makeglobal(el, mesh); 

[Jint]    = getMask(mesh,gl,el,mat); % masking
saveas(gcf,[Dir.results_dir '/' num2str(Dir.scan_ID) '.' ...
        num2str(Dir.Trail_number) '_selection.png']); close all

[el,mesh,gl,Jint]=StrainInegration(mat,mesh,el,Jint,gl);% STRAIN INTEGRATION
if Operation == 'Str'; saveas(gcf,fullfile(Dir.results_dir,[num2str(Dir.scan_ID) '.' ...
        num2str(Dir.Trail_number) '_integrated_displacements.png'])); close all; end
    
[Results,Jint]  = J_IntegralCacl(Jint,mat,mesh,el); % J-INTEGRAL CALCULATION
[Results,error] = WestergaardCorrection(Results,mat,Operation,mesh.x.xi,mesh.x.length); % correction
[Jtrue,Jdiv,Ktrue,Kdiv] = PlotAll(Results,Dir,error,mesh,el,gl,Jint); % Plotting 