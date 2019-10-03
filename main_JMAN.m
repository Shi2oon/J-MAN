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
% FE => mesh.winFE(x&y,point)-window size of image considered for the FE analysis
% Jo => Jint.nout(x&y,point)-outer Jintegral path     % Ji=>Jint.nin(x&y,point)-inner Jintegral path
%VARIABLES:             % gl.   - Global interpolation of FE analysis from Gauss points
% mat.  - material and element properties    % Jint. - J-integral elements, variables and results
% el.   - FE elements, variables and results % mesh. - FE mesh 

% close all;      clear;      colordef white;     clc;        warning off; 
addpath(genpath(pwd));              printCredit;        
set(0,'defaultAxesFontSize',20);    set(0,'DefaultLineMarkerSize',12)   
 
%% Setting the scene
Dir.fullpath     = 'P:\XiaoCa'; % filr Directory
Dir.fillname     = 'forJMAN';   % all in .mat w/ data in alldata var except DVC in .dat format
Dir.Trail_number = 3;           % file number %try to do at least 4
Dir.pixel_size   = 1;           % if DIC/DVC values are in pixel, 1 if in physical units;
Dir.input_unit   = 'mm';        % meter (m) or milmeter (mm) or micrometer(um);
Operation        = 'Dis';       % calculation mode\n Dis = DIC, Str = Strain, xED = xEBSD, DVC
Dir.results      = fullfile(Dir.fullpath ,'J-MAN');          mkdir(Dir.results);

%% INPUT MATERIAL PROPERTIES AND DATA
% Poisson's ratio,     Young's Modulus [Pa],   Yield Stress [Pa] 
  Dir.nu    = 0.2;     Dir.E  = 10.5E9;         Dir.yield = 820E6;
      
%% Load, MESHING & STRAIN/STRESSES CALCULATIONS
% loadUD output is in meter
[mesh,Dir] = loadUd(Dir.pixel_size,Dir.fullpath, Dir.fillname, Dir.input_unit,Operation,Dir);
[mat]      = matprop(Dir.E,Dir.nu,Dir.yield,Dir.stressstat,Dir.input_unit); % material props.
[el,mesh]  = meshDIC(mesh);
%runs isoparematric FE anlysis to slove for strain and stress at Gauss points
[el]      = FEanalysisStrains(el,mat,mesh);
%Interpolates FE results into a global domain at nodal co-ordinates
[gl]      = makeglobal(el, mesh); 

%% Checking if the crack is on the good axis
Pla_reg = isplastic(el,mat.yield);
H1 = prePlot(mesh,gl,0,Pla_reg,el);
[mesh,mat,el,gl] = Crack_align(mesh,mat,el,gl);
H1 = prePlot2(mesh,gl,0,Pla_reg,el);
close(H1);
%% J_intergal calc.
[Jint]    = getMask(mesh,gl,el,mat); % masking
saveas(gcf,[Dir.results '/' Dir.fillname '.' ...
        num2str(Dir.Trail_number) '_selection.png']); close all
[el,mesh,gl,Jint]=StrainInegration(mat,mesh,el,Jint,gl);% STRAIN INTEGRATION, gl.dy & dx
if mesh.Operation == 'Str'; saveas(gcf,fullfile(Dir.results,[Dir.fillname '.' ...
        num2str(Dir.Trail_number) '_integrated_displacements.fig'])); close all; end
    
[Results,Jint]  = J_IntegralCacl(Jint,mat,mesh,el); % J-INTEGRAL CALCULATION
[Results,error] = WestergaardCorrection(Results,mat,mesh.Operation,mesh.x.xi,mesh.x.length);
[Jtrue,Jdiv,Ktrue,Kdiv] = PlotAll(Results,Dir,1,mesh,el,gl,Jint); % Plotting 
% in results [J K J_corrected K_corrected]
save([Dir.results '\' Dir.fillname '.' num2str(Dir.Trail_number) '.mat']);