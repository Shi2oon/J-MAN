%% FE code to evaluate J-integral from HR-EBSD data,
% by Thorsten Becker - b.thorsten@gmail.com
% and Selim Barhli - selim.barhli@ensiacet.fr
% with modifications made by Phil Earp - philip.earp@materials.ox.ac.uk

%DIC IMAGE:         <--------------mesh.winsize(1)---------->
%                  ...........................................
%             |    :  FE----------------------------------   :                                            .
%             |    :  |  Jo----------------------------  |   :
%             |    :  |  |  Ji----------------------  |  |   :
% mesh.winsize(2)  :========> crack in pos x dir.  |  |  |   :
%             |    :  |  |  ----------------------Ji  |  |   :
%             |    :  |  ----------------------------Jo  |   :
%             |    :  ----------------------------------FE   :                                            .
%                  :.........................................:
% FE => mesh.winFE(x&y,point)-window size of image considered for the FE analysis
% Jo => Jint.nout(x&y,point) -outer Jintegral path
% Ji => Jint.nin(x&y,point)  -inner Jintegral path
%VARIABLES:
% mat.  - material and element properties
% mesh. - FE mesh
% el.   - FE elements, variables and results
% Jint. - J-integral elements, variables and results
% gl.   - Global interpolation of FE analysis from Gauss points

%% set the scene
% JMAN_EBSD_Dashboard

% clear all; 
colordef white; clc;  close all; warning('off'); clearvars -except input I Ii ij

% input.fillname='13';
% input.Trail_number=i; %try to do at least 4
% input.Dir='P:\Abdo\EBSD Data\Alex xEBSD\001_berko_3rd_face_XEBSD';
input.fullpath = fullfile(input.Dir,[num2str(input.fillname) '.mat']);
% input.results=fullfile(input.Dir,'J-MAN 3'); mkdir(input.results);
input.unit = 'mm'; %meter (m) or milmeter (mm) or micrometer(\mum) ;
input.StressUnits = 'MPa';
input.rotate180 = false;
input.isXEBSD = true;

%%
addpath([pwd '\functions']);
addpath([pwd '\functions\LoadEBSD']);

fprintf(1,'FE code to evaluate J-integral from HR-EBSD Data \nby Phil Earp - philip.earp@materials.ox.ac.uk\n');
fprintf(1,'S.Barhli and T.Becker\n\n');
fprintf('Trial number %d\n\n',input.Trail_number);

%% LOAD, ROTATE & CROP EBSD DATA
fprintf('\nLoading EBSD data...\n');
[EBSDdata,input] = LoadEBSD( input );
fprintf('...done\n')
%% INPUT MATERIAL PROPERTIES AND DATA
fprintf('\nLoading data...');
%Material properties
mat.nu = 0.36; % Poisson's ratio
mat.E  = 110E9; % Young's Modulus [Pa]
mat.yield = 1E9; % Yield Stress [Pa]
fprintf ('Poissons ratio %0.1i\nYoungs modulus %0.1i GPa\nTensile Strength %0.1i MPa\n\n'...
          ,mat.nu,mat.E/1e9,mat.yield/1e6);
      
mat.stressstate = 'plane_stress'; %set to'plane_strain' or 'plane_stress'
[mat] = matprop(mat);                %sets material peoperties: E,nu,etc...

%% INPUT MATERIAL PROPERTIES AND DATA
fprintf('\nPreparing data...');
[mat] = matprop();                %sets material peoperties: E,nu,etc...
[mesh] = loadUd(EBSDdata);
fprintf('done\n')
fprintf(1,'data points (gauss points) defined:%6.0f \n',length(mesh.dispgraddata));
fprintf(1,'nodes created:%6.0f \n',length(mesh.UDIC));
close all

%% MESHING & STRAIN/STRESSES CALCULATIONS
fprintf('\nMeshing...');
[el,mesh] = meshDIC(mesh);
fprintf('done\n')
fprintf(1,'Number of elements:%6.0f (%3.0f x%3.0f)\n',mesh.elFE,mesh.winFE(2,2)-mesh.winFE(1,2),mesh.winFE(2,1)-mesh.winFE(1,1));
% FEM (calculate element stress/strain & assemble into global system)
fprintf('\nFEM analysis');
[el] = FEanalysisStrains(el,mat,mesh);%runs isoparematric FE anlysis to slove for strain and stress at Gauss points
[gl] = makeglobal(el, mesh);         %Interpolates FE results into a global domain at nodal co-ordinates
fprintf(': done\n');

% save('tmp.mat')
%% CONTOUR DETERMINATION
% load('tmp.mat')

isDispOk = 'N';H2 = [];
while (isDispOk ~= 'Y')
    %Printing screen for manual choice of J-Integral contour, and mask function
    close(H2);
    H1=prePlot(mesh,gl,1,el);
    
    title('Select Crack Tip')
    [estcrktip(1), estcrktip(2)] = ginput(1);
    plot(estcrktip(1),estcrktip(2),'r.','MarkerSize',20); drawnow;
    
    title('Select (Multiple) Masks - Press ENTER when done')
    t=0;id=1;
    while (t<1)
        [xmask, ymask] = ginput(2);
        rectangle('Position',[min(xmask(1),xmask(2)),min(ymask(1),ymask(2)),abs(xmask(2)-xmask(1)),abs(ymask(2)-ymask(1))],'FaceColor','k'); drawnow;
        [xmask,ymask] = convCoord(xmask,ymask,mesh.winFE,mesh.UFE);
        Jint.mask(:,:,id) = [ymask(1) ymask(2);xmask(1) xmask(2)]; Jint.mask(Jint.mask<1)=1;
        id = id+1;
        t = waitforbuttonpress();
    end
    
    title('Select Outer J-Contour')
    [xo,yo] = ginput(2);
    rectangle('Position',[min(xo(1),xo(2)),min(yo(1),yo(2)),abs(xo(2)-xo(1)),abs(yo(2)-yo(1))],'EdgeColor','W'); drawnow;
    [xo,yo] = convCoord(xo,yo,mesh.winFE,mesh.UFE);
    
    title('Select Inner J-Contour')
    [xi,yi] = ginput(2);
    rectangle('Position',[min(xi(1),xi(2)),min(yi(1),yi(2)),abs(xi(2)-xi(1)),abs(yi(2)-yi(1))],'EdgeColor','W');drawnow;
    [xi,yi] = convCoord(xi,yi,mesh.winFE,mesh.UFE);
    
    Jint.nout = [max(yo) min(yo); min(xo) max(xo)];
    Jint.nin  = [max(yi) min(yi); min(xi) max(xi)];
    
    % prepare J-integral mask
    [Jint] = findJintmask(mat, mesh, el, Jint);
    isDispOk = 'Y';
end


%% J-INTEGRAL CALCULATION


[Ctrs,Nb_ct] = CtrCalc(Jint);
Results = zeros(Nb_ct,2);
fprintf('\nJ-Integrating...\n');
for i=1:Nb_ct
    Jint.nin(1,:)=Ctrs(i,1:2);
    Jint.nin(2,:)=Ctrs(i,3:4);
    Jint.nout(1,:)=Ctrs(i,5:6);
    Jint.nout(2,:)=Ctrs(i,7:8);
    %         Jint.mask = [ymask(1) ymask(2);xmask(1) xmask(2)]; Jint.mask(Jint.mask<1)=1;
    [Jint] = findJintelem(mat, mesh, el, Jint);  %finds elements within defined (Jint)ergral area, Jint.nout/nin
    [Jint] = getq(el,Jint, mesh);       %creates smooth virtual crack extension function q over Jint.nout/nin
    [Jint] = Jcalc(el,Jint, mat);
    
    Results(i,1)=Jint.J.*1000000; % [J/mm^2] --> [J/m^2]
    Results(i,2)=Jint.K; % [MPa sqrt(m)]
    
    if(i>Nb_ct)
        clear Jint;
    end
    %           plotDIC(mesh,el,mat,Jint);
end

    fprintf(1,'Average J is %3.4f J/m^2\n', mean(Results(:,1)));
    fprintf(1,'Average K is %3.4f MPa sqrt(m)\n', 1E-6 * mean(Results(:,2)));
    Jdiv=std(Results(:,1));
    fprintf(1,'Std Dev for J is %3f\n', Jdiv);
    fprintf(1,'Std Dev for K is %3f\n\n', 1E-6 * std(Results(:,2)));
    
    fprintf(1,'\nMean K is %3.4f MPa sqrt(m) over last 50pc\n', 1E-6 * mean(Results(round(Nb_ct/2):end, 2)));
    Kdiv=1e-6 * std(Results(:,2));
    fprintf(1,'Std Dev for K is %3.4f MPa sqrt(m) over last 50pc\n', Kdiv);
    Jtrue=((mean(Results(:,1))+max(Results(:,1)))/2);
    fprintf(1,'True J is %3.2f J/m^2\n', Jtrue);
    Ktrue=1E-6.*((mean(Results(:,2))+max(Results(:,2)))/2);
    fprintf(1,'True K is %3.2f MPa sqrt(m)\n', Ktrue);



%% Plot
fprintf('\nPlotting results...');
figure; plot(Results(:,1));
xlabel('Contour Number'); ylabel('J [J m^-^2]');
input.Save = fullfile(input.results,[ num2str(input.fillname) '.' ...
    num2str(input.Trail_number) '_Results J.png']);
saveas(gcf,input.Save); close all

plotFE(mesh,el,mat,Jint,gl,input);
input.Save = fullfile(input.results,[ num2str(input.fillname) '.' ...
    num2str(input.Trail_number) '_output.png']);
saveas(gcf,input.Save); close all

plotFEe(mesh,el,mat,Jint,gl,input)
input.Save = fullfile(input.results,[ num2str(input.fillname) '.' ...
    num2str(input.Trail_number) '_FE.png']);
saveas(gcf,input.Save); close all

alldata = table(Results(:,1), Results(:,2),'VariableNames',{ 'J', 'K' } );
input.Save = fullfile(input.results,[ num2str(input.fillname) '.' num2str(input.Trail_number) '_Result.mat']);
save(input.Save, 'alldata','Ktrue','Kdiv','Jtrue','Jdiv');

fprintf('trial %d is done\n',input.Trail_number);
