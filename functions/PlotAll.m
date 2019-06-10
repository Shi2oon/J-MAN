function [Jtrue,Jdiv,Ktrue,Kdiv] = PlotAll(Results,Dir,error,mesh,el,gl,Jint)
%     fprintf(1,'Average J is %3.2f J/m^2\n', mean(Results(:,3)));
%     fprintf(1,'Average K is %3.2f MPa sqrt(m)\n', 1E-6 * mean(Results(:,4)));
   
%     fprintf(1,'\nMean K is %3.2f MPa sqrt(m) over last 50pc\n', 1E-6 *...
%         mean(Results(round(length(Results)/2):end, 4)));
%     Kdiv=1e-6 * std(Results(:,4));
%     fprintf(1,'Std Dev for K is %3.2f MPa sqrt(m) over last 50pc\n\n', Kdiv);
clc; printCredit
    Jtrue = ((mean(Results(:,3))+max(Results(:,3)))/2);
    Jdiv  = std(Results(:,3))+Jtrue*abs(1-error);
    fprintf(1,'\nTrue J is %3.2f J/m^2 with ±%3.3f Std Dev\n', Jtrue, Jdiv);
    Ktrue = 1E-6.*((mean(Results(:,4))+max(Results(:,4)))/2);
    Kdiv  = 1E-6 * std(Results(:,4))+Ktrue*abs(1-error);
    fprintf(1,'True K is %3.2f MPa sqrt(m) with ±%3.3f Std Dev\n', Ktrue, Kdiv);
    
set(0,'defaultAxesFontSize',15)
set(0,'DefaultLineMarkerSize',12)    
fprintf('\nPlotting results ...'); close all

errorbar(1E-6* Results(:,4),1E-6* Results(:,4).*abs((error-1)/2),'k--o',...
    'HandleVisibility','off','LineWidth',1); hold on
plot(1E-6* Results(:,4),'r--o','MarkerEdgeColor','r','LineWidth',1.5);
set(gcf,'position',[600,100,950,650])
xlabel('Contour Number'); 
title (['Stress Intensity Factor (K_I) = ' num2str(round(Ktrue,2)) ' ± ' ...
    num2str(round(Kdiv,3)) ' MPa\surdm' ]);
ylabel('K (MPa m^{1/2})')
pyxe_D_path = fullfile(Dir.results_dir,[ num2str(Dir.scan_ID) '.'...
    num2str(Dir.Trail_number) '_ResultsK.png']);
saveas(gcf,pyxe_D_path); 

% plot(1E-6* Results(:,2),'--sb','MarkerEdgeColor','k','MarkerFaceColor',...
%     'b','LineWidth',1.5); hold off
% legend('Corrected','Original'); 
% pyxe_D_path = fullfile(Dir.results_dir,[ num2str(Dir.scan_ID) '.'...
%     num2str(Dir.Trail_number) '_Results.png']);saveas(gcf,pyxe_D_path); 
close all

plotFE(mesh,el,Jint,gl)
pyxe_D_path = fullfile(Dir.results_dir,[ num2str(Dir.scan_ID) '.'...
    num2str(Dir.Trail_number) '_In.png']);
saveas(gcf,pyxe_D_path); close all

plotFEe(mesh,el,gl)
pyxe_D_path = fullfile(Dir.results_dir,[ num2str(Dir.scan_ID) '.'...
    num2str(Dir.Trail_number) '_Out.png']);
saveas(gcf,pyxe_D_path); close all

%% save
alldata = table(Results(:,3),Results(:,4),Results(:,1),Results(:,2),...
    'VariableNames',{ 'J', 'K','J_Raw', 'K_Raw' } );
pyxe_D_path = fullfile(Dir.results_dir,[ num2str(Dir.scan_ID) '.'...
    num2str(Dir.Trail_number) '_Result.mat']);
error=abs(1-error);
save(pyxe_D_path, 'alldata','Ktrue','Jtrue','error');

fprintf(' done\n'); 
