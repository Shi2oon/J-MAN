function [Jtrue,Jdiv,Ktrue,Kdiv] = PlotAll(Results,Dir,ER,mesh,el,gl,Jint)
if ER == 0;      Results(:,3:4) = Results(:,1:2);        end

%% calculate J and K values at the last 1/5 contours
contrs = length(Results);           contrs = contrs - round(contrs/4);
clc; printCredit
    Jtrue = ((mean(Results(contrs:end,3))+max(Results(contrs:end,3)))/2);
    Jdiv  = std([std(Results(contrs:end,1)) std(Results(contrs:end,3))]);
    fprintf(1,'\nTrue J is %3.2f J/m^2 with ±%3.3f Std Dev\n', Jtrue, Jdiv);
    Ktrue = 1E-6.*((mean(Results(contrs:end,4))+max(Results(contrs:end,4)))/2);
    Kdiv  = 1E-6 * std([std(Results(contrs:end,2)) std(Results(contrs:end,4))]);
    fprintf(1,'True K is %3.2f MPa sqrt(m) with ±%3.3f Std Dev\n', Ktrue, Kdiv);
    fprintf('\nPlotting results ...'); close all

% errorbar(1E-6* Results(:,4),0.5E-6*abs(Results(:,2)-Results(:,4)),'k--o',...
%     'HandleVisibility','off','LineWidth',1); hold on
plot(1E-6* Results(:,4),'r--o','MarkerEdgeColor','r','LineWidth',1.5,'MarkerFaceColor','r');
set(gcf,'position',[600,100,950,650])
xlabel('Contour Number'); 
title ({['Stress Intensity Factor (K_I) = ' num2str(round(Ktrue,2)) ' ± ' ...
    num2str(round(Kdiv,3)) ' MPa\surdm' ];''});
ylabel('K (MPa m^{1/2})')
pyxe_D_path = fullfile(Dir.results,[ Dir.fillname '.'...
    num2str(Dir.Trail_number) '_ResultsK.png']);
saveas(gcf,pyxe_D_path); 
pyxe_D_path = fullfile(Dir.results,[ Dir.fillname '.'...
    num2str(Dir.Trail_number) '_ResultsK.fig']);
saveas(gcf,pyxe_D_path); close all;

% plot(1E-6* Results(:,2),'--sb','MarkerEdgeColor','k','MarkerFaceColor','b','LineWidth',1.5); 
% hold off;       legend('Corrected','Original'); 
% pyxe_D_path = fullfile(Dir.results_dir,[ Dir.scan_ID '.'...
%     num2str(Dir.Trail_number) '_Results.png']);saveas(gcf,pyxe_D_path); close all

%% plot Images
plotFE(mesh,el,Jint,gl,Dir) %
pyxe_D_path = fullfile(Dir.results,[ Dir.fillname '.' num2str(Dir.Trail_number) '_In.png']);
saveas(gcf,pyxe_D_path);
pyxe_D_path = fullfile(Dir.results,[ Dir.fillname '.' num2str(Dir.Trail_number) '_In.fig']);
saveas(gcf,pyxe_D_path); close all

plotFEe(mesh,el,gl,Dir) %uy displacement field
pyxe_D_path = fullfile(Dir.results,[ Dir.fillname '.' num2str(Dir.Trail_number) '_Out.png']);
saveas(gcf,pyxe_D_path); 
pyxe_D_path = fullfile(Dir.results,[ Dir.fillname '.' num2str(Dir.Trail_number) '_Out.fig']);
saveas(gcf,pyxe_D_path); close all

if  any(ismember(fields(mesh),'OpS'))==1
    if sum(strfind(mesh.OpS,'xED'))~=0
        CroppedPlot2(Dir,gl,el,Jint,20)
pyxe_D_path = fullfile(Dir.results,[ Dir.fillname '.' num2str(Dir.Trail_number) '_All.png']);
        saveas(gcf,pyxe_D_path); 
pyxe_D_path = fullfile(Dir.results,[ Dir.fillname '.' num2str(Dir.Trail_number) '_All.fig']);
        saveas(gcf,pyxe_D_path); close all
    end
end

%% save
alldata = table(Results(:,3),Results(:,4),Results(:,1),Results(:,2),...
    'VariableNames',{ 'J', 'K','J_Raw', 'K_Raw' } );
pyxe_D_path = fullfile(Dir.results,[ Dir.fillname '.'...
    num2str(Dir.Trail_number) '_Result.mat']);
ER=abs(1-ER);
save(pyxe_D_path, 'alldata','Ktrue','Jtrue','ER');

fprintf(' done\n'); 
