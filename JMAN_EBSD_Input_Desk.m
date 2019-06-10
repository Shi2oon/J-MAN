colordef white; clc;  close all; warning('off'); clear

input.source='P:\Abdo\EBSD Data\Alex xEBSD';
I.II1='P:\Abdo\EBSD Data\Alex xEBSD\001_berko_3rd_face_XEBSD'; I.ii1=13;
I.II6='P:\Abdo\EBSD Data\Alex xEBSD\Si_1000nm_unirr_32_XEBSD'; I.ii6=100032;
I.II2='P:\Abdo\EBSD Data\Alex xEBSD\1000nm unirradiated silicon'; I.ii2=10002;
I.II3='P:\Abdo\EBSD Data\Alex xEBSD\1700x2000nmBerkoindent_XEBSD'; I.ii3=17002;
I.II4='P:\Abdo\EBSD Data\Alex xEBSD\berko_2nd_face_XEBSD'; I.ii4=2;
I.II5='P:\Abdo\EBSD Data\Alex xEBSD\indent_25_XEBSD'; I.ii5=25;
%I.II7='P:\Abdo\EBSD Data\Alex xEBSD\Si_unirr_indent3-2_XEBSD'; I.ii7=32;
%I.II8='P:\Abdo\EBSD Data\Alex xEBSD\6H_750_irr_1000nm_8-2-19_XEBSD'; I.ii2=6750;

input.filesCount=6;
input.avJ=zeros(input.filesCount*3,3); input.avK=zeros(input.filesCount*3,3);
input.Cracks=3; %number for cracks to look into

for Ii=1:input.filesCount
    eval(sprintf('input.fillname=I.ii%d;',Ii));
    eval(sprintf('input.Dir=I.II%d;',Ii));
    input.results=fullfile(input.Dir,'J-MAN 2'); mkdir(input.results);
%     
%     for ij=1:9 % 3 for each tip
%         input.Trail_number=ij;   JMAN_EBSD
%     end
%     
    jii=0;  ij=0;
    for kk=1:input.Cracks
        ij=input.Cracks*kk; jii=ij-input.Cracks+1;
        Dir.operation='J';          KidealEBSD
        Dir.operation='K';          KidealEBSD
    end
end

[yaxis,~]=ind2sub(size(input.avJ),find(input.avJ==0));
input.avJ(min(yaxis):end,:) = [];
plot (input.avJ,'--o'); hold on;
set(gcf,'position',[400,100,1400,750])
legend (num2str((1:Ii)'),'Location','NorthEastOutside')
title(['J-intergal Values for ' num2str(input.Cracks) ' Cracks/Corners'])
plot (trimmean(input.avJ,50,2),':>','DisplayName','Avg.'); hold off;
xlabel('Contour Number'); ylabel('J [J m^-^2]'); 
Dir.pyxe_D_path = fullfile(input.source,[num2str(ij) ' Js 2.png']);
saveas(gcf,Dir.pyxe_D_path)

[yaxis,~]=ind2sub(size(input.avK),find(input.avK==0));
input.avK(min(yaxis):end,:) = [];
plot (input.avK,'--o'); hold on;
legend (num2str((1:Ii)'),'Location','NorthEastOutside')
title(['K Values for ' num2str(input.Cracks) ' Cracks/Corners'])
plot (trimmean(input.avK,50,2),':>','DisplayName','Avg.'); hold off;
xlabel('Contour Number'); ylabel('K / MPa m^{1/2}'); 
Dir.pyxe_D_path = fullfile(input.source,[num2str(ij) ' Ks 2.png']);
saveas(gcf,Dir.pyxe_D_path)

Table=table(input.avK(:,1),input.avEK(:,1),'VariableNames',{'K_in_MPa','K_Error'});
for i=2:input.Cracks
Table = [Table table(input.avK(:,i),input.avEK(:,i),'VariableNames',...
    {['K_in_MPa_' num2str(i)],['K_Error_' num2str(i)]})];
end
for i=1:input.Cracks
Table = [Table table(input.avJ(:,1),input.avEJ(:,1),'VariableNames',...
    {['J_in_MPa_' num2str(i)],['J_Error_' num2str(i)]})];
end

Dir.pyxe_D_path = fullfile(input.source,[num2str(Ii) ' Data sets with '...
    num2str(input.Cracks) ' Cracks.xlsx']);
save(Dir.pyxe_D_path, 'Table');