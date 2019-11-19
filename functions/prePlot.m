function DIC = prePlot(mesh,gl,y,Pla_reg,el)
%%  
 %Figure
 if(y)
  DIC = figure('Position',[50 100 1000 800]);
 elseif(y==0)
  DIC = figure('Position',[50 100 500 400]);
 end
  axesDIC = axes('Parent',DIC);
  hold(axesDIC,'all');
 
 %Displacment field (DIC)
  hdisp = imagesc(mesh.UDIC(1,1:mesh.winDIC(2):end), mesh.UDIC(2,1:mesh.winDIC(2)),...
          reshape(sqrt(mesh.dDIC(1,:).^2+mesh.dDIC(2,:).^2),mesh.winDIC)');
          shading flat;
 %DIC grid
  hgrid = gridxy(mesh.UDIC(1,1:mesh.winDIC(2):end), mesh.UDIC(2,1:mesh.winDIC(2)),...
          'LineStyle','-','Color',[0.5 0.5 0.5]); 
 
 %plot strain contours
  hdisp2 = pcolor(gl.Uxdef,gl.Uydef,gl.smises);
  
  %plasticity zone
  hPla = fill([ el.Uxdef(Pla_reg,1)';el.Uxdef(Pla_reg,2)';el.Uxdef(Pla_reg,3)';...
                el.Uxdef(Pla_reg,4)';el.Uxdef(Pla_reg,1)'], ...
                [el.Uydef(Pla_reg,1)';el.Uydef(Pla_reg,2)';el.Uydef(Pla_reg,3)';...
                el.Uydef(Pla_reg,4)';el.Uydef(Pla_reg,1)'], ... 
                'm','FaceAlpha',0.25);
 
 %Row and Colunm numbering
   text(mesh.UDIC(1,1),mesh.UDIC(2,1),'1',...
        'VerticalAlignment','top','HorizontalAlignment','left',...
        'FontSize',6,'Color',[1 1 1]);
  text(mesh.UDIC(1,2:mesh.winDIC(2)),mesh.UDIC(2,2:mesh.winDIC(2)),...
      num2str((2:mesh.winDIC(2))'),... %Row numbering
       'VerticalAlignment','top','HorizontalAlignment','left',...
       'FontSize',6,'Color',[1 1 1]);
  text(mesh.UDIC(1,(mesh.winDIC(2)+1):mesh.winDIC(2):end),...%Column numbering
      mesh.UDIC(2,(mesh.winDIC(2)+1):mesh.winDIC(2):end),num2str((2:mesh.winDIC(1))'),... 
       'VerticalAlignment','middle','HorizontalAlignment','left',...
       'FontSize',6,'Color',[1 1 1]);
  
 
 %FE displacments 
%   hquiv = quiver(mesh.UFE(1,(mesh.dFE(1,:)~=0)),mesh.UFE(2,(mesh.dFE(2,:)~=0)),...
%    mesh.dFE(1,(mesh.dFE(1,:)~=0)),mesh.dFE(2,(mesh.dFE(2,:)~=0)),...
%      'Color',[1 1 1],'AutoScaleFactor',1,'Parent',axesDIC);
 
  
 %Create xlabel
  xlabel('along specimen length (m)');
 %Create ylabel
  ylabel('along specimen width (m)');
 %Create title
  title({'DIC/Stain field showing FE elements and J integral elements';''});
 
 %Create colorbar
colorbar('peer',axesDIC);   
set(gcf,'position',[500,100,1200,800])
 % axis([min(mesh.UDIC(1,:)) max(mesh.UDIC(1,:)) min(mesh.UDIC(2,:)) max(mesh.UDIC(2,:))]);
  axis equal;
  colormap(jet)
axis([min(mesh.UDIC(1,:)) max(mesh.UDIC(1,:)) min(mesh.UDIC(2,:)) max(mesh.UDIC(2,:))]);

