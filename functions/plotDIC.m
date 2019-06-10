function plotDIC(mesh,el,mat,Jint)
%%  
 %Figure 
  DIC = figure('Position',[50 100 1000 800]);
  axesDIC = axes('Parent',DIC);
  hold(axesDIC,'all');
 
 %Displacment field (DIC)
  hdisp = imagesc(mesh.UDIC(1,1:mesh.winDIC(1)),mesh.UDIC(2,1:mesh.winDIC(2):end),...
          reshape(sqrt(mesh.dDIC(1,:).^2+mesh.dDIC(2,:).^2),mesh.winDIC)');
          shading flat;
 %DIC grid
  hgrid = gridxy(mesh.UDIC(1,1:mesh.winDIC(1)), mesh.UDIC(2,1:mesh.winDIC(1):end),...
          'LineStyle',':','Color',[0.5 0.5 0.5]); 
 
 %plot strain contours
%   hdisp2 = pcolor(gl.Uxdef,gl.Uydef,gl.smises);
  
 %Row and Colunm numbering
   text(mesh.UDIC(1,1),mesh.UDIC(2,1),'1',...
        'VerticalAlignment','top','HorizontalAlignment','left',...
        'FontSize',6,'Color',[1 1 1]);
  text(mesh.UDIC(1,2:mesh.winDIC(1)),mesh.UDIC(2,2:mesh.winDIC(1)),num2str((2:mesh.winDIC(1))'),...
       'VerticalAlignment','top','HorizontalAlignment','center',...
       'FontSize',6,'Color',[1 1 1]);
  text(mesh.UDIC(1,(mesh.winDIC(1)+1):mesh.winDIC(1):end),mesh.UDIC(2,(mesh.winDIC(1)+1):mesh.winDIC(1):end),num2str((2:mesh.winDIC(2))'),...
       'VerticalAlignment','middle','HorizontalAlignment','left',...
       'FontSize',6,'Color',[1 1 1]);
 %FE elements and nodes 
  switch mat.eltype 
      case 'Q4'
          hFEel = plot([el.Ux(:,1)';el.Ux(:,2)';el.Ux(:,3)';el.Ux(:,4)'], ...
               [el.Uy(:,1)';el.Uy(:,2)';el.Uy(:,3)';el.Uy(:,4)'], ...
               '.-w','MarkerSize',10);
      case 'Q8'
          hFEel = plot([el.Ux(:,1)';el.Ux(:,5)';el.Ux(:,2)';el.Ux(:,6)';el.Ux(:,3)';el.Ux(:,7)';el.Ux(:,4)';el.Ux(:,8)';el.Ux(:,1)'], ...
                       [el.Uy(:,1)';el.Uy(:,5)';el.Uy(:,2)';el.Uy(:,6)';el.Uy(:,3)';el.Uy(:,7)';el.Uy(:,4)';el.Uy(:,8)';el.Uy(:,1)'], ...
                       '.-W','MarkerSize',10);
      case 'Q9'
          plot([el.Ux(:,1)';el.Ux(:,2)';el.Ux(:,3)';el.Ux(:,4)';el.Ux(:,5)';el.Ux(:,6)';el.Ux(:,7)';el.Ux(:,8)';el.Ux(:,9)'], ...
               [el.Uy(:,1)';el.Uy(:,2)';el.Uy(:,3)';el.Uy(:,4)';el.Uy(:,5)';el.Uy(:,6)';el.Uy(:,7)';el.Uy(:,8)';el.Uy(:,9)'], ...
               '.w','MarkerSize',10);
          hFEel = plot([el.Ux(:,1)';el.Ux(:,5)';el.Ux(:,2)';el.Ux(:,6)';el.Ux(:,3)';el.Ux(:,7)';el.Ux(:,4)';el.Ux(:,8)';el.Ux(:,1)'], ...
               [el.Uy(:,1)';el.Uy(:,5)';el.Uy(:,2)';el.Uy(:,6)';el.Uy(:,3)';el.Uy(:,7)';el.Uy(:,4)';el.Uy(:,8)';el.Uy(:,1)'], ...
               '.-w','MarkerSize',10);
  end
  
 %FE displacments 
  hquiv = quiver(mesh.UFE(1,(mesh.dFE(1,:)~=0)),mesh.UFE(2,(mesh.dFE(2,:)~=0)),mesh.dFE(1,(mesh.dFE(1,:)~=0)),mesh.dFE(2,(mesh.dFE(2,:)~=0)),...
     'Color',[1 1 1],'AutoScaleFactor',1,'Parent',axesDIC);
 
 %Jintegral Elements 
  hJint = fill([el.Ux(Jint.el,1)';el.Ux(Jint.el,2)';el.Ux(Jint.el,3)';el.Ux(Jint.el,4)';el.Ux(Jint.el,1)'], ...
          [el.Uy(Jint.el,1)';el.Uy(Jint.el,2)';el.Uy(Jint.el,3)';el.Uy(Jint.el,4)';el.Uy(Jint.el,1)'], ... 
         'k','EdgeColor','none','FaceAlpha',0.4);
 
  
 %Create xlabel
  xlabel('along specimen length (mm)');
 %Create ylabel
  ylabel('along specimen width (mm)');
 %Create title
  title('DIC displacment field showing FE elements and J integral elements');
 %Create legend
  legend1 = legend(axesDIC,[hgrid hFEel(1) hquiv hJint(1)], ...
            'DIC data grid','FE elements','FE nodal displacments','J integral elements');
  set(legend1,'Orientation','horizontal','Location','SouthOutside','Color',[0.5 0.5 0.5]);
  legend boxoff;
 %Create colorbar
  colorbar('peer',axesDIC);   
    
  axis([min(mesh.UDIC(1,:)) max(mesh.UDIC(1,:)) min(mesh.UDIC(2,:)) max(mesh.UDIC(2,:))]);
  axis square;
  
  
