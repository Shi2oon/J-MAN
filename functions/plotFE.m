function plotFE(mesh,el,Jint,gl)
if mesh.Operation=='Dis'
warning off
 %Figure 
  FE = figure('Position',[50 100 1000 800]);
  axesFE = axes('Parent',FE);
  hold(axesFE,'all');
 
 %plot strain contours
  colormap(jet)
  hdisp = pcolor(gl.Uxdef,gl.Uydef,gl.dy,'Parent',axesFE);
           shading interp
hcont = contour(gl.Uxdef,gl.Uydef,gl.dy,10,'LineWidth',1,'LineColor',...
    [0 0 0],'Parent',axesFE);
%   hdisp = pcolor(gl.Uxdef,gl.Uydef,gl.smises,'Parent',axesFE);
%           shading interp
%   hcont = contour(gl.Uxdef,gl.Uydef,gl.smises,10,'LineWidth',1,...
%       'LineColor',[0 0 0],'Parent',axesFE);
%   
%   hdisp = pcolor(gl.Uxdef,gl.Uydef,gl.dy,'Parent',axesFE);
%           shading interp
%           set(gca,'YDir','normal');
%  hdisp = scatter(gl.Uxdef(:),gl.Uydef(:),3,gl.dy(:));
   
%   hcont = contour(gl.Uxdef,gl.Uydef,gl.smises,10,'LineWidth',1,'LineColor',[0 0 0],'Parent',axesFE);
  
%  %FE elements and nodes 
%   switch mat.eltype 
%       case 'Q4'
%          hFEel = plot([el.Uxdef(1:20,1)';el.Uxdef(1:20,2)';el.Uxdef(1:20,3)';el.Uxdef(1:20,4)'], ...
%                       [el.Uydef(1:20,1)';el.Uydef(1:20,2)';el.Uydef(1:20,3)';el.Uydef(1:20,4)'], ...
%                        '.w','MarkerSize',10);
%       case 'Q8'
%           hFEel = plot([el.Uxdef(:,1)';el.Uxdef(:,5)';el.Uxdef(:,2)';el.Uxdef(:,6)';el.Uxdef(:,3)';el.Uxdef(:,7)';el.Uxdef(:,4)';el.Uxdef(:,8)';el.Uxdef(:,1)'], ...
%                        [el.Uydef(:,1)';el.Uydef(:,5)';el.Uydef(:,2)';el.Uydef(:,6)';el.Uydef(:,3)';el.Uydef(:,7)';el.Uydef(:,4)';el.Uydef(:,8)';el.Uydef(:,1)'], ...
%                         '.:W','MarkerSize',10);
%       case 'Q9'
%           plot([el.Uxdef(:,1)';el.Uxdef(:,2)';el.Uxdef(:,3)';el.Uxdef(:,4)';el.Uxdef(:,5)';el.Uxdef(:,6)';el.Uxdef(:,7)';el.Uxdef(:,8)';el.Uxdef(:,9)'], ...
%                [el.Uydef(:,1)';el.Uydef(:,2)';el.Uydef(:,3)';el.Uydef(:,4)';el.Uydef(:,5)';el.Uydef(:,6)';el.Uydef(:,7)';el.Uydef(:,8)';el.Uydef(:,9)'], ...
%                '.w','MarkerSize',10);
%           hFEel = plot([el.Uxdef(:,1)';el.Uxdef(:,5)';el.Uxdef(:,2)';el.Uxdef(:,6)';el.Uxdef(:,3)';el.Uxdef(:,7)';el.Uxdef(:,4)';el.Uxdef(:,8)';el.Uxdef(:,1)'], ...
%                        [el.Uydef(:,1)';el.Uydef(:,5)';el.Uydef(:,2)';el.Uydef(:,6)';el.Uydef(:,3)';el.Uydef(:,7)';el.Uydef(:,4)';el.Uydef(:,8)';el.Uydef(:,1)'], ...
%                         '.:w','MarkerSize',10);
%   end
  
 %Jintegral Elements 
   hJint = fill([el.Ux(Jint.el,1)';el.Ux(Jint.el,2)';el.Ux(Jint.el,3)';el.Ux(Jint.el,4)';el.Ux(Jint.el,1)'], ...
               [el.Uy(Jint.el,1)';el.Uy(Jint.el,2)';el.Uy(Jint.el,3)';el.Uy(Jint.el,4)';el.Uy(Jint.el,1)'], ... 
                'k','EdgeColor','none','FaceAlpha',0.3);
  
 %Create xlabel
  xlabel('Along Specimen Length (mm)');
 %Create ylabel
  ylabel('Along Specimen Width (mm)');
 %Create title
  title('DIC displacment field showing FE elements and J integral elements');
 %Create legend
 legend2 = legend(axesFE,[hdisp hJint(end)], ...
           'FE elements','J integral elements');
 set(legend2,'Orientation','horizontal','Location','SouthOutside','Color',[0.5 0.5 0.5]);
 legend boxoff;
 %Create colorbar
  colorbar('peer',axesFE);   
    
  axis([min(min(el.Uxdef)) max(max(el.Uxdef)) min(min(el.Uydef)) max(max(el.Uydef))]);
  axis square;
warning on
 
elseif mesh.Operation =='Str'
    %%
warning off
 %Figure 
  FE = figure('Position',[50 100 1000 800]);
  axesFE = axes('Parent',FE);
  hold(axesFE,'all');
 
 %plot strain contours
  colormap(jet)
  hdisp = pcolor(gl.Uxdef,gl.Uydef,gl.smises,'Parent',axesFE);
          shading interp
  hcont = contour(gl.Uxdef,gl.Uydef,gl.smises,10,'LineWidth',1,'LineColor',...
      [0 0 0],'Parent',axesFE);
  
%  %FE elements and nodes 
%   switch mat.eltype 
%       case 'Q4'
%          hFEel = plot([el.Uxdef(1:20,1)';el.Uxdef(1:20,2)';el.Uxdef(1:20,3)';el.Uxdef(1:20,4)'], ...
%                       [el.Uydef(1:20,1)';el.Uydef(1:20,2)';el.Uydef(1:20,3)';el.Uydef(1:20,4)'], ...
%                        '.w','MarkerSize',10);
%       case 'Q8'
%           hFEel = plot([el.Uxdef(:,1)';el.Uxdef(:,5)';el.Uxdef(:,2)';el.Uxdef(:,6)';el.Uxdef(:,3)';el.Uxdef(:,7)';el.Uxdef(:,4)';el.Uxdef(:,8)';el.Uxdef(:,1)'], ...
%                        [el.Uydef(:,1)';el.Uydef(:,5)';el.Uydef(:,2)';el.Uydef(:,6)';el.Uydef(:,3)';el.Uydef(:,7)';el.Uydef(:,4)';el.Uydef(:,8)';el.Uydef(:,1)'], ...
%                         '.:W','MarkerSize',10);
%       case 'Q9'
%           plot([el.Uxdef(:,1)';el.Uxdef(:,2)';el.Uxdef(:,3)';el.Uxdef(:,4)';el.Uxdef(:,5)';el.Uxdef(:,6)';el.Uxdef(:,7)';el.Uxdef(:,8)';el.Uxdef(:,9)'], ...
%                [el.Uydef(:,1)';el.Uydef(:,2)';el.Uydef(:,3)';el.Uydef(:,4)';el.Uydef(:,5)';el.Uydef(:,6)';el.Uydef(:,7)';el.Uydef(:,8)';el.Uydef(:,9)'], ...
%                '.w','MarkerSize',10);
%           hFEel = plot([el.Uxdef(:,1)';el.Uxdef(:,5)';el.Uxdef(:,2)';el.Uxdef(:,6)';el.Uxdef(:,3)';el.Uxdef(:,7)';el.Uxdef(:,4)';el.Uxdef(:,8)';el.Uxdef(:,1)'], ...
%                        [el.Uydef(:,1)';el.Uydef(:,5)';el.Uydef(:,2)';el.Uydef(:,6)';el.Uydef(:,3)';el.Uydef(:,7)';el.Uydef(:,4)';el.Uydef(:,8)';el.Uydef(:,1)'], ...
%                         '.:w','MarkerSize',10);
%   end
  
 %Jintegral Elements 
hJint = fill([el.Ux(Jint.el,1)';el.Ux(Jint.el,2)';el.Ux(Jint.el,3)';el.Ux(Jint.el,4)';el.Ux(Jint.el,1)'], ...
        [el.Uy(Jint.el,1)';el.Uy(Jint.el,2)';el.Uy(Jint.el,3)';el.Uy(Jint.el,4)';el.Uy(Jint.el,1)'], ... 
                'k','EdgeColor','none','FaceAlpha',0.3);
  
 %Create xlabel
  xlabel('along specimen length (mm)');
 %Create ylabel
  ylabel('along specimen width (mm)');
 %Create title
  title('DIC displacment field showing FE elements and J integral elements');
 %Create legend
 legend2 = legend(axesFE,[hdisp hJint(end)], ...
           'FE elements','J integral elements');
 set(legend2,'Orientation','horizontal','Location','SouthOutside','Color',[0.5 0.5 0.5]);
 legend boxoff;
 %Create colorbar
  colorbar('peer',axesFE);   
    
  axis([min(min(el.Uxdef)) max(max(el.Uxdef)) min(min(el.Uydef)) max(max(el.Uydef))]);
  axis equal;
warning on
end