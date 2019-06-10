function plotFEe(mesh,el,gl)
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
% hcont = contour(gl.Uxdef,gl.Uydef,gl.dy,10,'LineWidth',1,'LineColor',...
%     [0 0 0],'Parent',axesFE);
%    hdisp = scatter(gl.Uxdef(:),gl.Uydef(:),3,gl.dy(:));
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

  
 %Create xlabel
 %Create xlabel
  xlabel('Along Specimen Length (mm)');
 %Create ylabel
  ylabel('Along Specimen Width (mm)');
 %Create title
  title('Displacement Field (U_y) showing FE elements');
 %Create colorbar
  colorbar('peer',axesFE);   
    
  axis([min(min(el.Uxdef)) max(max(el.Uxdef)) min(min(el.Uydef)) max(max(el.Uydef))]);
warning on
elseif mesh.Operation =='Str'
    warning off
 %Figure 
  FE = figure('Position',[50 100 1000 800]);
  axesFE = axes('Parent',FE);
  hold(axesFE,'all');
 
 %plot strain contours
  colormap(jet)
%   imagesc(gl.dy)
  hdisp = pcolor(gl.Uxdef*1000,gl.Uydef*1000,gl.dy,'Parent',axesFE);
          shading interp
  hcont = contour(gl.Uxdef*1000,gl.Uydef*1000,gl.dy,10,'LineWidth',1,'LineColor',...
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

  
 %Create xlabel
  xlabel('Along Specimen Length (mm)');
 %Create ylabel
  ylabel('Along Specimen Width (mm)');
 %Create title
  title('Displacement Field (U_y) showing FE elements');
 %Create colorbar
  colorbar('peer',axesFE);   
    
  axis([min(min(el.Uxdef*1000)) max(max(el.Uxdef*1000)) ...
      min(min(el.Uydef*1000)) max(max(el.Uydef*1000))]);
warning on
end
 