function [Jint]=getMask(mesh,gl,el,mat)

Pla_reg = isplastic(el,mat.yield);  %finding plasticity regions
if  mesh.Operation =='Dis'
    H1 = prePlot2(mesh,gl,1,Pla_reg,el);
else
    H1 = prePlot(mesh,gl,0,Pla_reg,el); 
end

if mesh.operation=='Norm'
    title({'Mask the Crack area';''})
[xmask, ymask] = ginput(2);
rectangle('Position',[min(xmask(1),xmask(2)),min(ymask(1),ymask(2)),...
    abs(xmask(2)-xmask(1)),abs(ymask(2)-ymask(1))],'FaceColor','k'); drawnow;
[xmask,ymask] = convCoord(xmask,ymask,mesh.winFE,mesh.UFE);

title({'Select Outer J-Contour';''})
[xo,yo] = ginput(2);
rectangle('Position',[min(xo(1),xo(2)),min(yo(1),yo(2)),abs(xo(2)-xo(1)),...
    abs(yo(2)-yo(1))],'EdgeColor','W'); drawnow;
[xo,yo] = convCoord(xo,yo,mesh.winFE,mesh.UFE);

 title({'Select Inner J-Contour';''})
[xi,yi] = ginput(2);
rectangle('Position',[min(xi(1),xi(2)),min(yi(1),yi(2)),abs(xi(2)-xi(1)),...
    abs(yi(2)-yi(1))],'EdgeColor','k');drawnow;
[xi,yi] = convCoord(xi,yi,mesh.winFE,mesh.UFE);

elseif mesh.operation=='West'
    RangeT=0.05*mesh.maxGrid; %crack tolerance
    xmask=[RangeT; -mesh.maxGrid]; ymask=[RangeT; -RangeT];
        rectangle('Position',[min(xmask(1),xmask(2)),min(ymask(1),ymask(2)),...
            abs(xmask(2)-xmask(1)),abs(ymask(2)-ymask(1))],'FaceColor','k'); drawnow;
        [xmask,ymask] = convCoord(xmask,ymask,mesh.winFE,mesh.UFE);


xo=[mesh.maxGrid; -mesh.maxGrid]; yo=[mesh.maxGrid; -mesh.maxGrid];
    rectangle('Position',[min(xo(1),xo(2)),min(yo(1),yo(2)),abs(xo(2)-xo(1)),...
        abs(yo(2)-yo(1))],'EdgeColor','W'); drawnow;
    [xo,yo] = convCoord(xo,yo,mesh.winFE,mesh.UFE);
    
xi=[2*RangeT; -2*RangeT]; yi=[2*RangeT; -2*RangeT];
    rectangle('Position',[min(xi(1),xi(2)),min(yi(1),yi(2)),abs(xi(2)-xi(1)),...
        abs(yi(2)-yi(1))],'EdgeColor','W');drawnow;
    [xi,yi] = convCoord(xi,yi,mesh.winFE,mesh.UFE);
    
end

Jint.nout = [max(yo) min(yo); min(xo) max(xo)];
Jint.nin  = [max(yi) min(yi); min(xi) max(xi)];

Jint.mask(:,:,1) = [ymask(1) ymask(2);xmask(1) xmask(2)]; 
Jint.mask(Jint.mask<1)=1;

