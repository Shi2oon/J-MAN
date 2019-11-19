function [el,mesh,gl,Jint]=StrainInegration(mat,mesh,el,Jint,gl)
close all
if mesh.Operation =='Str'
    [Jint]    = findJintmask(mat, mesh, el, Jint);
    [el,mesh] = FEintegrateStrains(el,mat,mesh,Jint);
    [el,mesh] = RotRemoval('true',mesh,el);
    [gl]      = makeglobal(el, mesh);
    figure;    surf(gl.dy);    view([0 90]);   colormap(jet);
end