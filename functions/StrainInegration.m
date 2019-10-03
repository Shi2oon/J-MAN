function [el,mesh,gl,Jint]=StrainInegration(mat,mesh,el,Jint,gl)
if mesh.Operation =='Str'
    [Jint]    = findJintmask(mat, mesh, el, Jint);
    [el,mesh] = FEintegrateStrains(el,mat,mesh,Jint);
    [el,mesh] = RotRemoval('true',mesh,el);
    [gl]      = makeglobal(el, mesh);
    H2        =  figure;    surf(gl.dy);    view([3 42]);   colormap(jet);
end