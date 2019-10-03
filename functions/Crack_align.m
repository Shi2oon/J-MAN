function [mesh,mat,el,gl] = Crack_align(mesh,mat,el,gl)

reply = input('Is the crack on the x axis? Y/N [Y]: ', 's');
if isempty(reply)
    reply = 'Y';
end
if(reply=='N')
    mesh.winDIC = fliplr(mesh.winDIC);
    mesh.winFE = fliplr(mesh.winFE);
    Tri(1:2,:) = flipud(mesh.UDIC);
    Tri(3:4,:) = flipud(mesh.dDIC);

    Tri = sortrows(Tri',[1 2])';
    mesh.UDIC = Tri(1:2,:);
    mesh.dDIC = Tri(3:4,:);
%     mesh.strainU = flipud(mesh.strainU);

if mesh.Operation=='Str'
    mesh.straindata(:,1:2) = fliplr(mesh.straindata(:,1:2));
    mesh.straindata(:,2) = -mesh.straindata(:,2);
    mesh.straindata(:,3:4) = fliplr(mesh.straindata(:,3:4));
    mesh.straindata(:,5) = -mesh.straindata(:,5);
end
    
    [el,mesh] = meshDIC(mesh);
    [el] = FEanalysisStrains(el,mat,mesh);
    [gl] = makeglobal(el, mesh);
end