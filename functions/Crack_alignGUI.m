function [mesh,mat,el,gl] = Crack_alignGUI(mesh,mat,el,gl)

    mesh.winDIC = fliplr(mesh.winDIC);
    mesh.winFE = fliplr(mesh.winFE);
    Tri(1:2,:) = flipud(mesh.UDIC);
    Tri(3:4,:) = flipud(mesh.dDIC);
    
    for ii=1:length(Tri)-1
        for i=1:length(Tri)-1
            if(Tri(2,i)>Tri(2,i+1))
                    buf=Tri(:,i);
                    Tri(:,i)=Tri(:,i+1);
                    Tri(:,i+1)=buf;
            end
        end
    end
    
    mesh.UDIC = Tri(1:2,:);
    mesh.dDIC = Tri(3:4,:);
    
    [el,mesh] = meshDIC(mat,mesh);
    [el] = FEanalysis(el,mat);
    [gl] = makeglobal(el,mesh);
