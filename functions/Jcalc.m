function [Jint] = Jcalc(el,Jint, mat,mesh)
%% J-integral calculation
[gaussloc,gaussweight] = gaussquad(mat);
Jint.gp=zeros([length(Jint.n) 1]);

for elnum = 1:size(Jint.n,1)
    Ue = zeros(1,size(Jint.n,2)*2);
    de = zeros(1,size(Jint.n,2)*2);
    
    xel = Jint.Ux(elnum,:);
    yel = Jint.Uy(elnum,:);
    dxel = Jint.dx(elnum,:);
    dyel = Jint.dy(elnum,:);
    
    Ue(1:2:size(Jint.n,2)*2) = [xel(4) xel(1) xel(2) xel(3)];
    Ue(2:2:size(Jint.n,2)*2) = [yel(4) yel(1) yel(2) yel(3)];
    de(1:2:size(Jint.n,2)*2) = [dxel(4) dxel(1) dxel(2) dxel(3)];
    de(2:2:size(Jint.n,2)*2) = [dyel(4) dyel(1) dyel(2) dyel(3)];
      
    eta = gaussloc(1);
    xi  = gaussloc(2);
    %forming [sxx sxy; sxy syy] matrix
    sgpm = [Jint.gpsxx(elnum) Jint.gpsxy(elnum); Jint.gpsxy(elnum) Jint.gpsyy(elnum)];
    
    %finding differentiation of (d)ispalcement and (q)crack ext for J-integral calculation
    [H,dHex] = shapefunct(eta, xi, mat);
    [~,J] = BJmat(dHex, Jint, Ue);
    Jint.gpq(elnum) = H*Jint.nQ(elnum,:)';
    dHxy = J\dHex;
    
    %dx is [Exx dUy/dx]    
    dx = (dHxy(1,:)*[de(1:2:end);de(2:2:end)]')';

    if mesh.Operation == 'Str'; dx = [Jint.gpexx(elnum) dx(2)]; end
    
    %dq/dy and dq/dx
    dq = roundn(dHxy*Jint.nQ(elnum,:)',-3);
    dq=flipud(dq);
    Jint.dq(elnum,:)=dq;
    
    %calculate Jint
    lin1 = (sgpm(1,1)*dx(1)+sgpm(1,2)*dx(2)-Jint.gpW(elnum))*-dq(1);
    lin2 = (sgpm(1,2)*dx(1)+sgpm(2,2)*dx(2)-0)*dq(2);
    Jtest = (lin1+lin2)*(det(J)*gaussweight);
    
    if(isnan(Jtest))
        disp 'Error';
        return
    end
    Jint.gp(elnum) = Jtest;
end
Jint.J = (sum(Jint.gp));

if strcmp(mat.stressstate,'plane_stress')
    Jint.K = sqrt(abs(Jint.J)*mat.E);
else
    Jint.K = sqrt(abs(Jint.J)*mat.E)/sqrt((1-mat.nu^2));
end
