function [Results,Jint]=J_IntegralCacl(Jint,mat,mesh,el)
    [Ctrs,Nb_ct] = CtrCalc(Jint);
    Results = zeros(Nb_ct,2);
    fprintf('\nJ-Integraling...\n');
for i=1:Nb_ct
        Jint.nin(1,:)=Ctrs(i,1:2);
        Jint.nin(2,:)=Ctrs(i,3:4);
        Jint.nout(1,:)=Ctrs(i,5:6);
        Jint.nout(2,:)=Ctrs(i,7:8);
        %finds elements within defined (Jint)ergral area, Jint.nout/nin
        [Jint] = findJintelem(mat, mesh, el, Jint); 
        %creates smooth virtual crack extension function q over Jint.nout/nin
        [Jint] = getq(el,Jint, mesh);       
        [Jint] = Jcalc(el,Jint, mat,mesh);
        Results(i,:)=[Jint.J Jint.K];
        if(i>Nb_ct)
            clear Jint;
        end
end
