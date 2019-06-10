function [Ctrs Nb] = CtrCalc(Jint)


Diff = abs(Jint.nin-Jint.nout);
Nb = max(Diff(:))-1;
fprintf(1,'Calculating J on %3.0f contours\n', Nb);

Ctrs = zeros(Nb,8);

for i=1:Nb
    if i<Diff(1,1)
        Jint.nout(1,1)=Jint.nin(1,1)+i;
    end
    if i<Diff(1,2)
        Jint.nout(1,2)=Jint.nin(1,2)-i;
    end
    if i<Diff(2,1)
        Jint.nout(2,1)=Jint.nin(2,1)-i;
    end
    if i<Diff(2,2)
        Jint.nout(2,2)=Jint.nin(2,2)+i;
    end
    Ctrs(i,:)= [Jint.nin(1,:) Jint.nin(2,:) Jint.nout(1,:) Jint.nout(2,:)];
end
