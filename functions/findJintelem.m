function [Jint] = findJintelem(~, mesh, el, Jint)

if mesh.Operation=='Dis'

nodes = reshape(1:length(mesh.UDIC),mesh.winDIC(2),mesh.winDIC(1));

Jint.nout(Jint.nout>max(mesh.winDIC)) = max(mesh.winDIC);

Jintnodes = zeros(size(nodes));
Jintnodes(min(Jint.nout(1,:)):max(Jint.nout(1,:)),min(Jint.nout(2,:)):max(Jint.nout(2,:)))=1;
Jintnodes((min(Jint.nin(1,:))+1):(max(Jint.nin(1,:))-1),(min(Jint.nin(2,:))+1):(max(Jint.nin(2,:))-1)) = 0;
Jintnodes(min(Jint.mask(1,:)):max(Jint.mask(1,:)),min(Jint.mask(2,:)):max(Jint.mask(2,:)))=0;

Jintnodes = Jintnodes.*nodes;
Jintnodes = unique(Jintnodes);
Jintnodes(Jintnodes==0)=[];

elnum =1;
for c = 1:length(el.n)
    if(sum(ismember(el.n(c,:),Jintnodes))==4)
        JDelem(elnum) = c;
        elnum = elnum+1;
    end
end

Jint.el = JDelem';

elseif mesh.Operation=='Str'
nodes = reshape(1:length(mesh.UDIC),mesh.winDIC(2),mesh.winDIC(1));

Jintnodes = zeros(size(nodes));
Jintnodes(min(Jint.nout(1,:)):max(Jint.nout(1,:)),min(Jint.nout(2,:)):max(Jint.nout(2,:)))=1;
Jintnodes((min(Jint.nin(1,:))+1):(max(Jint.nin(1,:))-1),(min(Jint.nin(2,:))+1):(max(Jint.nin(2,:))-1)) = 0;

% Jintnodes(min(Jint.mask(1,:)):max(Jint.mask(1,:)),min(Jint.mask(2,:)):max(Jint.mask(2,:)))=0;
for i=1:size(Jint.mask,3)
    Jintnodes(min(Jint.mask(1,:,i)):max(Jint.mask(1,:,i)),min(Jint.mask(2,:,i)):max(Jint.mask(2,:,i)))=0;
end

Jintnodes = Jintnodes.*nodes;
Jintnodes = unique(Jintnodes);
Jintnodes(Jintnodes==0)=[];

elnum =1;
for c = 1:length(el.n)
    if(sum(ismember(el.n(c,:),Jintnodes))==4)
        JDelem(elnum) = c;
        elnum = elnum+1;
    end
end

Jint.el = JDelem';
end
