
function [nelem,xelem,yelem,nedges,xedge,yedge,ledges,elemedge,BCs,Bcxelem,Bcyelem,LBF,PeriG,Vertsx,Vertsy] = readgri(grifile)

% Read mesh (gri file)

fid = fopen(grifile, 'r');
% Read in nodes
A = fscanf(fid,'%d', 3);
nnode    = A(1);
nelem = A(2);
dim      = A(3);
Verts = zeros(nnode, dim);
for inode = 1:nnode
  A = fscanf(fid, '%lf', 2);
  Verts(inode,:) = A(1:2)';
end
Vertsx = Verts(:,1);
Vertsy = Verts(:,2);
% Read through boundary info
A = fscanf(fid, '%d', 1);
nbfgrp = A(1);
alines = readlines(grifile);
BC1line = contains(alines,'BGroup1');
BC1linenum = find(BC1line);

BC2line = contains(alines,'BGroup2');
BC2linenum = find(BC2line);

BC3line = contains(alines,'BGroup3');
BC3linenum = find(BC3line);

BC4line = contains(alines,'BGroup4');
BC4linenum = find(BC4line);

BC5line = contains(alines,'BGroup5');
BC5linenum = find(BC5line);

BC6line = contains(alines,'BGroup6');
BC6linenum = find(BC6line);

BC7line = contains(alines,'BGroup7');
BC7linenum = find(BC7line);

BC8line = contains(alines,'BGroup8');
BC8linenum = find(BC8line);

BCs = ["BGroup1" "BGroup2" "BGroup3" "BGroup4" "BGroup5" "BGroup6" "BGroup7" "BGroup8"];
PeriG = [1 0 2 0 2 0 1 0];
BClinenum = [BC1linenum BC2linenum BC3linenum BC4linenum BC5linenum BC6linenum BC7linenum BC8linenum];
gg = 1;
for bc = 1:nbfgrp
    BCL = sscanf(alines(BClinenum(bc)),'%d');
    LBF(bc) = BCL(1);
    for ee = 1:BCL(1)
        bcnodes = str2num(alines(BClinenum(bc)+ee));
        Bcxelem(gg,1) = Verts(bcnodes(1),1);
        Bcyelem(gg,1) = Verts(bcnodes(1),2);
        Bcxelem(gg,2) = Verts(bcnodes(2),1);
        Bcyelem(gg,2) = Verts(bcnodes(2),2);
        gg = gg+1;
    end
end

lineelem = contains(alines,'TriLagrange');
linenum = find(lineelem);
xelem = zeros(nelem,3);
yelem = zeros(nelem,3);
kk = 1; nedges = 0;
ledges = zeros(nelem,1); edgeconnect = zeros(nelem,2); elemedge = zeros(nelem,3);
for ii = 1:nelem
    elemc = str2num(alines(linenum+ii));
    xelem(ii,1) = Verts(elemc(1),1);
    yelem(ii,1) = Verts(elemc(1),2);
    xelem(ii,2) = Verts(elemc(2),1);
    yelem(ii,2) = Verts(elemc(2),2);
    xelem(ii,3) = Verts(elemc(3),1);
    yelem(ii,3) = Verts(elemc(3),2);
    if ii == 1
        nedges = nedges+3;
        ledges(kk) = sqrt((xelem(ii,1)-xelem(ii,2))^2+(yelem(ii,1)-yelem(ii,2))^2);
        ledges(kk+1) = sqrt((xelem(ii,3)-xelem(ii,2))^2+(yelem(ii,3)-yelem(ii,2))^2);
        ledges(kk+2) = sqrt((xelem(ii,1)-xelem(ii,3))^2+(yelem(ii,1)-yelem(ii,3))^2);
        edgeconnect(kk,1) = elemc(1);
        edgeconnect(kk,2) = elemc(2);
        edgeconnect(kk+1,1) = elemc(2);
        edgeconnect(kk+1,2) = elemc(3);
        edgeconnect(kk+2,1) = elemc(3);
        edgeconnect(kk+2,2) = elemc(1);
        xedge(kk,1) = Verts(elemc(1),1);
        xedge(kk,2) = Verts(elemc(2),1);
        yedge(kk,1) = Verts(elemc(1),2);
        yedge(kk,2) = Verts(elemc(2),2);
        xedge(kk+1,1) = Verts(elemc(2),1);
        xedge(kk+1,2) = Verts(elemc(3),1);
        yedge(kk+1,1) = Verts(elemc(2),2);
        yedge(kk+1,2) = Verts(elemc(3),2);
        xedge(kk+2,1) = Verts(elemc(3),1);
        xedge(kk+2,2) = Verts(elemc(1),1);
        yedge(kk+2,1) = Verts(elemc(3),2);
        yedge(kk+2,2) = Verts(elemc(1),2);
        kk = kk+3;
        elemedge(ii,1) = 1;
        elemedge(ii,2) = 2;
        elemedge(ii,3) = 3;
    else
        nedges1 = nedges;
        edgel(1,1) = elemc(1);
        edgel(1,2) = elemc(2);
        edgel(2,1) = elemc(2);
        edgel(2,2) = elemc(3);
        edgel(3,1) = elemc(3);
        edgel(3,2) = elemc(1);
        for aa = 1:3
            fl = 0;
            for jj = 1:nedges1
                if (edgel(aa,1)==edgeconnect(jj,1)&&edgel(aa,2)==edgeconnect(jj,2)) || (edgel(aa,2)==edgeconnect(jj,1)&&edgel(aa,1)==edgeconnect(jj,2))
                    fl = 1;
                    elemedge(ii,aa) = jj;
                    break
                end
            end
            if fl == 0
                nedges = nedges+1;
                edgeconnect(kk,1) = edgel(aa,1);
                edgeconnect(kk,2) = edgel(aa,2);
                xedge(kk,1) = Verts(edgel(aa,1),1);
                xedge(kk,2) = Verts(edgel(aa,2),1);
                yedge(kk,1) = Verts(edgel(aa,1),2);
                yedge(kk,2) = Verts(edgel(aa,2),2);
                ledges(kk) = sqrt((Verts(edgel(aa,1),1)-Verts(edgel(aa,2),1))^2+(Verts(edgel(aa,1),2)-Verts(edgel(aa,2),2))^2);
                elemedge(ii,aa) = kk;
                kk=kk+1;
            end
        end
         
    end

end

end
