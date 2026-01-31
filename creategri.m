
function creategri(refelem,nnodes,xnode,ynode,nodecon,nbfgrp,LBF,BCs,Bcxelem,Bcyelem,filename)
fid = fopen(filename,'w');
fprintf(fid,'%d %d 2\n',[nnodes refelem]);
for ii = 1:nnodes
    fprintf(fid,'%d %d\n',[xnode(ii) ynode(ii)]);
end

fprintf(fid,'%d\n',nbfgrp); db = 1; 
for ii = 1:nbfgrp
    rr = 1;
    fprintf(fid,'%d 2 %s\n',LBF(ii),BCs(ii));
    for jj = db:db+LBF(ii)-1
        fl = [0 0];
        for kk = 1:nnodes
            if Bcxelem(jj,1)==xnode(kk)&&Bcyelem(jj,1)==ynode(kk)
                conBC1 = kk;
                fl(1) = 1;
            end
            if Bcxelem(jj,2)==xnode(kk)&&Bcyelem(jj,2)==ynode(kk)
                conBC2 = kk;
                fl(2) = 1;
            end
            if fl(1) == 1 && fl(2) == 1
                break
            end
        end
        fprintf(fid,'%d %d\n',[conBC1 conBC2]);
        if ii == 1
            p1(rr,1) = conBC1;
            p1(rr,2) = conBC2;
            rr = rr+1;
        end
        if ii == 3
            p3(rr,1) = conBC1;
            p3(rr,2) = conBC2;
            rr = rr+1;
        end
        if ii == 5
            p5(rr,1) = conBC1;
            p5(rr,2) = conBC2;
            rr = rr+1;
        end
        if ii == 7
            p7(rr,1) = conBC1;
            p7(rr,2) = conBC2;
            rr = rr+1;
        end
    end
    db = db+LBF(ii);
end
tt = 1;
for jj = 1:length(p1)
    if jj == 1
        p1N(tt) = p1(jj,1);
        p1N(tt+1) = p1(jj,2);
        tt = tt+2;
    else
        f = 0;
        for kk = 1:length(p1N)
            if p1(jj,1) == p1N(kk)
                f = 1;
                break
            end
        end
        if f == 0
            p1N(tt) = p1(jj,1);
            tt = tt+1;
        end
        f = 0;
        for kk = 1:length(p1N)
            if p1(jj,2) == p1N(kk)
                f = 1;
                break
            end
        end
        if f == 0
            p1N(tt) = p1(jj,2);
            tt = tt+1;
        end
    end
end

tt = 1;
for jj = 1:length(p3)
    if jj == 1
        p3N(tt) = p3(jj,1);
        p3N(tt+1) = p3(jj,2);
        tt = tt+2;
    else
        f = 0;
        for kk = 1:length(p3N)
            if p3(jj,1) == p3N(kk)
                f = 1;
                break
            end
        end
        if f == 0
            p3N(tt) = p3(jj,1);
            tt = tt+1;
        end
        f = 0;
        for kk = 1:length(p3N)
            if p3(jj,2) == p3N(kk)
                f = 1;
                break
            end
        end
        if f == 0
            p3N(tt) = p3(jj,2);
            tt = tt+1;
        end
    end
end

tt = 1;
for jj = 1:length(p5)
    if jj == 1
        p5N(tt) = p5(jj,1);
        p5N(tt+1) = p5(jj,2);
        tt = tt+2;
    else
        f = 0;
        for kk = 1:length(p5N)
            if p5(jj,1) == p5N(kk)
                f = 1;
                break
            end
        end
        if f == 0
            p5N(tt) = p5(jj,1);
            tt = tt+1;
        end
        f = 0;
        for kk = 1:length(p5N)
            if p5(jj,2) == p5N(kk)
                f = 1;
                break
            end
        end
        if f == 0
            p5N(tt) = p5(jj,2);
            tt = tt+1;
        end
    end
end

tt = 1;
for jj = 1:length(p7)
    if jj == 1
        p7N(tt) = p7(jj,1);
        p7N(tt+1) = p7(jj,2);
        tt = tt+2;
    else
        f = 0;
        for kk = 1:length(p7N)
            if p7(jj,1) == p7N(kk)
                f = 1;
                break
            end
        end
        if f == 0
            p7N(tt) = p7(jj,1);
            tt = tt+1;
        end
        f = 0;
        for kk = 1:length(p7N)
            if p7(jj,2) == p7N(kk)
                f = 1;
                break
            end
        end
        if f == 0
            p7N(tt) = p7(jj,2);
            tt = tt+1;
        end
    end
end

fprintf(fid,'%d 1 TriLagrange\n',refelem);

for jj = 1:refelem
    fprintf(fid,'%d %d %d\n',nodecon(jj,:));
end

fprintf(fid,'2 PeriodicGroup\n');


fprintf(fid,'%d Translational\n',LBF(1)+1);

for ii = 1:length(p1N)
    for jj = 1:length(p7N)
        if xnode(p1N(ii))==xnode(p7N(jj))
            fprintf(fid,'%d %d\n',[p1N(ii) p7N(jj)]);
            break
        end
    end
end
fprintf(fid,'%d Translational\n',LBF(3)+1);
for ii = 1:length(p3N)
    for jj = 1:length(p5N)
        if xnode(p3N(ii))==xnode(p5N(jj))
            fprintf(fid,'%d %d\n',[p3N(ii) p5N(jj)]);
            break
        end
    end
end

fclose(fid);
