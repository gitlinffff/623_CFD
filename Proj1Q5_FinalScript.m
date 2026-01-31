rl = readmatrix('h.dat');
hdat = rl(:,4);
xch = rl(:,2);
ych = rl(:,3);
InterpFunc = scatteredInterpolant(xch,ych,hdat,'linear','nearest');

dat = readmatrix('blade.txt');
xd = dat(:,1);
yd = dat(:,2);
xsplinel = [xd(2:71);];
ysplinel = [yd(2:71);];
ysplinel = ysplinel-18;
xsplineu = [xd(73:139)];
ysplineu = [yd(73:139)];
sp = 1;
for ii = 2:70
    xsl(sp:sp+9) = 0.0001+linspace(xsplinel(ii-1),xsplinel(ii),10);
    ysl(sp:sp+9) = 0.0001+linspace(ysplinel(ii-1),ysplinel(ii),10);
    sp = sp+10;
end
sp = 1;
for ii = 2:67
    xsu(sp:sp+9) = 0.0001+linspace(xsplineu(ii-1),xsplineu(ii),10);
    ysu(sp:sp+9) = 0.00001+linspace(ysplineu(ii-1),ysplineu(ii),10);
    sp = sp+10;
end

%[nelem,xelem,yelem,nedges,xedge,yedge,ledges,elemedge,BCs,Bcxelem,Bcyelem,LBF,PeriG,xinodes,yinodes] = readgri('initial_mesh.gri');
connectv = [1 2;2 3;3 1];
%BCf = zeros(nedges,1); LBRF=LBF; 
w = 0.2; 

for qq = 1:3
    if qq == 1
        filename = 'initial_mesh1.gri';
   
    end
    [nelem,xelem,yelem,nedges,xedge,yedge,ledges,elemedge,BCs,Bcxelem,Bcyelem,LBF,PeriG,xinodes,yinodes] = readgri(filename);
    BCf = zeros(nedges,1); LBRF=LBF;
    fedge = false(nedges,1);
    rt = 1;
    for ii = 1:nedges
        
        db = 1;
        for tv = 1:length(LBF)
            for bc = db:db+LBF(tv)-1
                if (xedge(ii,1)==Bcxelem(bc,1)&&yedge(ii,1)==Bcyelem(bc,1)&&xedge(ii,2)==Bcxelem(bc,2)&&yedge(ii,2)==Bcyelem(bc,2))||(xedge(ii,1)==Bcxelem(bc,2)&&yedge(ii,1)==Bcyelem(bc,2)&&xedge(ii,2)==Bcxelem(bc,1)&&yedge(ii,2)==Bcyelem(bc,1))
                    BCf(ii) = tv;
                    break
                end
            end
            if BCf(ii) == tv
                break
            end
            db = db+LBF(tv);
        end
        xmed = (xedge(ii,1)+xedge(ii,2))/2;
        ymed = (yedge(ii,1)+yedge(ii,2))/2;
        h = InterpFunc(xmed,ymed);
        if ledges(ii)>h
            fedge(ii) = true;
            if BCf(ii)>0
                LBRF(BCf(ii)) = LBRF(BCf(ii))+1;
                Bcrxelem(rt,1) = xedge(ii,1);
                Bcrxelem(rt,2) = xmed;
                Bcryelem(rt,1) = yedge(ii,1);
                Bcryelem(rt,2) = ymed;
                Bcrxelem(rt+1,2) = xedge(ii,2);
                Bcrxelem(rt+1,1) = xmed;
                Bcryelem(rt+1,2) = yedge(ii,2);
                Bcryelem(rt+1,1) = ymed;
                BCcon(rt) = BCf(ii);
                BCcon(rt+1) = BCf(ii);
                rt = rt+2;
            end
        else
            if BCf(ii)>0
                Bcrxelem(rt,1) = xedge(ii,1);
                Bcrxelem(rt,2) = xedge(ii,2);
                Bcryelem(rt,1) = yedge(ii,1);
                Bcryelem(rt,2) = yedge(ii,2);
                BCcon(rt) = BCf(ii);
                rt = rt+1;
            end
        end
    end
    %BCorderedx=Bcrxelem; BCorderedy=Bcryelem; 
    dd = 1; ll = 1;
    for ii = 1:length(LBF)
        for jj = 1:length(BCcon)
            if BCcon(jj) == ii
                BCorderedx(ll,1) = Bcrxelem(jj,1);
                BCorderedx(ll,2) = Bcrxelem(jj,2);
                BCorderedy(ll,1) = Bcryelem(jj,1);
                BCorderedy(ll,2) = Bcryelem(jj,2);
                BCOrdCon(ll) = ii;
                ll = ll+1;
            end
        end
    end
    refelem = nelem;
    kk = 1;
    xrelem = zeros(nelem,3);
    yrelem = zeros(nelem,3);
    ee = 1; 
    for ii = 1:nelem
        if fedge(elemedge(ii,1)) == true && fedge(elemedge(ii,2)) == true && fedge(elemedge(ii,3)) == true
            refelem = refelem+3;
            xrelem(kk,1) = xelem(ii,1);
            xrelem(kk,2) = (xelem(ii,1)+xelem(ii,2))/2;
            xrelem(kk,3) = (xelem(ii,1)+xelem(ii,3))/2;
            xrelem(kk+1,1) = (xelem(ii,1)+xelem(ii,2))/2;
            xrelem(kk+1,2) = xelem(ii,2);
            xrelem(kk+1,3) = (xelem(ii,2)+xelem(ii,3))/2;
            xrelem(kk+2,1) = (xelem(ii,1)+xelem(ii,3))/2;
            xrelem(kk+2,2) = (xelem(ii,2)+xelem(ii,3))/2;
            xrelem(kk+2,3) = xelem(ii,3);
            xrelem(kk+3,1) = (xelem(ii,1)+xelem(ii,3))/2;
            xrelem(kk+3,2) = (xelem(ii,2)+xelem(ii,3))/2;
            xrelem(kk+3,3) = (xelem(ii,1)+xelem(ii,2))/2;
            yrelem(kk,1) = yelem(ii,1);
            yrelem(kk,2) = (yelem(ii,1)+yelem(ii,2))/2;
            yrelem(kk,3) = (yelem(ii,1)+yelem(ii,3))/2;
            yrelem(kk+1,1) = (yelem(ii,1)+yelem(ii,2))/2;
            yrelem(kk+1,2) = yelem(ii,2);
            yrelem(kk+1,3) = (yelem(ii,2)+yelem(ii,3))/2;
            yrelem(kk+2,1) = (yelem(ii,1)+yelem(ii,3))/2;
            yrelem(kk+2,2) = (yelem(ii,2)+yelem(ii,3))/2;
            yrelem(kk+2,3) = yelem(ii,3);
            yrelem(kk+3,1) = (yelem(ii,1)+yelem(ii,3))/2;
            yrelem(kk+3,2) = (yelem(ii,2)+yelem(ii,3))/2;
            yrelem(kk+3,3) = (yelem(ii,1)+yelem(ii,2))/2;
            kk = kk+4;
            % ee = ee+8;
        elseif fedge(elemedge(ii,1)) == true && fedge(elemedge(ii,2)) == false && fedge(elemedge(ii,3)) == false
            refelem = refelem+1;
            xrelem(kk,1) = xelem(ii,1);
            xrelem(kk,2) = (xelem(ii,1)+xelem(ii,2))/2;
            xrelem(kk,3) = xelem(ii,3);
            xrelem(kk+1,1) = (xelem(ii,1)+xelem(ii,2))/2;
            xrelem(kk+1,2) = xelem(ii,2);
            xrelem(kk+1,3) = xelem(ii,3);
            yrelem(kk,1) = yelem(ii,1);
            yrelem(kk,2) = (yelem(ii,1)+yelem(ii,2))/2;
            yrelem(kk,3) = yelem(ii,3);
            yrelem(kk+1,1) = (yelem(ii,1)+yelem(ii,2))/2;
            yrelem(kk+1,2) = yelem(ii,2);
            yrelem(kk+1,3) = yelem(ii,3);
            kk = kk+2;
        elseif fedge(elemedge(ii,1)) == false && fedge(elemedge(ii,2)) == true && fedge(elemedge(ii,3)) == false
            refelem = refelem+1;
            xrelem(kk,1) = xelem(ii,1);
            xrelem(kk,2) = xelem(ii,2);
            xrelem(kk,3) = (xelem(ii,2)+xelem(ii,3))/2;
            xrelem(kk+1,1) = xelem(ii,1);
            xrelem(kk+1,2) = xelem(ii,3);
            xrelem(kk+1,3) = (xelem(ii,2)+xelem(ii,3))/2;
            yrelem(kk,1) = yelem(ii,1);
            yrelem(kk,2) = yelem(ii,2);
            yrelem(kk,3) = (yelem(ii,2)+yelem(ii,3))/2;
            yrelem(kk+1,1) = yelem(ii,1);
            yrelem(kk+1,2) = yelem(ii,3);
            yrelem(kk+1,3) = (yelem(ii,2)+yelem(ii,3))/2;
            kk = kk+2;
        elseif fedge(elemedge(ii,1)) == false && fedge(elemedge(ii,2)) == false && fedge(elemedge(ii,3)) == true
            refelem = refelem+1;
            xrelem(kk,1) = xelem(ii,1);
            xrelem(kk,2) = xelem(ii,2);
            xrelem(kk,3) = (xelem(ii,1)+xelem(ii,3))/2;
            xrelem(kk+1,1) = xelem(ii,2);
            xrelem(kk+1,2) = xelem(ii,3);
            xrelem(kk+1,3) = (xelem(ii,1)+xelem(ii,3))/2;
            yrelem(kk,1) = yelem(ii,1);
            yrelem(kk,2) = yelem(ii,2);
            yrelem(kk,3) = (yelem(ii,1)+yelem(ii,3))/2;
            yrelem(kk+1,1) = yelem(ii,2);
            yrelem(kk+1,2) = yelem(ii,3);
            yrelem(kk+1,3) = (yelem(ii,1)+yelem(ii,3))/2;
            kk = kk+2;
        elseif fedge(elemedge(ii,1)) == false && fedge(elemedge(ii,2)) == true && fedge(elemedge(ii,3)) == true
            refelem = refelem+2;
            a = sqrt((xelem(ii,1)-xelem(ii,3))^2+(yelem(ii,1)-yelem(ii,3))^2);
            b = sqrt((xelem(ii,2)-xelem(ii,3))^2+(yelem(ii,2)-yelem(ii,3))^2);
            c = sqrt((xelem(ii,1)-xelem(ii,2))^2+(yelem(ii,1)-yelem(ii,2))^2);
            t = (a^2-b^2+c^2)/(2*c);
            hei = sqrt(a^2-t^2);
            an1 = 180*asin(hei/a);
            an2 = 180*asin(hei/b);
            if an2>=an1
                xrelem(kk,1) = xelem(ii,1);
                xrelem(kk,2) = xelem(ii,2);
                xrelem(kk,3) = (xelem(ii,1)+xelem(ii,3))/2;
                xrelem(kk+1,1) = xelem(ii,2);
                xrelem(kk+1,2) = (xelem(ii,1)+xelem(ii,3))/2;
                xrelem(kk+1,3) = (xelem(ii,2)+xelem(ii,3))/2;
                xrelem(kk+2,1) = (xelem(ii,1)+xelem(ii,3))/2;
                xrelem(kk+2,2) = (xelem(ii,2)+xelem(ii,3))/2;
                xrelem(kk+2,3) = xelem(ii,3);
                yrelem(kk,1) = yelem(ii,1);
                yrelem(kk,2) = yelem(ii,2);
                yrelem(kk,3) = (yelem(ii,1)+yelem(ii,3))/2;
                yrelem(kk+1,1) = yelem(ii,2);
                yrelem(kk+1,2) = (yelem(ii,1)+yelem(ii,3))/2;
                yrelem(kk+1,3) = (yelem(ii,2)+yelem(ii,3))/2;
                yrelem(kk+2,1) = (yelem(ii,1)+yelem(ii,3))/2;
                yrelem(kk+2,2) = (yelem(ii,2)+yelem(ii,3))/2;
                yrelem(kk+2,3) = yelem(ii,3);
            else
                xrelem(kk,1) = xelem(ii,1);
                xrelem(kk,2) = xelem(ii,2);
                xrelem(kk,3) = (xelem(ii,2)+xelem(ii,3))/2;
                xrelem(kk+1,1) = xelem(ii,1);
                xrelem(kk+1,2) = (xelem(ii,2)+xelem(ii,3))/2;
                xrelem(kk+1,3) = (xelem(ii,1)+xelem(ii,3))/2;
                xrelem(kk+2,1) = (xelem(ii,1)+xelem(ii,3))/2;
                xrelem(kk+2,2) = (xelem(ii,2)+xelem(ii,3))/2;
                xrelem(kk+2,3) = xelem(ii,3);
                yrelem(kk,1) = yelem(ii,1);
                yrelem(kk,2) = yelem(ii,2);
                yrelem(kk,3) = (yelem(ii,2)+yelem(ii,3))/2;
                yrelem(kk+1,1) = yelem(ii,1);
                yrelem(kk+1,2) = (yelem(ii,2)+yelem(ii,3))/2;
                yrelem(kk+1,3) = (yelem(ii,1)+yelem(ii,3))/2;
                yrelem(kk+2,1) = (yelem(ii,1)+yelem(ii,3))/2;
                yrelem(kk+2,2) = (yelem(ii,2)+yelem(ii,3))/2;
                yrelem(kk+2,3) = yelem(ii,3);
            end
            kk = kk+3;
        elseif fedge(elemedge(ii,1)) == true && fedge(elemedge(ii,2)) == false && fedge(elemedge(ii,3)) == true
            refelem = refelem+2;
            b = sqrt((xelem(ii,1)-xelem(ii,3))^2+(yelem(ii,1)-yelem(ii,3))^2);
            c = sqrt((xelem(ii,2)-xelem(ii,3))^2+(yelem(ii,2)-yelem(ii,3))^2);
            a = sqrt((xelem(ii,1)-xelem(ii,2))^2+(yelem(ii,1)-yelem(ii,2))^2);
            t = (a^2-b^2+c^2)/(2*c);
            hei = sqrt(a^2-t^2);
            an1 = 180*asin(hei/a);
            an2 = 180*asin(hei/b);
            if an2>=an1
                xrelem(kk,1) = (xelem(ii,1)+xelem(ii,2))/2;
                xrelem(kk,2) = xelem(ii,2);
                xrelem(kk,3) = xelem(ii,3);
                xrelem(kk+1,1) = (xelem(ii,1)+xelem(ii,2))/2;
                xrelem(kk+1,2) = xelem(ii,3);
                xrelem(kk+1,3) = (xelem(ii,1)+xelem(ii,3))/2;
                xrelem(kk+2,1) = xelem(ii,1);
                xrelem(kk+2,2) = (xelem(ii,2)+xelem(ii,1))/2;
                xrelem(kk+2,3) = (xelem(ii,1)+xelem(ii,3))/2;
                yrelem(kk,1) = (yelem(ii,1)+yelem(ii,2))/2;
                yrelem(kk,2) = yelem(ii,2);
                yrelem(kk,3) = yelem(ii,3);
                yrelem(kk+1,1) = (yelem(ii,1)+yelem(ii,2))/2;
                yrelem(kk+1,2) = yelem(ii,3);
                yrelem(kk+1,3) = (yelem(ii,1)+yelem(ii,3))/2;
                yrelem(kk+2,1) = yelem(ii,1);
                yrelem(kk+2,2) = (yelem(ii,2)+yelem(ii,1))/2;
                yrelem(kk+2,3) = (yelem(ii,1)+yelem(ii,3))/2;
            else
                xrelem(kk,1) = xelem(ii,2);
                xrelem(kk,2) = xelem(ii,3);
                xrelem(kk,3) = (xelem(ii,1)+xelem(ii,3))/2;
                xrelem(kk+1,1) = xelem(ii,2);
                xrelem(kk+1,2) = (xelem(ii,1)+xelem(ii,3))/2;
                xrelem(kk+1,3) = (xelem(ii,1)+xelem(ii,2))/2;
                xrelem(kk+2,1) = xelem(ii,1);
                xrelem(kk+2,2) = (xelem(ii,2)+xelem(ii,1))/2;
                xrelem(kk+2,3) = (xelem(ii,1)+xelem(ii,3))/2;
                yrelem(kk,1) = yelem(ii,2);
                yrelem(kk,2) = yelem(ii,3);
                yrelem(kk,3) = (yelem(ii,1)+yelem(ii,3))/2;
                yrelem(kk+1,1) = yelem(ii,2);
                yrelem(kk+1,2) = (yelem(ii,1)+yelem(ii,3))/2;
                yrelem(kk+1,3) = (yelem(ii,1)+yelem(ii,2))/2;
                yrelem(kk+2,1) = yelem(ii,1);
                yrelem(kk+2,2) = (yelem(ii,2)+yelem(ii,1))/2;
                yrelem(kk+2,3) = (yelem(ii,1)+yelem(ii,3))/2;
            end
            kk = kk+3;
        elseif fedge(elemedge(ii,1)) == true && fedge(elemedge(ii,2)) == true && fedge(elemedge(ii,3)) == false
            refelem = refelem+2;
            c = sqrt((xelem(ii,1)-xelem(ii,3))^2+(yelem(ii,1)-yelem(ii,3))^2);
            b = sqrt((xelem(ii,2)-xelem(ii,3))^2+(yelem(ii,2)-yelem(ii,3))^2);
            a = sqrt((xelem(ii,1)-xelem(ii,2))^2+(yelem(ii,1)-yelem(ii,2))^2);
            t = (a^2-b^2+c^2)/(2*c);
            hei = sqrt(a^2-t^2);
            an1 = 180*asin(hei/a);
            an2 = 180*asin(hei/b);
            if an2>=an1
                xrelem(kk,1) = (xelem(ii,1)+xelem(ii,2))/2;
                xrelem(kk,2) = xelem(ii,1);
                xrelem(kk,3) = xelem(ii,3);
                xrelem(kk+1,1) = (xelem(ii,1)+xelem(ii,2))/2;
                xrelem(kk+1,2) = xelem(ii,3);
                xrelem(kk+1,3) = (xelem(ii,2)+xelem(ii,3))/2;
                xrelem(kk+2,1) = xelem(ii,2);
                xrelem(kk+2,2) = (xelem(ii,2)+xelem(ii,1))/2;
                xrelem(kk+2,3) = (xelem(ii,2)+xelem(ii,3))/2;
                yrelem(kk,1) = (yelem(ii,1)+yelem(ii,2))/2;
                yrelem(kk,2) = yelem(ii,1);
                yrelem(kk,3) = yelem(ii,3);
                yrelem(kk+1,1) = (yelem(ii,1)+yelem(ii,2))/2;
                yrelem(kk+1,2) = yelem(ii,3);
                yrelem(kk+1,3) = (yelem(ii,2)+yelem(ii,3))/2;
                yrelem(kk+2,1) = yelem(ii,2);
                yrelem(kk+2,2) = (yelem(ii,2)+yelem(ii,1))/2;
                yrelem(kk+2,3) = (yelem(ii,2)+yelem(ii,3))/2;
            else
                xrelem(kk,1) = xelem(ii,1);
                xrelem(kk,2) = xelem(ii,3);
                xrelem(kk,3) = (xelem(ii,2)+xelem(ii,3))/2;
                xrelem(kk+1,1) = xelem(ii,1);
                xrelem(kk+1,2) = (xelem(ii,2)+xelem(ii,3))/2;
                xrelem(kk+1,3) = (xelem(ii,1)+xelem(ii,2))/2;
                xrelem(kk+2,1) = xelem(ii,2);
                xrelem(kk+2,2) = (xelem(ii,2)+xelem(ii,1))/2;
                xrelem(kk+2,3) = (xelem(ii,2)+xelem(ii,3))/2;
                yrelem(kk,1) = yelem(ii,1);
                yrelem(kk,2) = yelem(ii,3);
                yrelem(kk,3) = (yelem(ii,2)+yelem(ii,3))/2;
                yrelem(kk+1,1) = yelem(ii,1);
                yrelem(kk+1,2) = (yelem(ii,2)+yelem(ii,3))/2;
                yrelem(kk+1,3) = (yelem(ii,1)+yelem(ii,2))/2;
                yrelem(kk+2,1) = yelem(ii,2);
                yrelem(kk+2,2) = (yelem(ii,2)+yelem(ii,1))/2;
                yrelem(kk+2,3) = (yelem(ii,2)+yelem(ii,3))/2;
            end
            kk = kk+3;
        else
            xrelem(kk,1) = xelem(ii,1);
            xrelem(kk,2) = xelem(ii,2);
            xrelem(kk,3) = xelem(ii,3);
            yrelem(kk,1) = yelem(ii,1);
            yrelem(kk,2) = yelem(ii,2);
            yrelem(kk,3) = yelem(ii,3);
            kk = kk+1;
        end
    end
    nodenum = zeros(refelem,1);
    xnode = zeros(refelem,1);
    ynode = zeros(refelem,1);
    cc = 1;
    nnodes = 0;
    nodecon = zeros(refelem,3);
    for bb = 1:refelem
        if bb == 1
            nnodes = nnodes+3;
            nodenum(cc) = 1;
            nodenum(cc+1) = 2;
            nodenum(cc+2) = 3;
            nodecon(bb,1) = nodenum(cc);
            nodecon(bb,2) = nodenum(cc+1);
            nodecon(bb,3) = nodenum(cc+2);
            xnode(cc) = xrelem(bb,1);
            ynode(cc) = yrelem(bb,1);
            xnode(cc+1) = xrelem(bb,2);
            ynode(cc+1) = yrelem(bb,2);
            xnode(cc+2) = xrelem(bb,3);
            ynode(cc+2) = yrelem(bb,3);
            cc = cc+3;
        else
            nnodes1 = nnodes;
            for dd = 1:3
                fl = 0;
                for ss = 1:nnodes1
                    if xrelem(bb,dd)==xnode(ss)&&yrelem(bb,dd)==ynode(ss)
                        nodecon(bb,dd) = ss;
                        fl = 1;
                        break
                    end
                end
                if fl == 0
                    nnodes = nnodes+1;
                    nodenum(cc) = cc;
                    nodecon(bb,dd) = nodenum(cc);
                    xnode(cc) = xrelem(bb,dd);
                    ynode(cc) = yrelem(bb,dd);
                    cc=cc+1;
                end
            end
        end
    end
    nbfgrp = length(BCs);
    filename = sprintf('refMesh%d.gri',qq);
    creategri(refelem,nnodes,xnode,ynode,nodecon,nbfgrp,LBRF,BCs,BCorderedx,BCorderedy,filename);
    %plot_fast(filename)
    [nelem,xelem,yelem,nedges,xedge,yedge,ledges,elemedge,BCs,Bcxelem,Bcyelem,LBF,PeriG,xnode,ynode] = readgri(filename);
    % 
    % BCf = zeros(nedges,1);
    % 
    % for ii = 1:nedges
    %     db = 1;
    %     for tv = 1:length(LBF)
    %         for bc = db:db+LBF(tv)-1
    %             if (xedge(ii,1)==Bcxelem(bc,1)&&yedge(ii,1)==Bcyelem(bc,1)&&xedge(ii,2)==Bcxelem(bc,2)&&yedge(ii,2)==Bcyelem(bc,2))||(xedge(ii,1)==Bcxelem(bc,2)&&yedge(ii,1)==Bcyelem(bc,2)&&xedge(ii,2)==Bcxelem(bc,1)&&yedge(ii,2)==Bcyelem(bc,1))
    %                 BCf(ii) = tv;
    %                 break
    %             end
    %         end
    %         if BCf(ii) == tv
    %             break
    %         end
    %         db = db+LBF(tv);
    %     end
    % end

    % Neighx = zeros(3,1);
    % Neighy = zeros(3,1);
    xsnode = xnode;
    ysnode = ynode;
    
    for tt = 1:nnodes
        
        fl = 0;
        [mv,ind]=mink(sqrt((xsnode(tt)-xsnode).^2+(ysnode(tt)-ysnode).^2),4);
        Neighx = xsnode(ind(2:4));
        Neighy = ysnode(ind(2:4));
        for kk = 1:length(Bcxelem(:,1))
            if (xsnode(tt) == Bcxelem(kk,1) && ysnode(tt) == Bcyelem(kk,1))||(xsnode(tt) == Bcxelem(kk,2) && ysnode(tt) == Bcyelem(kk,2))
                fl=1;
                % bound = BCOrdCon(kk);
                % if bound == 2
                %     [v,ind] = min(sqrt((xsnode(tt)-xsu).^2+(ysnode(tt)-ysu).^2));
                %     xsnode(tt) = xsu(ind(1));
                %     ysnode(tt) = ysu(ind(1));
                %     if xsnode(tt) == Bcxelem(kk,1) && ysnode(tt) == Bcyelem(kk,1)
                %         Bcxelem(kk,1) = xsu(ind(1));
                %         Bcyelem(kk,1) = ysu(ind(1));
                %     else
                %         Bcxelem(kk,2) = xsu(ind(1));
                %         Bcyelem(kk,2) = ysu(ind(1));
                %     end
                % elseif bound == 6
                %     [v,ind] = min(sqrt((xsnode(tt)-xsl).^2+(ysnode(tt)-ysl).^2));
                %     xsnode(tt) = xsl(ind(1));
                %     ysnode(tt) = ysl(ind(1));
                %     if xsnode(tt) == Bcxelem(kk,1) && ysnode(tt) == Bcyelem(kk,1)
                %         Bcxelem(kk,1) = xsl(ind(1));
                %         Bcyelem(kk,1) = ysl(ind(1));
                %     else
                %         Bcxelem(kk,2) = xsl(ind(1));
                %         Bcyelem(kk,2) = ysl(ind(1));
                %     end
                % end
                break
            end
        end

        if fl == 0
            xsnode(tt) = (1-w)*xnode(tt)+(w/3)*sum(Neighx);
            ysnode(tt) = (1-w)*ynode(tt)+(w/3)*sum(Neighy);
        end
    end
    filename = sprintf('refsmooth%d.gri',qq);
    creategri(nelem,nnodes,xsnode,ysnode,nodecon,nbfgrp,LBF,BCs,Bcxelem,Bcyelem,filename)
    
end

plot_fast(filename)
