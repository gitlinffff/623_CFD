[nelem,xelem,yelem,nedges,xedge,yedge,ledges,elemedge,BCs,Bcxelem,Bcyelem,LBF,PeriG] = readgri('Onerefined.gri');

BCf = zeros(nedges,1); LBRF=LBF; 

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
        
        if BCf(ii)>0
                Bcrxelem(rt,1) = xedge(ii,1);
                Bcrxelem(rt,2) = xmed;
                Bcryelem(rt,1) = yedge(ii,1);
                Bcryelem(rt,2) = ymed;
                Bcrxelem(rt+1,2) = xedge(ii,2);
                Bcrxelem(rt+1,1) = xmed;
                Bcryelem(rt+1,2) = yedge(ii,2);
                Bcryelem(rt+1,1) = ymed;
                LBRF(BCf(ii)) = LBRF(BCf(ii))+1;
                rt = rt+2;
        end
        
    end
    refelem = nelem;
    kk = 1;
    xrelem = zeros(nelem,3);
    yrelem = zeros(nelem,3);
    
    for ii = 1:nelem
        
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
nbfgrp = 8;
creategri(refelem,nnodes,xnode,ynode,nodecon,nbfgrp,LBF,BCs,Bcxelem,Bcyelem,'uniformRefined.gri')
plotgri('uniformRefined.gri')

% [nelem,xelem,yelem,nedges,xedge,yedge,ledges,elemedge,BCs,Bcxelem,Bcyelem,LBF,PeriG] = readgri('initialTest.gri');

% BCf = zeros(nedges,1);
% for ii = 1:nedges
%         db = 1;
%         for tv = 1:length(LBF)
%             for bc = db:db+LBF(tv)-1
%                 if (xedge(ii,1)==Bcxelem(bc,1)&&yedge(ii,1)==Bcyelem(bc,1)&&xedge(ii,2)==Bcxelem(bc,2)&&yedge(ii,2)==Bcyelem(bc,2))||(xedge(ii,1)==Bcxelem(bc,2)&&yedge(ii,1)==Bcyelem(bc,2)&&xedge(ii,2)==Bcxelem(bc,1)&&yedge(ii,2)==Bcyelem(bc,1))
%                     BCf(ii) = tv;
%                     break
%                 end
%             end
%             if BCf(ii) == tv
%                 break
%             end
%             db = db+LBF(tv);
%         end
% 
% 
% 
% end