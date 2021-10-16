function [fBeads]=mnl_RemoveBeadsOnEdge(Beads,dim)
%NB Remember that I've already switched the X and Y axes so dim(2) is the
%image X and dim(1) is the image Y
szB=size(Beads,2);
n=1;
r=1;
for i=1:szB
    Xpos=Beads(i).BeadCentre(1);
    Ypos=Beads(i).BeadCentre(2);
    Zpos=Beads(i).BeadCentre(3);
    Xsize=Beads(i).BeadSize(1);
    Ysize=Beads(i).BeadSize(2);
    Zsize=Beads(i).BeadSize(3);
    %Is it within the XY(MATLAB) and Z axis limit? I've added 200% the
    %size for z now to account for edges in Z as this is where the
    %chromatic aberration is much bigger
    if Xpos+Xsize<=dim(2) && Ypos+Ysize<=dim(1) && Zpos+(Zsize*2)<=dim(4) &&...
            Xpos-Xsize>=1 && Ypos-Ysize>=1 && Zpos-(Zsize*2)>=1 %Is the bead outside of the maximums or minimums 
        %if Xpos+Xsize<=dim(2) && Ypos+Ysize<=dim(1) && Zpos+(Zsize*1.5)<=dim(4) &&...
            %Xpos-Xsize>=1 && Ypos-Ysize>=1 && Zpos-(Zsize*1.5)>=1 %Is the bead outside of the maximums or minimums 
        fBeads(n)=Beads(i);
        n=n+1;
    else
        rBeads(r)=Beads(i);
        RemovedBead(r)=i;
        r=r+1;
    end
end
sprintf('%s%d','Number of Beads removed...',szB-n+1)
end