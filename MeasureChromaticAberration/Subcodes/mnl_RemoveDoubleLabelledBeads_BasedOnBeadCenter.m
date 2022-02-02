function [fBeads]=mnl_RemoveDoubleLabelledBeads_BasedOnBeadCenter(Beads,ImDim)
sz=size(Beads,2);
n=1;
Overlap=[];
for i=1:sz
    %BeadOneImLimits=Beads(i).ImageLimits;
    BeadOneBeadCentre=Beads(i).BeadCentre;
    BeadOneSize=Beads(i).BeadSize;
    %Add the diameter to the centre (i.e. radius times 2) to make sure there is a sufficient gap between the beads
    Xrange1=BeadOneBeadCentre(1)-(BeadOneSize(1)):BeadOneBeadCentre(1)+(BeadOneSize(1));
    Yrange1=BeadOneBeadCentre(2)-(BeadOneSize(2)):BeadOneBeadCentre(2)+(BeadOneSize(2));
    Zrange1=BeadOneBeadCentre(3)-(BeadOneSize(3)):BeadOneBeadCentre(3)+(BeadOneSize(3));
    for j=i+1:sz
        %BeadTwoExtremes=Beads(j).BeadExtremes; 
        BeadTwoBeadCentre=Beads(j).BeadCentre;
        BeadTwoSize=Beads(j).BeadSize;
        %Bead Two we just use the extremes (i.e. size of the bead)
        Xrange2=BeadTwoBeadCentre(1)-(BeadTwoSize(1)/2):BeadTwoBeadCentre(1)+(BeadTwoSize(1)/2);
        Yrange2=BeadTwoBeadCentre(2)-(BeadTwoSize(2)/2):BeadTwoBeadCentre(2)+(BeadTwoSize(2)/2);
        Zrange2=BeadTwoBeadCentre(3)-(BeadTwoSize(3)/2):BeadTwoBeadCentre(3)+(BeadTwoSize(3)/2);
        
       %Check if there is an overlap
        Xoverlap=ismember(Xrange1,Xrange2);
        Yoverlap=ismember(Yrange1,Yrange2);
        Zoverlap=ismember(Zrange1,Zrange2);
        if sum(Xoverlap)>0 && sum(Yoverlap)>0 && sum(Zoverlap)>0
            Overlap(n,:)=[i j];
            n=n+1;
        end
    end
end
% Now filter it
f=1;
if isempty(Overlap==1)
    fBeads=Beads;
else
    %Remove the overlapping Beads
    for i=1:sz
        a=Overlap==i; %Is this bead in the overlap?
        if sum(a(:))==0
            fBeads(f)=Beads(i);
            f=f+1;
        end
    end
end
end