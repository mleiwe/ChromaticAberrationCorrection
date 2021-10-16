function [fBeads]=mnl_RemoveDoubleLabelledBeads(Beads,ImDim)
sz=size(Beads,2);
n=1;
Overlap=[];
for i=1:sz
    BeadOneImLimits=Beads(i).ImageLimits;
    for dim=1:3
        %If it is smaller than 1
        if BeadOneImLimits(dim,1)<1
            BeadOneImLimits(dim,1)=1;
        end
        %If it extends beyond the dimension limits
        if BeadOneImLimits(dim,2)>ImDim(dim)
            BeadOneImLimits(dim,2)=ImDim(dim);
        end
    end
    Xrange1=BeadOneImLimits(1,1):BeadOneImLimits(1,2);
    Yrange1=BeadOneImLimits(2,1):BeadOneImLimits(2,2);
    Zrange1=BeadOneImLimits(3,1):BeadOneImLimits(3,2);
    for j=i+1:sz
        BeadTwoExtremes=Beads(j).BeadExtremes;
        for dim=1:3
            %If it is smaller than 1
            if BeadTwoExtremes(dim,1)<1
                BeadTwoExtremes(dim,1)=1;
            end
            %If it extends beyond the dimension limits
            if BeadTwoExtremes(dim,2)>ImDim(dim)
                BeadTwoExtremes(dim,2)=ImDim(dim);
            end
        end
        Xrange2=BeadTwoExtremes(1,1):BeadTwoExtremes(1,2);
        Yrange2=BeadTwoExtremes(2,1):BeadTwoExtremes(2,2);
        Zrange2=BeadTwoExtremes(3,1):BeadTwoExtremes(3,2);
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