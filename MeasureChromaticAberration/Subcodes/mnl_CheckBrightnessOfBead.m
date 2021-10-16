function [Include]=mnl_CheckBrightnessOfBead(Im,th)
szI=size(Im);
for c=1:szI(3)
    tempChan(:,:,:)=Im(:,:,c,:);
    t=isnan(tempChan);
    if sum(tempChan(:))==0 %if it is all zeros
        Chk(c)=1;
    elseif sum(t)>=1 %If there are NaNs
        Chk(c)=1;      
    else
        %Because the image size may be big, mark the background values (<Mean + 10SD) as zero
        tempChan=double(tempChan);
        Mean=mean(tempChan(:));
        SD=std(tempChan(:));
        Thresh=Mean+(th*SD);
        idx=tempChan<Thresh;
        tempChan(idx)=0;
        tempChan=medfilt3(tempChan,[1 1 3]);
        %If all the tempChan are zeros
        if sum(tempChan(:))<1
            Chk(c)=1;
        else
            Chk(c)=0;
        end
    end
end
Include=sum(Chk);
end