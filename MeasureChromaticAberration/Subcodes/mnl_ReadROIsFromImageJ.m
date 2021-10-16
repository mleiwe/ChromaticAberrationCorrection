function [Beads]=mnl_ReadROIsFromImageJ(fn,Data)
%Function that uses the code written by Dylan Muir Dylan Muir (2021).
%ReadImageJROI (https://github.com/DylanMuir/ReadImageJROI), GitHub.
%Retrieved June 3, 2021. And then ammends the values to fit the required
%format for mnl_Pipeline_MeasurePSFandCA
% Inputs
% fn - the filename
% Data - the full image (x*y*c*z)
%
% Beads - structure matrix for beads
%% Basic Info
ImDim=size(Data);

%% ReadImageJROI
[sROI]=ReadImageJROI(fn);
nROI=size(sROI,2);
%% Now convert into required bead format
Beads=struct('BeadLocation',[],'BeadSize',[],'BeadCentre',[],'BeadExtremes',[],'ImageLimits',[]);
nB=1;
for i=1:nROI
    RectanglePos=sROI{1,i}.vnRectBounds;
    %Top Left Corner
    Ytl=RectanglePos(1);
    Xtl=RectanglePos(2);
    if Ytl==0
        Ytl=1;
    end
    if Xtl==0
        Xtl=1;
    end
    %Bottom Right Corner
    Ybr=RectanglePos(3);
    Xbr=RectanglePos(4);
    if Ytl>ImDim(1)
        Ytl=ImDim(1);
    end
    if Xtl>ImDim(2)
        Xtl=ImDim(2);
    end
    %Find the Image Limits
    ImageLimits=[Xtl Xbr;Ytl Ybr;1 ImDim(4)];
    %Get Find the Bead for each channel (final positions will be based on channel 1)
    Zextremes=nan(ImDim(3),2);
    for j=1:ImDim(3)
        TrimData(:,:,:)=Data(Ytl:Ybr,Xtl:Xbr,j,:);
        [tBeadLocation,tBeadSize,tBeadCentre,tBeadExtremes]=mnl_SubFunctionFindBead(TrimData);
        if j==1
            BeadLocation=tBeadLocation;
            BeadSize=tBeadSize;
            BeadCentre=tBeadCentre;
            BeadExtremes=tBeadExtremes;
        end
        TempStruct(j).BeadLocation=tBeadLocation;
        TempStruct(j).BeadSize=tBeadSize;
        TempStruct(j).BeadCentre=tBeadCentre;
        TempStruct(j).BeadExtremes=tBeadExtremes;
        Zextremes(j,:)=tBeadExtremes(3,:);
        clear TrimData
    end
    BeadLocation=[BeadLocation(1,1)+Xtl BeadLocation(1,2)+Ytl BeadLocation(1,3)]; %Update the bead location to be for the whole image
    BeadCentre=[BeadCentre(1,1)+Xtl BeadCentre(1,2)+Ytl BeadCentre(1,3)]; %Update the bead centre to be for the whole image
    %Check to make sure the bead extremes don't hit the edge in z
    minZex=min(Zextremes(:,1));
    maxZex=max(Zextremes(:,2));
    if minZex>1 && maxZex<ImDim(4) 
        Beads(nB).BeadLocation=BeadLocation;
        Beads(nB).BeadSize=BeadSize;
        Beads(nB).BeadCentre=BeadCentre;
        Beads(nB).BeadExtremes=BeadExtremes;
        Beads(nB).ImageLimits=ImageLimits;
        nB=nB+1;
    end
    clear TrimData
end
end
%% Subfunction
function [BeadLocation,BeadSize,BeadCentre,BeadExtremes]=mnl_SubFunctionFindBead(TrimData)
%% Compress each dimension
szI=size(TrimData);
% Compress X
for i=1:szI(1)
    for j=1:szI(2)
        XY(i,j)=max(TrimData(i,j,:));
    end
    X(i)=max(XY(i,:));
end
% Compress Y
for i=1:szI(2)
    for j=1:szI(3)
        YZ(i,j)=max(TrimData(:,i,j));
    end
    Y(i)=max(YZ(i,:));
end
% Compress Z
for i=1:szI(3)
    for j=1:szI(1)
        ZX(i,j)=max(TrimData(j,:,i));
    end
    Z(i)=max(ZX(i,:));
end
%% Find the location of the bead by the brightest point along each dimension
%For X
[MxXLoc,xSize,Xextremes]=mnl_FindPeakAndFWHH(X);
%For Y
[MxYLoc,ySize,Yextremes]=mnl_FindPeakAndFWHH(Y);
%For Z
[MxZLoc,zSize,Zextremes]=mnl_FindPeakAndFWHH(Z);
%Calculate the centre of mass - remove all data below mean + 3SD
tempChan=double(TrimData);
Mean=mean(tempChan(:));
SD=std(tempChan(:));
Thresh=Mean+(10*SD);
idx=tempChan<Thresh;
tempChan(idx)=0;
CentrePoints=mnl_3DcentreofMass(tempChan);
% Now the outputs
BeadLocation=[MxXLoc MxYLoc MxZLoc];
BeadSize=[xSize ySize zSize];
BeadCentre=round(CentrePoints)-1;
BeadExtremes=[Xextremes(1) Xextremes(2);Yextremes(1) Yextremes(2);Zextremes(1) Zextremes(2)];
end
function [Location,FWHH,extremes]=mnl_FindPeakAndFWHH(X)
%First remove background signal by using the median
medVal=median(X);
X=X-medVal;
%Then smoothen it
[X]=mnl_CreateRollingAverage(X,3);%smooth by 3px either side
%Now calculate the centre and size
[mxX,Loc]=max(X);
Location=Loc+1;
thresh=mxX/2; %Half Max
index=X>=thresh;
%Find left side
Left_Idx=index(1:Loc);
szLI=size(Left_Idx,2);
nL=0;
for i=1:szLI
    pos=szLI-(i-1);
    if Left_Idx(pos)==0
        nL=nL+1;
        BelowThresh(nL)=pos;
    end
end
if nL>0
    Lpos=max(BelowThresh)+1;
    extremes(1,1)=Lpos;
else
    Lpos=1;
    extremes(1,1)=Lpos;
end
clear BelowThresh
%Find right side
Right_Idx=index(Location:end);
szRI=size(Right_Idx,2);
nR=0;
for i=1:szRI
    Rpos=i+Location-1;
    if Right_Idx(i)==0
        nR=nR+1;
        BelowThresh(nR)=Rpos;
    end
end
if nR>0
    Rpos=min(BelowThresh)-1;
    extremes(1,2)=Rpos;
else
    Rpos=size(index,2);
    extremes(1,2)=Rpos;
end
FWHH=Rpos-Lpos+1;
end
function [mX]=mnl_CreateRollingAverage(X,r)
%function that averages by r points in both directions to create a smoother
%trace
szX=size(X,2);
mX=nan(1,szX);
%Rolling Average
for i=1:szX
    St=i-r; 
    Ed=i+r;
    if St<=0
        St=1;
    end
    if Ed>szX
        Ed=szX;
    end
    mX(i)=mean(X(St:Ed));
end
end