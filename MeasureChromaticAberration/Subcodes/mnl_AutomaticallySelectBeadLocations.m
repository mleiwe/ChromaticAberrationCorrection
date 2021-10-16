function [Beads]=mnl_AutomaticallySelectBeadLocations(Data)
%Input
% Data - xycz format

dim=size(Data);
nColours=dim(3);
nRows=round(sqrt(nColours));
nColumns=ceil(sqrt(nColours));
%% Create the single channel MIPs with a median filter
r=3;
[mData]=mnl_MedianFilterEachChannel_2D(Data,r);
MIPs=zeros(dim(1),dim(2),nColours);
for i=1:nColours
    tColours(:,:,:)=mData(:,:,i,:);
    tMIP(:,:)=max(tColours,[],3);
    tn=sprintf('%s%d%s','MIP_Channel ',i,'.tiff');
    saveastiff(tMIP,tn);
    filenames{i}=tn;
    MIPs(:,:,i)=tMIP;
end
prompt='Please select the channel to detect the beads in';
ChosenChan=input(prompt);
%% For the chosen channel get the minimum value
nSD=10; %Intensity Threshold
cMIP=MIPs(:,:,ChosenChan);
mVal=mean(cMIP(:));
stVal=std(cMIP(:));
IntThresh=mVal+(nSD*stVal); %Signal must be mean +nSD standard deviations
%% Now find the biggest spots in 3D
% Binarise the image based on the minimum value
BW=zeros(dim(1),dim(2),dim(4));
cData(:,:,:)=mData(:,:,ChosenChan,:);
%idx=cData>window_min;
idx=cData>IntThresh;
BW(idx)=1;
%Get the regions
[L,n]=bwlabeln(BW);
nVxThresh=100; %Volum Threshold - the area must be at least this big
[L,n]=mnl_RemoveSmallROIs(L,n,nVxThresh);
stats=regionprops3(L,'Volume','Centroid','BoundingBox');
%% Store as bead structure
nBeads=n;
Beads=struct('BeadLocation',[],'BeadSize',[],'BeadCentre',[],'BeadExtremes',[],'ImageLimits',[]);
for i=1:nBeads
    Centroid=stats.Centroid(i,:);
    BoundBox=stats.BoundingBox(i,:);
    %Calculate the size
    Xedge=[BoundBox(1) BoundBox(1)+BoundBox(4)];
    Yedge=[BoundBox(2) BoundBox(2)+BoundBox(5)];
    Zedge=[BoundBox(3) BoundBox(3)+BoundBox(6)];
    Xsize=Xedge(2)-Xedge(1);
    Ysize=Yedge(2)-Yedge(1);
    Zsize=Zedge(2)-Zedge(1);
    Beads(i).BeadSize=[Xsize Ysize Zsize];
    %Centre is same as location
    Beads(i).BeadCentre=Centroid;
    %Bead extremes
    Beads(i).BeadExtremes=[Xedge(1) Xedge(2);Yedge(1) Yedge(2);Zedge(1) Zedge(2)];
    Beads(i).BeadLocation=[BoundBox(1) BoundBox(2) BoundBox(3)];%This is the topleft corner
    %Calculate the image limits
    s=2; %safety margin - how many sizes from the bead should the image limits be
    X_ImLim=[Centroid(1)-(Xsize*s) Centroid(1)+(Xsize*2)];
    Y_ImLim=[Centroid(2)-(Ysize*s) Centroid(2)+(Ysize*2)];
    Z_ImLim=[Centroid(3)-(Zsize*s*1.5) Centroid(3)+(Zsize*s*1.5)]; %NB Z is given an extra 50% as a safety margin
    Beads(i).ImageLimits=[X_ImLim(1) X_ImLim(2);Y_ImLim(1) Y_ImLim(2);Z_ImLim(1) Z_ImLim(2)];
end
%% Visual Confirmation
figure('Name','Detected Beads')
subplot(1,2,1)
imagesc(cMIP)
title('Autoscaled Median Filtered MIP')
subplot(1,2,2)
L_MIP=max(L,[],3);
imagesc(L_MIP)
title('Auto-Detected Beads')
end
%% Subfunctions
function [L,n]=mnl_RemoveSmallROIs(L,n,nVxThresh)
newID=0;
for i=1:n
    idx=L==i;
    nVx=sum(idx(:));
    if nVx<nVxThresh
        L(idx)=0;
    else
        newID=newID+1;
        L(idx)=newID;        
    end
    mnl_InsertProgressTrackerInLoops(i,n)
end
n=newID;
end