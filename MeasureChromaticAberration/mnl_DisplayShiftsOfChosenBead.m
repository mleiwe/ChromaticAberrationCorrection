function mnl_DisplayShiftsOfChosenBead(f2Beads,ChosenBead,Scale,Zrange_um)
%Function designed to plot the effects of chromatic aberration in each
%channel as RGB components
%Initial Set up
Bead=f2Beads(ChosenBead);
nChan=size(Bead.ChromaticAberration,2);
nplots=2+nChan;
cmap=colormap(lines(nChan));
BeadLocation=Bead.BeadLocation;
ImageLimits=Bead.ImageLimits;
UpperCorners=ImageLimits(:,1);
UpperCorners=UpperCorners';
BeadImageLocation=Bead.ChromaticAberration(1).BasePosition;
ZdepthPx_C1=BeadLocation(3);
ZRangePx=Zrange_um/Scale(3);
ChosenZlim=[ZdepthPx_C1-ZRangePx ZdepthPx_C1+ZRangePx];
Centroid_cmap=colormap(lines(nChan));
PxShifts=Bead.ChromaticAberration(1).DistToChannel_px;
%Calculate the centroid positions
CentroidPos=nan(nChan,3);
for i=1:nChan
    CentroidPos(i,:)=BeadImageLocation+PxShifts(i,:);
end
%Legend Names
for i=1:nChan
    LegNames{i}=sprintf('%s%d%s','Channel ',i,' Channel');
end
%% Now the figure
%Subplot 1 - the xy mip (merge)
subplot(nChan+1,nplots,1)
for c=1:nChan
    tXY(:,:,c)=Bead.MIPs.Channel(c).XY;
end
[RGB_data]=mnl_ConvertToRGB(tXY);
dim=size(RGB_data);
image(RGB_data)
hold on
axis equal
%Now plot the centroid positions
for i=1:nChan
    scatter(CentroidPos(i,2),CentroidPos(i,1),200,Centroid_cmap(i,:),'x','LineWidth',2.5)
end
title('XY MIP merge')
xlim([1 dim(2)]);
ylim([1 dim(1)]);
mnl_ScaleBar(Scale(1),1,'northeast','1 um')
%Subplot 2 - the xy mip single channels
for i=1:nChan
    subplot(nChan+1,nplots,(nplots*i)+1)
    XY(:,:)=tXY(:,:,i);
    XY_RGB=mnl_ConvertToRGB_SingleChannel(XY,i,nChan);
    image(XY_RGB)
    hold on
    axis equal
    scatter(CentroidPos(i,2),CentroidPos(i,1),200,Centroid_cmap(i,:),'x','LineWidth',2.5)
    xlim([1 dim(2)]);
    ylim([1 dim(1)]);
end
%Then YZ plots for each channel
for c=1:nChan
    %Image first
    tYZ(:,:,c)=Bead.MIPs.Channel(c).YZ';
    YZ(:,:)=tYZ(:,:,c);
    ZY_RGB=mnl_ConvertToRGB_SingleChannel(YZ,c,nChan);
    subplot(1,nplots,c+1)
    image(ZY_RGB)
    hold on
    %Now the centre position
    scatter(CentroidPos(c,2),CentroidPos(c,3),200,Centroid_cmap(c,:),'x','LineWidth',2.5)   
    tn=sprintf('%s%d','XZ MIP Channel ',c);
    title(tn)
    axis equal
    dim=size(ZY_RGB);
    xlim([1 dim(2)]);
    ylim([1 dim(1)]);
    %ylim(ChosenZlim);
end
%Then a merge
subplot(1,nplots,nplots)
[RGB_data]=mnl_ConvertToRGB(tYZ);
image(RGB_data)
hold on
for c=1:nChan
    scatter(CentroidPos(c,2),CentroidPos(c,3),200,Centroid_cmap(c,:),'x','LineWidth',2.5)  
end
title('Merge')
axis equal
xlim([1 dim(2)]);
ylim([1 dim(1)]);
%ylim(ChosenZlim);
mnl_ScaleBar(Scale(1),1,'northeast','1 um')
legend(LegNames)
end
%% Subfunctions
function [Image_cmap]=mnl_CalculateColourmaps(n)
if n==1
    Image_cmap=[1 0 0];
elseif n==2
    Image_cmap=[0 1 0;1 0 1];
elseif n==3
    Image_cmap=[0 0 1;0 1 0;1 0 1];
elseif n==4
    Image_cmap=[0.5 0 1;0 0 1;0 1 0;1 0 0];
elseif n==5
    Image_cmap=[0.5 0 1;0 0 1;0 1 1;0 1 0;1 0 0];
elseif n==6
    Image_cmap=[0.5 0 1;0 0 1;0 1 1;0 1 0;1 1 0;1 0 0];
elseif n==7
    Image_cmap=[0.5 0 1;0 0 1;0 1 1;0 1 0;1 1 0;1 0.5 0;1 0 0];
else
    Image_cmap=colormap(jet(n));
end
end
function [RGB_data]=mnl_ConvertToRGB(Image)
%Makes an RGB image 
dim=size(Image);
nChan=dim(3);
%Normalise each channel to the 99.9th percentile
for i=1:nChan
    Temp(:,:)=Image(:,:,i);
    mxT=prctile(Temp(:),99.9);
    nTemp=Temp./mxT;
    nImage(:,:,i)=nTemp;
    clear nTemp
end
%fix to make above 1 equal to 1
idx=nImage>1;
nImage(idx)=1;
%Decide what colours to pick
[ChanColour]=mnl_CalculateColourmaps(nChan);
%Apply these colours to the nImage
RGB_data=zeros(dim(1),dim(2),3);
for i=1:nChan
    %Go through RGB
    for j=1:3
        if ChanColour(i,j)>0
            temp(:,:)=nImage(:,:,i);
            tRGB(:,:)=temp.*ChanColour(i,j);            
            RGB_data(:,:,j)=RGB_data(:,:,j)+tRGB;
        end
    end
end
%If RGB data is >1 make it 1
idx=RGB_data(:)>1;
RGB_data(idx)=1;
end
function [RGB_data]=mnl_ConvertToRGB_SingleChannel(Image,ChanNum,TotalChan)
%Makes an RGB image 
dim=size(Image);
%Normalise channel
Temp(:,:)=Image(:,:);
mxT=prctile(Temp(:),99.9); %Norm to 99.9th prcentile
nTemp=Temp./mxT;
idx=nTemp>1;
nTemp(idx)=1;
nImage(:,:)=nTemp;
clear nTemp Temp
%Decide what colours to pick
[ChanColour]=mnl_CalculateColourmaps(TotalChan);
%Apply these colours to the nImage
RGB_data=zeros(dim(1),dim(2),3);
%Go through RGB
for j=1:3
    if ChanColour(ChanNum,j)>0
        temp(:,:)=nImage(:,:);
        tRGB(:,:)=temp.*ChanColour(ChanNum,j);
        RGB_data(:,:,j)=RGB_data(:,:,j)+tRGB;
    end
end
end