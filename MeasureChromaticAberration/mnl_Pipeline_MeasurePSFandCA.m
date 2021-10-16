function [f2Beads,Scale]=mnl_Pipeline_MeasurePSFandCA(fname)
%% Now Load in the Image
[Data,Scale,dim]=mnl_Load4Dimage;
Mdim=[dim(2) dim(1) dim(3) dim(4)]; %The MATLAB way of storing dimensions
%% Load in the Beads
if exist('fname','var')==1
    ID=fname(end-3:end);
else
    ID=[];
end
if strcmp(ID,'.xls')==1
    % Extract Bead Information from the ICY spot detection excell sheet
    %[Beads]=mnl_ExtractBeadInfoFromICY(fname); %NB The x and y have been switched to follow the MATLAB conventions
    [Beads]=mnl_ExtractBeadInfoFromICY_temp(fname); 
    % Filter out beads that are too close
    [fBeads]=mnl_RemoveDoubleLabelledBeads(Beads,dim);
    % Filter out beads that are partly outside the image
    [f2Beads]=mnl_RemoveBeadsOnEdge(fBeads,dim);
elseif strcmp(ID,'.zip')==1
    [Beads]=mnl_ReadROIsFromImageJ(fname,Data);
    % Filter out beads that are too close
    [fBeads]=mnl_RemoveDoubleLabelledBeads_BasedOnBeadCenter(Beads,dim);
    % Filter out beads that are partly outside the image
    [f2Beads]=mnl_RemoveBeadsOnEdge(fBeads,dim);
else
    prompt='Could not recognise input from either ICY or ImageJ, do you want to automatically detect the beads? y/n';
    Auto=input(prompt,'s');
    if strcmp(Auto,'y')==1
        [Beads]=mnl_AutomaticallySelectBeadLocations(Data);
        % Filter out beads that are too close
        [fBeads]=mnl_RemoveDoubleLabelledBeads_BasedOnBeadCenter(Beads,dim);
        % Filter out beads that are partly outside the image
        [f2Beads]=mnl_RemoveBeadsOnEdge(fBeads,dim);
        % Check the other channels
    else
        error('Please go back and detect ROIs first')
    end
end
%% Now find the PSF and Chromatic Aberration of each bead
prompt='Please input the Starting Z depth of the image (relative to the start of the block, the units of the image (usually microns)';
Zstart=input(prompt);
PSFfoldername='PSF figures';
mkdir(PSFfoldername);
CAfoldername='ChromaticAberration';
mkdir(CAfoldername);
szBeads=size(f2Beads,2); %The number of Beads
%Select 5 beads at random
PossibleBeads=randi(szBeads,5,1);
uPB=unique(PossibleBeads);%if by chance the same bead is selected it will use fewer
for i=1:szBeads
    [MkFig]=find(uPB==i);
    if isempty(MkFig)==0
        figures='y';
    else
        figures='n';
    end
    %DistanceFromCentre
    Xcent=f2Beads(i).BeadCentre(1);
    Ycent=f2Beads(i).BeadCentre(2);
    XCentIM=dim(1)/2;
    YCentIM=dim(2)/2;
    XDistFromCentre=(Xcent-XCentIM)*Scale(1);
    YDistFromCentre=(Ycent-YCentIM)*Scale(2);
    f2Beads(i).DistFromCentre=sqrt((XDistFromCentre^2)+(YDistFromCentre^2));
    f2Beads(i).XDistFromCentre=XDistFromCentre;
    f2Beads(i).YDistFromCentre=YDistFromCentre;
    %Z position
    f2Beads(i).Zdepth=Zstart+(f2Beads(i).BeadCentre(3)*Scale(3));
    %Generate the small image of the bead
    Data2=[];
    fData2=[];
    tData=[];
    Xmin=f2Beads(i).ImageLimits(1,1);
    if Xmin<1
        Xmin=1;
    end
    Xmax=f2Beads(i).ImageLimits(1,2);
    if Xmax>Mdim(1)
        Xmax=Mdim(1);
    end
    Ymin=f2Beads(i).ImageLimits(2,1);
    if Ymin<1
        Ymin=1;
    end
    Ymax=f2Beads(i).ImageLimits(2,2);
    if Ymax>Mdim(2)
        Ymax=Mdim(2);
    end
    Zmin=f2Beads(i).ImageLimits(3,1);
    if Zmin<1
        Zmin=1;
    end
    Zmax=f2Beads(i).ImageLimits(3,2);
    if Zmax>Mdim(4)
        Zmax=Mdim(4);
    end
    if strcmp(ID,'.zip')==1
        Data2=Data(Ymin:Ymax,Xmin:Xmax,:,Zmin:Zmax);
    else
        Data2=Data(Xmin:Xmax,Ymin:Ymax,:,Zmin:Zmax);
    end
    szData2=size(Data2);
    for j=1:dim(3)
        tData(:,:,:)=Data2(:,:,j,:);
        medVal(j)=median(tData(:));
        fData2(:,:,j,:)=Data2(:,:,j,:)-medVal(j);
    end
    %Check if the seed is bright enough
    th=3;
    [Include]=mnl_CheckBrightnessOfBead(fData2,th);
    if Include==0
        %Calculate the PSF
        for c=1:szData2(3)
            clear Im
            Im(:,:,:)=fData2(:,:,c,:);
            [tPSF,tPSFleft,tPSFright,XY,YZ,ZX]=mnl_CalculatePSF_3D(Im,Scale,i,c,PSFfoldername,figures);
            f2Beads(i).PSF.Channel(c).PSF=tPSF;
            f2Beads(i).PSF.Channel(c).PSFleft=tPSFleft;
            f2Beads(i).PSF.Channel(c).PSFright=tPSFright;
            f2Beads(i).MIPs.Channel(c).XY=XY;
            f2Beads(i).MIPs.Channel(c).YZ=YZ;
            f2Beads(i).MIPs.Channel(c).ZX=ZX;
        end
        %Calculate the chromatic aberration
        Bead=f2Beads(i);
        [tCA]=mnl_MeasureChromaticAberration(Bead,fData2,Scale,i,CAfoldername,figures);
        f2Beads(i).ChromaticAberration=tCA;
    end
end
%% Filter to remove the poor beads
n=1;
for i=1:szBeads
    if isempty(f2Beads(i).PSF)==0
        f3Beads(n)=f2Beads(i);
        n=n+1;
    end
end
f2Beads=f3Beads;
szBeads=size(f2Beads,2);
%% Measure PSF For this Stack
% For the X,Y, and Z spread
Cent_xPSF=nan(szBeads,Mdim(3),2); %Pre-allocation
Cent_yPSF=nan(szBeads,Mdim(3),2); %Pre-allocation
Cent_zPSF=nan(szBeads,Mdim(3),2); %Pre-allocation
Depth_xPSF=nan(szBeads,Mdim(3),2); %Pre-allocation
Depth_yPSF=nan(szBeads,Mdim(3),2); %Pre-allocation
Depth_zPSF=nan(szBeads,Mdim(3),2); %Pre-allocation
for i=1:szBeads %per bead
    for j=1:Mdim(3) %per channel
        %To the Centre
        Cent_xPSF(i,j,1)=f2Beads(i).DistFromCentre;
        Cent_xPSF(i,j,2)=f2Beads(i).XDistFromCentre;
        Cent_xPSF(i,j,3)=f2Beads(i).YDistFromCentre;
        Cent_xPSF(i,j,4)=f2Beads(i).PSF.Channel(j).PSF(1);
        Cent_yPSF(i,j,1)=f2Beads(i).DistFromCentre;
        Cent_yPSF(i,j,2)=f2Beads(i).XDistFromCentre;
        Cent_yPSF(i,j,3)=f2Beads(i).YDistFromCentre;
        Cent_yPSF(i,j,4)=f2Beads(i).PSF.Channel(j).PSF(2);
        Cent_zPSF(i,j,1)=f2Beads(i).DistFromCentre;        
        Cent_zPSF(i,j,2)=f2Beads(i).XDistFromCentre;
        Cent_zPSF(i,j,3)=f2Beads(i).YDistFromCentre;
        Cent_zPSF(i,j,4)=f2Beads(i).PSF.Channel(j).PSF(3);
        %To the Zdepth
        Depth_xPSF(i,j,1)=f2Beads(i).Zdepth;
        Depth_xPSF(i,j,2)=f2Beads(i).PSF.Channel(j).PSF(1);
        Depth_yPSF(i,j,1)=f2Beads(i).Zdepth;
        Depth_yPSF(i,j,2)=f2Beads(i).PSF.Channel(j).PSF(2);
        Depth_zPSF(i,j,1)=f2Beads(i).Zdepth;
        Depth_zPSF(i,j,2)=f2Beads(i).PSF.Channel(j).PSF(3);
    end        
end
%% Now to measure the chromatic aberration
figure('Name','Aberration to each channel')
for c=1:Mdim(3)
    for i=1:szBeads
        %Establish the Dist and dept values
        ToChan(i,1:Mdim(3),1)=f2Beads(i).DistFromCentre;
        ToChan(i,1:Mdim(3),2)=f2Beads(i).XDistFromCentre;
        ToChan(i,1:Mdim(3),3)=f2Beads(i).YDistFromCentre;
        ToChan(i,1:Mdim(3),4)=f2Beads(i).Zdepth;    
        %Now Shifts per channel
        for j=1:Mdim(3)
            ToChan(i,j,5)=f2Beads(i).ChromaticAberration(1).DistToChannel(j,1);
            ToChan(i,j,6)=f2Beads(i).ChromaticAberration(1).DistToChannel(j,2);
            ToChan(i,j,7)=f2Beads(i).ChromaticAberration(1).DistToChannel(j,3);
        end
    end
    subplot(2,Mdim(3),c)
    mnl_PlotXYZchromatic_ToCent(ToChan);
    tn=sprintf('%s%d','To Channel ',c);
    title(tn)
    subplot(2,Mdim(3),c+Mdim(3))
    mnl_PlotXYZchromatic_ToDepth(ToChan)
    if c==1
        ChosenChan=ToChan;
    end
end
%% Chromatic aberration to Chan 1 (X, Y and Z separate
figure('Name','Aberration to channel 1')
%Xdist vs Xshift
subp=1;
mnl_PlotXYZchromatic_ToXDist(ChosenChan,subp)
%Ydist
subp=4;
mnl_PlotXYZchromatic_ToYDist(ChosenChan,subp)
%Zdist
subp=7;
mnl_PlotXYZchromatic_ToZDist(ChosenChan,subp)
end

function mnl_PlotXYZchromatic_ToCent(ToChan1)
sz=size(ToChan1,2);
cmap=colormap(jet(sz));
ln=1;%legend counter
for i=1:sz%Loop for X drifts
    plot(ToChan1(:,i,1),ToChan1(:,i,5),'.','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Xdrift');
    ln=ln+1;
end
for i=1:sz%Loop for Y drifts
    plot(ToChan1(:,i,1),ToChan1(:,i,6),'x','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Ydrift');
    ln=ln+1;
end
for i=1:sz%Loop for Z drifts
    plot(ToChan1(:,i,1),ToChan1(:,i,7),'o','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Zdrift');
    ln=ln+1;
end
legend(LegName,'Location','northeastoutside')
xlabel('Distance to Centre of Image')
ylabel('Distance Between Bead Channels')
end
function mnl_PlotXYZchromatic_ToDepth(ToChan1)
sz=size(ToChan1,2);
cmap=colormap(jet(sz));
ln=1;%legend counter
for i=1:sz%Loop for X drifts
    plot(ToChan1(:,i,4),ToChan1(:,i,5),'.','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Xdrift');
    ln=ln+1;
end
for i=1:sz%Loop for Y drifts
    plot(ToChan1(:,i,4),ToChan1(:,i,6),'x','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Ydrift');
    ln=ln+1;
end
for i=1:sz%Loop for Z drifts
    plot(ToChan1(:,i,4),ToChan1(:,i,6),'o','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Zdrift');
    ln=ln+1;
end
legend(LegName,'Location','northeastoutside')
xlabel('Z depth of Bead')
ylabel('Distance Between Bead Channels')
end
function mnl_PlotXYZchromatic_ToXDist(ToChan1,subp)
sz=size(ToChan1,2);
cmap=colormap(jet(sz));
ln=1;%legend counter
subplot(3,3,subp)
for i=1:sz%Loop for X drifts
    plot(ToChan1(:,i,2),ToChan1(:,i,5),'.','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Xdrift');
    ln=ln+1;
end
xlabel('X Distance to Centre of Image (um)')
ylabel('X Distance Between Bead Channels (um)')
subplot(3,3,subp+1)
for i=1:sz%Loop for Y drifts
    plot(ToChan1(:,i,2),ToChan1(:,i,6),'.','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Ydrift');
    ln=ln+1;
end
xlabel('X Distance to Centre of Image (um)')
ylabel('Y Distance Between Bead Channels (um)')
subplot(3,3,subp+2)
for i=1:sz%Loop for Z drifts
    plot(ToChan1(:,i,2),ToChan1(:,i,7),'.','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Zdrift');
    ln=ln+1;
end
legend(LegName,'Location','northeastoutside')
xlabel('X Distance to Centre of Image (um)')
ylabel('Z Distance Between Bead Channels (um)')
end
function mnl_PlotXYZchromatic_ToYDist(ToChan1,subp)
sz=size(ToChan1,2);
cmap=colormap(jet(sz));
ln=1;%legend counter
subplot(3,3,subp)
for i=1:sz%Loop for X drifts
    plot(ToChan1(:,i,3),ToChan1(:,i,5),'.','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Xdrift');
    ln=ln+1;
end
xlabel('Y Distance to Centre of Image (um)')
ylabel('X Distance Between Bead Channels (um)')
subplot(3,3,subp+1)
for i=1:sz%Loop for Y drifts
    plot(ToChan1(:,i,3),ToChan1(:,i,6),'.','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Ydrift');
    ln=ln+1;
end
xlabel('Y Distance to Centre of Image (um)')
ylabel('Y Distance Between Bead Channels (um)')
subplot(3,3,subp+2)
for i=1:sz%Loop for Z drifts
    plot(ToChan1(:,i,3),ToChan1(:,i,7),'.','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Zdrift');
    ln=ln+1;
end
legend(LegName,'Location','northeastoutside')
xlabel('Y Distance to Centre of Image (um)')
ylabel('Z Distance Between Bead Channels (um)')
end
function mnl_PlotXYZchromatic_ToZDist(ToChan1,subp)
sz=size(ToChan1,2);
cmap=colormap(jet(sz));
ln=1;%legend counter
subplot(3,3,subp)
for i=1:sz%Loop for X drifts
    plot(ToChan1(:,i,4),ToChan1(:,i,5),'.','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Xdrift');
    ln=ln+1;
end
xlabel('Z Distance toSurface (um)')
ylabel('X Distance Between Bead Channels (um)')
subplot(3,3,subp+1)
for i=1:sz%Loop for Y drifts
    plot(ToChan1(:,i,4),ToChan1(:,i,6),'.','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Ydrift');
    ln=ln+1;
end
xlabel('Z Distance toSurface (um)')
ylabel('Y Distance Between Bead Channels (um)')
subplot(3,3,subp+2)
for i=1:sz%Loop for Z drifts
    plot(ToChan1(:,i,4),ToChan1(:,i,7),'.','Color',cmap(i,:))
    hold on
    LegName{ln}=sprintf('%s%d%s','Channel ',i,' Zdrift');
    ln=ln+1;
end
legend(LegName,'Location','northeastoutside')
xlabel('Z Distance toSurface (um)')
ylabel('Z Distance Between Bead Channels (um)')
end