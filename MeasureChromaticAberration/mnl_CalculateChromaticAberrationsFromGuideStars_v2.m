function [ChromaticCorrections,ROIs]=mnl_CalculateChromaticAberrationsFromGuideStars_v2
%Same as original but now with coloured figures

%% Input Channel Data
prompt='How many channels are used for imaging?';
nChan=input(prompt);
LaserEx=nan(nChan,1);
for i=1:nChan
    LaserQ=sprintf('%s%d%s','Please enter the excitation laser wavelength for channel ',i,' (nm)');
    prompt=LaserQ;
    LaserEx(i)=input(prompt);
end
%% Loop for each mini image drawn
n=1;
while 1
    %Load the image
    disp('Please load in the image... Make sure the scale is already in the metadata!')
    [Im,Scale,dim]=mnl_Load4Dimage;
    %Input the depths
    prompt='Please input the starting frame number - i.e. which point from the whole image';
    zStart=input(prompt);
    ROIs(n).Zdepth=(zStart-1)*Scale(3);
    %Calculate via cross correlation
    %[CA]=mnl_CalculateTheCrossCorrelation(Im,Scale);
    %Calculate the shift via centroid selection
    [CA]=mnl_CalculateTheAberrationWithCentroids(Im,Scale);
    ROIs(n).ChromaticAberration=CA;
    %Save the image then close it
    fign=sprintf('%s%d%s','ROI ',n,' Centroids');
    savefig(fign)
    h=gcf;
    mnl_ExportEPSdense(h,fign)
    close all
    %Check is this the last image?
    prompt='Any more images to load? y/n';
    OneMore=input(prompt,'s');
    if strcmp(OneMore,'n')
        break
    else
        disp('Input not "n", so loading another image')
        n=n+1;
    end
end
nROIs=n;
%% Now plot and calculate the regression for each channel
cmap=colormap(jet(nChan));
for i=1:nChan
    fign=sprintf('%s%d%s','Shifts relative to ',LaserEx(i),' nm');
    figure('Name',fign)
    %Scatter plot
    for j=1:nChan
        %Now get the XY values
        x=nan(1,nROIs);
        y=nan(1,nROIs);
        for k=1:nROIs
            x(k)=ROIs(k).Zdepth;
            y(k)=ROIs(k).ChromaticAberration(i).DistToChannel(j);
        end
        scatter(x,y,10,cmap(j,:),'filled')
        hold on
        [m(j),c(j),pval(j),r2(j),F(j),bint(j).Matrix,r(j,:),rint(j).Matrix]=mnl_CalculateRegression(x',y');
    end
    %Store Chromatic Aberrations
    ChromaticCorrections(i).ToWhichLaser=LaserEx(i);
    ChromaticCorrections(i).ForWhichLaser=LaserEx;
    ChromaticCorrections(i).mValues=m;
    ChromaticCorrections(i).cValues=c;
    ChromaticCorrections(i).pValues=pval;
    ChromaticCorrections(i).r2Values=r2;
    ChromaticCorrections(i).FValues=F;
    ChromaticCorrections(i).bint=bint;
    ChromaticCorrections(i).r=r;
    ChromaticCorrections(i).rint=rint;
    %Plot the regression lines
    xlimvals=xlim;
    Xrange=0:xlimvals(2);
    for j=1:nChan
        Yrange=(Xrange.*m(j))+c(j);
        plot(Xrange,Yrange,'Color',cmap(j,:))
    end
end

end
%% Sub-functions
function [CA]=mnl_CalculateTheAberrationWithCentroids(Im,Scale)
%Function to calculate the chromatic aberration by calculating the centre
%of mass for each image
szIm=size(Im);
Image_cmap=mnl_CalculateColourmaps(szIm(3)); %Calculate the colourmaps to apply
%Maximum intensity projection along XZ for each channel
figure('Name','Please select your dendritic "guidestars"...')
colormap(gray)
CombinedMIP=zeros(szIm(4),szIm(2),szIm(3));
for i=1:szIm(3)
    tIm(1:szIm(1),1:szIm(2),1:szIm(4))=Im(:,:,i,:); %Isolate each channel
    %NB remember y is the first because it is the row
    tMIP_XZ(:,:)=max(tIm,[],1);
    MIP_XZ=tMIP_XZ';
    CombinedMIP(:,:,i)=MIP_XZ;
    subplot(1,szIm(3)+1,i)
    %Covert to RGB
    MIP_XZ=double(MIP_XZ);
    rMIP_XZ=mnl_ConvertToRGB_SingleChannel(MIP_XZ,i,szIm(3));
    %Resize to make z the same size as x
    image(rMIP_XZ)
    axis equal
    axis off
    ylim([1 szIm(4)])
    hold on
    % Calculate the centre of mass
    Centroid(i,:)=centerOfMass(tIm);
    ROIs(i).xCent=Centroid(i,2);
    ROIs(i).yCent=Centroid(i,3);  
    scatter(Centroid(i,2),Centroid(i,3),250,'xk')
    axis off
    % Add in scale bar if it is the first image
    if i==1
        mnl_ScaleBar(Scale(2),1,'northwest','1 \mum')
    end
    %Add title
    tn=sprintf('%s%d','Channel ',i);
    title(tn)
end
% Now plot the merge
[MergeRGB]=mnl_ConvertToRGB(CombinedMIP);
subplot(1,szIm(3)+1,szIm(3)+1)
image(MergeRGB)
axis equal
axis off
ylim([1 szIm(4)])
hold on
[Image_cmap]=mnl_CalculateColourmaps(szIm(3));
for i=1:szIm(3)
    scatter(Centroid(i,2),Centroid(i,3),500,Image_cmap(i,:),'x')
end
title('Merge')  
%Now for each reference channel calculate z(y?) distance between the images
for i=1:szIm(3)
    RefZ=ROIs(i).yCent; %Reference Z
    for j=1:szIm(3)
        AltZ=ROIs(j).yCent; %Alternate Z
        Distance=AltZ-RefZ;
        %Store in the Structure
        CA(i).DistToChannel_px(j)=Distance; %-Negative values mean the channel is higher than the reference
        CA(i).DistToChannel(j)=(Distance).*Scale(3);
    end
end

end
function [Image_cmap]=mnl_CalculateColourmaps(n)
if n==1
    Image_cmap=[1 0 0];
elseif n==2
    Image_cmap=[0 1 0;1 0 1];
elseif n==3
    Image_cmap=[0 0 1;0 1 0;1 0 0];
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
function [m,c,pval,r2,F,bint,r,rint]=mnl_CalculateRegression(x,y)
%function to calculate the line of best fit. 
% Inputs
% x - list of X values
% y - list of Y values
%
% Outputs
% values for the equation y=mx+c
sz=size(x,1);
X=[ones(sz,1),x(:,1)]; %attach ones
[b,bint,r,rint,stats]=regress(y,X);
if stats(3)>0.05 %Justified via p value. Not looking to explain all the variance, just how good the model is. Hence why I'm not using the R (stats(1))
    sprintf('%s','Warning this is not a good fit')
end
m=b(2);c=b(1); %for y=mx+c
pval=stats(3);
r2=stats(1);
F=stats(2);
end