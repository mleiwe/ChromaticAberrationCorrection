function [PSF,PSFleft,PSFright,XY,YZ,ZX]=mnl_CalculatePSF_3D(Im,Scale,BeadID,ColourID,foldername,figures)
szI=size(Im);
%% Compress X
for i=1:szI(1)
    for j=1:szI(2)
        XY(i,j)=max(Im(i,j,:));
    end
    X(i)=max(XY(i,:));
end
%% Compress Y
for i=1:szI(2)
    for j=1:szI(3)
        YZ(i,j)=max(Im(:,i,j));
    end
    Y(i)=max(YZ(i,:));
end
%% Compress Z
for i=1:szI(3)
    for j=1:szI(1)
        ZX(i,j)=max(Im(j,:,i));
    end
    Z(i)=max(ZX(i,:));
end
%% Visualise the PSF
if strcmp(figures,'y')==1
    figure('Name','X and Y PSF')
    subplot(2,2,3)
    imagesc(XY)
    set(gca,'YDir','normal')
    xlabel('Y')
    ylabel('X')
    subplot(2,2,1)
    plot(1:1:szI(2),Y)
    xlim([1 szI(2)])
    xlabel('Distribution Along Y')
    subplot(2,2,4)
    plot(X,1:1:szI(1))
    ylim([1 szI(1)])
    ylabel('Distribution Along X')
    fn=sprintf('%s%s%d%s%d%s',foldername,'\Bead_',BeadID,'_Chan',ColourID,'_XY-PSF');
    savefig(fn)
    figure('Name','X and Z PSF')
    subplot(2,2,3)
    imagesc(ZX)
    xlabel('X')
    ylabel('Z')
    set(gca,'YDir','normal')
    subplot(2,2,1)
    plot(1:1:szI(1),X)
    xlim([1 szI(1)])
    xlabel('Distribution Along X')
    subplot(2,2,4)
    plot(Z,1:1:szI(3))
    ylim([1 szI(3)])
    ylabel('Distribution Along Z')
    fn=sprintf('%s%s%d%s%d%s',foldername,'\Bead_',BeadID,'_Chan',ColourID,'_XZ-PSF');
    savefig(fn)
end
%% Calculate the PSF
%For X
mxX=max(X);
thresh=mxX/2; %Half Max
index=find(X>=thresh);
MxXLoc=find(X==mxX);
xPSF=max(index)-min(index);
xleft_PSF=MxXLoc-min(index);
xright_PSF=max(index)-MxXLoc;
%For Y
mxY=max(Y);
thresh=mxY/2; %Half Max
index=find(Y>=thresh);
MxYLoc=find(Y==mxY);
yPSF=max(index)-min(index);
yleft_PSF=MxYLoc-min(index);
yright_PSF=max(index)-MxYLoc;
%For Z
mxZ=max(Z);
thresh=mxZ/2; %Half Max
index=find(Z>=thresh);
MxZLoc=find(Z==mxZ);
zPSF=max(index)-min(index);
zleft_PSF=MxZLoc-min(index);
zright_PSF=max(index)-MxZLoc;
% Collate to Key Values
PSF=[xPSF*Scale(1) yPSF*Scale(2) zPSF*Scale(3)];
PSFleft=[xleft_PSF*Scale(1) yleft_PSF*Scale(2) zleft_PSF*Scale(3)];
PSFright=[xright_PSF*Scale(1) yright_PSF*Scale(2) zright_PSF*Scale(3)];
close all
end