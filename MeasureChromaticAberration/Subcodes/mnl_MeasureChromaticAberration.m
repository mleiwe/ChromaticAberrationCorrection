function [CA]=mnl_MeasureChromaticAberration(Bead,Im,Scale,BeadID,foldername,figures)
%Function to measure the chromatic aberration of beads. NB this uses the
%code "centerOfMass" written by Jered R Wells. This is added as a
%subfunction below.

%% Calculate the centre of mass for each channel
if strcmp(figures,'y')==1
figure('WindowState','maximized')
end
szI=size(Im);
for c=1:szI(3)
    tempChan(:,:,:)=Im(:,:,c,:);
    t=isnan(tempChan);
    if sum(tempChan(:))==0 %if it is all zeros
        CentrePoints(c,:)=nan(1,3);
    elseif sum(t(:))>1 %if there is a single NaN in the image
        CentrePoints(c,:)=nan(1,3);
    else
        %Because the image size may be big, mark the background values (<Mean + 10SD) as zero
        tempChan=double(tempChan);
        Mean=mean(tempChan(:));
        SD=std(tempChan(:));
        %Thresh=Mean+(10*SD);
        Thresh=Mean+(5*SD);
        idx=tempChan<Thresh;
        tempChan(idx)=0;
        %Now add a median filter to remove single hot points
        tempChan=medfilt3(tempChan,[1 1 3]);
        if ~isnan(centerOfMass(tempChan))
            CentrePoints(c,:)=centerOfMass(tempChan);
        else
            CentrePoints(c,:)=nan(1,3);
        end
    end
    
    if strcmp(figures,'y')==1
        %For XY
        subplot(szI(3),3,((c-1)*3)+1)
        imagesc(Bead.MIPs.Channel(c).XY)
        set(gca,'YDir','normal')
        hold on
        scatter(CentrePoints(c,2),CentrePoints(c,1),'x','MarkerFaceColor',[1 1 0],'LineWidth',2)
        tn=sprintf('%s%d%s%d%s','Bead-',BeadID,' Chan-',c,' XY');
        title(tn)
        xlabel('Y')
        ylabel('X')
        %For YZ
        subplot(szI(3),3,((c-1)*3)+2)
        imagesc(Bead.MIPs.Channel(c).YZ)
        set(gca,'YDir','normal')
        hold on
        scatter(CentrePoints(c,3),CentrePoints(c,2),'x','MarkerFaceColor',[1 1 0],'LineWidth',2)
        tn=sprintf('%s%d%s%d%s','Bead-',BeadID,' Chan-',c,' YZ');
        title(tn)
        xlabel('Z')
        ylabel('Y')
        %For ZX
        subplot(szI(3),3,((c-1)*3)+3)
        imagesc(Bead.MIPs.Channel(c).ZX)
        set(gca,'YDir','normal')
        hold on
        scatter(CentrePoints(c,1),CentrePoints(c,3),'x','MarkerFaceColor',[1 1 0],'LineWidth',2)
        tn=sprintf('%s%d%s%d%s','Bead-',BeadID,' Chan-',c,' ZX');
        title(tn)
        xlabel('X')
        ylabel('Z')
        if c==szI(3)
            fn=sprintf('%s%s%d%s',foldername,'\Bead_',BeadID,'-EachChannelSeperately');
            savefig(fn)
        end
    end
end

%% Visualise the Centre Points
if strcmp(figures,'y')==1
    figure('WindowState','maximized')
    cmap=colormap(jet(szI(3)));
    % Make XY MIP
    subplot(2,2,1)
    for c=1:szI(3)
        tXY(:,:,c)=Bead.MIPs.Channel(c).XY;
    end
    [RGB_data]=mnl_ConvertToRGB(tXY);
    image(RGB_data)
    set(gca,'YDir','normal')
    hold on
    for c=1:szI(3)
        scatter(CentrePoints(c,2),CentrePoints(c,1),'x','MarkerFaceColor',cmap(c,:),'LineWidth',2)
        LegName{c}=sprintf('%s%d','Channel ',c);
    end
    title('XY MIP')
    xlabel('Y axis')
    ylabel('X axis')
    legend(LegName,'Location','northeastoutside')
    clear RGB_data
    % Make YZ MIP
    subplot(2,2,2)
    for c=1:szI(3)
        tYZ(:,:,c)=Bead.MIPs.Channel(c).YZ;
    end
    [RGB_data]=mnl_ConvertToRGB(tYZ);
    image(RGB_data)
    hold on
    for c=1:szI(3)
        scatter(CentrePoints(c,3),CentrePoints(c,2),'x','MarkerFaceColor',cmap(c,:),'LineWidth',2)
        LegName{c}=sprintf('%s%d','Channel ',c);
    end
    title('YZ MIP')
    xlabel('Z axis')
    ylabel('Y axis')
    legend(LegName,'Location','northeastoutside')
    clear RGB_data
    % Make ZX MIP
    subplot(2,2,3)
    for c=1:szI(3)
        tZX(:,:,c)=Bead.MIPs.Channel(c).ZX;
    end
    [RGB_data]=mnl_ConvertToRGB(tZX);
    image(RGB_data)
    hold on
    for c=1:szI(3)
        scatter(CentrePoints(c,1),CentrePoints(c,3),'x','MarkerFaceColor',cmap(c,:),'LineWidth',2)
        LegName{c}=sprintf('%s%d','Channel ',c);
    end
    title('ZX MIP')
    xlabel('X axis')
    ylabel('Z axis')
    legend(LegName,'Location','northeastoutside')
    fn=sprintf('%s%s%d',foldername,'\Bead_',BeadID);
    clear RGB_data
    savefig(fn)
end
%% Measure the distance to the points of each channel
for c=1:szI(3)
    Base=CentrePoints(c,:);
    CA(c).BasePosition=Base;
    for c2=1:szI(3)
        Points=CentrePoints(c2,:);
        CA(c).DistToChannel_px(c2,:)=Points-Base;
        CA(c).DistToChannel(c2,:)=(Points-Base).*Scale;
    end
end
end
%% Sub functions
function [RGB_data]=mnl_ConvertToRGB(Image)
%Makes an RGB image for upto 4 channels
dim=size(Image);
nChan=dim(3);
%Normalise each channel
for i=1:nChan
    Temp(:,:)=Image(:,:,i);
    mxT=max(Temp(:));
    nTemp=Temp./mxT;
    nImage(:,:,i)=nTemp;
    clear nTemp
end
%Decide what colours to pick
if nChan==2 %Green then Red
    ChanColour(1,:)=[0 1 0];
    ChanColour(2,:)=[1 0 0];
elseif  nChan==3 %Blue then Green then Red
    ChanColour(1,:)=[0 0 1];
    ChanColour(2,:)=[0 1 0];
    ChanColour(3,:)=[1 0 0];
elseif nChan==4 %Purple then Blue then Green then Red
    ChanColour(1,:)=[0.5 0 1];
    ChanColour(2,:)=[0 0 1];
    ChanColour(3,:)=[0 1 0];
    ChanColour(4,:)=[1 0 0];
end
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
end
function varargout = centerOfMass(A,varargin)
% CENTEROFMASS finds the center of mass of the N-dimensional input array
%
%   CENTEROFMASS(A) finds the gray-level-weighted center of mass of the
%   N-dimensional numerical array A. A must be real and finite. A warning
%   is issued if A contains any negative values. Any NaN elements of A will
%   automatically be ignored. CENTEROFMASS produces center of mass
%   coordinates in units of pixels. An empty array is returned if the
%   center of mass is undefined.
%
%   The center of mass is reported under the assumption that the first
%   pixel in each array dimension is centered at 1.
%
%   Also note that numerical arrays other than DOUBLE and SINGLE are
%   converted to SINGLE in order to prevent numerical roundoff error.
%
%   Examples:
%       A = rgb2gray(imread('saturn.png'));
%       C = centerOfMass(A);
%
%       figure; imagesc(A); colormap gray; axis image
%       hold on; plot(C(2),C(1),'rx')
%
%   See also: 
%
%

%
%   Jered R Wells
%   2013/05/07
%   jered [dot] wells [at] gmail [dot] com
%
%   v1.0
%
%   UPDATES
%       YYYY/MM/DD - jrw - v1.1
%
%

%% INPUT CHECK
narginchk(0,1);
nargoutchk(0,1);
fname = 'centerOfMass';

% Checked required inputs
validateattributes(A,{'numeric'},{'real','finite'},fname,'A',1);

%% INITIALIZE VARIABLES
A(isnan(A)) = 0;
if ~(strcmpi(class(A),'double') || strcmpi(class(A),'single'))
    A = single(A);
end
if any(A(:)<0)
    warning('MATLAB:centerOfMass:neg','Array A contains negative values.');
end

%% PROCESS
sz = size(A);
nd = ndims(A);
M = sum(A(:));
C = zeros(1,nd);
if M==0
    C = [];
else
    for ii = 1:nd
        shp = ones(1,nd);
        shp(ii) = sz(ii);
        rep = sz;
        rep(ii) = 1;
        ind = repmat(reshape(1:sz(ii),shp),rep);
        C(ii) = sum(ind(:).*A(:))./M;
    end
end

% Assemble the VARARGOUT cell array
varargout = {C};

end % MAIN