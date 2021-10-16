function [AllBeads,ChromaticCorrections,Scale]=mnl_MergeBeadData 
%% Input Channel Data
prompt='How many channels are used for imaging?';
nChan=input(prompt);
LaserEx=nan(nChan,1);
for i=1:nChan
    LaserQ=sprintf('%s%d%s','Please enter the excitation laser wavelength for channel ',i,' (nm)');
    prompt=LaserQ;
    LaserEx(i)=input(prompt);
end

%% Merge fData
[Wkspaces]=uipickfiles;
sz=size(Wkspaces,2);
% Create a matrix - 1st
% Col=CentDist,2nd=Zdist,3rd=XAb_Ch1,4th=YAb_Ch1,5th=ZAb_Ch1,6th=XAb_Ch2,7th=YAb_Ch2,8th=ZAb_Ch2,9th=XAb_Ch3,10th=YAb_Ch3,11th=ZAb_Ch3,12th=XAb_Ch3,13th=YAb_Ch3,14th=ZAb_Ch3,15th=XAb_Ch4,16th=YAb_Ch4,17th=ZAb_Ch4
n=1;
for i=1:sz
    load(Wkspaces{1,i});
    NumBeads=size(f2Beads,2);
    n3=n+NumBeads;
    n2=n3-1;
    AllBeads(n:n2)=f2Beads;
    n=n3;
end
TotNumBeads=n2;
nComps=nChan*(nChan-1);
for k=1:nChan
    Matrix=nan(n2,4+nComps);%Matrix2=nan(n2,16);Matrix3=nan(n2,16);Matrix4=nan(n2,16);
    % To Channel k
    for i=1:n2
        Matrix(i,1)=AllBeads(i).DistFromCentre;
        Matrix(i,2)=AllBeads(i).XDistFromCentre;
        Matrix(i,3)=AllBeads(i).YDistFromCentre;
        Matrix(i,4)=AllBeads(i).Zdepth;
        for j=1:nChan
            %Chan j
            Matrix(i,4+(((j-1)*3)+1))=AllBeads(i).ChromaticAberration(k).DistToChannel(j,1);
            Matrix(i,4+(((j-1)*3)+2))=AllBeads(i).ChromaticAberration(k).DistToChannel(j,2);
            Matrix(i,4+(((j-1)*3)+3))=AllBeads(i).ChromaticAberration(k).DistToChannel(j,3);
        end
    end
    if k==1
        fign=sprintf('%s%d','Chromatic Aberration To Ch',k);
        figure('Name',fign)
        cmap=colormap(jet(nChan));
        subplot(1,3,1) %X Aberration
        for j=1:nChan
            scatter3(Matrix(:,1),Matrix(:,4),Matrix(:,4+(((j-1)*3)+1)),'MarkerFaceColor',cmap(j,:)) %Chj
            hold on
            legnames{j}=sprintf('%s%d','Channel ',j);
        end
        legend(legnames)
        xlabel('Distance From Centre')
        ylabel('Z depth')
        zlabel('Shift')
        title('X shift')
        subplot(1,3,2) %Y Aberration
        for j=1:nChan
            scatter3(Matrix(:,1),Matrix(:,4),Matrix(:,(((j-1)*3)+2)),'MarkerFaceColor',cmap(j,:)) %Ch1
            hold on
        end
        legend(legnames)
        xlabel('Distance From Centre')
        ylabel('Z depth')
        zlabel('Shift')
        title('Y shift')
        subplot(1,3,3) %Z Aberration
        for j=1:nChan
            scatter3(Matrix(:,1),Matrix(:,4),Matrix(:,(((j-1)*3)+3)),'MarkerFaceColor',cmap(1,:)) %Ch1
            hold on
        end
        legend(legnames)
        xlabel('Distance From Centre')
        ylabel('Z depth')
        zlabel('Shift')
        title('Z shift')
    end
    %Now plot each xyz distance separately
    fign2=sprintf('%s%d','To Channel ',k);
    figure('Name',fign2)
    mnl_PlotAllChromaticAberrationsXYZ(Matrix,LaserEx)
    %% Calculate the linear regression along Z to Channel k
    xRange=Matrix(:,4);
    for j=1:nChan
        yRange=Matrix(:,4+(((j-1)*3)+3));%Ch1
        [m(j),c(j)]=mnl_CalculateRegression(xRange,yRange);        
        WhichLaser(j)=LaserEx(j);
    end
    ChromaticCorrections(k).ToWhichLaser=LaserEx(k);
    ChromaticCorrections(k).ForWhichLaser=WhichLaser;
    ChromaticCorrections(k).mValues=m;
    ChromaticCorrections(k).cValues=c;
    %% Now just Zdrift and dist to centre as 2D plots
    fign3=sprintf('%s%d','To Laser ',k);
    figure('Name',fign3)
    mValues=ChromaticCorrections(k).mValues;
    cValues=ChromaticCorrections(k).cValues;
    mnl_PlotAllChromaticAberrations(Matrix,LaserEx,mValues,cValues)
    mnl_PlotAllChromaticAberrationsPerChanPerZdepth(Matrix,LaserEx,k,100,25)
end
end
%% Subfunctions
function mnl_PlotAllChromaticAberrationsXYZ(Matrix,LaserEx)
szL=size(LaserEx,1);
cmap=colormap(jet(szL));
for i=1:szL
    LegNames{i}=sprintf('%d%s',LaserEx(i),'nm');
end
tylims=[0 0];
%Row Xshift of bead
xvals=Matrix(:,2);
%Xdist and Xshift
subplot(3,3,1)
for i=1:szL
    yvalues1=Matrix(:,4+((i-1)*3)+1);
    plot(xvals,yvalues1,'.','color',cmap(i,:))
    hold on
    clear yvalues1
end
xlabel('X distance to the centre (um)')
ylabel('X shift of bead (um)')
%Now check the min and max scales of subplots
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end
%Xdist and Yshift
subplot(3,3,2)
for i=1:szL
    yvalues=Matrix(:,4+((i-1)*3)+2);
    plot(xvals,yvalues,'.','color',cmap(i,:))
    hold on
    clear yvalues1
end
xlabel('X distance to the centre (um)')
ylabel('Y shift of bead (um)')
%Now check the min and max scales of subplots
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end
%Xdist and Zshift
subplot(3,3,3)
for i=1:szL
    yvalues1=Matrix(:,4+((i-1)*3)+3);
    plot(xvals,yvalues1,'.','color',cmap(1,:))
    hold on
    clear yvalues1
end
set(gca, 'YDir','reverse')
xlabel('X distance to the centre (um)')
ylabel('Z shift of bead (um)')
%Now check the min and max scales of subplots
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end

%Row Yshift of bead
xvals=Matrix(:,3);
%Ydist and Xshift
subplot(3,3,4)
for i=1:szL
    yvalues1=Matrix(:,4+((i-1)*3)+1);
    plot(xvals,yvalues1,'.','color',cmap(i,:))
    hold on
    clear yvalues1
end
xlabel('Y distance to the centre (um)')
ylabel('X shift of bead (um)')
%Now check the min and max scales of subplots
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end
%Ydist and Yshift
subplot(3,3,5)
for i=1:szL
    yvalues1=Matrix(:,4+((i-1)*3)+2);
    plot(xvals,yvalues1,'.','color',cmap(i,:))
    hold on
    clear yvalues1
end
xlabel('Y distance to the centre (um)')
ylabel('Y shift of bead (um)')
%Now check the min and max scales of subplots
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end
%Ydist and Zshift
subplot(3,3,6)
for i=1:szL
     yvalues1=Matrix(:,4+((i-1)*3)+3);
    plot(xvals,yvalues1,'.','color',cmap(i,:))
    hold on
    clear yvalues1
end
set(gca, 'YDir','reverse')
xlabel('Y distance to the centre (um)')
ylabel('Z shift of bead (um)')
%Now check the min and max scales of subplots
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end

%Row Zshift of bead
xvals=Matrix(:,4);
%Ydist and Xshift
subplot(3,3,7)
for i=1:szL
    yvalues1=Matrix(:,4+((i-1)*3)+1);
    plot(xvals,yvalues1,'.','color',cmap(i,:))
    hold on
    clear yvalues1
end
xlabel('Z distance to the surface (um)')
ylabel('X shift of bead (um)')
%Now check the min and max scales of subplots
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end
%Ydist and Yshift
subplot(3,3,8)
for i=1:szL
    yvalues1=Matrix(:,4+((i-1)*3)+2);
    plot(xvals,yvalues1,'.','color',cmap(i,:))
    hold on
    clear yvalues1
end
xlabel('Z distance to the surface (um)')
ylabel('Y shift of bead (um)')
%Now check the min and max scales of subplots
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end
%Ydist and Zshift
subplot(3,3,9)
for i=1:szL
    yvalues1=Matrix(:,4+((i-1)*3)+3);
    plot(xvals,yvalues1,'.','color',cmap(i,:))
    hold on
    clear yvalues1
end
set(gca, 'YDir','reverse')
xlabel('Z distance to the surface (um)')
ylabel('Z shift of bead (um)')
%Now check the min and max scales of subplots
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end
legend(LegNames)
%% Now rescale all y axis
ylimdist=tylims(2)-tylims(1);
ylimrange=ylimdist*1.05; %Add 5% to the range
ylimmid=tylims(1)+(0.5*ylimdist);
fylims=[ylimmid-(0.5*ylimrange) ylimmid+(0.5*ylimrange)];
for i=1:9
    subplot(3,3,i)
    ylim(fylims)
end
end
function mnl_PlotAllChromaticAberrations(Matrix,LaserEx,mValues,cValues)
szL=size(LaserEx,1);
cmap=colormap(jet(szL));
for i=1:szL
    LegNames{((i-1)*2)+1}=sprintf('%d%s',LaserEx(i),'nm');
    LegNames{((i-1)*2)+2}=sprintf('%d%s',LaserEx(i),'nm - Linear fit');
end
tylims=[0 0];
%Row distance to centre
xvals=Matrix(:,1);
%XYdist and XY shift
subplot(2,2,1)
%Get XYshift values
for i=1:szL
    %To Chan1
    tXval=Matrix(:,4+((i-1)*3)+1);
    tYval=Matrix(:,4+((i-1)*3)+2);
    yvalues1=sqrt((tXval.^2)+(tYval.^2));
    % %Now plot it
    plot(xvals,yvalues1,'.','color',cmap(i,:))
    hold on
    clear yvalues1
end
xlabel('XY distance to the centre of image (um)')
ylabel('XY shift of bead (um)')
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end

%XYdist and Zshift
subplot(2,2,2)
for i=1:szL
    yvalues1=Matrix(:,4+((i-1)*3)+3);
    plot(xvals,yvalues1,'.','color',cmap(i,:))
    hold on
    clear yvalues1
end
set(gca, 'YDir','reverse')
xlabel('XY distance to the centre of image (um)')
ylabel('Z shift of bead (um)')
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end

%Row Z depth
xvals=Matrix(:,4);
%Zdist and XYshift
subplot(2,2,3)
for i=1:szL
    %Get XYshift values
    %To Chan1
    tXval=Matrix(:,4+((i-1)*3)+1);
    tYval=Matrix(:,4+((i-1)*3)+2);
    yvalues1=sqrt((tXval.^2)+(tYval.^2));
    %Now plot it
    plot(xvals,yvalues1,'.','color',cmap(i,:))
    hold on
    clear yvalues1
end
xlabel('Z distance to the surface (um)')
ylabel('XY shift of bead (um)')
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end

%Zdist and Zshift
subplot(2,2,4)
for i=1:szL
    yvalues1=Matrix(:,4+((i-1)*3)+3);
    plot(xvals,yvalues1,'.','color',cmap(i,:))
    hold on
    clear yvalues1
    %Now fit regression lines    
    tx=xlim;
    fit_xvals=tx(1):tx(2);
    yvalues1=(mValues(i)*fit_xvals)+cValues(i);
    plot(fit_xvals,yvalues1,'-','color',cmap(i,:))
end
set(gca, 'YDir','reverse')
xlabel('Z distance to the surface (um)')
ylabel('Z shift of bead (um)')
tmplim=ylim;
if tmplim(1)<tylims(1)
    tylims(1)=tmplim(1);
end
if tmplim(2)>tylims(2)
    tylims(2)=tmplim(2);
end
legend(LegNames)
%% Now rescale all y axis
ylimdist=tylims(2)-tylims(1);
ylimrange=ylimdist*1.05; %Add 5% to the range
ylimmid=tylims(1)+(0.5*ylimdist);
fylims=[ylimmid-(0.5*ylimrange) ylimmid+(0.5*ylimrange)];
for i=1:4
    subplot(2,2,i)
    ylim(fylims)
end
end
function mnl_PlotAllChromaticAberrationsPerChanPerZdepth(Matrix,LaserEx,RefChan,ZintSize,XYintSize)
%%Initial Settings
szL=size(LaserEx,1);
%Caluclate the z intervals
AllDepths(:,1)=Matrix(:,4);
MaxDepth=max(AllDepths);
nInts=ceil(MaxDepth/ZintSize);
BiggestLimit=nInts*ZintSize;
Zints=0:ZintSize:BiggestLimit;
%Calculate the xy intervals
AllXYdists=Matrix(:,1);
MaxXY=max(AllXYdists);
nXYints=ceil(MaxXY/XYintSize);
BiggestLimit=nXYints*XYintSize;
XYints=0:XYintSize:BiggestLimit;

%% Now the graphs
for h=1:szL %For each other channel to the reference
    if h~=RefChan
        fign=sprintf('%s%d%s%d%s','Reference LaserEx ',LaserEx(RefChan),'Distance to LaserEx',LaserEx(h),'nm');
        figure('Name',fign);
        tylims=[0 0];
        %Per Zint
        cmap=colormap(jet(nInts));
        for i=1:nInts
            LegNames{i}=sprintf('%s%d%s%d%s','Depth ',Zints(i),'um to ',Zints(i+1),'um');
        end
        
        for j=1:nInts
            idx=Matrix(:,4)>=Zints(j) & Matrix(:,4)<Zints(i+1);
            tMatrix=Matrix(idx,:);
            %Row distance to centre
            xvals=tMatrix(:,1);
            %XYdist and XY shift
            subplot(2,2,1)
            %Get XYshift values
            %To Chan h
            tXval=tMatrix(:,4+((h-1)*3)+1);
            tYval=tMatrix(:,4+((h-1)*3)+2);
            yvalues1=sqrt((tXval.^2)+(tYval.^2));
            % %Now plot it
            plot(xvals,yvalues1,'.','color',cmap(j,:),'MarkerSize',10)
            hold on
            clear yvalues1
            if j==nInts
                %Label axis
                xlabel('XY distance to the centre of image (um)')
                ylabel('XY shift of bead (um)')
                %Find the axes range
                tmplim=ylim;
                if tmplim(1)<tylims(1)
                    tylims(1)=tmplim(1);
                end
                if tmplim(2)>tylims(2)
                    tylims(2)=tmplim(2);
                end
            end
            
            
            %XYdist and Zshift
            subplot(2,2,2)
            yvalues1=tMatrix(:,4+((h-1)*3)+3);
            plot(xvals,yvalues1,'.','color',cmap(j,:),'MarkerSize',10)
            hold on
            clear yvalues1
            if j==nInts
                set(gca, 'YDir','reverse')
                xlabel('XY distance to the centre of image (um)')
                ylabel('Z shift of bead (um)')
                tmplim=ylim;
                if tmplim(1)<tylims(1)
                    tylims(1)=tmplim(1);
                end
                if tmplim(2)>tylims(2)
                    tylims(2)=tmplim(2);
                end
                legend(LegNames)
            end
        end
        %Row Z depth
        %Per XY int
        cmap=colormap(jet(nXYints));
        for j=1:nXYints
            LegNames{j}=sprintf('%s%d%s%d%s','Dist to centre ',XYints(j),'um to ',XYints(j+1),'um');
            idx=Matrix(:,1)>=XYints(j) & Matrix(:,1)<XYints(j+1);
            tMatrix=Matrix(idx,:);
            xvals=tMatrix(:,4);
            %Zdist and XYshift
            subplot(2,2,3)
            %Get XYshift values
            %To Chan1
            tXval=tMatrix(:,4+((h-1)*3)+1);
            tYval=tMatrix(:,4+((h-1)*3)+2);
            yvalues1=sqrt((tXval.^2)+(tYval.^2));
            %Now plot it
            plot(xvals,yvalues1,'.','color',cmap(j,:),'MarkerSize',10)
            hold on
            clear yvalues1
            if j==nXYints
                xlabel('Z distance to the surface (um)')
                ylabel('XY shift of bead (um)')
                tmplim=ylim;
                if tmplim(1)<tylims(1)
                    tylims(1)=tmplim(1);
                end
                if tmplim(2)>tylims(2)
                    tylims(2)=tmplim(2);
                end
            end
            
            %Zdist and Zshift
            subplot(2,2,4)
            yvalues1=tMatrix(:,4+((h-1)*3)+3);
            plot(xvals,yvalues1,'.','color',cmap(j,:),'MarkerSize',15)
            hold on
            clear yvalues1
            if j==nXYints
                set(gca, 'YDir','reverse')
                xlabel('Z distance to the surface (um)')
                ylabel('Z shift of bead (um)')
                tmplim=ylim;
                if tmplim(1)<tylims(1)
                    tylims(1)=tmplim(1);
                end
                if tmplim(2)>tylims(2)
                    tylims(2)=tmplim(2);
                end
            end
        end
        legend(LegNames)
        %% Now rescale all y axis
        ylimdist=tylims(2)-tylims(1);
        ylimrange=ylimdist*1.05; %Add 5% to the range
        ylimmid=tylims(1)+(0.5*ylimdist);
        fylims=[ylimmid-(0.5*ylimrange) ylimmid+(0.5*ylimrange)];
        for i=1:4
            subplot(2,2,i)
            ylim(fylims)
        end
    end
end
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
end