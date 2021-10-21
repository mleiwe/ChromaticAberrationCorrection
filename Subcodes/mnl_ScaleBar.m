function []=mnl_ScaleBar(scale,dist,location,Label)
% draws a scale bar on the current figure
% Inputs
% scale - distance per pixel, if blank default is 1
% dist - how big do you want the scale bar to be
% location - what position on the figure do you want the scale bar (e.g.
% 'northeast','northwest','southeast','southwest')
% Label - text to go underneath or above the scale bar

%% Establish Parameters
%Scale
chk=isempty(scale);
if chk==1
    scale=1;
end
chk=isempty(Label);
if chk==1
    Label=' ';
end
%Distance
TruDist=dist/scale;

%%   GET IMAGE AND AXES DATA
hAxes = gca;
axeslims = [get(hAxes,'xlim')' get(hAxes,'ylim')'];
axesdir=[1,1];

%% Determine the location of vector for the scale bar
if ischar(location)
    switch location
        case 'northeast'
            anchor = [axeslims(2,1) - axesdir(1)*range(axeslims(:,1))*0.05, ... %x axis position 5% from edge
                axeslims(2,2) - axesdir(2)*range(axeslims(:,2))*0.075]; %y axis position 5% from edge
            anchor2 = [axeslims(2,1) - axesdir(1)*range(axeslims(:,1))*0.1, ... %x axis position 10% from edge
                axeslims(2,2) - axesdir(2)*range(axeslims(:,2))*0.1]; %y axis position 10% from edge
            % Now Calculate Positions
            Yvec=[anchor(2),anchor(2)];
            Xvec=[anchor(1)-TruDist;anchor(1)];
            Yvec2=[anchor2(2),anchor2(2)+30];
            Xvec2=[anchor2(1)-TruDist;anchor2(1)];
        case 'northwest'
            anchor = [axeslims(1,1) + axesdir(1)*range(axeslims(:,1))*0.05, ...  %x axis position 5% from edge
                axeslims(2,2) - axesdir(2)*range(axeslims(:,2))*0.075]; %y axis position 5% from edge
            anchor2 = [axeslims(1,1) + axesdir(1)*range(axeslims(:,1))*0.1, ...  %x axis position 10% from edge
                axeslims(2,2) - axesdir(2)*range(axeslims(:,2))*0.1]; %y axis position 10% from edge
            % Now Calculate Positions
            Yvec=[anchor(2),anchor(2)];
            Xvec=[anchor(1),anchor(1)+TruDist];
            Yvec2=[anchor2(2),anchor2(2)+30];
            Xvec2=[anchor(1),anchor(1)];
        case 'southwest'
            anchor = [axeslims(1,1) + axesdir(1)*range(axeslims(:,1))*0.05, ... %x axis position 5% from edge
                axeslims(1,2) + axesdir(2)*range(axeslims(:,2))*0.075]; %y axis position 5% from edge
            anchor2 = [axeslims(1,1) + axesdir(1)*range(axeslims(:,1))*0.10, ... %x axis position 10% from edge
                axeslims(1,2) + axesdir(2)*range(axeslims(:,2))*0.10]; %y axis position 10% from edge
            % Now Calculate Positions
            Xvec=[anchor(1),anchor(1)+TruDist];
            Yvec=[anchor(2),anchor(2)];
            Xvec2=[anchor2(1),anchor2(1)-TruDist];
            Yvec2=[anchor2(2),anchor2(2)-30];
        case 'southeast'
            anchor = [axeslims(2,1) - axesdir(1)*range(axeslims(:,1))*0.05, ... %x axis position 5% from edge
                axeslims(1,2) + axesdir(2)*range(axeslims(:,2))*0.075]; %y axis position 5% from edge
            anchor2 = [axeslims(2,1) - axesdir(1)*range(axeslims(:,1))*0.1, ... %x axis position 5% from edge
                axeslims(1,2) + axesdir(2)*range(axeslims(:,2))*0.1]; %y axis position 5% from edge
            % Now Calculate Positions
            Yvec=[anchor(2),anchor(2)];
            Xvec=[anchor(1)-TruDist;anchor(1)];
            Yvec2=[anchor2(2),anchor2(2)-30];
            Xvec2=[anchor2(1)+TruDist;anchor2(1)];
    end
else
    anchor = location;
    if location
        dirToCentre = min(axeslims)+range(axeslims)/2 - location.*axesdir;
        direction = directions{ceil((-1*atan2(dirToCentre(2),dirToCentre(1))+pi)/(2*pi)*4)};
    end
end

line(Xvec,Yvec,'Color','w','LineWidth',2)
text(mean(Xvec2),mean(Yvec2),Label,'Color','w','FontSize',14)

