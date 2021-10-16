function [Beads]=mnl_ExtractBeadInfoFromICY_temp(fname)
% Function to extract the bead information 
%% Pre-allocate the structure
Beads=struct('BeadLocation',[],'BeadSize',[],'BeadCentre',[],'BeadExtremes',[],'ImageLimits',[]);
%% Info from excel sheet
[~,txt,~]=xlsread(fname);
NumSpots=size(txt,1)-1;
for i=1:NumSpots
    %Pixel Locations- NB Switch the X and Y axes around for MATLAB, this
    %will be the opposite of what the excell file says
    Ypos(i)=str2num(txt{i+1,3});
    Xpos(i)=str2num(txt{i+1,4});
    Zpos(i)=str2num(txt{i+1,5});
    %Pixel Size
    Ysize(i)=str2num(txt{i+1,8});
    Xsize(i)=str2num(txt{i+1,9});
    Zsize(i)=str2num(txt{i+1,10});
    % Bead Extremes
    Xextremes(i,:)=[Xpos(i) Xpos(i)+Xsize(i)];
    Yextremes(i,:)=[Ypos(i) Ypos(i)+Ysize(i)];
    Zextremes(i,:)=[Zpos(i) Zpos(i)+Zsize(i)];
    %CentrePoints
    Ycent(i)=round(str2num(txt{i+1,13}));
    Xcent(i)=round(str2num(txt{i+1,14}));
    Zcent(i)=round(str2num(txt{i+1,15}));
    %ImageLimits
    X_ImLim(i,:)=[Xcent(i)-Xsize(i) Xcent(i)+Xsize(i)];
    Y_ImLim(i,:)=[Ycent(i)-Ysize(i) Ycent(i)+Ysize(i)];
    Z_ImLim(i,:)=[Zcent(i)-(Zsize(i)*2) Zcent(i)+(Zsize(i)*2)]; % Add an extra 150% to control for the increased chromatic abberation
end

for i=1:NumSpots
    Beads(i).BeadLocation=[Xpos(i) Ypos(i) Zpos(i)];
    Beads(i).BeadSize=[Xsize(i) Ysize(i) Zsize(i)];
    Beads(i).BeadCentre=[Xcent(i) Ycent(i) Zcent(i)];
    Beads(i).BeadExtremes=[Xextremes(i,1) Xextremes(i,2);Yextremes(i,1) Yextremes(i,2);Zextremes(i,1) Zextremes(i,2)];
    Beads(i).ImageLimits=[X_ImLim(i,1) X_ImLim(i,2);Y_ImLim(i,1) Y_ImLim(i,2);Z_ImLim(i,1) Z_ImLim(i,2)];
end
end