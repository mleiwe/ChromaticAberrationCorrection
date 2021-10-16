function [Data2]=mnl_ChromaticAbberationCorrectionZ4(Data,Scale,mValues,cValues)
%Function to correct for the chromatic abberation, based on bead data
%extracted from the codes 'mnl_Pipeline_MeasurePSFandCA' and
%'mnl_MergeBeadData'
%
% Inputs
% Data - Original Data (x*y*c*z)
% Scale - x*y*z Scale of the image
% WhichLasers - 1*numChannels of which excitation laser was used
% mValues - 1*numChannels of the corresponding m values for the equation
% y=mx+c
% cValues - 1*numChannels of the corresponding c values for the equation
% y=mx+c
%
% Outputs
% Data2 - corrected image
sz=size(Data);
%% Step 1 - Calculate the Z-range
prompt='Please input the starting depth of the image relative to the surface (this code assumes the lens will go deeper)';
Zstart=input(prompt);
Zrange_um=linspace(Zstart,Zstart+(sz(4)-1)*Scale(3),sz(4));
%% Step 2 - Calculate the chromatic abberation equations for each channel
Zshifts=nan(sz(3),sz(4));
for i=1:sz(3) %For each channel
    Zshifts(i,:)=(Zrange_um.*mValues(i))+cValues(i)+Zrange_um;
end
minZrange=min(Zshifts(:));
maxZrange=max(Zshifts(:));
newRangeSize=maxZrange-minZrange;
Range2=Zstart:Scale(3):Zstart+newRangeSize;
for i=1:sz(3) %For each channel
    Zshifts2(i,:)=(Range2.*mValues(i))+cValues(i);
end
%% Step 3 - Translate these into Zshifts
Nplanes=size(Range2,2);
Zshifts_vx=round(Zshifts2./Scale(3));
%Now re-code to the one that has to be shifted up the most
StartPoint=max(Zshifts_vx(:,1)); %Biggest shift up
EndPoint=min(Zshifts_vx(:,Nplanes));
cZshifts_vx=Zshifts_vx-StartPoint;
%% Now place in the corrections
disp('Adjusting the Image...')
Data2=Data.*0; %Just for preallocation purposes
for i=1:Nplanes
    for j=1:sz(3)
        Zplane=i+cZshifts_vx(j,i);
        if Zplane>=1 && Zplane<=sz(4)
            Data2(:,:,j,i)=Data(:,:,j,Zplane);
        end
    end
    mnl_InsertProgressTrackerInLoops(i,Nplanes);
end

end