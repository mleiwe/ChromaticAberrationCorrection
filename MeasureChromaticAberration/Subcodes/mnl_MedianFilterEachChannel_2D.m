function [mData]=mnl_MedianFilterEachChannel_2D(Data,r)
%Function to perform a 2D median filter on a 3D stack of each channel (i.e.
%a 4D image)

%Re-arrange to split channels into individual Z stacks
if isempty(r)==1
    r=3;
end
dim=size(Data);
DataSize_pts=dim(1)*dim(2)*dim(3)*dim(4);
if isa(Data,'uint16')==1
    DataSizeBits=DataSize_pts*16;
    DataSizeBytes=DataSizeBits/8;
    DataSizeGB=DataSizeBytes/1073741824;
else
    DataSizeBytes=DataSize_pts; %Default assumption is it is 8 bit
    DataSizeGB=DataSizeBytes/1073741824;
end
if DataSizeGB>5
    display('Warning this image is large it may take some time to process')
end
fprintf('%s\n','Splitting Each Z stack')
Chan=struct('Zstack',[]);
for i=1:dim(3)
    temp(:,:,:)=Data(:,:,i,:);
    Chan(i).Zstack=double(temp);
end
clear temp Data
% Median Filter
fprintf('%s\n','Median Filtering')
mData=nan(dim);
for i=1:dim(3)
    a=single(Chan(i).Zstack);
    a2=medfilt3(a,[r r 1]);
%     for j=1:dim(4)
%         A=a(:,:,j);
%         a2(:,:,j)=nanmedfilt2(A,r);
%         %[a2(:,:,j)]=mnl_NanMedFilt2D(a(:,:,j),r);
%     end
    Chan(i).MedZstack=a2;
    clear a a2
end
%Now Make the mData file
fprintf('%s\n','Putting the Image Back Together')
for i=1:dim(3)
    mData(:,:,i,:)=Chan(i).MedZstack;
end
end