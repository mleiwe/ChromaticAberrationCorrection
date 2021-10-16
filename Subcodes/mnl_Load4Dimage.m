function [Data,Scale,dim,metadata]=mnl_Load4Dimage
%Function to load the image, and relevant OME metadata using bfopen
% Outputs
% Data - The image as a matrix (x*y*c*z)
% Scale - the size of a voxel (x*y*z)
% dim - the dimensions of the image
% metadata - the OME metadata

file=bfopen;
series1=file{1, 1};
metadata=file{1,4};
% ID image dimensions
x_num=metadata.getPixelsSizeX(0).getValue();
y_num=metadata.getPixelsSizeY(0).getValue();
z_num=metadata.getPixelsSizeZ(0).getValue();
c_num=metadata.getPixelsSizeC(0).getValue();
dim=[x_num y_num c_num z_num];
%Scale
x=double(metadata.getPixelsPhysicalSizeX(0).value());
y=double(metadata.getPixelsPhysicalSizeY(0).value());
z=double(metadata.getPixelsPhysicalSizeZ(0).value());
Scale=[x y z];
% Extract Data
%Data=nan(y_num,x_num,c_num,z_num);
t=1;
for i=1:z_num
    for i2=1:c_num
        a=series1{t,1};
        Data(:,:,i2,i)=a;
        t=t+1;
    end
end
end