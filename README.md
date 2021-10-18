# ChromaticAberrationCorrection
Codes to measure and then correct chromatic aberration. 
These codes are written for use in MATLAB. NB You will need the Statistics and Machine Learning tooldbox for the regress function. But with some searching I believe you can find some others who have written the same function
Other credits... I've used the following codes from the MATLAB file exchange
UIpickfiles - by Douglas Schwarz to open multiple files in MATLAB https://jp.mathworks.com/matlabcentral/fileexchange/10867-uipickfiles-uigetfile-on-steroids

centerOfMass - a faster way of calculating the centre of mass than my own codes, by Jered Wells https://jp.mathworks.com/matlabcentral/fileexchange/41675-center-of-mass

bfopen - by the Open Microscope Environment team (bioformats) to open and read 4D images and obtain the metadata https://www.openmicroscopy.org/bio-formats/downloads/

ReadImageJROI - by Dylan Muir to read in ImageJ/FIJI roi selections for bead detections. https://jp.mathworks.com/matlabcentral/fileexchange/32479-readimagejroi


So how do you use it?
There are two main options
1) Create a template using bead data
- Image fluorescently labelled beads under the same settings that you wish to use for biological samples. If the z size (depth) of your block is to big you can create multiple sub-volumes instead (e.g. 0-25um, 50-75um, 100-125um,...,250-275um).
- Then run the code mnl_Pipeline_MeasurePSFandCA, you can either automatically detect the beads here or use the imports from the "Spot Detector" in ICY or ROIs from ImageJ. Remember to save the workspace (especially if you have multiple subvolumes).
- Then calculate the linear regression in Z with mnl_MergeBeadData - outputs are...
    AllBeads - structure containing the relevant PSF and chromatic aberration information for all beads
    ChromaticCorrections - the structure containing the regression data used to correct future images.
    Scale - The original scale of the volume, NB you need to delete this from the workspace before running mnl_CorrectForChromaticAberration otherwise the scale information will be applied to the image to be corrected
- Then to create a corrected image, run mnl_CorrectForChromaticAberration, you will be prompted to load both the image and the ChromaticCorrections structure. Each channel will be produced as a separate 3D tiff stack saved to your current directory.


2) Use guide stars
- Firstly in select your guide stars and crop them out as individual images (in the paper we used FIJI). Make sure you have the correct scale loaded in the metadata, and also record the z range that you have obtained the cropped sections from
- Then use mnl_CalculateChromaticAberrationsFromGuideStars (v2 - saves a figure showing the centre of mass for each guidestar as an EPS and .fig file). Load in each guide star and then you will be asked to input the starting z plane of each crop. 
  Outputs are...
  ChromaticCorrections - the structure containing the regression data used to correct future images. There are also metrics measuring the quality of the regressioni
  ROIs-The chromatic aberration of each guide star and it's Z depth
  - Remember to save the workspace
- Finally create a corrected image, run mnl_CorrectForChromaticAberration, you will be prompted to load both the image and the ChromaticCorrections structure. Each channel will be produced as a separate 3D tiff stack saved to your current directory.
  
Any more questions? Feel free to let us know!
