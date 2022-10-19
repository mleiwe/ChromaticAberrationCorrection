## Python version

A simple Python port of the original MATLAB chromatic aberration correction codes. This version is not intended to replace the MATLAB one, rather it's an additional option. Recommended to use MATLAB codes as its more optimized and more "plug and play" and has more options. 

## Publication:
![Screenshot 2022-10-19 160745](https://user-images.githubusercontent.com/29883365/196620792-3c50156e-451d-4cd7-baf7-b9c475638cb3.png)

DOI: [10.3389/fnana.2021.760063](https://www.frontiersin.org/articles/10.3389/fnana.2021.760063/full)


## Dependencies

Following the python libraries needed. All packages are installed using **pip** from pypi website.

    matplotlib==3.5.1
	napari==0.4.12
	numpy==1.21.4
	pandas==1.3.5
	scikit_image==0.19.0
	scikit_learn==1.0.2
	scipy==1.7.3
	seaborn==0.11.2
	skimage==0.0
	tifffile

## Usage

Code is divided into two parts: `Chromatic_aberration_calculation.py` and `Chromatic_aberration_correction.py` .  Currently, it will only work with ICY spot detection results.xls

 - `Chromatic_abberation_calculation.py` -> Read the image, read the ICY spot detector .xls files, calculate the regressions, draw the regression lines, save the regression values as a csv.
 - `Chromatic_abberation_correction.py` -> Read the regression values csv file, read the image to be corrected, make corrections, save the image as tiff(32 bit) `tifffile_imwrite` loves to save as a 32bit image. Saved image is ImageJ compatible and metadata can also be embedded while saving the image.

*Optional napari codes also have been added for quick visualization using the powerful n-D napari viewer.*

***For more detailed explanation please refer to the original paper*** 

## Contributing:
Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are  **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement". Don't forget to give the project a star! Thanks again!

1.  Fork the Project
2.  Create your Feature Branch 
3.  Commit your Changes
4.  Push to the Branch 
5.  Open a Pull Request

## License:
Distributed under GPL-3.0 license.

## Acknowledgements:
Credits to the original authors of the paper:
1. [Marcus Lewie](https://loop.frontiersin.org/people/1503524/overview)
2. [Satoshi Fujimoto](https://loop.frontiersin.org/people/1445674/overview)
3. [Takeshi Imai](https://loop.frontiersin.org/people/1036459/overview)

*Department of Developmental Neurophysiology, Graduate School of Medical Sciences, Kyushu University, Fukuoka, Japan*

### How to cite:
If you have used this repository please cite the original paper as: 

*Leiwe MN, Fujimoto S and Imai T (2021) _Post hoc_ Correction of Chromatic Aberrations in Large-Scale Volumetric Images in Confocal Microscopy. _Front. Neuroanat._ 15:760063. doi: 10.3389/fnana.2021.760063*

