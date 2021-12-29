## Python version of the Chromatic abberation codes

A simple Python version of the MATLAB codes. Python has some disadvantages and so the MATLAB codes takes precedence but the Python codes provide a starting point to implement this making it more open source friendly.

## Depedencies:
All the packages can be installed with pip or pip3(for Python3) [pip](https://pip.pypa.io/en/stable/)
Included in the requirements.txt
```bash
matplotlib==3.5.1
napari==0.4.12
numpy==1.21.4
pandas==1.3.5
scikit_image==0.19.0
scikit_learn==1.0.2
scipy==1.7.3
seaborn==0.11.2
skimage==0.0
```

## Usage
The folder consists of two scripts `chromatic_abb_measure.py` and `abberation_correction.py`. The first one does all the image processing and regression calculation and the 2nd one does the corrections based on those calculations. Make sure to keep both of them in the same directory. All the codes are efficiently commented.

The code uses the current distribution of [napari](https://napari.org/), highly interactive image viewer built with Python. 
```python
pip install napari[all]

```
<img src="https://user-images.githubusercontent.com/29883365/147657273-b30f9746-30d4-4bca-ad6f-a85fc767d784.png" width="45%"></img> 

The viewer looks like this, and is highly interactive in nature and support n-dim viewing. The images shown are multichannel bead images, and the colors are connected-component labelling (`skimage.measure.label`) obtained using `chromatic_abb_measure.py` script.

## Comments:
This code as of now, works perfectly fine but there definitely room for improvement. As of now, the code is easy to edit, and image processing functions can be tuned based on the user's actions. Also, it can be bit slower if the image file is too big. Also, ***be careful of the image axis notations, the input image will have the order (Z,Y,X,C) and the numpy indexing starts from 0.***

## Contributing:
Good Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change. No hostile forks. 

## License
[MIT](https://choosealicense.com/licenses/mit/)
