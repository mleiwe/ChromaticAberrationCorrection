#!/usr/bin/env python
# coding: utf-8




import numpy as np
import matplotlib.pyplot as plt
from skimage.io import*
from scipy import ndimage as ndi
import napari
from skimage.morphology import*
from skimage.segmentation import*
from skimage.measure import*
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.metrics import*
from matplotlib.figure import*
import pandas as pd
import os
from skimage.filters import median
from skimage.color import label2rgb
import matplotlib.patches as mpatches
from sklearn.metrics import r2_score
import pandas as pd
import warnings; warnings.simplefilter('ignore')
sns.set_style('ticks')
import tifffile as tiff





def chromatic_abberation_correction(img,regression_file_path,zscale,zstart):
    """
    This function corrects chromatic abberation in the image.
    Input:
    img: image to be corrected
    regression_file: outfile of chromatic abberation calculation
    zscale: pixel size along z axis
    zstart: starting z position of the image(assuming the lens go deeper)
    Output:
    img_corrected: image with corrected chromatic abberation

    Change the laser powers based on your imaging conditions
    """

    #load the regression file
    reg_vals = pd.read_csv(regression_file_path)
    zrange_um = np.linspace(zstart,zstart+(img.shape[0]-1)*zscale,img.shape[0])
    base_laser = (input('Enter the base laser wavelength (in nm): '))

    ### correct chromatic abberation based on the base laser value###
    if base_laser == 405:   #take channel 0 as the base
        zshifts = []
        for i in range(img.shape[3]):                                   ##iterating all the 4 channels of the image
            zshifts.append((zrange_um*reg_vals['m0'][i])+reg_vals['c0'][i]+zrange_um)

        zshifts = np.array(zshifts)
        minrange = np.amin(zshifts)
        maxrange = np.amax(zshifts)
        newrange = maxrange-minrange

        range2 = np.arange(zstart,zstart+newrange,zscale)
        zshifts2 = []
        for i in range(img.shape[3]):
            zshifts2.append((range2*reg_vals['m0'][i])+reg_vals['c0'][i])
        zshifts2 = np.array(zshifts2)
        nplanes = np.size(range2)
        zshifts_vx = np.round(zshifts2/zscale)
        startpoint = np.amax(zshifts_vx)
        endpoint = np.amin(zshifts_vx)

        czshifts = zshifts_vx-startpoint
        czshifts = czshifts.astype(int)

        img2 = np.zeros((nplanes,img.shape[1],img.shape[2],img.shape[3]))

        for i in range(nplanes):
            for j in range(img.shape[3]):
                zplanes = i+czshifts[j,i]
                if zplanes >0 and zplanes <img.shape[0]:
                    print(zplanes)
                    img2[i,:,:,j] = img[zplanes,:,:,j]
        #tiff.imwrite('corr.tif',img2,metadata={'axes':'ZYXC'})
        return img2
    
    if base_laser == 488:
        zshifts = []
        for i in range(img.shape[3]):                                   ##iterating all the 4 channels of the image
            zshifts.append((zrange_um*reg_vals['m1'][i])+reg_vals['c1'][i]+zrange_um)

        zshifts = np.array(zshifts)
        minrange = np.amin(zshifts)
        maxrange = np.amax(zshifts)
        newrange = maxrange-minrange

        range2 = np.arange(zstart,zstart+newrange,zscale)
        zshifts2 = []
        for i in range(img.shape[3]):
            zshifts2.append((range2*reg_vals['m1'][i])+reg_vals['c1'][i])
        zshifts2 = np.array(zshifts2)
        nplanes = np.size(range2)
        zshifts_vx = np.round(zshifts2/zscale)
        startpoint = np.amax(zshifts_vx)
        endpoint = np.amin(zshifts_vx)

        czshifts = zshifts_vx-startpoint
        czshifts = czshifts.astype(int)

        img2 = np.zeros((nplanes,img.shape[1],img.shape[2],img.shape[3]))

        for i in range(nplanes):
            for j in range(img.shape[3]):
                zplanes = i+czshifts[j,i]
                if zplanes >0 and zplanes <img.shape[0]:
                    print(zplanes)
                    img2[i,:,:,j] = img[zplanes,:,:,j]
        tiff.imwrite('corr.tif',img2,metadata={'axes':'ZYXC'})
        return img2
    
    if base_laser == 552:
        zshifts = []
        for i in range(img.shape[3]):                                   ##iterating all the 4 channels of the image
            zshifts.append((zrange_um*reg_vals['m2'][i])+reg_vals['c2'][i]+zrange_um)

        zshifts = np.array(zshifts)
        minrange = np.amin(zshifts)
        maxrange = np.amax(zshifts)
        newrange = maxrange-minrange

        range2 = np.arange(zstart,zstart+newrange,zscale)
        zshifts2 = []
        for i in range(img.shape[3]):
            zshifts2.append((range2*reg_vals['m2'][i])+reg_vals['c2'][i])
        zshifts2 = np.array(zshifts2)
        nplanes = np.size(range2)
        zshifts_vx = np.round(zshifts2/zscale)
        startpoint = np.amax(zshifts_vx)
        endpoint = np.amin(zshifts_vx)

        czshifts = zshifts_vx-startpoint
        czshifts = czshifts.astype(int)

        img2 = np.zeros((nplanes,img.shape[1],img.shape[2],img.shape[3]))

        for i in range(nplanes):
            for j in range(img.shape[3]):
                zplanes = i+czshifts[j,i]
                if zplanes >0 and zplanes <img.shape[0]:
                    print(zplanes)
                    img2[i,:,:,j] = img[zplanes,:,:,j]
        tiff.imwrite('corr.tif',img2,metadata={'axes':'ZYXC'})
        return img2
    
    if base_laser==638:
        zshifts = []
        for i in range(img.shape[3]):                                   ##iterating all the 4 channels of the image
            zshifts.append((zrange_um*reg_vals['m3'][i])+reg_vals['c3'][i]+zrange_um)

        zshifts = np.array(zshifts)
        minrange = np.amin(zshifts)
        maxrange = np.amax(zshifts)
        newrange = maxrange-minrange

        range2 = np.arange(zstart,zstart+newrange,zscale)
        zshifts2 = []
        for i in range(img.shape[3]):
            zshifts2.append((range2*reg_vals['m3'][i])+reg_vals['c3'][i])
        zshifts2 = np.array(zshifts2)
        nplanes = np.size(range2)
        zshifts_vx = np.round(zshifts2/zscale)
        startpoint = np.amax(zshifts_vx)
        endpoint = np.amin(zshifts_vx)

        czshifts = zshifts_vx-startpoint
        czshifts = czshifts.astype(int)

        img2 = np.zeros((nplanes,img.shape[1],img.shape[2],img.shape[3]))

        for i in range(nplanes):
            for j in range(img.shape[3]):
                zplanes = i+czshifts[j,i]
                if zplanes >0 and zplanes <img.shape[0]:
                    print(zplanes)
                    img2[i,:,:,j] = img[zplanes,:,:,j]
        tiff.imwrite('corr.tif',img2,metadata={'axes':'ZYXC'})
    
    
   

