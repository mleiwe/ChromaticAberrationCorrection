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





def chromatic_abberation_calculation(img,output_dir):
    '''' 
    Calculates the chromatic abberation from regression equations and return a dataframe with the m and c values.

    Parameters:
    -----------
    img: n-D numpy array
    channel_data: xls file from ICY spot detector results, usually named 'results.xls'

    Returns:
    --------
    df: pandas dataframe containing regression values for each channel
    '''''
    ## Do edit this based on the number of channels the image will have
    ch0 = read_channel_data('ch0.xls')
    ch1 = read_channel_data('ch1.xls')
    ch2 = read_channel_data('ch2.xls')
    ch3 = read_channel_data('ch3.xls')
    
    ## get centroids for each channel along the x, y and z axis
    centroids_0_z = get_centroids(ch0, 'Z')
    centroids_0_y = get_centroids(ch0, 'Y')
    centroids_0_x = get_centroids(ch0, 'X')

    centroids_1_z = get_centroids(ch1, 'Z')
    centroids_1_y = get_centroids(ch1, 'Y')
    centroids_1_x = get_centroids(ch1, 'X')

    centroids_2_z = get_centroids(ch2, 'Z')
    centroids_2_y = get_centroids(ch2, 'Y')
    centroids_2_x = get_centroids(ch2, 'X')

    centroids_3_z = get_centroids(ch3, 'Z')
    centroids_3_y = get_centroids(ch3, 'Y')
    centroids_3_x = get_centroids(ch3, 'X')

    ## calculate the difference between the centroids for each channel wrt to each channel in microns
    centroid_diff_00 = centroid_diff(centroids_0_z, centroids_0_z)
    centroid_diff_01 = centroid_diff(centroids_0_z, centroids_1_z)
    centroid_diff_02 = centroid_diff(centroids_0_z, centroids_2_z)
    centroid_diff_03 = centroid_diff(centroids_0_z, centroids_3_z)

    centroid_diff_10 = centroid_diff(centroids_1_z, centroids_0_z)
    centroid_diff_11 = centroid_diff(centroids_1_z, centroids_1_z)
    centroid_diff_12 = centroid_diff(centroids_1_z, centroids_2_z)
    centroid_diff_13 = centroid_diff(centroids_1_z, centroids_3_z)

    centroid_diff_20 = centroid_diff(centroids_2_z, centroids_0_z)
    centroid_diff_21 = centroid_diff(centroids_2_z, centroids_1_z)
    centroid_diff_22 = centroid_diff(centroids_2_z, centroids_2_z)
    centroid_diff_23 = centroid_diff(centroids_2_z, centroids_3_z)

    centroid_diff_30 = centroid_diff(centroids_3_z, centroids_0_z)
    centroid_diff_31 = centroid_diff(centroids_3_z, centroids_1_z)
    centroid_diff_32 = centroid_diff(centroids_3_z, centroids_2_z)
    centroid_diff_33 = centroid_diff(centroids_3_z, centroids_3_z)

    ##calculate regression for each channel wrt each channel
    reg00,m00,c00 = calc_regressions(centroids_0_z, centroid_diff_00)
    reg01,m01,c01 = calc_regressions(centroids_0_z, centroid_diff_01)
    reg02,m02,c02 = calc_regressions(centroids_0_z, centroid_diff_02)
    reg03,m03,c03 = calc_regressions(centroids_0_z, centroid_diff_03)

    reg10,m10,c10 = calc_regressions(centroids_1_z, centroid_diff_10)
    reg11,m11,c11 = calc_regressions(centroids_1_z, centroid_diff_11)
    reg12,m12,c12 = calc_regressions(centroids_1_z, centroid_diff_12)
    reg13,m13,c13 = calc_regressions(centroids_1_z, centroid_diff_13)

    reg20,m20,c20 = calc_regressions(centroids_2_z, centroid_diff_20)
    reg21,m21,c21 = calc_regressions(centroids_2_z, centroid_diff_21)
    reg22,m22,c22 = calc_regressions(centroids_2_z, centroid_diff_22)
    reg23,m23,c23 = calc_regressions(centroids_2_z, centroid_diff_23)

    reg30,m30,c30 = calc_regressions(centroids_3_z, centroid_diff_30)
    reg31,m31,c31 = calc_regressions(centroids_3_z, centroid_diff_31)
    reg32,m32,c32 = calc_regressions(centroids_3_z, centroid_diff_32)
    reg33,m33,c33 = calc_regressions(centroids_3_z, centroid_diff_33)


    ##Create a list of m and c values
    m0  = m00.tolist() + m01.tolist() + m02.tolist() + m03.tolist()
    c0 = [c00,c01,c02,c03]
    m1 = m10.tolist() + m11.tolist() + m12.tolist() + m13.tolist()
    c1 = [c10,c11,c12,c13]
    m2 = m20.tolist() + m21.tolist() + m22.tolist() + m23.tolist()
    c2 = [c20,c21,c22,c23]
    m3 = m30.tolist() + m31.tolist() + m32.tolist() + m33.tolist()
    c3 = [c30,c31,c32,c33]


    reg_vals = pd.DataFrame(list(zip(m0, c0, m1, c1, m2, c2, m3, c3)),columns=['m0','c0','m1','c1','m2','c2','m3','c3'])

    reg_vals.to_csv(os.path.join(output_dir, 'regression_vals.csv'))

    return reg_vals


def read_channel_data(channel_data_file_path):
    #check if the filename has .xls extension if not throw an error
    if not channel_data_file_path.endswith('.xls'):
        raise ValueError('File must be an excel file')
    #read the excel file
    channel_data = pd.read_excel(channel_data_file_path)
    return channel_data

def get_centroids(channel_data, axis):  #get the center x, y and z separately for each channel from the data table
    centroids = []
    for i in range(len(channel_data)):
        centroids.append(channel_data['Center'+' '+str(axis)][i])
    centroids = np.array(centroids)
    return centroids
def centroid_diff(centroids_1, centroids_2,voxel_size=0.24):        # do change the voxel size for the image
    diff = centroids_1 - centroids_2
    diff = diff*voxel_size
    return diff

def calc_regressions(centroids,centroid_diff):
    centroids_reshape = centroids.reshape(-1,1)
    reg = LinearRegression().fit(centroids_reshape,centroid_diff)
    mvals = reg.coef_
    cvals = reg.intercept_
    return reg,mvals,cvals

