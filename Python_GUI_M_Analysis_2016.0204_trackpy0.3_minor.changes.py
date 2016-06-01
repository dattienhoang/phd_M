# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 17:43:42 2014

@author: 
"""


import scipy.optimize as opt
import numpy as np
from scipy.optimize import curve_fit

import pylab as plt
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib
import matplotlib.animation
import matplotlib.gridspec as gridspec
#import matplotlib.patches as mpatches

#import pandas as pd
import trackpy as tp
import pims
import os
import time
import math

import wx
from sys import exit

#from pandas import DataFrame, Series  # for convenience


wildcard = "Tiff file (*.tif; *.tiff)|*.tif;*.tiff|" \
         "All files (*.*)|*.*"
 
wildcardDataFile = "txt file (*.txt)|*.txt|" \
         "All files (*.*)|*.*"


StartingPanel = 1           # 1: Single channel M calculation, 2: Two channel, 3: Data analysis




#######################################################
ThrethholdIntensity = 5000         #Typically 50000 in feature finding
PixelSize = 7 # typically 7         #in feature finding


ColorScaleMin = 1100
ColorScaleMax = 2200

FilePrefixInput = ''

#######################################################




#######################################################
# M analysis default values  ##########################
FitPeriodGuess = 90   # in frames

ThresholdAveIntensity = 0

MinRsquareFitGoodness = 0.7 # in %



# tracing rotation
trace_rotation_x1 = 168 
trace_rotation_y1 = 209
trace_rotation_F2N = 91
trace_rotation_x2 = 153
trace_rotation_y2 = 198


#######################################################





#######################################################
### 2 channel analysis ################################

# A pair x,y in ImageJ
ReferencePairx1, ReferencePairy1 = 110 , 67
ReferencePairx2, ReferencePairy2 = 354 , 68


# boundary limits in finding features.
left_Ch_x_left, left_Ch_x_right = 25, 243
right_Ch_x_left, right_Ch_x_right = 269, 488

#upper_y, bottom_y = 



#######################################################





manual_openFile = 'Open Image File: only for grayscale tiff file' + \
'\n\nChoose Frame #s : enter two frame numbers for finding features. The two frames are averaged. ' +\
'\n\nAverage # frames: It increases the number of frames to be averaged. \nEg) if frame #0 and # 10 are chosen with 3 Average # frames, frames 0,1,3, 10,11,12 are averaged' + \
'\n\nParameters of Crocker-Weeks algorithm\nFeature Size: approximate size of features in pixels \nMin Intensity: threshold integrated brightness ' +\
'\n\nMore(x,y) in finding features: Add (x,y) values from the Finding Features window to include additional features'



manual_showM = 'Show Raw Images: feature images are shown without any calculation.' + \
'\n\nIntensity Method: feature and background intensities are calculated based on the method chosen.' + \
'\nFeature intensity: Within center 7x7 pixels, max 5 ave means that average of the brightest 5 pixels is used as a feature intensity' +\
'\nBackground intensity: edge of center 9x9 pixels: so 32 edge pixels are used to calculate background intensity' + \
'\n\nFit initial guess for period in frames: \nEg) If the polarizer is rotated at 10 deg/s and movie is recorded every 0.2 s, \nit takes 18 s for a half rotation of the polarizer. This is one period for the fit.' + \
'\nThe 18 s corresponds to 90 frames. So enter 90' + \
'\n\nThreshold average intensity: Some measured data are excluded in the anisotropy calculation,\n  if the average intensity of each trajectory is below this value' + \
'\n\nInclude Fit Graphs: Fit graphs are shown together with the intensity trajectories' + \
'\n\nTrace Rotation: If images are rotating when a polarizer is rotating, check this option' + \
'\nTwo reference points need to be entered to calculate angular velocity, initial angle, rotation radius\nChoose a feature in the first frame in ImageJ, enter the x, y in ImageJ' + \
'\nand then, trace the feature until the frame is at 180deg of the polarizer, Enter the frame # and x, y in ImageJ'

manual_saveData = 'Data Files Prefix: Letters to be added at the beginning of the data files to be saved'





def dataplot(title_text1, title_text2, M_before_bleaching, M_after_bleaching, Intensity_before, Intensity_after):
    print '\n\nstarting data plot'
   

    
    '''
    scatterSize = -0.0002*len(M_before_bleaching) + 10
    if scatterSize <1:
        scatterSize = 1
    elif scatterSize > 40:
        scatterSize = 40
        
    scatterAlpha = -0.00001 * len(M_before_bleaching) + 0.4
    if scatterAlpha < 0.02:
        scatterAlpha = 0.02
    elif scatterAlpha > 0.8:
        scatterAlpha = 0.8
    '''
    
    if len(M_before_bleaching) < 1000:
        scatterSize = 20
        scatterAlpha = 0.5
    elif 1000<=len(M_before_bleaching) < 10000:
        scatterSize = 7
        scatterAlpha = 0.2
    elif 10000<=len(M_before_bleaching) < 50000:
        scatterSize = 5
        scatterAlpha = 0.15
    elif 50000<=len(M_before_bleaching) < 70000:
        scatterSize = 3
        scatterAlpha = 0.1
    elif 70000<=len(M_before_bleaching) < 100000:
        scatterSize = 2
        scatterAlpha = 0.05
    elif 100000<=len(M_before_bleaching) < 500000:
        scatterSize = 1
        scatterAlpha = 0.04
    elif 500000<=len(M_before_bleaching):
        scatterSize = 1
        scatterAlpha = 0.02
    
    
    
    
    print 'len(M_before_bleaching) ', len(M_before_bleaching)
    print 'scatterSize ', scatterSize
    print 'scatterAlpha ', scatterAlpha
    
    
    
    
    
    
    Intensity_diff = np.array(Intensity_after) - np.array(Intensity_before)
        
    Intensity_diff_per = 100.*Intensity_diff/np.array(Intensity_before)
    

    
    
    
    
    
    
    
    
    
    
    bin_width = 0.1
    
    fig_3 = plt.figure('DataPlot', figsize=(20, 15), facecolor='white')
    fig_3.subplots_adjust(top = 0.88, bottom = 0.1, left = 0.09, right = 0.97, wspace = 0.4, hspace = 0.7)
    #plt.plot(x, y, '-ro')
    
    plt.figtext(0.02, 0.96, 'Experiment', color='red', size = 16)
    plt.figtext(0.02, 0.92, 'Total #: ' + str(len(M_before_bleaching)), color='b', size = 12)
    
    plt.figtext(0.1, 0.95, title_text1, size = 11)
    plt.figtext(0.1, 0.91, title_text2, size = 10)
    
    
    
    plt.subplot(4,5,1)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels along the bottom edge are off
    
    plt.hist(M_before_bleaching, bins = np.arange(0,1.3,bin_width) , width=bin_width, alpha=1, label = 'dd', color='c')
    plt.title('before, Median: ' + str(np.around(np.median(M_before_bleaching),2)), size=12)
    plt.xlabel('M')
    plt.xlim(0,1.2)
    
    
    
    
    
    
    plt.subplot(4,5,2)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels along the bottom edge are off
    
    plt.hist(M_after_bleaching, bins = np.arange(0,1.3,bin_width) , width=bin_width, alpha=1, label = 'dd', color='c')
    plt.title('after, Median: ' + str(np.around(np.median(M_after_bleaching),2)), size=12)
    plt.xlabel('M')
    plt.xlim(0,1.2)
    
    
    
    
    
    
    
    
    plt.subplot(4,5,3)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels along the bottom edge are off
    
    
    plt.hist(M_before_bleaching, bins = np.arange(0,1.3,bin_width) , width=bin_width, alpha=0.6, label = ' ', color='c')
    plt.hist(M_after_bleaching, bins = np.arange(0,1.3,bin_width) , width=bin_width, alpha=0.6, label = ' ', color='m')
    
    plt.xlabel('M')
    
    plt.xlim(0,1.2)
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    x_temp = np.array([-0.2, 1.3])
    initial_guess = [0, -50]
    
    popt, pcov = opt.curve_fit(Linear_func, M_before_bleaching, Intensity_diff_per, p0 = initial_guess)
    #print '\nFIONA popt: ', popt
    #print '\nFIONA pcov: \n', pcov
    Afit = popt[0]
    Bfit = popt[1]
    AfitErr = np.sqrt(pcov[0,0])                    
    #BfitErr = np.sqrt(pcov[1,1])
    
    
    
    
    plt.subplot(4,5,6)
    plt.title('% intensity change\nAve: ' + str(np.around(np.mean(Intensity_diff_per),1)) + ',  Median: ' + str(np.around(np.median(Intensity_diff_per),1)), size=12)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='on',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='on') # labels along the bottom edge are off
    
    plt.scatter(M_before_bleaching, Intensity_diff_per, s=scatterSize, marker='o', edgecolor='black', linewidth='0', facecolor='red', alpha=scatterAlpha)
    plt.plot(x_temp, Afit*x_temp + Bfit, 'c', label='Fit: slope='+str(np.around(Afit,2)) + ' Err: ' + str(np.around(AfitErr,2)), linewidth = 2)
    plt.rcParams.update({'mathtext.default': 'regular' })
    plt.xlabel('$M_{before}$', size=14)
    plt.ylabel('%', size=14)
    plt.xticks([0.0, 0.5, 1.0])
    plt.yticks([-100, -50, 0, 50])
    plt.xlim(0,1.2)
    plt.ylim(-100,50)
    plt.legend(fontsize=11, loc='upper right',frameon=False)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    median_before = np.around(np.median(M_before_bleaching),3) 
    median_after = np.around(np.median(M_after_bleaching),3)
    
    median_text = 'Median:  before = ' + str(median_before)  + ' ,  after = ' + str(median_after)
    
    
    
    plt.subplot(4,5,7)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels along the bottom edge are off
    
    plt.title(median_text, size = 11)
    plt.hist(M_before_bleaching, bins = np.arange(0,1.3,bin_width) , width=bin_width, alpha=0.6, label = 'Before', color='c')
    plt.hist(M_after_bleaching, bins = np.arange(0,1.3,bin_width) , width=bin_width, alpha=0.6, label = 'After ', color='m')
    plt.xlabel('M')
    plt.xlim(0,1.2)
    
    plt.legend(fontsize=10, loc='upper right',frameon=False)
    
    
    
    
    
    
    
    
    
    
    plt.subplot(4,5,8, aspect='equal')
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='on',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='on') # labels along the bottom edge are off
    
    plt.title('', size=11)
    plt.plot([0,2],[0,2], 'r', linewidth=2)
    plt.scatter(M_before_bleaching, M_after_bleaching, s=scatterSize, marker='o', edgecolor='black', linewidth='0', facecolor='blue', alpha=scatterAlpha)
    plt.rcParams.update({'mathtext.default': 'regular' })
    plt.xlabel('$M_{before}$')
    plt.ylabel('$M_{after}$')
    
    plt.xticks([0.0, 0.5, 1.0])
    plt.yticks([0.0, 0.5, 1.0])
    #plt.ylim(0, 1.1)
    
    plt.xlim(0,1.2)
    plt.ylim(0,1.2)
    
    plt.legend(fontsize=12, loc='upper right',frameon=False)
    
    
    
    
    




    
    
    
    
    plt.subplot(4,5,9, aspect='equal')
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='on',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='on') # labels along the bottom edge are off
    
    plt.title('', size=11)
    plt.plot([0,2],[0,2], 'r', linewidth=2)
    colors = np.array(Intensity_diff_per)
    plt.scatter(M_before_bleaching, M_after_bleaching, c=colors, s=30, marker='o', linewidth='0', alpha=0.5 )
    plt.colorbar()
    plt.rcParams.update({'mathtext.default': 'regular' })
    plt.xlabel('$M_{before}$')
    plt.ylabel('$M_{after}$')
    
    plt.xticks([0.0, 0.5, 1.0])
    plt.yticks([0.0, 0.5, 1.0])
    #plt.ylim(0, 1.1)
    plt.xlim(0,1.2)
    plt.ylim(0,1.2)
    plt.legend(fontsize=12, loc='upper right',frameon=False)
    
    
    
    

    
    
    
    
    
    
    
    M_1 = M_before_bleaching
    M_diff = np.array(M_after_bleaching) - np.array(M_before_bleaching)
        
        
    M_1_bins = np.arange(0.0, 1.2, 0.1)        
    M_1_bins_elements = [[] for _ in xrange(len(M_1_bins))] 
    M_1_bins_numbers_of_elements = [0 for _ in xrange(len(M_1_bins))]     
    M_1_bins_numbers_of_elements_increased = [0 for _ in xrange(len(M_1_bins))]
    M_1_bins_numbers_of_elements_decreased = [0 for _ in xrange(len(M_1_bins))]
    M_1_bins_Ave_M_diff = [0 for _ in xrange(len(M_1_bins))]         
    
    
    for n in range(len(M_1)):
    
        M_temp = M_1[n]
        Mdiff_temp = M_diff[n]
        
        for m in range(len(M_1_bins)):
            if m*0.1 <= M_temp <= m*0.1 + 0.1:
                M_1_bins_elements[m].append(Mdiff_temp)
                
                if Mdiff_temp > 0:
                    M_1_bins_numbers_of_elements_increased[m] += 1
                    
                elif Mdiff_temp <0:
                    M_1_bins_numbers_of_elements_decreased[m] += 1
                    
                
                
    for n in range(len(M_1_bins)):
        M_1_bins_numbers_of_elements[n] = len(M_1_bins_elements[n])
        M_1_bins_Ave_M_diff[n] = np.mean(M_1_bins_elements[n])
        
        
        
    
            
    
    
    x_temp = np.array([-0.2, 1.3])
    initial_guess = [-0.5, 0.5]
    
    popt, pcov = opt.curve_fit(Linear_func, M_1, M_diff, p0 = initial_guess)
    #print '\nFIONA popt: ', popt
    #print '\nFIONA pcov: \n', pcov
    Afit = popt[0]
    Bfit = popt[1]
    AfitErr = np.sqrt(pcov[0,0])                    
    #BfitErr = np.sqrt(pcov[1,1])
    
    
    
    
    plt.subplot(4,5,11)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='on',      # ticks along the bottom edge are off
    right='on',         # ticks along the top edge are off
    labelleft='on') # labels along the bottom edge are off
    
    plt.title('M change vs M ', size=11)
    plt.scatter(M_1, M_diff, s = scatterSize, marker='o', edgecolor='black', linewidth='0', facecolor='red', alpha=scatterAlpha)
    plt.plot([0,2],[0,0], color='black')
    plt.plot(x_temp, Afit*x_temp + Bfit, 'c', label='Fit: slope='+str(np.around(Afit,2)) + ' Err: ' + str(np.around(AfitErr,2)), linewidth = 2)
    plt.xlim(0, 1.2)        
    plt.ylabel('M change', size=14)
    plt.xlabel('$M_{before}$', size=14)
    plt.legend(fontsize=9, loc='upper right',frameon=False)
    plt.locator_params(nbins=5)
    
    
    
    plt.subplot(4,5,12)
    
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='on',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='on') # labels along the bottom edge are off
    
    plt.title('Ave M change vs M ', size=12)
    plt.bar(M_1_bins, M_1_bins_Ave_M_diff, 0.1)
    
    plt.xlabel('$M_{before}$', size=14)
    plt.xlim(0,1.2)
    plt.legend(fontsize=12, loc='upper right',frameon=False)
    plt.locator_params(nbins=5)
    
    
    
    
    
    plt.subplot(4,5,13)
    
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='on',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='on') # labels along the bottom edge are off
    
    plt.title('# vs M ', size=11)
    plt.bar(M_1_bins, M_1_bins_numbers_of_elements, 0.1)
    plt.xlim(0, 1.2)
    plt.ylabel('#', size=14)
    plt.xlabel('$M_{before}$', size=14)
    plt.legend(fontsize=9, loc='upper right')
    plt.locator_params(nbins=5)
    
            
           
    
    
    plt.subplot(4,5,14)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='on',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='on') # labels along the bottom edge are off
    
    plt.title('# vs M ', size=11)
    plt.bar(M_1_bins, M_1_bins_numbers_of_elements_increased, 0.1, label='increased', color='m', alpha=0.4)
    plt.bar(M_1_bins, M_1_bins_numbers_of_elements_decreased, 0.1, label='decreased', color='c', alpha=0.4)
    plt.xlim(0, 1.2)
    plt.ylabel('#', size=14)
    plt.xlabel('$M_{before}$', size=14)
    plt.legend(fontsize=9, loc='upper right',frameon=False)
    plt.locator_params(nbins=5)




    M_1_bins_net_numbers_of_elements = np.array(M_1_bins_numbers_of_elements_increased) - np.array(M_1_bins_numbers_of_elements_decreased)
    
    plt.subplot(4,5,15)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='on',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='on') # labels along the bottom edge are off
    
    plt.title('Net # vs M ', size=11)
    plt.bar(M_1_bins, M_1_bins_net_numbers_of_elements, 0.1, label='', color='b', alpha=1.0)
    plt.xlim(0, 1.2)
    plt.ylabel('#', size=14)
    plt.xlabel('$M_{before}$', size=14)
    plt.legend(fontsize=9, loc='upper right',frameon=False)
    plt.locator_params(nbins=5)




    
    
    font = {'family' : 'normal','weight' : 'normal', 'size'   : 14}
    matplotlib.rc('font', **font)
    majortick_size = 4
    majortick_width = 1.5
    plt.rc('axes', linewidth = 1.5)
    #matplotlib.rcParams['lines.linewidth'] = 4
    matplotlib.rcParams['xtick.major.size'] = majortick_size
    matplotlib.rcParams['xtick.major.width'] = majortick_width
    matplotlib.rcParams['xtick.minor.size'] = 5
    matplotlib.rcParams['xtick.minor.width'] = 4
    matplotlib.rcParams['ytick.major.size'] = majortick_size
    matplotlib.rcParams['ytick.major.width'] = majortick_width
    matplotlib.rcParams['ytick.minor.size'] = 5
    matplotlib.rcParams['ytick.minor.width'] = 4
    matplotlib.rcParams.update({'font.size': 14, 'family' : 'normal', 'weight' : 'normal',  'size': 8})
    
           
    


    
    
      
    
    plt.subplot(4,5,17, aspect=1.8)
    
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='on',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
    
    plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='on',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='on') # labels along the bottom edge are off
    
    #plt.title('Ave M change vs M ', size=12)
    plt.plot([0,2],[0,0], 'k--', linewidth=1.4)
    plt.bar(M_1_bins, M_1_bins_Ave_M_diff, 0.1, edgecolor='c', linewidth=2, alpha=0.4, facecolor='c')
    
    plt.xlabel('$M_{before}$', size=14)
    plt.xlim(0,1.2)
    plt.legend(fontsize=12, loc='upper right',frameon=False)
    plt.locator_params(nbins=5)
      
    plt.ylabel(r'$\Delta M \ $', size=20)
    plt.yticks([-0.1, 0, 0.1, 0.2, 0.3])
    plt.ylim(-0.12, 0.32)
    
    
    
    
    #todayDate = time.strftime("%Y%m%d_%Hh%Mm%Ss")
    #fig_3.savefig('exp_data_' + todayDate + '.png', dpi=600)
    
    
    #exit()
    

    #plt.ylabel(r'$\Delta I_{ave} \, (a.u.)$', size=20)
















def MessageBox(self, msg, title):
    """ Display Info Dialog Message """
    font = wx.Font(14, wx.MODERN, wx.NORMAL, wx.NORMAL)
    style = wx.OK | wx.ICON_INFORMATION | wx.STAY_ON_TOP
    dialog = wx.MessageDialog(self, msg, title, style)
    dialog.CenterOnParent()
    dialog.SetFont(font)
    result = dialog.ShowModal()
    if result == wx.ID_OK:
        print dialog.GetFont().GetFaceName()
    dialog.Destroy()
    return



def xytohist(xdata, ydata, binstep, Nhist):
    ' making histogram from x, y data set '
    histscaling = float(Nhist)
    Nxdata = len(xdata)
    
    areaydata = np.trapz(ydata, xdata)
    
    
    
    xdatamin = np.min(xdata)
    #xdatamax = np.max(xdata)
    histdata = []
    histdataEach = []
    
    for n in range(Nxdata):  # to find starting position and end position of the list tausincILT for step integration
        
        xi =  xdatamin + binstep * n
        xiIndex = findindex(xi, xdata)
        xf = xdatamin + binstep*(n+1)
        xfIndex = findindex(xf, xdata)
        
        areatemp = np.trapz(ydata[xiIndex:xfIndex], xdata[xiIndex:xfIndex])
        
        arearatio = areatemp / areaydata
        
        xcentervalue = (xi + xf)/2.0    
        print 'xcentervalue = ', xcentervalue
        histtempvalue = [xcentervalue] * int(np.around(arearatio * histscaling))
        histdata.extend(histtempvalue)
        histdataEach.extend([xcentervalue])
        print histdataEach
        
        if xf > xdata[Nxdata-1]:
            break
        
    #print histdata
    return histdata, histdataEach

def histtoxy2(histdata, binstart, binend, binstep, normalizationYN):
    'Converting histogram data to x-y data'
    #print histdata
    ydata = np.array(histdata)
    xresult = []
    yresult = []
    for n in np.arange(binstart, binend, binstep):
        #print '\nn', n, n+ binstep
        ytemp = np.where( (n<=ydata) & (ydata< (n+binstep)))
        Ntemp = len(ytemp[0])
        #print ydata
        #print ytemp[0]
        #print 'Ntemp', Ntemp
        xresult.append(n+0.5*binstep)
        yresult.append(float(Ntemp))
        
    
    if normalizationYN == 'Y':
        area = np.sum(yresult) * float(binstep)
        yresult = yresult/area
        yresult = yresult.tolist()
    
    return xresult, yresult
    #xyresult = [np.array(xresult), np.array(yresult)]
    #return xyresult
        


def R_SqureFitGoodness(dataArray, fitArray):  # non-weighted
    data_y = np.array(dataArray)    
    fit_y = np.array(fitArray)
    
    dataAve = np.mean(data_y)
    
    TSS = np.sum( (data_y - dataAve)**2 )
    RSS = np.sum( (data_y - fit_y)**2 )
    
    Rsquare = 1.0 - float(RSS)/float(TSS)
    
    return Rsquare
    
    






         
def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()
    
    


def Sine_fit_func(x, A, T, phi, y0):
    #print 'fitting sine func'
    return A*np.cos((2*np.pi/T)*x - phi) + y0

       


def ModulationDepth_M_func(x, A, M, T, phi):
    #print 'fitting sine func'
    return A*(1 + M*np.cos( (2*np.pi/T)*x - phi ) )

       



    
def Linear_func(x, A, B):
    return A*np.array(x)+B
    


    
def Noise_func(x, A, B):
    return A*np.sqrt(x)+B
    
    


def CreateDataFiles(filedirectoryOnly, filename, xpos_ImageJ, ypos_ImageJ, A, T, phi, y0, M, I_noBGS, I_BGS, I_BGS_fit, I_BG, I_BG_fit, CenterImages, cb_SaveImageJ_xyData_yn, cb_Complete_Data_yn, BGS_type):
    print '\n\nCreateDataFiles'
        

    todayDate = time.strftime("%Y%m%d_%Hh%Mm")
    print "\nToday's Date: ", time.strftime("%Y-%m-%d_%Hh%Mm")

    
    if cb_SaveImageJ_xyData_yn:

        ff = open(filedirectoryOnly + '\\' + filename +'_' + todayDate +  '_6_xy_A_T(frames)_phi(deg)_y0_M.txt','w')
        
        for i in range(len(xpos_ImageJ)):
            
            ff.write(str(xpos_ImageJ[i]) + ' ' + str(ypos_ImageJ[i]) 
            + ' ' + str(A[i]) + ' ' + str(T[i]) + ' ' + str(phi[i]) + ' ' + str(y0[i]) + ' ' + str(M[i]) + '\n'   )    
            
        ff.close()       

            

    if cb_Complete_Data_yn:
        
        ff = open(filedirectoryOnly + '\\' + filename +'_' + todayDate +  '_10_CompleteDataAll' + '.txt','w')
        
        ff.write('%% x y A T phi y0 M Background_Subtraction_Type: '+ BGS_type + ' - rawIntensity , intensityBGS , intensityBGS_fit , BG , BGfit , image0 , image1' + '\n')
        
        for n in range(len(xpos_ImageJ)):
       
            ff.write('#' + str(n) + '\n')
                
            ff.write(str(xpos_ImageJ[n]) + ' ' + str(ypos_ImageJ[n]) + ' ' + str(A[n])
            + ' ' + str(T[n]) + ' ' + str(phi[n]) 
            + ' ' + str(y0[n]) + ' ' + str(M[n]) + '\n'   )    
       
            for i in range(len(I_noBGS[n])):
                ff.write(str(int(I_noBGS[n][i])) + ' ')
            ff.write('\n')

            for i in range(len(I_BGS[n])):
                ff.write(str(int(I_BGS[n][i])) + ' ')
            ff.write('\n')
            
            for i in range(len(I_BGS_fit[n])):
                ff.write(str(int(I_BGS_fit[n][i])) + ' ')
            ff.write('\n')
       
            for i in range(len(I_BG[n])):
                ff.write(str(int(I_BG[n][i])) + ' ')
            ff.write('\n')

            for i in range(len(I_BG_fit[n])):
                ff.write(str(int(I_BG_fit[n][i])) + ' ')
            ff.write('\n')
     
            for i in range(len(CenterImages[n][0])):
                for k in range(len(CenterImages[n][0][i])):
                    ff.write(str(int(CenterImages[n][0][i][k])) + ' ')
                ff.write(' , ')
                
            ff.write('\n')
            
            for i in range(len(CenterImages[n][1])):
                for k in range(len(CenterImages[n][1][i])):
                    ff.write(str(int(CenterImages[n][1][i][k])) + ' ')
                ff.write(' , ')
                
            ff.write('\n\n')
                
        ff.close()       
            
        
  



    


    
def FIONA_FramesByFrame(filepath, xpos, ypos, frameStart, frameEnd, SvsB_Tolerance, maxEcc, maxPer, maxPx):   
    # using all features for all frames by analyzing frame by frame

    print '\n\n## FIONA _ Frame By Frame ##'
    print 'frameStart ', frameStart
    print 'frameEnd ', frameEnd
    print 'max Ecc = ', maxEcc
    print 'max Pencentage err = ', maxPer
    print 'max Pixel err = ', maxPx 
    print 'data file = ', filepath
    
    
    
    
    frames = tp.TiffStack(filepath)
    
    TotalNframes = len(frames)
    print 'TotalNframes: ', TotalNframes
    
    
    
    
    if frameEnd >= TotalNframes:
        frameEnd = TotalNframes - 1
        print '\nframeEnd changed to ', TotalNframes - 1, '\n'
    
    
    
    
    N_xPixel = frames.frame_shape[0]
    N_yPixel = frames.frame_shape[1]
    print 'N_xPixel: ', N_xPixel
    print 'N_yPixel: ', N_yPixel
    
    
    
    
    ################################################
    x_ImageJ, y_ImageJ = xpos, ypos
    
    xc = y_ImageJ
    yc = N_xPixel - x_ImageJ + 1
    ################################################

    





    NCenterSubPixel = 7  # number of pixels from the center
    
    
    
    #frameTemp0 = frame0[yc-NCenterSubPixel:yc+NCenterSubPixel + 1, xc-NCenterSubPixel:xc+NCenterSubPixel+1]
    x3 = np.linspace(0, 2 * NCenterSubPixel, 2* NCenterSubPixel + 1)
    y3 = np.linspace(0, 2 * NCenterSubPixel, 2*NCenterSubPixel+1)
    x3, y3 = np.meshgrid(x3, y3)  # for 2D gaussian fit
    
    
    BG_IntensityAve = [[] for _ in xrange(len(xc))] 
    centerIntensityAve = [[] for _ in xrange(len(xc))] 
    centerIntensityMax = [[] for _ in xrange(len(xc))] 
    centerImage = [[] for _ in xrange(len(xc))] # centerImage[Molecule#][Frame#]
    
    for frameTemp in np.arange(0, frameEnd+1): # from the 0th frame to the end frame.
        print 'calculating center intensity data for frame # ', frameTemp
        for k in range(len(xpos)):            
            
            frameLargeTemp = frames[frameTemp][yc[k]-4:yc[k]+5, xc[k]-4:xc[k]+5]
            frameCenterTemp = frames[frameTemp][yc[k]-2:yc[k]+3, xc[k]-2:xc[k]+3]
            
            if (frameTemp >= frameStart) and (frameTemp <= frameStart + 4):
                centerImage[k].append(frameLargeTemp) # it's not matching to the actual frame #, only first a few are saved
            
            centerIntensityAve[k].append(np.sum(frameCenterTemp)/25)
            centerIntensityMax[k].append(np.amax(frameCenterTemp))
            BG_IntensityAve[k].append( (np.sum(frameLargeTemp) - np.sum(frameCenterTemp))/56 )
            
    SBratioAve = np.float16(np.array(centerIntensityAve)) / np.array(BG_IntensityAve)    
    SBratioMax = np.float16(np.array(centerIntensityMax)) / np.array(BG_IntensityAve)    
    print 'SBratioAve = ', SBratioAve   
    print 'SBratioMax = ', SBratioMax
            
    
            
    
    #print 'centerIntensityAve = ', centerIntensityAve
    #print 'centerIntensityAve len = ', len(centerIntensityAve)
    
    BG_AveTotal = np.average(BG_IntensityAve)
    BG_AveEachFrame = np.around(np.average(BG_IntensityAve, axis = 0), 0) 
    print 'BG_AveTotal = ', BG_AveTotal
    print 'BG_AveEachFrame = ', BG_AveEachFrame
    #print 'frameCenterTemp ', frameCenterTemp.ravel()


    '''
    plt.ioff()
    plt.figure('centerAve')
    plt.plot(centerIntensityAve[0])
    plt.plot(centerIntensityAve[1])
    plt.figure('centerMax')
    plt.plot(centerIntensityMax[0])
    plt.plot(centerIntensityMax[1])
    plt.figure('BG_Ave')
    plt.plot(BG_IntensityAve[0])
    plt.plot(BG_IntensityAve[1])
    
    plt.figure('centerAve vs BG_Ave')

    #plt.plot(centerIntensityAve[0])
    #plt.plot(BG_IntensityAve[0])
    #plt.plot(centerIntensityAve[1])
    #plt.plot(BG_IntensityAve[1])

    plt.plot(centerIntensityAve[2])
    plt.plot(BG_IntensityAve[2])
    
    plt.show()
    '''
    
    xposfit_FIONA_0 = [[] for _ in xrange(len(xc))] # with S/B
    yposfit_FIONA_0 = [[] for _ in xrange(len(xc))] # with S/B
    
    xposfit_FIONA_1 = [[] for _ in xrange(len(xc))] # with S/B and max e
    yposfit_FIONA_1 = [[] for _ in xrange(len(xc))] # with S/B and max e
    
    xposfit_FIONA_2 = [[] for _ in xrange(len(xc))] # with S/B, max e, max AfitErr/Afit %, max px
    yposfit_FIONA_2 = [[] for _ in xrange(len(xc))] # with S/B, max e, max AfitErr/Afit %, max px
        
    
    
    eccentricityFor_FIONA = [[] for _ in xrange(len(xc))] 





    
    data_FIONA = [[] for _ in xrange(len(xc))] 
    

    for n in range(len(xc)):
        data_FIONA[n].append('M# ' + str(n) + ' (x,y): ' + str(x_ImageJ[n]) + ' ' + str(y_ImageJ[n])      )
        data_FIONA[n].append('frame#,  x0fitFionaImageJ, y0fitFionaImageJ, Afit, z0fit, sigma_x, sigma_y, theta, eccentricity, BG_IntensityAve, centerIntensityAve, centerIntensityMax, SBratioAve, SBratioMax ::: AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaEff') 
    
    
    
    for f1n in range(0, frameEnd + 1):
        
        #################################################################################
        ###  FIONA only part                              ############
        print '\n FIONA start for frame # ', f1n
        for MN_fiona in range(len(xpos)): #for each feature
        
            if SBratioAve[MN_fiona][f1n] <= SvsB_Tolerance :
                print 'Low SB, FIONA passSed for M# ', MN_fiona, '  frame # ', f1n
                
                pass

            else:           
                print '\nFIONA for M# ', MN_fiona, '  frame # ', f1n
                frameFIONA = frames[f1n][yc[MN_fiona]-NCenterSubPixel:yc[MN_fiona]+NCenterSubPixel+1, xc[MN_fiona]-NCenterSubPixel:xc[MN_fiona]+NCenterSubPixel+1]


                frameFIONA_raveled = frameFIONA.ravel()
            
            
                # FIONA fit the data  for all frames
                initial_guess = (6000, NCenterSubPixel+1,  NCenterSubPixel+1, 2, 2, 10, 1000)
 
                k = 1
               
                Afit = 0
                x0fit = 0
                y0fit = 0
                z0fit = 0
                
                try:
                    popt, pcov = opt.curve_fit(twoD_Gaussian, (x3, y3), frameFIONA_raveled, p0 = initial_guess)
                    #print '\nFIONA popt: ', popt
                    #print '\nFIONA pcov: \n', pcov
                    Afit = popt[0]
                    x0fit = popt[1]
                    y0fit = popt[2]
                    
                    x0fitImageJ = xpos[MN_fiona] - (NCenterSubPixel - y0fit)
                    y0fitImageJ = ypos[MN_fiona] - (NCenterSubPixel - x0fit)
                    
                    z0fit = popt[6]
                    
                    sigma_x = popt[3]
                    sigma_y = popt[4]
                    theta = popt[5]

                    
                    #FIONA AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaEff
                    AfitErr = np.sqrt(pcov[0,0])                    
                    x0fitImageJErr = np.sqrt(pcov[1,1])
                    y0fitImageJErr = np.sqrt(pcov[2,2]) 
                    z0fitErr = np.sqrt(pcov[6,6]) 
                    sigma_xErr = np.sqrt(pcov[3,3]) 
                    sigma_yErr = np.sqrt(pcov[4,4]) 
                    thetaErr = np.sqrt(pcov[5,5])
                                         

                    
                    # FIONA Percentage of the errors: AfitErr x0fitImageJErr y0fitImageJErr z0fitErr sigma_xErr sigma_yErr thetaEff
                    AfitErrP = np.around(100.0 * AfitErr / Afit, 0)                    
                    x0fitImageJErrPx = np.around(np.sqrt(pcov[1,1]),2)
                    y0fitImageJErrPx = np.around(np.sqrt(pcov[2,2]),2) 
                    z0fitErrP = np.around(100.0 * z0fitErr / z0fit, 0)      
                    sigma_xErrP = np.around(100.0 * sigma_xErr/sigma_x, 0)   
                    sigma_yErrP = np.around(100.0 * sigma_yErr/sigma_y, 0)  
                    thetaErrDeg = np.around(np.sqrt(pcov[5,5]),2)
                    
                    
                    
                    
                    ########### eccentricicity calculation  ############
                    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
                    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
                    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2) 
                    
                    #print 'a, b, c =', a,b,c
                    
                    
                    A = a
                    B = 2*b
                    C = c
                            
                    eccentricityTemp = np.sqrt(2*(np.sqrt((A-C)**2 + B**2)) / (A+C + np.sqrt(((A-C)**2 + B**2)))  )
                    eccentricityFor_FIONA[MN_fiona].append(eccentricityTemp)
                    #print 'eccentricityTemp = ', eccentricityTemp
                    
                    #'frame#,  x0fitFionaImageJ, y0fitFionaImageJ, Afit, z0fit, sigma_x, sigma_y, theta, eccentricity
                    #, BG_IntensityAve, centerIntensityAve, centerIntensityMax, SBratioAve, SBratioMax)'
                    datatempFIONA = (str(f1n) + ' ' + str(x0fitImageJ) + ' ' + str(y0fitImageJ) + ' ' 
                            + str(Afit) +  ' ' + str(z0fit) + ' ' + str(sigma_x) + ' ' + str(sigma_y) + ' ' +  str(theta) + ' ' + str(eccentricityTemp) + ' '
                            + str(BG_IntensityAve[MN_fiona][f1n]) + ' ' + str(centerIntensityAve[MN_fiona][f1n]) + ' ' 
                            + str(centerIntensityMax[MN_fiona][f1n]) + ' ' + str(SBratioAve[MN_fiona][f1n]) + ' ' + str(SBratioMax[MN_fiona][f1n]) + ' Err: '
                            + str(AfitErr) +  ' ' + str(x0fitImageJErr) + ' ' + str(y0fitImageJErr) + ' ' + str(z0fitErr) + ' ' 
                            + str(sigma_xErr) + ' ' + str(sigma_yErr) + ' ' +  str(thetaErr)
                            )
                    
                    
                    
                    data_FIONA[MN_fiona].append(datatempFIONA)

                    
                    #print 'M#, datatemp FIONA: ' , MN_fiona , '   ' ,  datatempFIONA
                    
                    xposfit_FIONA_0[MN_fiona].append(np.around(x0fitImageJ,2))
                    yposfit_FIONA_0[MN_fiona].append(np.around(y0fitImageJ,2))
                    
                    if (eccentricityTemp <= maxEcc):
                        xposfit_FIONA_1[MN_fiona].append(np.around(x0fitImageJ,2))
                        yposfit_FIONA_1[MN_fiona].append(np.around(y0fitImageJ,2))
                        print '\nFIONA_1 OK'
                        

                        if (AfitErrP <= maxPer) & (x0fitImageJErrPx <= maxPx) & (y0fitImageJErrPx <= maxPx):
                            xposfit_FIONA_2[MN_fiona].append(np.around(x0fitImageJ,2))
                            yposfit_FIONA_2[MN_fiona].append(np.around(y0fitImageJ,2))
                            print '\nFIONA_2 OK'
                            
                            
                        else:
                            print 'FINOA_2 failed'
                            print 'AfitErrP ', AfitErrP
                            print 'x0fitImageJErrPx ', x0fitImageJErrPx
                            print 'y0fitImageJErrPx ', y0fitImageJErrPx
                            
                    else:
                        print 'FIONA_1 failed'
                        print 'eccentricityTemp ', eccentricityTemp
                        

                

                 
                    ####################################################
                
                except:
                    print 'FIONA fitting failed'
                    eccentricityTemp = 1.2
                    datatempFIONA = (str(f1n) + ' ' + str(999) + ' ' + str(999) + ' ' 
                            + str(999) +  ' ' + str(999) + ' ' + str(999) + ' ' + str(999) + ' ' +  str(999) + ' ' + str(eccentricityTemp) + ' '
                            + str(BG_IntensityAve[MN_fiona][f1n]) + ' ' + str(centerIntensityAve[MN_fiona][f1n]) + ' ' 
                            + str(centerIntensityMax[MN_fiona][f1n]) + ' ' + str(SBratioAve[MN_fiona][f1n]) + ' ' + str(SBratioMax[MN_fiona][f1n]) + ' Err: '       
                            + str(999) +  ' ' + str(999) + ' ' + str(999) + ' ' + str(999) + ' ' 
                            + str(999) + ' ' + str(999) + ' ' +  str(999)
                            )
                     
                    data_FIONA[MN_fiona].append(datatempFIONA)
                    eccentricityFor_FIONA[MN_fiona].append(eccentricityTemp)

                    
                    #print 'M#, datatemp FIONA: ' , MN_fiona , '   ' ,  datatempFIONA
                    
        
        print 'FIONA end\n'
        ######################  FIONA only part                              ############
        #################################################################################
      
        



    return (centerImage, xposfit_FIONA_0, yposfit_FIONA_0, data_FIONA, xposfit_FIONA_1, yposfit_FIONA_1, xposfit_FIONA_2, yposfit_FIONA_2
            , eccentricityFor_FIONA)






      
      
def IntensityTrajectoryPlotsFuncRotForM(loadedframes, filename, xpos, ypos, radius, theta0, ref_dx, ref_dy, dtheta, IntensityMethod): # for rotating image movies
    print '\n\n Starting IntensityTrajectoryPlotsFuncRot for M' 
    
    print 'xpos ', xpos
    print 'ypos ', ypos
    
    
    print 'filename = ', filename
    frames = loadedframes
    
    Nframes = len(frames)
    print 'Nframes = ', Nframes
    
    
    ################################################
    x_ImageJ, y_ImageJ = np.array(xpos) + ref_dx , np.array(ypos) + ref_dy # making (x,y) as rotation centers
    
    #theta0 = np.arccos( (xpos[0] - x_ImageJ[0])/4.12 )
    
        
    #dtheta = (2*3.14159)/180
    
    FrameN = np.arange(0, Nframes)
    
    dx = (radius * np.cos(theta0 - dtheta * FrameN))
    dy = (radius * np.sin(theta0 - dtheta * FrameN))
    
    ImageJ_dx = dx.astype(int)
    ImageJ_dy = dy.astype(int)
    
    print 'ImageJ_dx ', ImageJ_dx
    print 'ImageJ_dy ', ImageJ_dy
    
    
    print 'ImageJ_x ' , x_ImageJ[0] + ImageJ_dx
    print 'ImageJ_y ' , y_ImageJ[0] + ImageJ_dy
    
    ################################################

    print 'x_ImageJ', x_ImageJ
    print 'y_ImageJ', y_ImageJ
    

    
    N_xPixel = frames.frame_shape[0]
    N_yPixel = frames.frame_shape[1]
    
    
    print N_xPixel
    print 'type(N_xPixel)', type(N_xPixel)    
    print N_yPixel
    print 'type(N_yPixel)', type(N_yPixel)
    
    #xc = y_ImageJ
    #yc = N_xPixel - x_ImageJ -1
    
    xc = x_ImageJ - 0 #trackpy 0.3
    yc = y_ImageJ  #trackpy 0.3
     
    
    #xc_dx = ImageJ_dy
    #yc_dy = - ImageJ_dx
    
    xc_dx = - ImageJ_dx
    yc_dy = ImageJ_dy
    
    
    
    print 'xc_dx ', xc_dx
    print 'yc_dy ', yc_dy
    
    
    
    frameStart = 0
    frameEnd = len(frames) - 1
    
    print '# of frames ', len(frames)    
    



    
    # calculating each frame for all features => 
    centerIntensityAveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))] 
    centerIntensityTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    backgroundIntensityAveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    backgroundIntensityTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    
    #print 'centerIntensityMaxTrajectory ', centerIntensityMaxTrajectory
    #print 'centerIntensityMaxTrajectory[0] len : ', len(centerIntensityMaxTrajectory[0])
    #int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
    
    
    try:
        edgeArray = frames[0][yc[0]-4:yc[0]+5, xc[0]-4:xc[0]+5] * 0
        print yc[0], xc[0]
        print 'edgeArray zeros try \n', edgeArray
        for n in range(9):
            for k in range(9):
                if n == 0 or n == 8 or k == 0 or k==8:
                    edgeArray[n][k] = 1
        
    except:
        edgeArray = frames[0][yc[1]-4:yc[1]+5, xc[1]-4:xc[1]+5] * 0
        print 'edgeArray zeros except \n', edgeArray
        for n in range(9):
            for k in range(9):
                if n == 0 or n == 8 or k == 0 or k==8:
                    edgeArray[n][k] = 1
                    
    print 'edgeArray \n ' , edgeArray
    
    
    
    
    edgeArrayTop = edgeArray * 0
    edgeArrayBottom = edgeArray * 0
    edgeArrayLeft = edgeArray * 0
    edgeArrayRight = edgeArray * 0
    
    for n in range(9):
        for k in range(9):
            if n == 0:
                edgeArrayTop[n][k] = 1
            if n == 8:
                edgeArrayBottom[n][k] = 1
            if k == 0:
                edgeArrayLeft[n][k] = 1
            if k == 8:
                edgeArrayRight[n][k] = 1
    
    print 'edgeArrayTop \n ' , edgeArrayTop
    print 'edgeArrayBottom \n ' , edgeArrayBottom
    print 'edgeArrayLeft \n ' , edgeArrayLeft
    print 'edgeArrayRight \n ' , edgeArrayRight
    
    
    for k in range(frameStart,frameEnd+1):
        print 'frame # ', k
        frametemp = frames[k]
        for n in range(len(xc)): # calculating background max5ave using a single side which has min std
            frameLargeTemp = frametemp[yc[n]+yc_dy[k] - 4 : yc[n]+yc_dy[k] + 5, xc[n]+xc_dx[k] - 4 : xc[n]+xc_dx[k] + 5]
            frameCenterTemp = frametemp[yc[n]+yc_dy[k] - 3 : yc[n]+yc_dy[k] + 4, xc[n]+xc_dx[k] - 3 : xc[n]+xc_dx[k] + 4]
            
            
            try:
                frameBackgroundTemp = frameLargeTemp * edgeArray
                #print 'frameBackgroundTempTop, frameBackgroundTempBottom, frameBackgroundTempLeft, frameBackgroundTempRight', frameBackgroundTempTop, frameBackgroundTempBottom, frameBackgroundTempLeft, frameBackgroundTempRight
                #print 'std top, bottom, left, right ', bgTop, bgBottom, bgLeft, bgRight
                #print 'frameBackgroundTemp ', frameBackgroundTemp
                
            except:
                frameBackgroundTemp = frameLargeTemp * 0
                
            
            frameBackgroundTemp_32pixels = np.sort(frameBackgroundTemp.ravel())[-32:]
            frameBackgroundTemp_49pixels = frameBackgroundTemp_32pixels.tolist() + np.random.choice(frameBackgroundTemp_32pixels, size=17).tolist()
            
            
            
            
            
            
            
            
            
            
            '''
            if IntensityMethod == 0:
                centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/49)
                centerIntensityTrajectory[n][k] = int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]))
                backgroundIntensityAveTrajectory[n][k] = int( (np.sum(frameLargeTemp) - np.sum(frameCenterTemp))/32 )
                backgroundIntensityTrajectory[n][k] = int(np.mean( np.sort(frameBackgroundTemp.ravel())[-5:])) 
                
            
            if IntensityMethod == 1:
                centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/49)
                centerIntensityTrajectory[n][k] = int(np.median( np.sort(frameCenterTemp.ravel())[-5:]) )
                backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameLargeTemp) - np.sum(frameCenterTemp))/32)
                backgroundIntensityTrajectory[n][k] = int(np.median( np.sort(frameBackgroundTemp.ravel())[-5:])) 
                
            
            if IntensityMethod == 2:
                centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/49)
                centerIntensityTrajectory[n][k] = int(np.amax( np.sort(frameCenterTemp.ravel())[-5:]))
                backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameLargeTemp) - np.sum(frameCenterTemp))/32)
                backgroundIntensityTrajectory[n][k] = int(np.amax(np.sort(frameBackgroundTemp.ravel())[-5:])) 
               
            '''
            


            if IntensityMethod == 0:
                centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/49)
                centerIntensityTrajectory[n][k] = int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
                backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameBackgroundTemp_32pixels))/32)
                backgroundIntensityTrajectory[n][k] = int(np.mean( np.sort(frameBackgroundTemp_49pixels)[-5:])) # max5mean 
                
            
            elif IntensityMethod == 1:
                centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/49)
                centerIntensityTrajectory[n][k] = int(np.median( np.sort(frameCenterTemp.ravel())[-5:]) )
                backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameBackgroundTemp_32pixels))/32)
                backgroundIntensityTrajectory[n][k] = int(np.median( np.sort(frameBackgroundTemp_49pixels)[-5:]) ) # max5median
                
            
            elif IntensityMethod == 2:
                centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/49)
                centerIntensityTrajectory[n][k] = int(np.amax( np.sort(frameCenterTemp.ravel())[-5:]) )
                backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameBackgroundTemp_32pixels))/32)
                backgroundIntensityTrajectory[n][k] = int(np.amax( np.sort(frameBackgroundTemp_49pixels)[-5:]) ) # max1 
                
                
            elif IntensityMethod == 3:
                centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/49)
                centerIntensityTrajectory[n][k] = int(np.sum(frameCenterTemp))
                backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameBackgroundTemp_32pixels))/32)
                backgroundIntensityTrajectory[n][k] = int(np.median(frameBackgroundTemp_32pixels)*49) # sum == median * 49
                
                
            
            
            
            
    
    
        
    print 'End of Intensity Trajectory func \n'
    
    
    ''' the centerMax intensity is average of center max 5 pixels   '''        
    return centerIntensityAveTrajectory, centerIntensityTrajectory, backgroundIntensityAveTrajectory, backgroundIntensityTrajectory
    

  
  
  
  
  
  
  
        

  
def IntensityTrajectoryPlotsFuncForM(loadedframes, filename, xpos, ypos, IntensityMethod, centerLength):
    print '\n Starting IntensityTrajectoryPlotsFunc for M' 
    ################################################
    x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
    
    ################################################

    print 'x_ImageJ', x_ImageJ
    print 'y_ImageJ', y_ImageJ
    print 'intensity method', IntensityMethod
    
    print '\ncenterLength', centerLength
    
    
    print 'filename = ', filename
    frames = loadedframes
    
    Nframes = len(frames)
    print 'Nframes = ', Nframes
    
    N_xPixel = frames.frame_shape[0]
    N_yPixel = frames.frame_shape[1]
    
    
    print N_xPixel
    print 'type(N_xPixel)', type(N_xPixel)    
    print N_yPixel
    print 'type(N_yPixel)', type(N_yPixel)
    
    
    #xc = y_ImageJ - 0
    #yc = N_xPixel - x_ImageJ - 1
    
    xc = x_ImageJ - 0 #trackpy 0.3
    yc = y_ImageJ  #trackpy 0.3
    
        
    
    frameStart = 0
    frameEnd = len(frames) - 1
    
    print '# of frames ', len(frames)    
    


    
    '''
    # calculating each feature for all frames => slow
    centerIntensityAveTrajectory = []
    centerIntensityMaxTrajectory = []
    for kTemp in range(len(xc)):
        print 'Intensity Trajectory feature # ', kTemp
        xctemp = xc[kTemp]
        yctemp = yc[kTemp]      
        
        centerIntensityAveTemp = []
        centerIntensityMaxTemp = []
        
        for nTemp in np.arange(frameStart, frameEnd+1):
            frameCenterTemp = frames[nTemp][yctemp-2:yctemp+3, xctemp-2:xctemp+3]
            centerIntensityAveTemp.append(np.sum(frameCenterTemp))
            centerIntensityMaxTemp.append(np.amax(frameCenterTemp))
            
        centerIntensityAveTrajectory.append(centerIntensityAveTemp)   
        centerIntensityMaxTrajectory.append(centerIntensityMaxTemp)    
    '''
    
    
    
    # calculating each frame for all features => 
    centerIntensityAveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))] 
    centerIntensityTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    backgroundIntensityAveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    backgroundIntensityTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    
    #print 'centerIntensityMaxTrajectory ', centerIntensityMaxTrajectory
    #print 'centerIntensityMaxTrajectory[0] len : ', len(centerIntensityMaxTrajectory[0])
    #int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
    
    
    
    CenterEdgeLengthHalf = int((centerLength-1)/2)
    BGEdgeLengthHalf = int((centerLength+1)/2)
    

    edgeArray = frames[0][yc[0]-BGEdgeLengthHalf:yc[0]+BGEdgeLengthHalf+1, xc[0]-BGEdgeLengthHalf:xc[0]+BGEdgeLengthHalf+1] * 0
    print 'edgeArray zeros try \n', edgeArray
    for n in range(centerLength+2):
        for k in range(centerLength+2):
            if n == 0 or n == BGEdgeLengthHalf*2 or k == 0 or k==BGEdgeLengthHalf*2:
                edgeArray[n][k] = 1
    
                    
    print 'edgeArray \n ' , edgeArray
    
    




    """

    if centerLength == 7:
        
        for k in range(frameStart,frameEnd+1):
            print 'frame # ', k
            frametemp = frames[k]
            for n in range(len(xc)):
                frameLargeTemp = frametemp[yc[n]-4:yc[n]+5, xc[n]-4:xc[n]+5]
                frameCenterTemp = frametemp[yc[n]-3:yc[n]+4, xc[n]-3:xc[n]+4]
                
                #   '''
                plt.figure()
                plt.subplot(121)
                plt.imshow(frameLargeTemp, interpolation = 'None')
                plt.subplot(122)
                plt.imshow(frameCenterTemp, interpolation = 'None')
                plt.ioff()
                plt.show()
                #'''
                

                frameBackgroundTemp = frameLargeTemp * edgeArray
                    
                
                #print 'frameBackgroundTemp \n', frameBackgroundTemp
                
                #print 'B_Max5Ave: ', int(np.mean( np.sort(frameBackgroundTemp.ravel())[-5:]) )
                
                frameBackgroundTemp_32pixels = np.sort(frameBackgroundTemp.ravel())[-32:]
    
                #print 'frameBackgroundTemp ', frameBackgroundTemp
                #print 'frameBackgroundTemp_32pixels ', frameBackgroundTemp_32pixels
                
                frameBackgroundTemp_49pixels = frameBackgroundTemp_32pixels.tolist() + np.random.choice(frameBackgroundTemp_32pixels, size=17).tolist()
                
                
                
                #print 'frameBackgroundTemp_49pixels ', frameBackgroundTemp_49pixels
                #print ' ', 
                
                
                
                if IntensityMethod == 0:
                    centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/49)
                    centerIntensityTrajectory[n][k] = int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
                    backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameBackgroundTemp_32pixels))/32)
                    backgroundIntensityTrajectory[n][k] = int(np.mean( np.sort(frameBackgroundTemp_49pixels)[-5:])) # max5mean 
                    
                
                elif IntensityMethod == 1:
                    centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/49)
                    centerIntensityTrajectory[n][k] = int(np.median( np.sort(frameCenterTemp.ravel())[-5:]) )
                    backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameBackgroundTemp_32pixels))/32)
                    backgroundIntensityTrajectory[n][k] = int(np.median( np.sort(frameBackgroundTemp_49pixels)[-5:]) ) # max5median
                    
                
                elif IntensityMethod == 2:
                    centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/49)
                    centerIntensityTrajectory[n][k] = int(np.amax( np.sort(frameCenterTemp.ravel())[-5:]) )
                    backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameBackgroundTemp_32pixels))/32)
                    backgroundIntensityTrajectory[n][k] = int(np.amax( np.sort(frameBackgroundTemp_49pixels)[-5:]) ) # max1 
                    
                    
                elif IntensityMethod == 3:
                    centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/49)
                    centerIntensityTrajectory[n][k] = int(np.sum(frameCenterTemp))
                    backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameBackgroundTemp_32pixels))/32)
                    backgroundIntensityTrajectory[n][k] = int(np.median(frameBackgroundTemp_32pixels)*49) # sum == median * 49
                    
                    
     #  """                
                    
                    


    N_centerPixel = centerLength**2
    N_BG_edgePixel = centerLength*2 + (centerLength+2)*2
    N_diff_center_edge = N_centerPixel - N_BG_edgePixel
    
    print 'N_centerPixel : ', N_centerPixel
    print 'N_BG_edgePixel : ', N_BG_edgePixel
    print 'N_diff_center_edge : ', N_diff_center_edge
 




   
    for k in range(frameStart,frameEnd+1):
        print 'frame # ', k
        frametemp = frames[k]
        for n in range(len(xc)):
            frameLargeTemp = frametemp[yc[n]-BGEdgeLengthHalf:yc[n]+BGEdgeLengthHalf+1, xc[n]-BGEdgeLengthHalf:xc[n]+BGEdgeLengthHalf+1]
            frameCenterTemp = frametemp[yc[n]-CenterEdgeLengthHalf:yc[n]+CenterEdgeLengthHalf+1, xc[n]-CenterEdgeLengthHalf:xc[n]+CenterEdgeLengthHalf+1]
            
            '''
            plt.figure()
            plt.subplot(121)
            plt.imshow(frameLargeTemp, interpolation = 'None')
            plt.subplot(122)
            plt.imshow(frameCenterTemp, interpolation = 'None')
            plt.ioff()
            plt.show()
            #'''
            

            frameBackgroundTemp = frameLargeTemp * edgeArray
                
            
            #print 'frameBackgroundTemp \n', frameBackgroundTemp
            
            #print 'B_Max5Ave: ', int(np.mean( np.sort(frameBackgroundTemp.ravel())[-5:]) )
            
            frameBackgroundTemp_EdgePixels = np.sort(frameBackgroundTemp.ravel())[-N_BG_edgePixel:]

            #print 'frameBackgroundTemp ', frameBackgroundTemp
            #print 'frameBackgroundTemp_32pixels ', frameBackgroundTemp_32pixels
            
            frameBackgroundTemp_EdgePixels_adjusted = frameBackgroundTemp_EdgePixels.tolist() + np.random.choice(frameBackgroundTemp_EdgePixels, size=N_diff_center_edge).tolist()
            
            
            
            #print 'frameBackgroundTemp_49pixels ', frameBackgroundTemp_49pixels
            #print ' ', 
            
            
            
            if IntensityMethod == 0:
                centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/N_centerPixel)
                centerIntensityTrajectory[n][k] = int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
                backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameBackgroundTemp_EdgePixels))/N_BG_edgePixel)
                backgroundIntensityTrajectory[n][k] = int(np.mean( np.sort(frameBackgroundTemp_EdgePixels_adjusted)[-5:])) # max5mean 
                
            
            elif IntensityMethod == 1:
                centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/N_centerPixel)
                centerIntensityTrajectory[n][k] = int(np.median( np.sort(frameCenterTemp.ravel())[-5:]) )
                backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameBackgroundTemp_EdgePixels))/N_BG_edgePixel)
                backgroundIntensityTrajectory[n][k] = int(np.median( np.sort(frameBackgroundTemp_EdgePixels_adjusted)[-5:]) ) # max5median
                
            
            elif IntensityMethod == 2:
                centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/N_centerPixel)
                centerIntensityTrajectory[n][k] = int(np.amax( np.sort(frameCenterTemp.ravel())[-5:]) )
                backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameBackgroundTemp_EdgePixels))/N_BG_edgePixel)
                backgroundIntensityTrajectory[n][k] = int(np.amax( np.sort(frameBackgroundTemp_EdgePixels_adjusted)[-5:]) ) # max1 
                
                
            elif IntensityMethod == 3:
                centerIntensityAveTrajectory[n][k] = int(np.sum(frameCenterTemp)/N_centerPixel)
                centerIntensityTrajectory[n][k] = int(np.sum(frameCenterTemp))
                backgroundIntensityAveTrajectory[n][k] = int((np.sum(frameBackgroundTemp_EdgePixels))/N_BG_edgePixel)
                backgroundIntensityTrajectory[n][k] = int(np.median(frameBackgroundTemp_EdgePixels)*N_centerPixel) # sum == median * 49
                
                
                    
                    
                    
                
                
                
             
                
                
                
                
                
                
                
                
    print 'End of Intensity Trajectory func \n'
    
    
    ''' the centerMax intensity is average of center max 5 pixels   '''        
    return centerIntensityAveTrajectory, centerIntensityTrajectory, backgroundIntensityAveTrajectory, backgroundIntensityTrajectory
    

      
      
      
      
      
      
      
      
      
      
      
      
      
def IntensityTrajectoryPlotsFuncRot(filename, xpos, ypos, radius, theta0, ref_dx, ref_dy, dtheta): # for rotating image movies
    print '\n\n Starting IntensityTrajectoryPlotsFuncRot' 
    
    print 'xpos ', xpos
    print 'ypos ', ypos
    
    
    print 'filename = ', filename
    frames = tp.TiffStack(filename)
    
    Nframes = len(frames)
    print 'Nframes = ', Nframes
    
    
    ################################################
    x_ImageJ, y_ImageJ = np.array(xpos) + ref_dx , np.array(ypos) + ref_dy
    
    #theta0 = np.arccos( (xpos[0] - x_ImageJ[0])/4.12 )
    
        
    #dtheta = (2*3.14159)/180
    
    FrameN = np.arange(0, Nframes)
    
    dx = (radius * np.cos(theta0 + dtheta * FrameN))
    dy = (-radius * np.sin(theta0 + dtheta * FrameN))
    
    ImageJ_dx = dx.astype(int)
    ImageJ_dy = dy.astype(int)
    
    print 'ImageJ_dx ', ImageJ_dx
    print 'ImageJ_dy ', ImageJ_dy
    
    
    print 'ImageJ_x ' , x_ImageJ[0] + ImageJ_dx
    print 'ImageJ_y ' , y_ImageJ[0] + ImageJ_dy
    
    ################################################

    print 'x_ImageJ', x_ImageJ
    print 'y_ImageJ', y_ImageJ
    

    
    N_xPixel = frames.frame_shape[0]
    N_yPixel = frames.frame_shape[1]
    
    
    print N_xPixel
    print 'type(N_xPixel)', type(N_xPixel)    
    print N_yPixel
    print 'type(N_yPixel)', type(N_yPixel)
    
    xc = y_ImageJ
    yc = N_xPixel - x_ImageJ
    
    xc_dx = ImageJ_dy
    yc_dy = - ImageJ_dx
    
    print 'xc_dx ', xc_dx
    print 'yc_dy ', yc_dy
    
    
    
    frameStart = 0
    frameEnd = len(frames) - 1
    
    print '# of frames ', len(frames)    
    


    
    '''
    # calculating each feature for all frames => slow
    centerIntensityAveTrajectory = []
    centerIntensityMax5AveTrajectory = []
    for kTemp in range(len(xc)):
        print 'Intensity Trajectory feature # ', kTemp
        xctemp = xc[kTemp]
        yctemp = yc[kTemp]      
        
        centerIntensityAveTemp = []
        centerIntensityMaxTemp = []
        
        for nTemp in np.arange(frameStart, frameEnd+1):
            frameCenterTemp = frames[nTemp][yctemp-2:yctemp+3, xctemp-2:xctemp+3]
            centerIntensityAveTemp.append(np.sum(frameCenterTemp))
            centerIntensityMaxTemp.append(np.amax(frameCenterTemp))
            
        centerIntensityAveTrajectory.append(centerIntensityAveTemp)   
        centerIntensityMax5AveTrajectory.append(centerIntensityMaxTemp)    
    '''
    
    # calculating each frame for all features => 
    centerIntensityAveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))] 
    centerIntensityMax5AveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    backgroundIntensityAveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    backgroundIntensityMax5AveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    
    #print 'centerIntensityMax5AveTrajectory ', centerIntensityMax5AveTrajectory
    #print 'centerIntensityMax5AveTrajectory[0] len : ', len(centerIntensityMax5AveTrajectory[0])
    #int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
    
    
    try:
        edgeArray = frames[0][yc[0]-4:yc[0]+5, xc[0]-4:xc[0]+5] * 0
        print yc[0], xc[0]
        print 'edgeArray zeros try \n', edgeArray
        for n in range(9):
            for k in range(9):
                if n == 0 or n == 8 or k == 0 or k==8:
                    edgeArray[n][k] = 1
        
    except:
        edgeArray = frames[0][yc[1]-4:yc[1]+5, xc[1]-4:xc[1]+5] * 0
        print 'edgeArray zeros except \n', edgeArray
        for n in range(9):
            for k in range(9):
                if n == 0 or n == 8 or k == 0 or k==8:
                    edgeArray[n][k] = 1
                    
    print 'edgeArray \n ' , edgeArray
    
    edgeArrayTop = edgeArray * 0
    edgeArrayBottom = edgeArray * 0
    edgeArrayLeft = edgeArray * 0
    edgeArrayRight = edgeArray * 0
    
    for n in range(9):
        for k in range(9):
            if n == 0:
                edgeArrayTop[n][k] = 1
            if n == 8:
                edgeArrayBottom[n][k] = 1
            if k == 0:
                edgeArrayLeft[n][k] = 1
            if k == 8:
                edgeArrayRight[n][k] = 1
    
    print 'edgeArrayTop \n ' , edgeArrayTop
    print 'edgeArrayBottom \n ' , edgeArrayBottom
    print 'edgeArrayLeft \n ' , edgeArrayLeft
    print 'edgeArrayRight \n ' , edgeArrayRight
    
    
    for k in range(frameStart,frameEnd+1):
        print 'frame # ', k
        frametemp = frames[k]
        for n in range(len(xc)): # calculating background max5ave using a single side which has min std
            frameLargeTemp = frametemp[yc[n]+yc_dy[k] - 4 : yc[n]+yc_dy[k] + 5, xc[n]+xc_dx[k] - 4 : xc[n]+xc_dx[k] + 5]
            frameCenterTemp = frametemp[yc[n]+yc_dy[k] - 3 : yc[n]+yc_dy[k] + 4, xc[n]+xc_dx[k] - 3 : xc[n]+xc_dx[k] + 4]
            centerIntensityAveTrajectory[n][k] = (np.sum(frameCenterTemp))
            
            try:
                frameBackgroundTempTop = frameLargeTemp * edgeArrayTop
                frameBackgroundTempBottom = frameLargeTemp * edgeArrayBottom
                frameBackgroundTempLeft = frameLargeTemp * edgeArrayLeft
                frameBackgroundTempRight = frameLargeTemp * edgeArrayRight
                
                bgTop = np.std( np.sort(frameBackgroundTempTop.ravel())[-5:])
                bgBottom = np.std( np.sort(frameBackgroundTempBottom.ravel())[-5:])
                bgLeft = np.std( np.sort(frameBackgroundTempLeft.ravel())[-5:])
                bgRight = np.std( np.sort(frameBackgroundTempRight.ravel())[-5:])
                
                bgIndex = np.argmin([ bgTop, bgBottom, bgLeft, bgRight])
                
                if bgIndex == 0:
                    frameBackgroundTemp = frameBackgroundTempTop
                elif bgIndex == 1:
                    frameBackgroundTemp = frameBackgroundTempBottom
                elif bgIndex == 2:
                    frameBackgroundTemp = frameBackgroundTempLeft
                elif bgIndex == 3:
                    frameBackgroundTemp = frameBackgroundTempRight
                        
                    
                #print 'frameBackgroundTempTop, frameBackgroundTempBottom, frameBackgroundTempLeft, frameBackgroundTempRight', frameBackgroundTempTop, frameBackgroundTempBottom, frameBackgroundTempLeft, frameBackgroundTempRight
                #print 'std top, bottom, left, right ', bgTop, bgBottom, bgLeft, bgRight
                #print 'frameBackgroundTemp ', frameBackgroundTemp




                
            except:
                frameBackgroundTemp = frameLargeTemp * 0
                
            
            #print 'frameBackgroundTemp \n', frameBackgroundTemp
            
            #print 'B_Max5Ave: ', int(np.mean( np.sort(frameBackgroundTemp.ravel())[-5:]) )
            
            
            
            
            #centerIntensityMax5AveTrajectory[n][k] = (np.amax(frameCenterTemp))      
            centerIntensityMax5AveTrajectory[n][k] = int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
            
            backgroundIntensityAveTrajectory[n][k] = ( (np.sum(frameLargeTemp) - np.sum(frameCenterTemp))/32.0 )
            
            #backgroundIntensityMax5AveTrajectory[n][k] = int(np.mean( np.sort(frameBackgroundTemp.ravel())[-1:]) ) # max1ave for a single side
            backgroundIntensityMax5AveTrajectory[n][k] = int(np.mean( np.sort(frameBackgroundTemp.ravel())[-5:]) ) # max5ave for all single sides
        
    
        
    print 'End of Intensity Trajectory func \n'
    
    
    ''' the centerMax intensity is average of center max 5 pixels   '''        
    return centerIntensityAveTrajectory, centerIntensityMax5AveTrajectory, backgroundIntensityAveTrajectory, backgroundIntensityMax5AveTrajectory
    

  
  
  
  
  
        

  
def IntensityTrajectoryPlotsFunc(filename, xpos, ypos):
    print '\n Starting IntensityTrajectoryPlotsFunc' 
    ################################################
    x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
    
    ################################################

    print 'x_ImageJ', x_ImageJ
    print 'y_ImageJ', y_ImageJ
    
    print 'filename = ', filename
    frames = tp.TiffStack(filename)
    
    Nframes = len(frames)
    print 'Nframes = ', Nframes
    
    N_xPixel = frames.frame_shape[0]
    N_yPixel = frames.frame_shape[1]
    
    
    print N_xPixel
    print 'type(N_xPixel)', type(N_xPixel)    
    print N_yPixel
    print 'type(N_yPixel)', type(N_yPixel)
    
    xc = y_ImageJ
    yc = N_xPixel - x_ImageJ + 1
    
    
    
    frameStart = 0
    frameEnd = len(frames) - 1
    
    print '# of frames ', len(frames)    
    


    
    '''
    # calculating each feature for all frames => slow
    centerIntensityAveTrajectory = []
    centerIntensityMax5AveTrajectory = []
    for kTemp in range(len(xc)):
        print 'Intensity Trajectory feature # ', kTemp
        xctemp = xc[kTemp]
        yctemp = yc[kTemp]      
        
        centerIntensityAveTemp = []
        centerIntensityMaxTemp = []
        
        for nTemp in np.arange(frameStart, frameEnd+1):
            frameCenterTemp = frames[nTemp][yctemp-2:yctemp+3, xctemp-2:xctemp+3]
            centerIntensityAveTemp.append(np.sum(frameCenterTemp))
            centerIntensityMaxTemp.append(np.amax(frameCenterTemp))
            
        centerIntensityAveTrajectory.append(centerIntensityAveTemp)   
        centerIntensityMax5AveTrajectory.append(centerIntensityMaxTemp)    
    '''
    
    # calculating each frame for all features => 
    centerIntensityAveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))] 
    centerIntensityMax5AveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    backgroundIntensityAveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    backgroundIntensityMax5AveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    
    #print 'centerIntensityMax5AveTrajectory ', centerIntensityMax5AveTrajectory
    #print 'centerIntensityMax5AveTrajectory[0] len : ', len(centerIntensityMax5AveTrajectory[0])
    #int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
    
    
    try:
        edgeArray = frames[0][yc[0]-4:yc[0]+5, xc[0]-4:xc[0]+5] * 0
        print 'edgeArray zeros try \n', edgeArray
        for n in range(9):
            for k in range(9):
                if n == 0 or n == 8 or k == 0 or k==8:
                    edgeArray[n][k] = 1
        
    except:
        edgeArray = frames[0][yc[1]-4:yc[1]+5, xc[1]-4:xc[1]+5] * 0
        print 'edgeArray zeros except \n', edgeArray
        for n in range(9):
            for k in range(9):
                if n == 0 or n == 8 or k == 0 or k==8:
                    edgeArray[n][k] = 1
                    
    print 'edgeArray \n ' , edgeArray
    
    edgeArrayTop = edgeArray * 0
    edgeArrayBottom = edgeArray * 0
    edgeArrayLeft = edgeArray * 0
    edgeArrayRight = edgeArray * 0
    
    for n in range(9):
        for k in range(9):
            if n == 0:
                edgeArrayTop[n][k] = 1
            if n == 8:
                edgeArrayBottom[n][k] = 1
            if k == 0:
                edgeArrayLeft[n][k] = 1
            if k == 8:
                edgeArrayRight[n][k] = 1
    
    print 'edgeArrayTop \n ' , edgeArrayTop
    print 'edgeArrayBottom \n ' , edgeArrayBottom
    print 'edgeArrayLeft \n ' , edgeArrayLeft
    print 'edgeArrayRight \n ' , edgeArrayRight
    
    
    for k in range(frameStart,frameEnd+1):
        print 'frame # ', k
        frametemp = frames[k]
        for n in range(len(xc)): # calculating background max5ave using a single side which has min std
            frameLargeTemp = frametemp[yc[n]-4:yc[n]+5, xc[n]-4:xc[n]+5]
            frameCenterTemp = frametemp[yc[n]-3:yc[n]+4, xc[n]-3:xc[n]+4]
            centerIntensityAveTrajectory[n][k] = (np.sum(frameCenterTemp))
            
            try:
                frameBackgroundTempTop = frameLargeTemp * edgeArrayTop
                frameBackgroundTempBottom = frameLargeTemp * edgeArrayBottom
                frameBackgroundTempLeft = frameLargeTemp * edgeArrayLeft
                frameBackgroundTempRight = frameLargeTemp * edgeArrayRight
                
                bgTop = np.std( np.sort(frameBackgroundTempTop.ravel())[-5:])
                bgBottom = np.std( np.sort(frameBackgroundTempBottom.ravel())[-5:])
                bgLeft = np.std( np.sort(frameBackgroundTempLeft.ravel())[-5:])
                bgRight = np.std( np.sort(frameBackgroundTempRight.ravel())[-5:])
                
                bgIndex = np.argmin([ bgTop, bgBottom, bgLeft, bgRight])
                
                if bgIndex == 0:
                    frameBackgroundTemp = frameBackgroundTempTop
                elif bgIndex == 1:
                    frameBackgroundTemp = frameBackgroundTempBottom
                elif bgIndex == 2:
                    frameBackgroundTemp = frameBackgroundTempLeft
                elif bgIndex == 3:
                    frameBackgroundTemp = frameBackgroundTempRight
                        
                    
                #print 'frameBackgroundTempTop, frameBackgroundTempBottom, frameBackgroundTempLeft, frameBackgroundTempRight', frameBackgroundTempTop, frameBackgroundTempBottom, frameBackgroundTempLeft, frameBackgroundTempRight
                #print 'std top, bottom, left, right ', bgTop, bgBottom, bgLeft, bgRight
                #print 'frameBackgroundTemp ', frameBackgroundTemp




                
            except:
                frameBackgroundTemp = frameLargeTemp * 0
                
            
            #print 'frameBackgroundTemp \n', frameBackgroundTemp
            
            #print 'B_Max5Ave: ', int(np.mean( np.sort(frameBackgroundTemp.ravel())[-5:]) )
            
            
            
            
            #centerIntensityMax5AveTrajectory[n][k] = (np.amax(frameCenterTemp))      
            centerIntensityMax5AveTrajectory[n][k] = int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) )
            
            backgroundIntensityAveTrajectory[n][k] = ( (np.sum(frameLargeTemp) - np.sum(frameCenterTemp))/32.0 )
            
            backgroundIntensityMax5AveTrajectory[n][k] = int(np.mean( np.sort(frameBackgroundTemp.ravel())[-1:]) ) # max1ave for a single side
        
    
        
    print 'End of Intensity Trajectory func \n'
    
    
    ''' the centerMax intensity is average of center max 5 pixels   '''        
    return centerIntensityAveTrajectory, centerIntensityMax5AveTrajectory, backgroundIntensityAveTrajectory, backgroundIntensityMax5AveTrajectory
    

  
  
  
  
class TwoChannelMainTracking(wx.Frame):
    #ImageFilePath = ''
 
    #----------------------------------------------------------------------
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "Tracking_Trackpy_for_SM1", size=(800, 1000), pos=(0,0))
                          
        self.MainPanel = wx.Panel(self, wx.ID_ANY)
        

        self.btn_SwitchTo1chAnalysis = wx.Button(self.MainPanel, pos=(600,10), label="Switch To 1 Ch Analysis")
        self.btn_SwitchTo1chAnalysis.Bind(wx.EVT_BUTTON, self.SwitchTo1chAnalysis)
        
        
        wx.StaticText(self.MainPanel, -1, 'Enter a pair (x,y) in ImageJ:       x1, y1\n For Two Channel Data', pos=(50, 20))
        wx.StaticText(self.MainPanel, -1, 'x2, y2', pos=(380, 20))
        self.ImageJPair_x1 = wx.TextCtrl(self.MainPanel, -1, str(ReferencePairx1), pos=(250, 17), size=(40,-1))
        self.ImageJPair_y1 = wx.TextCtrl(self.MainPanel, -1, str(ReferencePairy1), pos=(300, 17), size=(40,-1))
        self.ImageJPair_x2 = wx.TextCtrl(self.MainPanel, -1, str(ReferencePairx2), pos=(420, 17), size=(40,-1))
        self.ImageJPair_y2 = wx.TextCtrl(self.MainPanel, -1, str(ReferencePairy2), pos=(470, 17), size=(40,-1))
        
                
        btn2 = wx.Button(self.MainPanel, pos=(50,50), label="Open Image File")
        self.ImageFilePath = btn2.Bind(wx.EVT_BUTTON, self.onOpenImageFile)
        self.btn2text = wx.StaticText(self.MainPanel, -1, str(self.ImageFilePath), pos=(170, 50))
        wx.StaticText(self.MainPanel, -1, "Choose a Frame #:", pos=(170, 73))
        self.Nframe = wx.TextCtrl(self.MainPanel, -1, "0", pos=(280, 70))
        
        
        btn_FindingFeature = wx.Button(self.MainPanel, pos=(50,100), label="Find Features")
        btn_FindingFeature.Bind(wx.EVT_BUTTON, self.FindingFeatureFunc)
        
        
        wx.StaticText(self.MainPanel, -1, "Average # frames:", pos=(170, 103))
        self.NFramesForAveraging = wx.TextCtrl(self.MainPanel, -1, "10", pos=(280, 100))
        
        wx.StaticText(self.MainPanel, -1, "Feature Size:", pos=(400, 103))
        self.FeatureSize = wx.TextCtrl(self.MainPanel, -1, "7", pos=(480, 100))
        
        wx.StaticText(self.MainPanel, -1, "Min Intensity:", pos=(400, 133))
        self.MinIntensity = wx.TextCtrl(self.MainPanel, -1, str(ThrethholdIntensity), pos=(480, 130))
        
        wx.StaticText(self.MainPanel, -1, "Number of Features:", pos=(170, 133))
        self.NofFeatures = wx.StaticText(self.MainPanel, -1, "None", pos=(300, 133), style=5)
        
        
        
        wx.StaticText(self.MainPanel, -1, 'More (x, y) in Finding Features', pos=(50, 173))
        self.MorePairs_x = [[] for _ in xrange(5)] 
        self.MorePairs_y = [[] for _ in xrange(5)]        
        for n in range(len(self.MorePairs_x)):
            self.MorePairs_x[n] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(240 + n*90, 170), size=(30,-1))
            self.MorePairs_y[n] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(273+ n*90, 170), size=(30,-1))

        
        
        
        
        self.btn_FindAllPairs = wx.Button(self.MainPanel, pos=(50,260), label="Find All Pairs")
        self.btn_FindAllPairs.Bind(wx.EVT_BUTTON, self.FindingPairs)
        wx.StaticText(self.MainPanel, -1, "Total # of pairs:", pos=(170, 264))
        self.TotalNofPairs = wx.StaticText(self.MainPanel, -1, "None", pos=(270, 264))
        wx.StaticText(self.MainPanel, -1, ',      Color Scale       Min', pos=(350, 264))
        self.ColorScaleMin = wx.TextCtrl(self.MainPanel, -1, str(ColorScaleMin), pos=(470, 260), size=(50, -1))        
        wx.StaticText(self.MainPanel, -1, 'Max', pos=(550, 264))        
        self.ColorScaleMax = wx.TextCtrl(self.MainPanel, -1, str(ColorScaleMax), pos=(580, 260), size=(50, -1))





        
        self.btn_ShowAllPairsMovie = wx.Button(self.MainPanel, pos=(600,300), label="Show 4-Pair Movies")
        self.btn_ShowAllPairsMovie.Bind(wx.EVT_BUTTON, self.ShowAllPairsMovie)
        wx.StaticText(self.MainPanel, -1, "Four pair #: ", pos=(170, 304))
        
        self.MoviePairN = [None for _ in xrange(4)] 
        self.MoviePairN[0] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(250, 300), size=(30, -1))
        self.MoviePairN[1] = wx.TextCtrl(self.MainPanel, -1, "1", pos=(300, 300), size=(30, -1))
        self.MoviePairN[2] = wx.TextCtrl(self.MainPanel, -1, "2", pos=(350, 300), size=(30, -1))
        self.MoviePairN[3] = wx.TextCtrl(self.MainPanel, -1, "3", pos=(400, 300), size=(30, -1))
        
        
        


        self.btn_ShowAllPairsMovieSingleN = wx.Button(self.MainPanel, pos=(600,350), label="Show A Pair Movie")
        self.btn_ShowAllPairsMovieSingleN.Bind(wx.EVT_BUTTON, self.ShowAllPairsMovieSingleN)
        wx.StaticText(self.MainPanel, -1, "A pair #: ", pos=(170, 354))
        
        self.MoviePairNSingleN = [None for _ in xrange(1)] 
        self.MoviePairNSingleN[0] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(250, 350), size=(30, -1))





        self.btn_ShowSinglePairMovie = wx.Button(self.MainPanel, pos=(600,400), label="Show A Pair Movie in ImageJ")
        self.btn_ShowSinglePairMovie.Bind(wx.EVT_BUTTON, self.ShowSinglePairsMovie)
        wx.StaticText(self.MainPanel, -1, "For a pair in ImageJ :           x1, y1", pos=(100, 404))
        wx.StaticText(self.MainPanel, -1, "x2, y2 ", pos=(400, 404))

        self.SingleMovieImageJPair_x1 = wx.TextCtrl(self.MainPanel, -1, "0", pos=(300, 400), size=(40,-1))
        self.SingleMovieImageJPair_y1 = wx.TextCtrl(self.MainPanel, -1, "0", pos=(343, 400), size=(40,-1))
        
        self.SingleMovieImageJPair_x2 = wx.TextCtrl(self.MainPanel, -1, "0", pos=(440, 400), size=(40,-1))
        self.SingleMovieImageJPair_y2 = wx.TextCtrl(self.MainPanel, -1, "0", pos=(483, 400), size=(40,-1))






        
        
        self.btn_ShowLDTrajectory = wx.Button(self.MainPanel, pos=(50,450), label="Show LD Trajectory")
        self.btn_ShowLDTrajectory.Bind(wx.EVT_BUTTON, self.ShowLDTrajectory)
        wx.StaticText(self.MainPanel, -1, "LD, Intensity", pos=(200, 454))
        
        
        
        
        
        
        
        wx.StaticText(self.MainPanel, -1, "Start frame #", pos=(200, 512))
        self.frameStart = wx.TextCtrl(self.MainPanel, -1, "0", pos=(280, 510), size=(40,-1))
        wx.StaticText(self.MainPanel, -1, "End frame #", pos=(450, 512))
        self.frameEnd = wx.TextCtrl(self.MainPanel, -1, "4", pos=(520, 510), size=(40,-1))
        wx.StaticText(self.MainPanel, -1, "Min S(ave)/B(ave) ", pos=(180, 542))
        self.SvsB_Tolerance = wx.TextCtrl(self.MainPanel, -1, "1.1", pos=(280, 540), size=(40,-1))
        wx.StaticText(self.MainPanel, -1, "Max eccentricity", pos=(420, 542))
        self.MaxEccentricity = wx.TextCtrl(self.MainPanel, -1, "0.6", pos=(520, 540), size=(40,-1))
        wx.StaticText(self.MainPanel, -1, "Max fit err %: A", pos=(180, 572))
        self.MaxPercentage_Tolerance = wx.TextCtrl(self.MainPanel, -1, "10", pos=(280, 570), size=(40,-1))        
        wx.StaticText(self.MainPanel, -1, "Max fit err px: x, y", pos=(410, 572))
        self.MaxPixel_Tolerance = wx.TextCtrl(self.MainPanel, -1, "0.2", pos=(520, 570), size=(40,-1))        
        
                
                
        self.btn_ShowFIONAtracing = wx.Button(self.MainPanel, pos=(50,600), label="FIONA Tracing")
        self.btn_ShowFIONAtracing.Hide()
        self.btn_ShowFIONAtracing.Bind(wx.EVT_BUTTON, self.ShowFIONA_tracing)
        wx.StaticText(self.MainPanel, -1, "For pair #: ", pos=(200, 604))
        
        self.FIONA_pair_N = [[] for _ in xrange(10)]
        for n in range(len(self.FIONA_pair_N)):
            self.FIONA_pair_N[n] = wx.TextCtrl(self.MainPanel, -1, "-", pos=(270 + 40*n, 600), size=(30,-1))
            
        
        
        
        
        
        
        
        
        

        self.btn_ShowAnisotropyM = wx.Button(self.MainPanel, pos=(50,650), label="Show Excitation M")
        self.btn_ShowAnisotropyM.Bind(wx.EVT_BUTTON, self.ShowExcitationM)
        self.btn_ShowAnisotropyM.Hide()
        wx.StaticText(self.MainPanel, -1, "Excitation M & Intensity", pos=(200, 654))
        
        wx.StaticText(self.MainPanel, -1, 'Background Reference Pair #', pos=(400, 654))        
        self.BG_ReferncePair = wx.TextCtrl(self.MainPanel, -1, 'None', pos=(550, 650), size=(50, -1))
        



        
        
                
        wx.StaticText(self.MainPanel, -1, "Figure Files Prefix", pos=(50, 692))
        self.FilePrefix = wx.TextCtrl(self.MainPanel, -1, str(FilePrefixInput), pos=(160, 690))
        
        self.btn_Save_figures = wx.Button(self.MainPanel, pos=(50,720), label="Save Figures")  
        self.btn_Save_figures.Hide()
        self.btn_Save_figures.Bind(wx.EVT_BUTTON, self.SaveFigures)





     #----------------------------------------------------------------------
 
    def SwitchTo1chAnalysis(self, event):
        SingleChannelMain = SingleChannelMainTracking()
        SingleChannelMain.Show()
        self.Close()
        
 
 
    def onOpenFile1(self, event):
        plt.close("all")

        """
        Create and show the Open FileDialog
        """
        dlg = wx.FileDialog(
            self, message="Choose a file",
            defaultFile="",
            wildcard=wildcard,defaultDir=os.getcwd(),
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            print "You chose the following file(s):"
            for path in paths:
                print path
        print dlg
        print 'path = ', dlg.GetPath()
        print 'filename = ', dlg.GetFilename()
        #dlg.Destroy()
        
        self.Show()
        
        
    def onOpenImageFile(self, event):
        self.CloseAllFigures()
        self.btn_FindAllPairs.Hide()
        self.btn_ShowAllPairsMovie.Hide()
        self.btn_ShowAllPairsMovieSingleN.Hide()
        self.btn_ShowLDTrajectory.Hide()
        self.btn_Save_figures.Hide()
                
                
                
        """
        Create and show the Open FileDialog
        """
        plt.close()
        dlg2 = wx.FileDialog(
            self, message="Choose a file",
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.MULTIPLE) # | wx.CHANGE_DIR  )
            
        if dlg2.ShowModal() == wx.ID_OK:
            paths = dlg2.GetPaths()
            print "You chose the following file(s):"
            for path in paths:
                print path
        print ''
        print 'dlg2 = ', dlg2
        print 'path = ', dlg2.GetPath()
        print 'path type', type(dlg2.GetPath())
        print 'filename = ', dlg2.GetFilename()
        #dlg.Destroy()
        
        self.Show()
        
        #wx.StaticText(self.MainPanel, -1, str(dlg2.GetFilename()), pos=(50, 150))
        
        filenameTemp = dlg2.GetFilename()
        self.btn2text.SetLabel(str(filenameTemp))
        
        self.filepath = dlg2.GetPath()
        self.filenameOnly = dlg2.GetFilename()
        self.filedirectoryOnly = dlg2.GetDirectory()
        
        self.ShowAFewFrames()       
        

    def CloseAllFigures(self):
        try:
            self.fig_ShowAfewFrames.close()
            self.fig_FindingFeatures.close()
            for n in range(len(self.fig_PairImages)):
                self.fig_PairImages[n].close()
            
        except:
            pass
        
        
        
        try:
            for n in range(len(self.fig_PairImages)):
                self.fig_PairImages[n].close()
        except:
            pass


        
        try:
            self.FigShowAllPairsMovie.close()
        except:
            pass
        
        
        
        try:
            for n in range(len(self.fig_IntensityLDTrajectory)):
                self.fig_IntensityLDTrajectory[n].close()
        except:
            pass
        
                        
                        
    def ShowAFewFrames(self):
        frames = tp.TiffStack(self.filepath)
        self.TotalNframes = frames._count
        print 'TotalNframes = ', self.TotalNframes
        self.fig_ShowAfewFrames = plt.figure('frames')
        plt.subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.55, hspace = 0.2)
        plt.figtext(0.5, 0.95, str(os.path.basename(self.filepath))
            ,ha='center', color='black', weight='bold', size='small')
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(810,30,800, 800)
        
        if self.TotalNframes > 15:
            k = 16
        else:
            k = self.TotalNframes
            
        for n in range(k):
            print 'n = ', n
            plt.subplot(4,4,n+1)
            plt.title('frame # ' + str(n))
            
            
            try:
                print 'tp.__version__\n', tp.__version__
                tp_ver = tp.__version__[:3]
                if float(tp_ver[:3]) == 0.3:
                    plt.imshow(frames[n], origin='upper')  # for trackpy 3.0
                else:
                    plt.imshow(np.rot90(frames[n],3), origin='upper')
            except:
                plt.imshow(np.rot90(frames[n],3), origin='upper')
                    
                
                   
            
        
            
            #plt.imshow(np.rot90(frames[n],3), origin='upper')
            
            
        plt.ion()
        plt.show()
        

        
        
        
        
        
        
        

    def FindingFeatureFunc(self, event):
        
        
        plt.close('Finding Features')
        featureSize = int(self.FeatureSize.GetValue())
        minIntensity = int(self.MinIntensity.GetValue())
        frames = tp.TiffStack(self.filepath)
        NframeTemp = int(self.Nframe.GetValue())        
        NFramesForAveraging = int(self.NFramesForAveraging.GetValue())
        
        print 'Starting Frame # ', NframeTemp
        
        framesSummed = np.rot90(frames[0].T, 2) * 0
        for n in range(NFramesForAveraging):
            framesSummed += np.rot90(frames[n].T, 2)
        framesSummed /= NFramesForAveraging
        
        f = tp.locate(framesSummed, featureSize, minmass=minIntensity, invert=False)
        NFeatures = f.values.size/8
        print 'NFeatures = ', NFeatures
        self.NofFeatures.SetLabel(str(NFeatures))
        
        #plt.close()
        self.fig_FindingFeatures = plt.figure('Finding Features')
        plt.figtext(0.5, 0.95, str(os.path.basename(self.filepath)) +'\n# of frames for averaging: ' + str(NFramesForAveraging)
            ,ha='center', color='black', weight='bold', size=9)

        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(810,30,1000, 1000)
        #tp.annotate(f, np.rot90(frames[NframeTemp].T, 2))
        tp.annotate(f, framesSummed)
        
        self.xpos = []
        self.ypos = []
        for n in range(NFeatures):
            x_temp = f.values[n][0]
            y_temp = f.values[n][1]
            self.xpos.append(x_temp)
            self.ypos.append(y_temp)

        self.N_xPixel = frames.frame_shape[0]
        self.N_yPixel = frames.frame_shape[1]
        print 'self.N_xPixel: ', self.N_xPixel
        print 'self.N_yPixel: ', self.N_yPixel
         
        #self.xposImageJFound = self.N_yPixel - np.array(self.ypos)
        #self.yposImageJFound = np.array(self.xpos)
        
        self.xposImageJFound = np.array(self.xpos)
        self.yposImageJFound = self.N_xPixel - np.array(self.ypos)
        
        print 'xposImageJ', self.xposImageJFound, '\nyposImageJ', self.yposImageJFound

        print 'xpos' , self.xpos
        print 'ypos' , self.ypos
        
        
        #self.FindingPairs(self.xposImageJ, self.yposImageJ)
        
        
        
        plt.ion()
        plt.show()
        
        
        self.btn_FindAllPairs.Show()
        
        
    def FindingPairs(self, event):
        print '\n\nFindingPairs'

        moreXtemp = []
        moreYtemp = []
        for n in range(len(self.MorePairs_x)):
            if int(self.MorePairs_x[n].GetValue()) != 0:
                moreXtemp.append(int(self.MorePairs_x[n].GetValue()))
                moreYtemp.append(512 - int(self.MorePairs_y[n].GetValue()))
            
            
            
            
            
        xposImageJ = np.append(np.around(self.xposImageJFound, 0) , moreXtemp)
        yposImageJ = np.append(np.around(self.yposImageJFound, 0) , moreYtemp)
        
        print 'xposImageJ: ', xposImageJ
        print 'yposImageJ: ', yposImageJ
        
        

        dx = int(self.ImageJPair_x2.GetValue()) - int(self.ImageJPair_x1.GetValue())
        dy = int(self.ImageJPair_y2.GetValue()) - int(self.ImageJPair_y1.GetValue())
        print 'dx, dy = ', dx, dy
        
        AllPairX1 = []
        AllPairY1 = []
        AllPairX2 = []
        AllPairY2 = []
        for n in range(len(xposImageJ)):
            print 'xposImageJ[n]: ', xposImageJ[n]
            if left_Ch_x_left < xposImageJ[n] < left_Ch_x_right:
                if xposImageJ[n] in AllPairX1:
                    pass
                elif xposImageJ[n] - 1 in AllPairX1:
                    pass
                elif xposImageJ[n] + 1 in AllPairX1:
                    pass
                elif xposImageJ[n] - 2 in AllPairX1:
                    pass
                elif xposImageJ[n] + 2 in AllPairX1:
                    pass
                elif yposImageJ[n] - 1 in AllPairY1:
                    pass
                elif yposImageJ[n] + 1 in AllPairY1:
                    pass
                elif yposImageJ[n] - 2 in AllPairY1:
                    pass
                elif yposImageJ[n] + 2 in AllPairY1:
                    pass
                
                else:
                    AllPairX1.append(np.around(xposImageJ[n],0))
                    AllPairY1.append(np.around(yposImageJ[n],0))
                    AllPairX2.append(np.around(xposImageJ[n] + dx,0))
                    AllPairY2.append(np.around(yposImageJ[n] + dy,0))
                
            elif right_Ch_x_left < xposImageJ[n] < right_Ch_x_right:
                if xposImageJ[n] in AllPairX2:
                    pass
                elif xposImageJ[n] - 1 in AllPairX2:
                    pass
                elif xposImageJ[n] + 1 in AllPairX2:
                    pass
                elif xposImageJ[n] - 2 in AllPairX2:
                    pass
                elif xposImageJ[n] + 2 in AllPairX2:
                    pass
                elif yposImageJ[n] - 1 in AllPairY2:
                    pass
                elif yposImageJ[n] + 1 in AllPairY2:
                    pass
                elif yposImageJ[n] - 2 in AllPairY2:
                    pass
                elif yposImageJ[n] + 2 in AllPairY2:
                    pass
                else:
                    AllPairX1.append(np.around(xposImageJ[n] - dx,0))
                    AllPairY1.append(np.around(yposImageJ[n] - dy,0))
                    AllPairX2.append(np.around(xposImageJ[n],0))
                    AllPairY2.append(np.around(yposImageJ[n],0))
                
            else:
                pass
            
        print 'AllPairX1: ', AllPairX1
        print 'AllPairy1: ', AllPairY1
        print 'AllPairX2: ', AllPairX2
        print 'AllPairy2: ', AllPairY2
        
        self.ShowAllPairs(AllPairX1, AllPairY1, AllPairX2, AllPairY2)
        
        self.AllPairX1 = AllPairX1
        self.AllPairY1 = AllPairY1
        self.AllPairX2 = AllPairX2
        self.AllPairY2 = AllPairY2
        
        self.TotalNofPairs.SetLabel(str(len(self.AllPairX1)))
                
                
    def ShowAllPairs(self, AllPairX1, AllPairY1, AllPairX2, AllPairY2):
        colorScaleMin = int(self.ColorScaleMin.GetValue())
        colorScaleMax = int(self.ColorScaleMax.GetValue())
        
        DataFrameAll = tp.TiffStack(self.filepath)
        self.N_xPixel = DataFrameAll.frame_shape[0]
        self.N_yPixel = DataFrameAll.frame_shape[1]
        NframeTemp = int(self.Nframe.GetValue())
        DataFrameOne = DataFrameAll[NframeTemp]
        
        
        
        
        Nframes = DataFrameAll._count
        print 'Nframes = ', Nframes
        
        
        
   
        AnalysisConditions = ( 'Chosen f # ' + str(self.Nframe.GetValue()) + '   Ave # f: ' + str(self.NFramesForAveraging.GetValue()) + 
        '   f size: '+ str(self.FeatureSize.GetValue()) +  '   Min I: ' + str(self.MinIntensity.GetValue()) )


     
        NFeatures = len(AllPairX1)
        Nfig = int(  math.ceil((NFeatures/30.0)))
        
        TotalFeatureN_text = 'Total # ' + str(NFeatures)
        
        #self.fig_IntensityTrajectory = [[]] * Nfig
        self.fig_PairImages = [[] for _ in xrange(Nfig)]
  
        
        k = 0
        fn = 0
        for n in range(Nfig):

            print 'Pair image n = ', n
            self.fig_PairImages[n] = plt.figure('Pair_Images_'+ str(n), figsize = (18, 9))
            #fig_ShowAfewFrames.suptitle('test title', fontsize=20)
            self.fig_PairImages[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.25, hspace = 0.35)
            plt.figtext(0.5, 0.95, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath) + '\n' + AnalysisConditions   )
                ,ha='center', color='black', weight='normal', size='small')
            self.fig_PairImages[n].text(0.05, 0.96, TotalFeatureN_text, ha="left", va="bottom", size="medium",color="black", weight='bold')
            
         
            for m in range(30):

                xc1 = AllPairY1[m + fn]
                yc1 = self.N_xPixel - AllPairX1[m + fn]
                xc2 = AllPairY2[m + fn]
                yc2 = self.N_xPixel - AllPairX2[m + fn]
                
                
                PairImg_temp1 = DataFrameOne[yc1-10:yc1+11, xc1-10:xc1+11]
                PairImg_temp2 = DataFrameOne[yc2-10:yc2+11, xc2-10:xc2+11]
                
                #plt.hold(True)
                
                Pair1xy = '(' + str(AllPairX1[m + fn]) + ', ' + str(AllPairY1[m + fn]) + ')'
                Pair2xy = '(' + str(AllPairX2[m + fn]) + ', ' + str(AllPairY2[m + fn]) + ')'
       
                #print 'subplot # = ', m
               
                if m%2 == 0:
                    titleColor = 'red'
                else:
                    titleColor = 'green'
                
                plt.subplot(5,12, 2*m+1)
                plt.tick_params(labelsize=7)
                '''
                plt.title('Pair # ' + str(m + fn) + ',    (' + str(int(self.AllPairX1[m+fn])) + ', ' + str(int(self.AllPairY1[m+fn]))+')' 
                    +',  (' + str(int(self.AllPairX2[m+fn])) + ', ' + str(int(self.AllPairY2[m+fn]))+')' , fontsize=7
                    , color = titleColor, weight = 'bold')
                '''


                plt.title('Pair # ' + str(m + fn) + '\n' + Pair1xy, color = titleColor, size = 8, weight='bold')
                plt.imshow(np.rot90(PairImg_temp1,3),interpolation = 'None', origin='upper', vmin=colorScaleMin, vmax=colorScaleMax)
                
                plt.subplot(5,12,1 + 2*m + 1)
                plt.tick_params(labelsize=7)
                plt.title('Pair # ' + str(m + fn) + '\n' + Pair2xy, color = titleColor, size = 8, weight='bold')
                plt.imshow(np.rot90(PairImg_temp2,3),interpolation = 'None', origin='upper', vmin=colorScaleMin, vmax=colorScaleMax)
                
                k += 1
                if k%30 ==0:
                    fn += 30
                if k == NFeatures:
                    break
                
                


        plt.ion()
        plt.show()     

        self.btn_ShowAllPairsMovie.Show()
        self.btn_ShowAllPairsMovieSingleN.Show()
        self.btn_ShowLDTrajectory.Show()
        self.btn_ShowFIONAtracing.Show()        
    
    
    
    def ShowAllPairsMovie(self,event):
        
        MoviePairNumbers = []
        for n in range(4):
            MoviePairNumbers.append(int(self.MoviePairN[n].GetValue()))
            print 'MoviePairNumbers = ', MoviePairNumbers
            

        AllMoviePairX1 = []
        AllMoviePairY1 = []
        AllMoviePairX2 = []
        AllMoviePairY2 = []            
        for n in MoviePairNumbers:
            AllMoviePairX1.append(self.AllPairX1[n])
            AllMoviePairY1.append(self.AllPairY1[n])
            AllMoviePairX2.append(self.AllPairX2[n])
            AllMoviePairY2.append(self.AllPairY2[n])       
            
        self.aniRun = self.ShowAllPairsMovieAnimation(MoviePairNumbers, AllMoviePairX1, AllMoviePairY1, AllMoviePairX2, AllMoviePairY2)
        print 'aniRun type =', type(self.aniRun)
       
       
       
       
 

    def ShowAllPairsMovieSingleN(self,event):
        
        MoviePairNumbers = []
        for n in range(1):
            MoviePairNumbers.append(int(self.MoviePairNSingleN[n].GetValue()))
            print 'MoviePairNumbersSingleN = ', MoviePairNumbers
            

        AllMoviePairX1 = []
        AllMoviePairY1 = []
        AllMoviePairX2 = []
        AllMoviePairY2 = []            
        for n in MoviePairNumbers:
            AllMoviePairX1.append(self.AllPairX1[n])
            AllMoviePairY1.append(self.AllPairY1[n])
            AllMoviePairX2.append(self.AllPairX2[n])
            AllMoviePairY2.append(self.AllPairY2[n])       
            
        self.aniRun = self.ShowAllPairsMovieAnimation(MoviePairNumbers, AllMoviePairX1, AllMoviePairY1, AllMoviePairX2, AllMoviePairY2)
        print 'aniRun type =', type(self.aniRun)
        
 




       
        
    def ShowSinglePairsMovie(self,event):
        
        MoviePairNumbers = [0]
        

        AllMoviePairX1 = [int(self.SingleMovieImageJPair_x1.GetValue())]
        AllMoviePairY1 = [int(self.SingleMovieImageJPair_y1.GetValue())]
        AllMoviePairX2 = [int(self.SingleMovieImageJPair_x2.GetValue())]
        AllMoviePairY2 = [int(self.SingleMovieImageJPair_y2.GetValue())]            
    
            
        self.aniRunSingle = self.ShowAllPairsMovieAnimation(MoviePairNumbers, AllMoviePairX1, AllMoviePairY1, AllMoviePairX2, AllMoviePairY2)
        print 'aniRunSingle type =', type(self.aniRunSingle)
        
        
        
        
        
        
            
    def ShowAllPairsMovieAnimation(self, MoviePairNumbers, AllMoviePairX1, AllMoviePairY1, AllMoviePairX2, AllMoviePairY2):
        #plt.close('ShowAllPairs')
        '''
        AllMoviePairX1 = self.AllPairX1
        AllMoviePairY1 = self.AllPairY1
        AllMoviePairX2 = self.AllPairX2
        AllMoviePairY2 = self.AllPairY2
        '''
        
        colorScaleMin = int(self.ColorScaleMin.GetValue())
        colorScaleMax = int(self.ColorScaleMax.GetValue()) 
        
        vminValue = colorScaleMin
        vmaxValue = colorScaleMax
        
        
        print 'AllMoviePairX1 = ', AllMoviePairX1
        
        #colorScaleMin = 1000
        #colorScaleMax = 3000
        
        DataFrameAll = tp.TiffStack(self.filepath)
        self.N_xPixel = DataFrameAll.frame_shape[0]
        self.N_yPixel = DataFrameAll.frame_shape[1]
        frameNTemp = int(self.Nframe.GetValue())
        
       
        

        
        print 'vminValue ', vminValue
        print 'vmaxValue ', vmaxValue
        
                
        
        TotalNframes = DataFrameAll._count
        print 'Nframes = ', TotalNframes
        
        
        #plt.ion()
        self.FigShowAllPairsMovie = plt.figure('ShowAllPairsMovie')
        plt.figtext(0.5, 0.97, str(os.path.basename(self.filepath)),ha='center', color='black', weight='bold', size='medium')  
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(210,30,1600, 800)
        
        plt.subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.4, hspace = 0.5)
        
        

        #for m in np.arange(frameNTemp, TotalNframes, 1):
        ax = [None] * 32

        nPlusTemp = 0
        k = 0
        for n in range(len(AllMoviePairX1)):
            print 'Pair # = ', n
            print 'ax ', 2*n, 2*n +1
            #plt.hold(True)
            Pair1xy = '(' + str(AllMoviePairX1[n]) + ', ' + str(AllMoviePairY1[n]) + ')'
            Pair2xy = '(' + str(AllMoviePairX2[n]) + ', ' + str(AllMoviePairY2[n]) + ')'
            
            if n%2 == 0:
                titleColor = 'red'
            else:
                titleColor = 'green'
            
            
            if len(AllMoviePairX1) == 1:
                nPlusTemp = 4
                
            
            ax[2*n] = self.FigShowAllPairsMovie.add_subplot(4,8, 2*n + 1 + nPlusTemp)
            ax[2*n].set_title('Pair # ' + str(MoviePairNumbers[n]) + '\n' + Pair1xy, color = titleColor)
            ax[2*n + 1] = self.FigShowAllPairsMovie.add_subplot(4,8, 2*n + 2 + nPlusTemp)
            ax[2*n + 1 ].set_title('Pair # ' + str(MoviePairNumbers[n]) + '\n' + Pair2xy, color = titleColor)
            
            k += 1
            if k == 16:
                break


        movieDataAllFrame = []
        
        for m in np.arange(frameNTemp, TotalNframes, 1):
                    
            DataFrameOne = DataFrameAll[m]
            '''
            plt.ioff()
            plt.figure('test')
            plt.imshow(DataFrameOne)
            plt.show()
            '''
            print '\nframe # ', m

            k = 0
            movieDataTemp = []
            for n in range(len(AllMoviePairX1)):
                print 'Pair # = ', n
            
                                
                xc1 = AllMoviePairY1[n]
                yc1 = self.N_xPixel - AllMoviePairX1[n]
                xc2 = AllMoviePairY2[n]
                yc2 = self.N_xPixel - AllMoviePairX2[n]
                
                
                PairImg_temp1 = DataFrameOne[yc1-10:yc1+11, xc1-10:xc1+11]
                PairImg_temp2 = DataFrameOne[yc2-10:yc2+11, xc2-10:xc2+11]
                
                #print 'np.amax(PairImg_temp1.ravel())', np.amax(PairImg_temp1.ravel())
                #print 'np.amax(PairImg_temp1)', np.amax(PairImg_temp1)
                '''
                vmin1 = np.amin(PairImg_temp1)
                vmax1 = np.amax(PairImg_temp1)
                vmin2 = np.amin(PairImg_temp2)
                vmax2 = np.amax(PairImg_temp2)
                
                vminValue = np.amin([vmin1, vmin2])
                vmaxValue = np.amax([vmax1, vmax2])
                '''
                
                s1temp = ax[2*n].imshow(np.rot90(PairImg_temp1,3),interpolation = 'None', origin='upper', vmin= vminValue, vmax=vmaxValue)                
                #s1titleTemp = ax[2*n].text(0.02, 0.94, str(m), fontsize=12, color = 'white')
                #s1titleTemp = ax[2*n].annotate('F# '+str(m), xy=(1, 0), xycoords='axes fraction', fontsize=10, xytext=(0, -15), textcoords='offset points', ha='right', va='top')                
                s2temp = ax[2*n + 1].imshow(np.rot90(PairImg_temp2,3),interpolation = 'None', origin='upper', vmin= vminValue, vmax=vmaxValue)                                
                #s2titleTemp = ax[2*n + 1].text(0.02, 0.94, str(m), fontsize=12, color = 'white')
                s2titleTemp = ax[2*n + 1].annotate('Frame: '+str(m), xy=(1, 0), xycoords='axes fraction', fontsize=10, xytext=(0, -15), textcoords='offset points', ha='right', va='top')
                
                movieDataTemp.append(s1temp)
                #movieDataTemp.append(s1titleTemp)                
                movieDataTemp.append(s2temp)
                movieDataTemp.append(s2titleTemp)
                print 'movieDtaTmepLength = ', len(movieDataTemp)
                
                k += 1
                if k == 16:
                    break
            movieDataAllFrame.append(movieDataTemp)
            #print 'movieDataTemp = ', movieDataTemp
            #print 'movieDataTempLength = ', len(movieDataTemp)
            print 'frame # ', m, '\n'
            
        print 'movieDataTemp size = ', len(movieDataTemp)
            
        print 'movieDataAllFrame size = ', len(movieDataAllFrame)
        
        #print movieDataAllFrame
            
        ani = matplotlib.animation.ArtistAnimation(self.FigShowAllPairsMovie, movieDataAllFrame, interval=100, blit=False)
        ani.repeat = True
                
        plt.ion()                       
        plt.show()
            
        return ani   

        
        
    



    def ShowLDTrajectory(self, event):
        print 'running ShowLDTrajectory'
        
        Allxpos = self.AllPairX1 + self.AllPairX2
        Allypos = self.AllPairY1 + self.AllPairY2
        
        print 'self.AllPairX1', self.AllPairX1
        print 'type(Allxpos)', type(self.AllPairX1)
        
        print 'Allxpos', Allxpos
        print 'type(Allxpos)', type(Allxpos)



        ''' the centerMax intensity is average of center max 5 pixels   '''    
        centerIntensityAveTrajectory, centerIntensityMax5AveTrajectory, backgroundIntensityAveTrajectory = IntensityTrajectoryPlotsFunc(self.filepath, Allxpos, Allypos)

        
        
        LeftIntensity = np.array(centerIntensityMax5AveTrajectory[: len(Allxpos)/2])
        RightIntensity = np.array(centerIntensityMax5AveTrajectory[len(Allxpos)/2 :])
        
        LD = (LeftIntensity - RightIntensity)/(LeftIntensity + RightIntensity)
        
        #background subtracted LD
        
        LeftIntensityBS = np.array(centerIntensityMax5AveTrajectory[: len(Allxpos)/2]) - np.array(backgroundIntensityAveTrajectory[: len(Allxpos)/2])
        RightIntensityBS = np.array(centerIntensityMax5AveTrajectory[len(Allxpos)/2 :]) - np.array(backgroundIntensityAveTrajectory[: len(Allxpos)/2])
        
        LDBS = (LeftIntensityBS - RightIntensityBS)/(LeftIntensityBS + RightIntensityBS)
        LDBS_AveAll = np.average(LDBS, axis=1)
        LDBS_Ave10 = np.average(LDBS[:, 0:10], axis=1)
        LDBS_Ave1 = np.average(LDBS[:, 0:1], axis=1)
        



        AnalysisConditions = ( 'Chosen f # ' + str(self.Nframe.GetValue()) + '   Ave # f: ' + str(self.NFramesForAveraging.GetValue()) + 
        '   f size: '+ str(self.FeatureSize.GetValue()) +  '   Min I: ' + str(self.MinIntensity.GetValue()) 
        + '     Total # Frames:' + str(self.TotalNframes) )

        


        
        NFeatures = len(LeftIntensity)
        Nfig = int(  math.ceil((NFeatures/15.0)))
        
        TotalFeatureN_text = 'Total # ' + str(NFeatures)
        
        #self.fig_IntensityTrajectory = [[]] * Nfig
        self.fig_IntensityLDTrajectory = [[] for _ in xrange(Nfig)]
        
        k = 0
        fn = 0
        
        for n in range(Nfig):
            print 'Intensity Trajectory Nfig n = ', n
            self.fig_IntensityLDTrajectory[n] = plt.figure('IntensityLD_Trajectory_'+ str(n), figsize = (18, 9))
            #fig_ShowAfewFrames.suptitle('test title', fontsize=20)
            self.fig_IntensityLDTrajectory[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.25, hspace = 0.35)
            plt.figtext(0.5, 0.95, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath) + '\n' + AnalysisConditions   )
                ,ha='center', color='black', weight='normal', size='small')
            self.fig_IntensityLDTrajectory[n].text(0.05, 0.96, TotalFeatureN_text, ha="left", va="bottom", size="medium",color="black", weight='bold')
            
            
            for m in range(15):
                #print 'subplot # = ', m
                plt.tick_params(labelsize=9)
                
                if m%3 == 0:
                    titleColor = 'red'
                elif m%3 == 1:
                    titleColor = 'blue'
                else:
                    titleColor = 'green'
                
                plt.subplot(5,6, 2*m+1)
                plt.tick_params(labelsize=9)
                plt.title('Pair # ' + str(m + fn) + ',    (' + str(int(self.AllPairX1[m+fn])) + ', ' + str(int(self.AllPairY1[m+fn]))+')' 
                    +',  (' + str(int(self.AllPairX2[m+fn])) + ', ' + str(int(self.AllPairY2[m+fn]))+')' , fontsize=9
                    , color = titleColor, weight = 'bold')
                plt.plot(LeftIntensity[m + fn])
                plt.plot(RightIntensity[m + fn])
                plt.plot((LeftIntensity[m + fn] + RightIntensity[m + fn]), color = 'DeepPink')
                
                
                '''
                ax0 = plt.subplot(5,4, 2*m+1)
                ax0.set_title('Pair # ' + str(m + fn) + ',    (' + str(int(self.AllPairX1[m+fn])) + ', ' + str(int(self.AllPairY1[m+fn]))+')' 
                    +',  (' + str(int(self.AllPairX2[m+fn])) + ', ' + str(int(self.AllPairY2[m+fn]))+')' , fontsize=9)
                ax1 = ax0.twinx()
                ax0.plot(LeftIntensity[m + fn])
                ax0.plot(RightIntensity[m + fn])
                ax0.plot((LeftIntensity[m + fn] + RightIntensity[m + fn]), color = 'HotPink')
                ax1.plot(LD[m + fn], color = 'red')
                '''
                
                # Calculating AutoCorrelation
                LDacr = LD[m + fn]
                yunbiased = LDacr-np.mean(LDacr)
                ynorm = np.sum(yunbiased**2)
                acor = np.correlate(yunbiased, yunbiased, "full")/ynorm
                acor = acor[len(acor)/2:] # use only second half
                acor2 = np.append(acor[1:], [0])

                LDacrBS = LDBS[m + fn]
                yunbiasedBS = LDacrBS - np.mean(LDacrBS)
                ynormBS = np.sum(yunbiasedBS**2)
                acorBS = np.correlate(yunbiasedBS, yunbiasedBS, "full")/ynormBS
                acorBS = acorBS[len(acorBS)/2:] # use only second half
                acor2BS = np.append(acorBS[1:], [0])

                
                #print 'acor: \n', acor
                
                

                plt.subplot(5,6, 2*m+2)
                plt.tick_params(labelsize=9)
                plt.title('Pair # ' + str(m + fn) + ',    LD  ACOR', fontsize=9, color = titleColor, weight = 'bold')
                #plt.plot(LD[m + fn], color = 'orangered', label = 'LD')                
                #plt.plot(acor2, color = 'black', label = 'ACOR')
                plt.plot(LDBS[m + fn], color = 'orange', label = 'LDBS')
                plt.plot(acor2BS, color = 'maroon', label = 'ACORBS')
                plt.legend(loc='upper right',prop={'size':7}, frameon = False)
                plt.ylim(-1,1)
                
                
                k += 1
                if k%15 ==0:
                    fn += 15
                if k == NFeatures:
                    break
                
                
        self.fig_LDAveHist = plt.figure('LD average histogram', figsize = (12,6))
        self.fig_LDAveHist.text(0.05, 0.96, TotalFeatureN_text, ha="left", va="bottom", size="medium",color="black", weight='bold')
        plt.figtext(0.5, 0.92, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath) + '\n' + AnalysisConditions   )
                ,ha='center', color='black', weight='normal', size='small')
                
        plt.subplot(2,3,1)
        plt.hist(LDBS_Ave1)
        plt.tick_params(labelsize=9)
        plt.xlabel('1st Frame LD', fontsize = 11)
        plt.xlim(-1.0, 1.0)
        
        plt.subplot(2,3,2)
        plt.hist(LDBS_Ave10)
        plt.tick_params(labelsize=9)
        plt.xlabel('Ave. 10 Frames LD', fontsize = 11)
        plt.xlim(-1.0, 1.0)
        
        plt.subplot(2,3,3)
        plt.hist(LDBS_AveAll)
        plt.tick_params(labelsize=9)
        plt.xlabel('Ave. All Frames LD', fontsize = 11)
        plt.xlim(-1.0, 1.0)

        
        
        
        plt.show()

        
        
        

        print 'len(LD) ', len(LD)

        self.btn_Save_figures.Show()
        










    def ShowExcitationM(self, event):
        print 'running ShowExcitationM'
        
        Allxpos = self.AllPairX1 + self.AllPairX2
        Allypos = self.AllPairY1 + self.AllPairY2
        
        print 'self.AllPairX1', self.AllPairX1
        print 'type(Allxpos)', type(self.AllPairX1)
        
        print 'Allxpos', Allxpos
        print 'type(Allxpos)', type(Allxpos)



        ''' the center Max intensity is the average of center max 5 pixels   '''    
        centerIntensityAveTrajectory, centerIntensityMax5AveTrajectory, backgroundIntensityAveTrajectory, backgroundIntensityMax5AveTrajectory = IntensityTrajectoryPlotsFunc(self.filepath, Allxpos, Allypos)


        
        LeftIntensity = np.array(centerIntensityMax5AveTrajectory[: len(Allxpos)/2])
        RightIntensity = np.array(centerIntensityMax5AveTrajectory[len(Allxpos)/2 :])
        

        
        LeftIntensityBS = np.array(centerIntensityMax5AveTrajectory[: len(Allxpos)/2]) - np.array(backgroundIntensityAveTrajectory[: len(Allxpos)/2])
        RightIntensityBS = np.array(centerIntensityMax5AveTrajectory[len(Allxpos)/2 :]) - np.array(backgroundIntensityAveTrajectory[: len(Allxpos)/2])
        TotalIntensityBS = LeftIntensityBS + RightIntensityBS
        
        TotalIntensityBS_x = np.arange(len(TotalIntensityBS[0]))
        
        
        

        AnalysisConditions = ( 'Chosen f # ' + str(self.Nframe.GetValue()) + '   Ave # f: ' + str(self.NFramesForAveraging.GetValue()) + 
        '   f size: '+ str(self.FeatureSize.GetValue()) +  '   Min I: ' + str(self.MinIntensity.GetValue()) 
        + '     Total # Frames:' + str(self.TotalNframes) )

        




        if self.BG_ReferncePair.GetValue() != 'None':
            BG_ref_Pair_N = int(self.BG_ReferncePair.GetValue())
            popt, pcov = curve_fit(Sine_fit_func, TotalIntensityBS_x, TotalIntensityBS[BG_ref_Pair_N], p0=[100, 0.042, 0, 350])
            print 'popt\n', popt
            BG_ref_fit = Sine_fit_func(TotalIntensityBS_x, popt[0],popt[1],popt[2],popt[3])
            
            plt.figure()
            plt.plot(TotalIntensityBS_x, TotalIntensityBS[BG_ref_Pair_N])
            plt.plot(TotalIntensityBS_x, BG_ref_fit)
            plt.show()
            
            TotalIntensityBS = TotalIntensityBS - BG_ref_fit
            
            


        
        NFeatures = len(LeftIntensity)
        Nfig = int(  math.ceil((NFeatures/15.0)))
        
        TotalFeatureN_text = 'Total # ' + str(NFeatures)
        
        #self.fig_IntensityTrajectory = [[]] * Nfig
        self.fig_ExcitationM_Trajectory = [[] for _ in xrange(Nfig)]
        
        k = 0
        fn = 0
        
        for n in range(Nfig):
            print 'fig_ExcitationM_Trajectory Nfig n = ', n
            self.fig_ExcitationM_Trajectory[n] = plt.figure('ExcitationM_'+ str(n), figsize = (18, 9))
            #fig_ShowAfewFrames.suptitle('test title', fontsize=20)
            self.fig_ExcitationM_Trajectory[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.25, hspace = 0.35)
            plt.figtext(0.5, 0.95, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath) + '\n' + AnalysisConditions   )
                ,ha='center', color='black', weight='normal', size='small')
            self.fig_ExcitationM_Trajectory[n].text(0.05, 0.96, TotalFeatureN_text, ha="left", va="bottom", size="medium",color="black", weight='bold')
            
            
            for m in range(15):
                #print 'subplot # = ', m
                plt.tick_params(labelsize=9)
                
                if m%3 == 0:
                    titleColor = 'red'
                elif m%3 == 1:
                    titleColor = 'blue'
                else:
                    titleColor = 'green'
                

                
                
                plt.subplot(5,6, 2*m+1)
                plt.tick_params(labelsize=9)
                plt.title('Pair # ' + str(m + fn) + ',    (' + str(int(self.AllPairX1[m+fn])) + ', ' + str(int(self.AllPairY1[m+fn]))+')' 
                    +',  (' + str(int(self.AllPairX2[m+fn])) + ', ' + str(int(self.AllPairY2[m+fn]))+')' , fontsize=9
                    , color = titleColor, weight = 'bold')
                plt.plot(LeftIntensity[m + fn])
                plt.plot(RightIntensity[m + fn])
                plt.plot((LeftIntensity[m + fn] + RightIntensity[m + fn]), color = 'DeepPink')
                



                

                
                
                plt.subplot(5,6, 2*m+2)
                plt.tick_params(labelsize=9)
                plt.title('Pair # ' + str(m + fn) + ',  Total Intensity BS', fontsize=9, color = titleColor, weight = 'bold')

                plt.plot((TotalIntensityBS[m + fn]), color = 'OrangeRed')
                



                
                
                
                
                k += 1
                if k%15 ==0:
                    fn += 15
                if k == NFeatures:
                    break
                



        
        '''        
        self.fig_LDAveHist = plt.figure('LD average histogram', figsize = (12,6))
        self.fig_LDAveHist.text(0.05, 0.96, TotalFeatureN_text, ha="left", va="bottom", size="medium",color="black", weight='bold')
        plt.figtext(0.5, 0.92, str(self.filedirectoryOnly)+'\n' + str(os.path.basename(self.filepath) + '\n' + AnalysisConditions   )
                ,ha='center', color='black', weight='normal', size='small')
                
        plt.subplot(2,3,1)
        plt.hist(LDBS_Ave1)
        plt.tick_params(labelsize=9)
        plt.xlabel('1st Frame LD', fontsize = 11)
        plt.xlim(-1.0, 1.0)
        
        plt.subplot(2,3,2)
        plt.hist(LDBS_Ave10)
        plt.tick_params(labelsize=9)
        plt.xlabel('Ave. 10 Frames LD', fontsize = 11)
        plt.xlim(-1.0, 1.0)
        
        plt.subplot(2,3,3)
        plt.hist(LDBS_AveAll)
        plt.tick_params(labelsize=9)
        plt.xlabel('Ave. All Frames LD', fontsize = 11)
        plt.xlim(-1.0, 1.0)
        '''
        
        
        
        plt.show()



        self.btn_Save_figures.Show()
        







    def ShowFIONA_tracing(self, event):
        print '\nrunning ShowFIONA_trace'
        
        FIONA_pair_N_for_processing = []
        FIONA_pair_x_L = []
        FIONA_pair_y_L = []
        FIONA_pair_x_R = []
        FIONA_pair_y_R = []
        
        for n in range(len(self.FIONA_pair_N)):
            if self.FIONA_pair_N[n].GetValue() != '-':
                FIONA_pair_N_for_processing.append(int(self.FIONA_pair_N[n].GetValue()))
                FIONA_pair_x_L.append(self.AllPairX1[int(self.FIONA_pair_N[n].GetValue())])
                FIONA_pair_y_L.append(self.AllPairY1[int(self.FIONA_pair_N[n].GetValue())])
                FIONA_pair_x_R.append(self.AllPairX2[int(self.FIONA_pair_N[n].GetValue())])
                FIONA_pair_y_R.append(self.AllPairY2[int(self.FIONA_pair_N[n].GetValue())])
                
        xpos = np.append(FIONA_pair_x_L, FIONA_pair_x_R)
        ypos = np.append(FIONA_pair_y_L, FIONA_pair_y_R)
                
        

        

        CenterImages, fx0, fy0, self.data_FIONA, fx1,fy1, fx2, fy2, FeccentricityData = FIONA_FramesByFrame(self.filepath, xpos, ypos,
            int(self.frameStart.GetValue()), int(self.frameEnd.GetValue()) , float(self.SvsB_Tolerance.GetValue()), float(self.MaxEccentricity.GetValue())
            , float(self.MaxPercentage_Tolerance.GetValue()), float(self.MaxPixel_Tolerance.GetValue()) )

    

      
     
        
        CeterImage_Left = np.array(CenterImages[: len(xpos)/2])
        CeterImage_Right = np.array(CenterImages[len(xpos)/2 :])
        
        fx0_Left = np.array(fx0[: len(xpos)/2])
        fy0_Left = np.array(fy0[: len(xpos)/2])
        fx0_Right = np.array(fx0[len(xpos)/2 :])
        fy0_Right = np.array(fy0[len(xpos)/2 :])
        
        fx1_Left = np.array(fx1[: len(xpos)/2])
        fy1_Left = np.array(fy1[: len(xpos)/2])
        fx1_Right = np.array(fx1[len(xpos)/2 :])
        fy1_Right = np.array(fy1[len(xpos)/2 :])


        fx2_Left = np.array(fx2[: len(xpos)/2])
        fy2_Left = np.array(fy2[: len(xpos)/2])
        fx2_Right = np.array(fx2[len(xpos)/2 :])
        fy2_Right = np.array(fy2[len(xpos)/2 :])

        FionaEcc_Left = np.array(FeccentricityData[: len(xpos)/2])
        FionaEcc_Right = np.array(FeccentricityData[len(xpos)/2 :])





        AnalysisConditions = ( 'Chosen f # ' + str(self.Nframe.GetValue()) + '   Ave # f: ' + str(self.NFramesForAveraging.GetValue()) + 
        '   f size: '+ str(self.FeatureSize.GetValue()) +  '   Min I: ' + str(self.MinIntensity.GetValue()) 
        + '     Total # Frames:' + str(self.TotalNframes) )
        
        
        
        

        NFeatures = len(xpos)/2
      
        NCfig = int(  math.ceil((NFeatures/6.0)))
        self.fig_FIONA_2CH = [[] for _ in xrange(NCfig)]
        
        
        
        self.fig_FIONA_2CH[0] = plt.figure('FIONA_2CH', figsize = (18.6, 9))
        
        x3 = np.linspace(0, 9, 10)
        y3 = np.linspace(0, 9, 10)
        x3, y3 = np.meshgrid(x3, y3)
        
        
        gs = gridspec.GridSpec(3, 7)
        gs.update(left=0.05, right=0.95, wspace=0.18, hspace=0.42)
        plt.suptitle(str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly) + '\n' + AnalysisConditions
            , size=10)
 

        
        
        
        ax2_0 = plt.subplot(gs[0,0])
        plt.tick_params(labelsize=9)
        ax2_0.set_title('F# 0', size =10)
        ax2_0.imshow(CeterImage_Left[0][0].reshape(9, 9), interpolation = 'None', cmap=plt.cm.jet, origin='bottom',
                                   extent=(x3.min(), x3.max(), y3.min(), y3.max()) )

        ax2_10 = plt.subplot(gs[0,1])
        plt.tick_params(labelsize=9)
        ax2_10.set_title('F# 0', size =10)
        ax2_10.imshow(CeterImage_Right[0][0].reshape(9, 9), interpolation = 'None', cmap=plt.cm.jet, origin='bottom',
                                   extent=(x3.min(), x3.max(), y3.min(), y3.max()) )

        ax4 = plt.subplot(gs[1:2,0:1])
        plt.tick_params(labelsize=9)
        ax4.set_title('FIONA_0 Left', size =10)
        ax4.scatter(fx0_Left[0], fy0_Left[0], marker = 'x', c = range(len(fx0_Left[0])),  s=35, vmin=0, vmax= len(fx0_Left[0]))
        plt.xlim( (np.around(FIONA_pair_x_L[0]) - 2, np.around(FIONA_pair_x_L[0]) + 3 ) )
        plt.ylim( (np.around(FIONA_pair_y_L[0]) - 2, np.around(FIONA_pair_y_L[0]) + 3 ) )     

        ax5 = plt.subplot(gs[1:2,1:2])
        plt.tick_params(labelsize=9)
        ax5.set_title('FIONA_0 Right', size =10)
        ax5.scatter(fx0_Right[0], fy0_Right[0], marker = 'x', c = range(len(fx0_Right[0])),  s=35, vmin=0, vmax= len(fx0_Right[0]))
        plt.xlim( (np.around(FIONA_pair_x_R[0]) - 2, np.around(FIONA_pair_x_R[0]) + 3 ) )
        plt.ylim( (np.around(FIONA_pair_y_R[0]) - 2, np.around(FIONA_pair_y_R[0]) + 3 ) )     
        
        


        
        ax2_0 = plt.subplot(gs[0,2])
        plt.tick_params(labelsize=9)
        ax2_0.set_title('F# 0', size =10)
        ax2_0.imshow(CeterImage_Left[1][0].reshape(9, 9), interpolation = 'None', cmap=plt.cm.jet, origin='bottom',
                                   extent=(x3.min(), x3.max(), y3.min(), y3.max()) )
        
        ax2_10 = plt.subplot(gs[0,3])
        plt.tick_params(labelsize=9)
        ax2_10.set_title('F# 0', size =10)
        ax2_10.imshow(CeterImage_Right[1][0].reshape(9, 9), interpolation = 'None', cmap=plt.cm.jet, origin='bottom',
                                   extent=(x3.min(), x3.max(), y3.min(), y3.max()) )
                                   
        ax4 = plt.subplot(gs[1:2,2:3])
        plt.tick_params(labelsize=9)
        ax4.set_title('FIONA_0 Left', size =10)
        ax4.scatter(fx0_Left[1], fy0_Left[1], marker = 'x', c = range(len(fx0_Left[1])),  s=35, vmin=0, vmax= len(fx0_Left[1]))
        plt.xlim( (np.around(FIONA_pair_x_L[1]) - 2, np.around(FIONA_pair_x_L[1]) + 3 ) )
        plt.ylim( (np.around(FIONA_pair_y_L[1]) - 2, np.around(FIONA_pair_y_L[1]) + 3 ) )     
        
        ax5 = plt.subplot(gs[1:2,3:4])
        plt.tick_params(labelsize=9)
        ax5.set_title('FIONA_0 Right', size =10)
        ax5.scatter(fx0_Right[1], fy0_Right[1], marker = 'x', c = range(len(fx0_Right[1])),  s=35, vmin=0, vmax= len(fx0_Right[1]))
        plt.xlim( (np.around(FIONA_pair_x_R[1]) - 2, np.around(FIONA_pair_x_R[1]) + 3 ) )
        plt.ylim( (np.around(FIONA_pair_y_R[1]) - 2, np.around(FIONA_pair_y_R[1]) + 3 ) )   
        












           
                
                
      
    def SaveFigures(self, event):
        prefix = self.FilePrefix.GetValue()
        figFileName = prefix + self.filenameOnly[:40]
        
        
        todayDate = time.strftime("%Y%m%d")
        print "\nToday's Date: ", time.strftime("%Y-%m-%d")
        
        self.fig_ShowAfewFrames.savefig(self.filedirectoryOnly + '\\' + figFileName + '_' + todayDate + '_1_RawImageFrames.png')
        print 'Raw Frames figures are saved'
        
        self.fig_FindingFeatures.savefig(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_2_Molecules.png')
        print 'Feature finding figure is saved'
        
        save_fig_error = 0
        


                

        try:
            for n in range(len(self.fig_PairImages)):
                self.fig_PairImages[n].savefig(self.filedirectoryOnly + '\\' + figFileName + '_' + todayDate + '_3_Pairs_' + str(n)+ '.png')
                #plt.close(self.fig_PairImages[n])
            print 'Pair Image figures are saved'
        except:
            print 'Pair Image figures are NOT saved'      
            save_fig_error += 1
            


            
        try:        
            for n in range(len(self.fig_IntensityLDTrajectory)):    
                self.fig_IntensityLDTrajectory[n].savefig(self.filedirectoryOnly + '\\' + figFileName + '_' + todayDate + '_4_LD_' + str(n) + '.png')
                #plt.close(self.fig_IntensityLDTrajectory[n])
            self.fig_LDAveHist.savefig(self.filedirectoryOnly + '\\' + figFileName + '_' + todayDate + '_4_LD_Histogram.png')
            print 'LD figures are saved'
            
        except:
            print 'LD figures are NOT saved: error'  
            save_fig_error += 1
            
            

        
        '''
        try:        
            for n in range(len(self.fig_ExcitationM_Trajectory)):    
                self.fig_ExcitationM_Trajectory[n].savefig(self.filedirectoryOnly + '\\' + figFileName + '_' + todayDate + '_5_ExcitationM_' + str(n) + '.png')
                #plt.close(self.fig_ExcitationM_Trajectory[n])
            #self.fig_LDAveHist.savefig(self.filedirectoryOnly + '\\' + figFileName + '_' + todayDate + '_4_LD_Histogram.png')
            print 'fig_ExcitationM_Trajectory figures are saved'
            
        except:
            print 'fig_ExcitationM_Trajectory figures are NOT saved: error'  
            save_fig_error += 1
        '''





        
        
        if save_fig_error == 0 :
            save_fig_text = 'All figures are saved'
            print '\n', save_fig_text, '\n'
        else:
            save_fig_text = 'NOT All figures are saved'
            print '\n NOT All figures are saved\n # error: ', save_fig_error
            
        
        wx.StaticText(self.MainPanel, -1, save_fig_text, pos=(160, 723)) 
     
     
     



  
        

  
def IntensityCalculationForSingleNumpyFrame(loadedframe, xpos, ypos, IntensityMethod, centerLength):
    print '\n Starting IntensityCalculationForSingleNumpyFrame' 
    ################################################
    x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
    
    ################################################

    print 'x_ImageJ', x_ImageJ
    print 'y_ImageJ', y_ImageJ
    print 'intensity method', IntensityMethod
    
    print '\ncenterLength', centerLength
    
    
    frame = loadedframe
    
    print 'frame\n', frame
    

    N_xPixel = len(frame)
    N_yPixel = len(frame[0])
    
    
    print N_xPixel
    print 'type(N_xPixel)', type(N_xPixel)    
    print N_yPixel
    print 'type(N_yPixel)', type(N_yPixel)
    
    
    #xc = y_ImageJ - 0
    #yc = N_xPixel - x_ImageJ - 1
    
    xc = x_ImageJ - 0 #trackpy 0.3
    yc = y_ImageJ  #trackpy 0.3
    


    
    
    centerIntensity = []
    backgroundIntensity = []
    
    
    
    CenterEdgeLengthHalf = int((centerLength-1)/2)
    BGEdgeLengthHalf = int((centerLength+1)/2)
    

    edgeArray = frame[yc[0]-BGEdgeLengthHalf:yc[0]+BGEdgeLengthHalf+1, xc[0]-BGEdgeLengthHalf:xc[0]+BGEdgeLengthHalf+1] * 0
    print 'edgeArray zeros try \n', edgeArray
    for n in range(centerLength+2):
        for k in range(centerLength+2):
            if n == 0 or n == BGEdgeLengthHalf*2 or k == 0 or k==BGEdgeLengthHalf*2:
                edgeArray[n][k] = 1
    
                    
    print 'edgeArray \n ' , edgeArray
    
    



           


    N_centerPixel = centerLength**2
    N_BG_edgePixel = centerLength*2 + (centerLength+2)*2
    N_diff_center_edge = N_centerPixel - N_BG_edgePixel
    
    print 'N_centerPixel : ', N_centerPixel
    print 'N_BG_edgePixel : ', N_BG_edgePixel
    print 'N_diff_center_edge : ', N_diff_center_edge
 


    for n in range(len(xc)):
        frameLargeTemp = frame[yc[n]-BGEdgeLengthHalf:yc[n]+BGEdgeLengthHalf+1, xc[n]-BGEdgeLengthHalf:xc[n]+BGEdgeLengthHalf+1]
        frameCenterTemp = frame[yc[n]-CenterEdgeLengthHalf:yc[n]+CenterEdgeLengthHalf+1, xc[n]-CenterEdgeLengthHalf:xc[n]+CenterEdgeLengthHalf+1]
        
        '''
        plt.figure()
        plt.subplot(121)
        plt.imshow(frameLargeTemp, interpolation = 'None')
        plt.subplot(122)
        plt.imshow(frameCenterTemp, interpolation = 'None')
        plt.ioff()
        plt.show()
        #'''
        

        frameBackgroundTemp = frameLargeTemp * edgeArray
            
        
        #print 'frameBackgroundTemp \n', frameBackgroundTemp
        
        #print 'B_Max5Ave: ', int(np.mean( np.sort(frameBackgroundTemp.ravel())[-5:]) )
        
        frameBackgroundTemp_EdgePixels = np.sort(frameBackgroundTemp.ravel())[-N_BG_edgePixel:]

        #print 'frameBackgroundTemp ', frameBackgroundTemp
        #print 'frameBackgroundTemp_32pixels ', frameBackgroundTemp_32pixels
        
        frameBackgroundTemp_EdgePixels_adjusted = frameBackgroundTemp_EdgePixels.tolist() # + np.random.choice(frameBackgroundTemp_EdgePixels, size=N_diff_center_edge).tolist()
        
        
        
        #print 'frameBackgroundTemp_49pixels ', frameBackgroundTemp_49pixels
        #print ' ', 
        
        
        
        if IntensityMethod == 0: # max5mean 
            centerIntensity.append( int(np.mean( np.sort(frameCenterTemp.ravel())[-5:]) ))
            backgroundIntensity.append( int(np.mean( np.sort(frameBackgroundTemp_EdgePixels_adjusted)[-5:])) )
            
        
        elif IntensityMethod == 1: # max5median
            
            centerIntensity.append( int(np.median( np.sort(frameCenterTemp.ravel())[-5:]) ))
            backgroundIntensity.append( int(np.median( np.sort(frameBackgroundTemp_EdgePixels_adjusted)[-5:])) )
            
        
        elif IntensityMethod == 2: # max1 
            centerIntensity.append( int(np.amax( np.sort(frameCenterTemp.ravel())[-5:])) )
            backgroundIntensity.append( int(np.amax( np.sort(frameBackgroundTemp_EdgePixels_adjusted)[-5:])) )
            
            
        elif IntensityMethod == 3: # sum 
            
            centerIntensity.append( int(np.sum(frameCenterTemp)) )
            backgroundIntensity.append( int(np.median(frameBackgroundTemp_EdgePixels)*N_centerPixel) )
            
            
                
                
                
                
                
    print 'End of IntensityCalculationForSingleNumpyFrame \n'
    
    
    
    return centerIntensity, backgroundIntensity
    
    

      
      
      
      
      
      
      
     
     
     
     
     
'''     
###########################################################################################################################
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################
###########################################################################################################################
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
'''


     
class SingleChannelMainTracking(wx.Frame):
    #ImageFilePath = ''
 
    #----------------------------------------------------------------------
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, "M Measurements", size=(1000, 1000), pos=(0,0))
        self.MainPanel = wx.ScrolledWindow(self, -1)
        self.MainPanel.SetScrollbars(1, 1, 800, 900)
        self.MainPanel.Show()

        
        self.MainPanel.Bind(wx.EVT_PAINT, self.on_paint)
        
        
        



        self.btn_SwitchToDataAnalysis = wx.Button(self.MainPanel, pos=(250,5), label="Switch To Data Analysis")
        self.btn_SwitchToDataAnalysis.Bind(wx.EVT_BUTTON, self.SwitchToDataAnalysis)
        


        self.btn_SwitchTo2chAnalysis = wx.Button(self.MainPanel, pos=(500,5), label="Switch To 2 Ch Analysis")
        self.btn_SwitchTo2chAnalysis.Bind(wx.EVT_BUTTON, self.SwitchTo2chAnalysis)
        

        
        self.btn_Manual = wx.Button(self.MainPanel, pos=(650,5), label="Manual")
        self.btn_Manual.Bind(wx.EVT_BUTTON, self.ShowManual)
        self.btn_Manual.Show()
        
        




              
              
              
              


        wx.StaticText(self.MainPanel, -1, 'SM Single Channel Data Analysis', pos=(50, 10))
        

                
        btn2 = wx.Button(self.MainPanel, pos=(50,30), label="Open Image File")
        self.ImageFilePath = btn2.Bind(wx.EVT_BUTTON, self.onOpenImageFile)
        self.btn2text = wx.StaticText(self.MainPanel, -1, str(self.ImageFilePath), pos=(170, 33))
        
        wx.StaticText(self.MainPanel, -1, "Choose Frame #s:", pos=(170, 73))
        self.Nframe = wx.TextCtrl(self.MainPanel, -1, "0", pos=(280, 70), size=(40,-1))
        self.Nframe2 = wx.TextCtrl(self.MainPanel, -1, "44", pos=(340, 70), size=(40,-1))
        
        
        self.btn_FindingFeature = wx.Button(self.MainPanel, pos=(50,80), label="Find Features")
        self.btn_FindingFeature.Bind(wx.EVT_BUTTON, self.FindingFeatureFunc)
        self.btn_FindingFeature.Hide()

        
        wx.StaticText(self.MainPanel, -1, "Average # frames:", pos=(170, 103))
        self.NFramesForAveraging = wx.TextCtrl(self.MainPanel, -1, "2", pos=(280, 100), size=(40,-1))
        
        wx.StaticText(self.MainPanel, -1, "Feature Size:", pos=(400, 103))
        self.FeatureSize = wx.TextCtrl(self.MainPanel, -1, str(PixelSize), pos=(480, 100), size=(40,-1))
        
        wx.StaticText(self.MainPanel, -1, "Min Intensity:", pos=(600, 103))
        self.MinIntensity = wx.TextCtrl(self.MainPanel, -1, str(ThrethholdIntensity), pos=(680, 100), size=(60,-1))
        
        

        
        self.btn_SaveFindingFeatureData = wx.Button(self.MainPanel, pos=(800,100), label="Save Find-Features Data")
        self.btn_SaveFindingFeatureData.Bind(wx.EVT_BUTTON, self.SaveFindFeatureData)
        self.btn_SaveFindingFeatureData.Hide()
        


        
        
        wx.StaticText(self.MainPanel, -1, "   # of \nfeatures", pos=(5, 73))
        self.NofFeatures = wx.StaticText(self.MainPanel, -1, "_       ", pos=(15, 100), style=5)
        


        self.NofFeatures_intensity_filtered = wx.StaticText(self.MainPanel, -1, "_      ", pos=(15, 140), style=5)
        

        self.btn_FindingFeature_IntensityFilter = wx.Button(self.MainPanel, pos=(50,130), label="Apply Intensity Filter")
        self.btn_FindingFeature_IntensityFilter.Bind(wx.EVT_BUTTON, self.FindingFeature_ApplyIntensityFilter)
        self.btn_FindingFeature_IntensityFilter.Hide()

        
        wx.StaticText(self.MainPanel, -1, "Intensity Filter Threshold (by the methods below):", pos=(200, 133))
        self.IntensityFilterThreshold = wx.TextCtrl(self.MainPanel, -1, "50", pos=(480, 130), size=(40,-1))
        


        
        self.btn_SaveFindingFeatureData_filtered = wx.Button(self.MainPanel, pos=(800,130), label="Save Find-Features-Filtered Data")
        self.btn_SaveFindingFeatureData_filtered.Bind(wx.EVT_BUTTON, self.SaveFindFeatureData_filtered)
        self.btn_SaveFindingFeatureData_filtered.Hide()
        


        
        
        
        self.MorePairs_x = [[] for _ in xrange(10)] 
        self.MorePairs_y = [[] for _ in xrange(10)]
        
        wx.StaticText(self.MainPanel, -1, '        More (x, y)\n in Finding Features', pos=(70, 175))
        self.MorePairs_x[0] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(200, 170), size=(30,-1))
        self.MorePairs_y[0] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(233, 170), size=(30,-1))
        
        self.MorePairs_x[1] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(300, 170), size=(30,-1))
        self.MorePairs_y[1] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(333, 170), size=(30,-1))
        
        self.MorePairs_x[2] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(400, 170), size=(30,-1))
        self.MorePairs_y[2] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(433, 170), size=(30,-1))
        
        self.MorePairs_x[3] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(500, 170), size=(30,-1))
        self.MorePairs_y[3] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(533, 170), size=(30,-1))

        self.MorePairs_x[4] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(600, 170), size=(30,-1))        
        self.MorePairs_y[4] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(633, 170), size=(30,-1))
        
        

        self.MorePairs_x[5] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(200, 200), size=(30,-1))        
        self.MorePairs_y[5] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(233, 200), size=(30,-1))
                        
        self.MorePairs_x[6] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(300, 200), size=(30,-1))
        self.MorePairs_y[6] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(333, 200), size=(30,-1))

        self.MorePairs_x[7] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(400, 200), size=(30,-1))
        self.MorePairs_y[7] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(433, 200), size=(30,-1))

        self.MorePairs_x[8] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(500, 200), size=(30,-1))
        self.MorePairs_y[8] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(533, 200), size=(30,-1))

        self.MorePairs_x[9] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(600, 200), size=(30,-1))
        self.MorePairs_y[9] = wx.TextCtrl(self.MainPanel, -1, "0", pos=(633, 200), size=(30,-1))





        self.cb_apply_xyLimits_yn = wx.CheckBox(self.MainPanel, -1, '  Apply (x,y) Limits', (50, 227))
        self.cb_apply_xyLimits_yn.SetValue(False)
        self.cb_apply_xyLimits_yn.Show()
        
        wx.StaticText(self.MainPanel, -1, 'in Finding Features', pos=(70, 240))
        wx.StaticText(self.MainPanel, -1, 'x:', pos=(200, 233))
        wx.StaticText(self.MainPanel, -1, '-', pos=(260, 233))
        self.x_Limit_Min = wx.TextCtrl(self.MainPanel, -1, "100", pos=(220, 230), size=(30,-1))
        self.x_Limit_Max = wx.TextCtrl(self.MainPanel, -1, "420", pos=(270, 230), size=(30,-1))
        wx.StaticText(self.MainPanel, -1, 'y:', pos=(330, 233))
        wx.StaticText(self.MainPanel, -1, '-', pos=(390, 233))
        self.y_Limit_Min = wx.TextCtrl(self.MainPanel, -1, "100", pos=(350, 230), size=(30,-1))
        self.y_Limit_Max = wx.TextCtrl(self.MainPanel, -1, "420", pos=(400, 230), size=(30,-1))






        
        
        
        self.btn_LoadFeatureXY = wx.Button(self.MainPanel, pos=(50,263), label="Load (x,y)")
        self.btn_LoadFeatureXY.Bind(wx.EVT_BUTTON, self.LoadFeatureXY)
        self.btn_LoadFeatureXY.Hide()
        
        
        self.text_LoadFeatureXY = wx.StaticText(self.MainPanel, -1, '_      ', pos=(15, 265))
        self.text_LoadFeatureXY_filename = wx.StaticText(self.MainPanel, -1, ' ', pos=(140, 265))
        
        
        
        
        self.cb_use_loaded_xy_yn = wx.CheckBox(self.MainPanel, -1, 'Use loaded x,y', (50, 303))
        self.cb_use_loaded_xy_yn.SetValue(False)
        self.cb_use_loaded_xy_yn.Show()
        
        
        self.cb_xy_offset_yn = wx.CheckBox(self.MainPanel, -1, 'x,y offsets', (150, 303))
        self.cb_xy_offset_yn.SetValue(False)
        self.cb_xy_offset_yn.Show()
        

        
        #wx.StaticText(self.MainPanel, -1, 'Off sets:', pos=(190, 303))
        #wx.StaticText(self.MainPanel, -1, 'y ', pos=(320, 283))
        self.OffSet_y = wx.SpinCtrl(self.MainPanel, -1, pos=(230, 300), size=(50,-1))
        self.OffSet_y.SetRange(-100, 100)  
        self.OffSet_y.SetValue(0)        
        self.OffSet_x = wx.SpinCtrl(self.MainPanel, -1, pos=(330, 300), size=(50,-1))
        self.OffSet_x.SetRange(-100, 100)  
        self.OffSet_x.SetValue(0)
        wx.StaticText(self.MainPanel, -1, 'move up                    move right ', pos=(230, 323))
        #wx.StaticText(self.MainPanel, -1, 'move right -> ImageJ move up \nmove up -> ImageJ move right', pos=(430, 303))
        
        
        
        
        



        self.btn_ShowRawImages = wx.Button(self.MainPanel, pos=(50,350), label="Show Raw Images")
        self.btn_ShowRawImages.Bind(wx.EVT_BUTTON, self.ShowRawImages)
        self.btn_ShowRawImages.Hide()
        wx.StaticText(self.MainPanel, -1, "Raw Data Images F#:", pos=(180, 354))
        self.rawImageFN = wx.SpinCtrl(self.MainPanel, -1, pos=(300, 350), size=(50,-1))
        self.rawImageFN.SetRange(0, 100000)  
        self.rawImageFN.SetValue(0)
        
        
        


        self.BG_ref_x = [[] for _ in xrange(4)] 
        self.BG_ref_y = [[] for _ in xrange(4)]
        
        wx.StaticText(self.MainPanel, -1, 'Background reference (x, y)\n       in Finding Features', pos=(435, 360))
        self.BG_ref_x[0] = wx.TextCtrl(self.MainPanel, -1, "200", pos=(600, 350), size=(30,-1))
        self.BG_ref_y[0] = wx.TextCtrl(self.MainPanel, -1, "400", pos=(633, 350), size=(30,-1))
        
        self.BG_ref_x[1] = wx.TextCtrl(self.MainPanel, -1, "250", pos=(700, 350), size=(30,-1))
        self.BG_ref_y[1] = wx.TextCtrl(self.MainPanel, -1, "400", pos=(733, 350), size=(30,-1))
        
        self.BG_ref_x[2] = wx.TextCtrl(self.MainPanel, -1, "300", pos=(600, 370), size=(30,-1))
        self.BG_ref_y[2] = wx.TextCtrl(self.MainPanel, -1, "400", pos=(633, 370), size=(30,-1))
        
        self.BG_ref_x[3] = wx.TextCtrl(self.MainPanel, -1, "300", pos=(700, 370), size=(30,-1))
        self.BG_ref_y[3] = wx.TextCtrl(self.MainPanel, -1, "300", pos=(733, 370), size=(30,-1))



        self.radio_IntensityMethod = wx.RadioBox(self.MainPanel, -1, choices=['Max 5 Average', 'Max 5 Median', 'Max 1', 'Sum'], label='Intensity Method', pos=(50, 420), size=wx.Size(120, 115), style=wx.RA_SPECIFY_ROWS)
        self.radio_IntensityMethod.SetSelection(0)
        
        self.radio_BSMethod = wx.RadioBox(self.MainPanel, -1, choices=['frame by frame', 'Fit Ave'], label='Background Subtraction Method', pos=(250, 420), size=wx.Size(180, 98), style=wx.RA_SPECIFY_ROWS)
        self.radio_BSMethod.SetSelection(0)
        
        
                
        
        


        
        self.cb_Excitation_Pol_correction_yn = wx.CheckBox(self.MainPanel, -1, 'Excitation Polarization Correction by', (50, 450))
        self.cb_Excitation_Pol_correction_yn.SetValue(False)
        self.cb_Excitation_Pol_correction_yn.Hide()
        self.Excitation_Pol_correction_factor = wx.TextCtrl(self.MainPanel, -1, "0.1", pos=(280, 476), size=(30, -1))
        self.Excitation_Pol_correction_factor.Hide()


        
        
        
        self.radio_FitOptions = wx.RadioBox(self.MainPanel, -1, choices=['Modulation Depth M                                       , ', 'Median intensity y0 only'], label='Fit options', pos=(450, 420), size=wx.Size(510, 90), style=wx.RA_SPECIFY_COLS)
        self.radio_FitOptions.SetSelection(0)
                
                
        wx.StaticText(self.MainPanel, -1, 'Fit initial guess for period in frames:', pos=(460, 463))  
        self.fitInitialGuess_T = wx.TextCtrl(self.MainPanel, -1, str(FitPeriodGuess), pos=(650, 460), size=(30, -1))
        
        #wx.StaticText(self.MainPanel, -1, 'Threshold T Fit Goodness % ', pos=(460, 483))
        self.cb_MaxTFitGoodnessErr_yn = wx.CheckBox(self.MainPanel, -1, 'Threshold T Fit Difference %', (460, 483))
        self.cb_MaxTFitGoodnessErr_yn.SetValue(True)
        self.cb_MaxTFitGoodnessErr_yn.Show()
        self.MaxTFitGoodnessErr = wx.TextCtrl(self.MainPanel, -1, str(10), pos=(630, 480), size=(40, -1))
    


    
        #wx.StaticText(self.MainPanel, -1, 'Threshold average intensity ', pos=(460, 503))  
        self.cb_MinAveIntensity_yn = wx.CheckBox(self.MainPanel, -1, 'Threshold average intensity', (460, 523))
        self.cb_MinAveIntensity_yn.SetValue(True)
        self.cb_MinAveIntensity_yn.Show()
        self.MinAveIntensity = wx.TextCtrl(self.MainPanel, -1, str(ThresholdAveIntensity), pos=(630, 520), size=(40, -1))
        #wx.StaticText(self.MainPanel, -1, 'Threshold R^2 Fit Goodness ', pos=(460, 523))  
        #self.MinRsquareFitGoodness = wx.TextCtrl(self.MainPanel, -1, str(MinRsquareFitGoodness), pos=(620, 520), size=(50, -1))
        
        
        

        self.btn_ShowAnisotropyM = wx.Button(self.MainPanel, pos=(50,550), label="Show Anisotropy M")
        self.btn_ShowAnisotropyM.Bind(wx.EVT_BUTTON, self.ShowAnisotropyM)
        self.btn_ShowAnisotropyM.Hide()
        #wx.StaticText(self.MainPanel, -1, "Anisotropy M & Intensity", pos=(180, 554))

        self.text_NFeatures = wx.StaticText(self.MainPanel, -1, '_   ', pos=(20, 555)) 
        
                
        
        wx.StaticText(self.MainPanel, -1, 'F#s:', pos=(175, 553))  
        self.rawImageFN1 = wx.SpinCtrl(self.MainPanel, -1, pos=(200, 550), size=(50,-1))
        self.rawImageFN1.SetRange(0, 100000)  
        self.rawImageFN1.SetValue(0)
        
        self.rawImageFN2 = wx.SpinCtrl(self.MainPanel, -1, pos=(260, 550), size=(50,-1))
        self.rawImageFN2.SetRange(0, 100000)  
        self.rawImageFN2.SetValue(49)
        
        
        
        wx.StaticText(self.MainPanel, -1, 'Center pixels (7,9,11,13): ', pos=(370, 553))  
        self.centerPixelN = wx.SpinCtrl(self.MainPanel, -1, pos=(510, 550), size=(50,-1))
        self.centerPixelN.SetRange(7, 13)  
        self.centerPixelN.SetValue(9)
        
        
        
        
        
        
        self.cb_includeFitGraphs_yn = wx.CheckBox(self.MainPanel, -1, 'Include Fit Graphs', (50, 580))
        self.cb_includeFitGraphs_yn.SetValue(True)

        self.cb_BG_subtraction_yn = wx.CheckBox(self.MainPanel, -1, 'Background Subtraction', (50, 600))
        self.cb_BG_subtraction_yn.SetValue(True)

        

        
        
        self.cb_TraceRotation_yn = wx.CheckBox(self.MainPanel, -1, 'Trace Rotation (7x7 only):', (200, 580))
        self.cb_TraceRotation_yn.SetValue(False)
        #self.radio_rotationDirection = wx.RadioBox(self.MainPanel, -1, choices=['CW', 'CCW'], label='rotation', pos=(200, 600), size=wx.Size(130, 40), style=wx.RA_SPECIFY_COLS)
        self.radio_rotationDirection = wx.RadioBox(self.MainPanel, -1, choices=['CW', 'CCW'], pos=(200, 595), style=wx.RA_SPECIFY_COLS)
        self.radio_rotationDirection.SetSelection(1)
        
        
        wx.StaticText(self.MainPanel, -1, "___________________ Reference two points in ImageJ ___________________", pos=(330, 594))
        wx.StaticText(self.MainPanel, -1, "0 deg: F#,     x1,         y1", pos=(340, 614))
        wx.StaticText(self.MainPanel, -1, "180 deg: F#,       x2,        y2", pos=(530, 614))
        wx.StaticText(self.MainPanel, -1, "1", pos=(380, 613))        
        self.RefPoint_x1 = wx.TextCtrl(self.MainPanel, -1, str(trace_rotation_x1), pos=(400, 630), size=(40, -1))
        self.RefPoint_y1 = wx.TextCtrl(self.MainPanel, -1, str(trace_rotation_y1), pos=(450, 630), size=(40, -1))
        self.RefPoint_F2N = wx.TextCtrl(self.MainPanel, -1, str(trace_rotation_F2N), pos=(570, 630), size=(30, -1))
        self.RefPoint_x2 = wx.TextCtrl(self.MainPanel, -1, str(trace_rotation_x2), pos=(610, 630), size=(40, -1))
        self.RefPoint_y2 = wx.TextCtrl(self.MainPanel, -1, str(trace_rotation_y2), pos=(660, 630), size=(40, -1))
        
        
        
        
                
        wx.StaticText(self.MainPanel, -1, "Data Files Prefix", pos=(50, 702))
        self.FilePrefix = wx.TextCtrl(self.MainPanel, -1, str(FilePrefixInput), pos=(150, 700), size=(180, -1))
        
        self.btn_Save_data_figures = wx.Button(self.MainPanel, pos=(50,740), label="Save Data")  
        self.btn_Save_data_figures.Hide()
        self.btn_Save_data_figures.Bind(wx.EVT_BUTTON, self.SaveDataAndFigures)
        self.text_Save_figures = wx.StaticText(self.MainPanel, -1, '_    ', pos=(20, 790)) 
        
        
        self.cb_SaveFigs_yn = wx.CheckBox(self.MainPanel, -1, 'Figures', (350,722))
        self.cb_SaveFigs_yn.SetValue(True)
        self.cb_SaveImageJ_xyData_yn = wx.CheckBox(self.MainPanel, -1, 'ImageJ x,y data & Fit: A, T, phase, y0, M', (350, 742))
        self.cb_SaveImageJ_xyData_yn.SetValue(False)
        self.cb_SaveIntensityData_yn = wx.CheckBox(self.MainPanel, -1, 'Intensity Trajectory Data for all', (350, 762))
        self.cb_SaveIntensityData_yn.SetValue(False)
        self.cb_SaveIntensityFitData_yn = wx.CheckBox(self.MainPanel, -1, 'Fit Trajectory Data for all', (350, 782))
        self.cb_SaveIntensityFitData_yn.SetValue(False)


        self.cb_SaveImageJ_xyDataGoodFit_yn = wx.CheckBox(self.MainPanel, -1, 'ImageJ x,y data & Fit: A, T, phase, y0, M for Good T fit', (350, 822))
        self.cb_SaveImageJ_xyDataGoodFit_yn.SetValue(False)
        self.cb_Complete_Data_yn = wx.CheckBox(self.MainPanel, -1, 'Complete Data File', (350, 842))
        self.cb_Complete_Data_yn.SetValue(True)











        self.btn_SelectData = wx.Button(self.MainPanel, pos=(50,810), label="Select Data")
        self.btn_SelectData.Bind(wx.EVT_BUTTON, self.SelectData)
        self.btn_SelectData.Hide()
        
        
        
        




     #----------------------------------------------------------------------


        

    def on_paint(self, event):
        dc = wx.PaintDC(event.GetEventObject())
        dc.Clear()
        dc.SetPen(wx.Pen("BLACK", 4))
        dc.DrawLine(0, 340, 1000, 340)  
        dc.DrawLine(0, 690, 1000, 690)  
        
        
        
        
    

    

    def ShowManual(self, event):
        self.manual_window = wx.Frame(None, title='Manual ', size=(800, 1000), pos=(800,0))
        self.manual_window.Show()
        
        wx.StaticText(self.manual_window, -1, manual_openFile, pos=(50, 30)) 
        wx.StaticText(self.manual_window, -1, manual_showM, pos=(50, 350)) 
        wx.StaticText(self.manual_window, -1, manual_saveData, pos=(50, 700)) 
        
        

        

    def CloseAllfigures(self, event):
        try:
            plt.close('all')
            
        except:pass
    
 
    def SwitchTo2chAnalysis(self, event):
        TwoChannelMain = TwoChannelMainTracking()
        TwoChannelMain.Show()
        self.Close()
        
         
    def SwitchToDataAnalysis(self, enven):
        DataAnalysisMain = DataAnalysisMainClass()
        DataAnalysisMain.Show()
        self.Close()
        

    def resetting(self, event):
        
        print '\n\nResetting'
        
        try:
            plt.close()
            plt.close()
            plt.close()
            plt.close()
            
        except: pass
                
        try:
            for n in range(len(self.fig_ExcitationM_Trajectory)):
                print 'fig_ExcitationM_Trajectory n = ', n
                plt.close()
        except: pass
        
        try:
            for n in range(len(self.fig_AllDataCenterImages)):
                print 'fig_AllDataCenterImages n = ', n
                plt.close()
        except: pass
        
        self.Close()                      
        frame2 = SingleChannelMainTracking()
        frame2.Show()
          
        

        
        
    def onOpenImageFile(self, event):
        #self.btn_FindAllPairs.Hide()
        #self.btn_ShowAllPairsMovie.Hide()
        #self.btn_ShowAllPairsMovieSingleN.Hide()
    
    
        self.btn_FindingFeature.Hide()
        self.btn_LoadFeatureXY.Hide()

        self.btn_Save_data_figures.Hide()
        self.btn_ShowRawImages.Hide()
        self.text_Save_figures.SetLabel('      ')
        self.text_LoadFeatureXY.SetLabel('_      ')
        self.text_LoadFeatureXY_filename.SetLabel('     ')
        self.btn_ShowAnisotropyM.Hide()
        self.text_NFeatures.SetLabel('_      ')
        self.cb_use_loaded_xy_yn.SetValue(False)
        self.cb_xy_offset_yn.SetValue(False)
        self.btn_SelectData.Hide()
        self.cb_apply_xyLimits_yn.SetValue(False)
        self.btn_FindingFeature_IntensityFilter.Hide()
        self.btn_SaveFindingFeatureData.Hide()
        self.btn_SaveFindingFeatureData_filtered.Hide()
        
        self.NofFeatures.SetLabel('_     ')
        
        self.MaxTFitGoodnessErr.SetValue('10')
        self.OffSet_y.SetValue(0)        
        self.OffSet_x.SetValue(0)
        

        #tp_ver = tp.__version__[:3]
        #if float(tp_ver[:3]) != 0.3:
        #    MessageBox(self, 'The trackpy version is not 0.3.#\nUpdate trackpy', 'Error')




        try:
            plt.close(self.fig_BG_ref)
        except: pass
    
        
        try:
            plt.close(self.fig_ShowAfewFrames)
        except: pass
        
        try: plt.close(self.fig_ExcitationM_Histogram)
        except: pass
    
        try:        
            for n in range(len(self.fig_AllDataCenterImages)):    
                plt.close(self.fig_AllDataCenterImages[n])
        except: pass

        try:        
            for n in range(len(self.fig_AllDataCenterImagesOnly)):    
                plt.close(self.fig_AllDataCenterImagesOnly[n])
        except: pass


        try:
            self.SelectData_frame.Destroy()
        except: pass

                
                
        """
        Create and show the Open FileDialog
        """
        plt.close()
        dlg2 = wx.FileDialog(
            self, message="Choose a file",
            defaultFile="",
            wildcard=wildcard,
            style=wx.OPEN | wx.MULTIPLE) # | wx.CHANGE_DIR  )
            
        if dlg2.ShowModal() == wx.ID_OK:
            paths = dlg2.GetPaths()
            print "You chose the following file(s):"
            for path in paths:
                print path
        print ''
        print 'dlg2 = ', dlg2
        print 'path = ', dlg2.GetPath()
        print 'path type', type(dlg2.GetPath())
        print 'filename = ', dlg2.GetFilename()
        #dlg.Destroy()
        
        self.Show()
        
        #wx.StaticText(self.MainPanel, -1, str(dlg2.GetFilename()), pos=(50, 150))
        
        filenameTemp = dlg2.GetFilename()
        self.btn2text.SetLabel(str(filenameTemp))
        
        self.filepath = dlg2.GetPath()
        self.filenameOnly = dlg2.GetFilename()
        self.filedirectoryOnly = dlg2.GetDirectory()
        print 'self.filedirectoryOnly ', self.filedirectoryOnly
        
        self.IntensityFilter_Applied_yn = 'n'
        self.NofFeatures_intensity_filtered.SetLabel('_      ')
        
        self.ShowAFewFrames()



    
    def ShowAFewFrames(self):
        try: 
            frames = pims.TiffStack(self.filepath)
        except:
            print 'Error: no file was chosen'
            return

        
        self.TotalNframes = len(frames)
            
        
        self.framesloaded = frames



                
        
        
        print 'TotalNframes = ', self.TotalNframes
        self.fig_ShowAfewFrames = plt.figure('frames')
        plt.subplots_adjust(top = 0.9, bottom = 0.05, left = 0.05, right = 0.95, wspace = 0.55, hspace = 0.2)
        plt.figtext(0.5, 0.95, str(os.path.basename(self.filepath))
            ,ha='center', color='black', weight='bold', size='small')
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(810,30,800, 800)
        
        if self.TotalNframes > 15:
            k = 16
        else:
            k = self.TotalNframes
            
        for n in range(k):
            print 'n = ', n
            plt.subplot(4,4,n+1)
            plt.title('frame # ' + str(n))
            
            
            plt.imshow(frames[n], origin='upper')  # for trackpy 3.0
  
            #plt.imshow(np.rot90(frames[n],3), origin='upper')
            
        
        
        self.N_xPixel = frames.frame_shape[0]
        self.N_yPixel = frames.frame_shape[1]
        

        
        self.btn_FindingFeature.Show()
        self.btn_LoadFeatureXY.Show()
        
        
        if self.TotalNframes < 45:
            self.Nframe2.SetValue(str(self.TotalNframes-1))
            self.NFramesForAveraging.SetValue(str(1))
            self.rawImageFN2.SetValue(self.TotalNframes-1)
            self.radio_FitOptions.SetSelection(1)
            
                        
            

        
        plt.ion()
        plt.show()
                
        
        
        
        
        
        

    def FindingFeatureFunc(self, event):
        
        
        plt.close('Finding Features')
        
        self.IntensityFilter_Applied_yn = 'n'
        
        try:
            plt.close(self.fig_FindingFeatures_filtered)
        except: pass
    
        
        featureSize = int(self.FeatureSize.GetValue())
        minIntensity = int(self.MinIntensity.GetValue())
        
        frames = self.framesloaded

        NframeTemp = int(self.Nframe.GetValue())        
        NFramesForAveraging = int(self.NFramesForAveraging.GetValue())
        
        
        if featureSize%2 == 0:
            MessageBox(self, 'Feature size must be an odd number', 'Error')
            return
            
        
        
        

        try:
            for n in range(len(self.fig_AllDataCenterImagesOnly)):
                plt.close(self.fig_AllDataCenterImagesOnly[n])
        except: pass
            
        
        
        if self.TotalNframes  < int(self.Nframe2.GetValue()) + int(self.NFramesForAveraging.GetValue()):
            MessageBox(self, 'Total # of frames is ' + str(self.TotalNframes) + '\nCheck the number of the frames', 'Error')
            return
        
        
        
        
        
        print 'Starting Frame # ', NframeTemp
        
        
        
        
       
        
        framesSummed = frames[0] * 0 # for trackpy 3.0
        #framesSummed = np.rot90(frames[0].T, 2) * 0
        
        
        
        
        for n in range(NFramesForAveraging):
            
            #framesSummed += np.rot90(frames[NframeTemp + n].T, 1) # for trackpy 3.0
            framesSummed += frames[NframeTemp + n] # for trackpy 3.0
            #framesSummed += np.rot90(frames[NframeTemp + n].T, 2)
            
        framesSummed /= NFramesForAveraging
        
        
        if self.Nframe2.GetValue() != 'na':
            NframeTemp2 = int(self.Nframe2.GetValue())
            framesSummed2 = frames[0] * 0 # for trackpy 3.0
            #framesSummed2 = np.rot90(frames[0].T, 1) * 0 # for trackpy 3.0
            #framesSummed2 = np.rot90(frames[0].T, 2) * 0
            
            for n in range(NFramesForAveraging):
                
                #framesSummed2 += np.rot90(frames[NframeTemp2 + n].T, 1) # for trackpy 3.0
                framesSummed2 += frames[NframeTemp2 + n] # for trackpy 3.0
                #framesSummed2 += np.rot90(frames[NframeTemp2 + n].T, 2)
            framesSummed2 /= NFramesForAveraging
            framesSummed = (framesSummed + framesSummed2)/2
            

        self.AveragedImage = framesSummed
        
        
        f = tp.locate(framesSummed, featureSize, minmass=minIntensity, invert=False)
        #NFeatures = f.values.size/8
        NFeatures = len(f.index) #trackpy 0.3
        print 'NFeatures = ', NFeatures
        self.NofFeatures.SetLabel(str(NFeatures))
        
        
        
        
        if NFeatures == 0:
            MessageBox(self, '0 features. Adjust the Min Intensity', 'Error')
            return
            
            
            
            
            
        
        #plt.close()
        self.fig_FindingFeatures = plt.figure('Finding Features')
        plt.figtext(0.5, 0.95, str(os.path.basename(self.filepath)) + '\nChosen F#s: ' + str(self.Nframe.GetValue()) + ', ' + str(self.Nframe2.GetValue()) +' ,  # of frames for averaging: ' + str(NFramesForAveraging)
        + '\n trackpy min intensity threshold: ' + str(minIntensity) + ' ,  Feature size: ' + str(featureSize)  + ' ,   Total # of features: '+ str(NFeatures)    ,ha='center', color='black', weight='bold', size=9)

 
        print 'f ', f
        

        
        
        intensityData = []
        for n in range(NFeatures):
            intensity_temp = f.values[n][2]
            intensityData.append(intensity_temp)
        
        
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(510,30,1400, 700)
        #tp.annotate(f, frames[NframeTemp])
        
        plt.subplot(121)
        tp.annotate(f, framesSummed)

        plt.subplot(122)
        plt.hist(intensityData)
        plt.xlabel('Crocker Grier Feature Intensity (trackpy 0.3.0)')

        
   
        
        self.N_xPixel = frames.frame_shape[0]
        self.N_yPixel = frames.frame_shape[1]
        print 'self.N_xPixel: ', self.N_xPixel
        print 'self.N_yPixel: ', self.N_yPixel        
        
        self.xpos = []
        self.ypos = []
        pixel_dn = 30
        for n in range(NFeatures):
            x_temp = f.values[n][0]
            y_temp = f.values[n][1]
            
            print 'x_temp ', x_temp
            print 'y_temp ', y_temp
            if x_temp >=pixel_dn and x_temp <= self.N_xPixel - pixel_dn and y_temp >=pixel_dn and y_temp <= self.N_yPixel - pixel_dn:
                self.xpos.append(x_temp)
                self.ypos.append(y_temp)

        
         
        #self.xposImageJFound = self.N_yPixel - np.array(self.ypos)
        #self.yposImageJFound = np.array(self.xpos)
         
         
         
        
        self.xposImageJFound_temp = np.array(self.xpos) + 0                      # the last - 1 is for correction between ImageJ and the Crocker Grier (x,y)

        #self.yposImageJFound_temp = self.N_xPixel - np.array(self.ypos) - 1    # the last subtraction is for correction between ImageJ and Crock's Grier (x,y)
        self.yposImageJFound_temp = np.array(self.ypos) + 0    #trackpy 0.3
        
        self.FindingFeature_xposImageJ = self.xposImageJFound_temp
        self.FindingFeature_yposImageJ = self.yposImageJFound_temp
        self.FindingFeature_IntensityData = intensityData
        
        
     
        
        
        
        
        print '\nxpos' , self.xpos
        print '\nypos' , self.ypos
        
                
        print '\nxposImageJ_temp', self.xposImageJFound_temp, '\nyposImageJ_temp', self.yposImageJFound_temp

        
        
        #self.FindingPairs(self.xposImageJ, self.yposImageJ)
        
        
        
        plt.ion()
        plt.show()
        
        self.btn_ShowAnisotropyM.Show()
        self.btn_ShowRawImages.Show()
        self.text_NFeatures.SetLabel('_     ')
        self.NofFeatures_intensity_filtered.SetLabel('_     ')
        self.cb_use_loaded_xy_yn.SetValue(False)
        self.cb_xy_offset_yn.SetValue(False)
        self.btn_SelectData.Hide()
        
        self.btn_FindingFeature_IntensityFilter.Show()
        self.btn_SaveFindingFeatureData.Show()
        self.btn_SaveFindingFeatureData_filtered.Hide()
        
        self.cb_MaxTFitGoodnessErr_yn.SetValue(True)
        self.cb_MinAveIntensity_yn.SetValue(True)
        
                
        #self.btn_FindAllPairs.Show()
     
     
     
     
     
    def FindingFeature_ApplyIntensityFilter(self, event):
        print '\nApplyIntensityFilter'
        
        self.IntensityFilter_Applied_yn = 'y'
        

        IntensityMethod = self.radio_IntensityMethod.GetSelection()

        if IntensityMethod == 3:
            MessageBox(self, 'Summed intensity method is not supported', 'Error')
            return
        else: pass
            
        
        if self.cb_BG_subtraction_yn.GetValue():
            BG_ref_correction_text = 'Background Subtraction: Yes, '
        else:
            BG_ref_correction_text = 'Background Subtraction: No, '
            
        if IntensityMethod == 0:
            BG_ref_correction_text += ' Max 5 Ave '
        elif IntensityMethod == 1:
            BG_ref_correction_text += ' Max 5 Med '
        elif IntensityMethod == 2:
            BG_ref_correction_text += ' Max 1 '
        
        
        
            
            
        
        centerLength = int(self.centerPixelN.GetValue())
        if centerLength == 7 or centerLength == 9 or centerLength == 11 or centerLength == 13:
            pass
        else:
            MessageBox(self, 'Choose 7, 9, 11 or 13 for the center pixel #s', 'Error')
            return
        

        BG_ref_correction_text += ',  Center pixels: ' + str(centerLength) + 'x' + str(centerLength)
        

        
        centerIntensity,backgroundIntensity = IntensityCalculationForSingleNumpyFrame(self.AveragedImage, self.xposImageJFound_temp, self.yposImageJFound_temp, IntensityMethod, centerLength)
        
        centerIntensityBS = np.array(centerIntensity) - np.array(backgroundIntensity)
        print '\n centerIntensity',     centerIntensity
        print '\n backgroundIntensity',     backgroundIntensity
        print '\n centerIntensityBS',     centerIntensityBS

        print '\n ',     
        
        
        
        IntensityFilterThreshold = int( self.IntensityFilterThreshold.GetValue() )
        
        
        if self.cb_BG_subtraction_yn.GetValue():
            intensity_temp = centerIntensityBS
            BS_yn_text = 'Background Subtraction: Yes'
        else:
            intensity_temp = centerIntensity
            BS_yn_text = 'Background Subtraction: No'
            
        
        xTemp = []
        yTemp = []
        intensityData = []
        for n in range(len(self.xposImageJFound_temp)):
            if intensity_temp[n] >= IntensityFilterThreshold:
                xTemp.append(self.xposImageJFound_temp[n])
                yTemp.append(self.yposImageJFound_temp[n])
                intensityData.append(intensity_temp[n])
                
        self.xposImageJFound_temp_filtered = xTemp
        self.yposImageJFound_temp_filtered = yTemp
        self.intensity_temp_filtered = intensityData
        
        
        
        
        NFeatures_filtered = len(self.xposImageJFound_temp_filtered)
        NFramesForAveraging = int(self.NFramesForAveraging.GetValue())
        
        
        self.NofFeatures_intensity_filtered.SetLabel(str(NFeatures_filtered))
                
            
        self.fig_FindingFeatures_filtered = plt.figure('Finding Features filtered')
        self.fig_FindingFeatures_filtered.clf()
        
        plt.figtext(0.5, 0.95, str(os.path.basename(self.filepath)) + '\nChosen F#s: ' + str(self.Nframe.GetValue()) + ', ' + str(self.Nframe2.GetValue()) + ',  # of frames for averaging: ' + str(NFramesForAveraging) + ' ,  '+  BG_ref_correction_text
        + '\n Intensity threshold: ' + str(IntensityFilterThreshold) + ' ,    Total # of features: '+ str(NFeatures_filtered) ,ha='center', color='black', weight='bold', size=9)


 
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(510,230,1400, 700)
        #tp.annotate(f, frames[NframeTemp])
        
        plt.subplot(121)
            
        plt.imshow(self.AveragedImage, origin='upper', cmap = plt.cm.Greys_r)  # for trackpy 3.0
        plt.plot(self.xposImageJFound_temp_filtered, self.yposImageJFound_temp_filtered, 'o', markersize=15, markeredgewidth=2, markeredgecolor= 'c', markerfacecolor='None' )
        
        plt.xlim(0, len(self.AveragedImage))
        plt.ylim(len(self.AveragedImage[0]), 0)
        
          

        plt.subplot(122)
        plt.title(BS_yn_text)
        plt.hist(intensityData)
        plt.xlabel('Feature Intensity by the Intensity Methods')

        
        
        plt.ion()
        plt.show()
        
        
        
        
        self.btn_SaveFindingFeatureData_filtered.Show()
        
        






        
        
        
        



    def LoadFeatureXY(self, event):
        print '\n Starting LoadFeatureXY'
        
        self.IntensityFilter_Applied_yn = 'n'
        
        
        try:
            for n in range(len(self.fig_AllDataCenterImagesOnly)):
                plt.close(self.fig_AllDataCenterImagesOnly[n])
        except: pass
    
    
        dlg2 = wx.FileDialog(
            self, message="Choose a file",
            defaultFile="",
            wildcard = wildcardDataFile,
            style=wx.OPEN | wx.MULTIPLE | wx.CHANGE_DIR
            )
        if dlg2.ShowModal() == wx.ID_OK:
            paths = dlg2.GetPaths()
            print "The following file(s) are chosen:"
            for path in paths:
                print path
        
        
        self.filenameOnlyLoadeXY = dlg2.GetFilename()
        
        
        fdata = open(paths[0], 'r')
        linestemp = fdata.readlines()
        fdata.close()
        
        
        
        


        if (linestemp[0].split())[0] == '#0':
                        
            self.Nfeatures = len(linestemp)/7

            self.Loaded_ImageJxpos = []
            self.Loaded_ImageJypos = []
            
            for n in range(self.Nfeatures):
                p1 = linestemp[7*n + 1].split()
                self.Loaded_ImageJxpos.append(int(float(p1[0])) ) 
                self.Loaded_ImageJypos.append(int(float(p1[1])) )
                


        elif (linestemp[0].split())[0] == '##':
                        
            self.Nfeatures = len(linestemp)/8

            self.Loaded_ImageJxpos = []
            self.Loaded_ImageJypos = []
            
            for n in range(self.Nfeatures):
                p1 = linestemp[8*n + 2].split()
                self.Loaded_ImageJxpos.append(int(float(p1[0])) ) 
                self.Loaded_ImageJypos.append(int(float(p1[1])) )
                


        elif (linestemp[0].split())[0] == '%%':
                        
            self.Nfeatures = len(linestemp)/10
            print 'len(linestemp) ', len(linestemp)
            print  'self.Nfeatures ', self.Nfeatures
            
            self.Loaded_ImageJxpos = []
            self.Loaded_ImageJypos = []
            
            for n in range(self.Nfeatures):
                p1 = linestemp[10*n + 2].split()
                self.Loaded_ImageJxpos.append(int(float(p1[0])) ) 
                self.Loaded_ImageJypos.append(int(float(p1[1])) )
                


        else:
            self.Loaded_ImageJxpos = []
            self.Loaded_ImageJypos = []
            
            for line in linestemp:
                p = line.split()
                self.Loaded_ImageJxpos.append(int(float(p[0])) ) 
                self.Loaded_ImageJypos.append(int(float(p[1])) )
               


        

        print 'int(float(self.OffSet_x.GetValue()))   ', int(float(self.OffSet_x.GetValue()))
        
        self.xposImageJFound = self.Loaded_ImageJxpos
        self.yposImageJFound = self.Loaded_ImageJypos
        
        
        #self.MinAveIntensity.SetValue('-9999')
        #self.MaxTFitGoodnessErr.SetValue('9999')
        self.cb_MaxTFitGoodnessErr_yn.SetValue(False)
        self.cb_MinAveIntensity_yn.SetValue(False)
        
        
        self.text_LoadFeatureXY.SetLabel(str(len(self.Loaded_ImageJxpos)))
        self.text_LoadFeatureXY_filename.SetLabel(self.filenameOnlyLoadeXY)
        
        self.btn_ShowAnisotropyM.Show()
        self.btn_ShowRawImages.Show()
        
        self.text_NFeatures.SetLabel('_     ')
        self.cb_use_loaded_xy_yn.SetValue(True)
        self.cb_xy_offset_yn.SetValue(True)
        self.cb_xy_offset_yn.Show()
        
        self.cb_apply_xyLimits_yn.SetValue(False)
    






    def ShowRawImages(self, event):
        print '\nstarting ShowRawImages'
        
        if self.IntensityFilter_Applied_yn == 'y':
            self.xposImageJFound_temp = self.xposImageJFound_temp_filtered
            self.yposImageJFound_temp = self.yposImageJFound_temp_filtered
        else: pass
            
        
        
        
        

        
        if self.cb_use_loaded_xy_yn.GetValue():
            if self.cb_xy_offset_yn.GetValue():
                #self.xposImageJFound = [x + int(float(self.OffSet_x.GetValue())) for x in self.Loaded_ImageJxpos]
                #self.yposImageJFound = [y - int(float(self.OffSet_y.GetValue())) for y in self.Loaded_ImageJypos]
                
                self.xposImageJFound = [x - int(float(self.OffSet_x.GetValue())) for x in self.Loaded_ImageJxpos] #trackpy 0.3
                self.yposImageJFound = [y + int(float(self.OffSet_y.GetValue())) for y in self.Loaded_ImageJypos] #trackpy 0.3
                
            else:
                self.xposImageJFound = self.Loaded_ImageJxpos
                self.yposImageJFound = self.Loaded_ImageJypos
           
        else:
            if self.cb_xy_offset_yn.GetValue():
                #self.xposImageJFound = [x + int(float(self.OffSet_x.GetValue())) for x in self.xposImageJFound_temp]
                #self.yposImageJFound = [y - int(float(self.OffSet_y.GetValue())) for y in self.yposImageJFound_temp]
                
                self.xposImageJFound = [x - int(float(self.OffSet_x.GetValue())) for x in self.xposImageJFound_temp] #trackpy 0.3
                self.yposImageJFound = [y + int(float(self.OffSet_y.GetValue())) for y in self.yposImageJFound_temp] #trackpy 0.3
                
            else:
                self.xposImageJFound = self.xposImageJFound_temp
                self.yposImageJFound = self.yposImageJFound_temp
                
            


        moreXtemp = []
        moreYtemp = []
        for n in range(len(self.MorePairs_x)):
            if int(self.MorePairs_x[n].GetValue()) != 0:
                #moreXtemp.append(int(self.MorePairs_x[n].GetValue()) + 0)
                #moreYtemp.append(self.N_xPixel - int(self.MorePairs_y[n].GetValue()) -1 )
            
                moreXtemp.append(int(self.MorePairs_x[n].GetValue()) - 0) #trackpy 0.3
                moreYtemp.append(int(self.MorePairs_y[n].GetValue()) + 0) #trackpy 0.3
                
            
        
        
        xposTemp = np.append(np.around(self.xposImageJFound, 0) , moreXtemp)
        yposTemp = np.append(np.around(self.yposImageJFound, 0) , moreYtemp)
        




        if self.cb_apply_xyLimits_yn.GetValue():
                
            x_limit_min = int(float(self.x_Limit_Min.GetValue()))
            x_limit_max = int(float(self.x_Limit_Max.GetValue()))
            y_limit_min = self.N_xPixel - int(float(self.y_Limit_Max.GetValue()))
            y_limit_max = self.N_xPixel - int(float(self.y_Limit_Min.GetValue()))
            
    
            print 'xposTemp ', xposTemp
            print 'x_limit_min ', x_limit_min
            print 'x_limit_max ', x_limit_max
            print 'y_limit_min ', y_limit_min
            print 'y_limit_max ', y_limit_max
    
            
    
            
            xpos = []
            ypos = []
            
           
            
            for n in range(len(xposTemp)):
                if x_limit_min <= xposTemp[n] <= x_limit_max and y_limit_min <= yposTemp[n] <= y_limit_max:
                    xpos.append(xposTemp[n])
                    ypos.append(yposTemp[n])
                    
                            
            xpos = np.array(xpos)                        
            ypos = np.array(ypos)
            
        else:
            xpos = xposTemp
            ypos = yposTemp


        
        print 'len(xpos) ', len(xpos)
        
        
        
        
        
        ################################################
        x_ImageJ, y_ImageJ = xpos, ypos
        
        
        #xc = y_ImageJ
        #yc = self.N_xPixel - x_ImageJ - 1
        xc = x_ImageJ - 0 # for trackpy 3.0
        yc = y_ImageJ # for trackpy 3.0
            
        ################################################
        
        
        
        filepath = self.filepath
        frameStart = int(self.Nframe.GetValue())
        frameEnd = frameStart 
        
        
        
        
        print 'frameStart ', frameStart
        print 'frameEnd ', frameEnd
        print 'data file = ', filepath
        
    
        
        frames = self.framesloaded
        
        TotalNframes = len(frames)
        print 'TotalNframes: ', TotalNframes
        
      
        

        NCenterSubPixel = 7  # number of pixels from the center
        
        
        
        #frameTemp0 = frame0[yc-NCenterSubPixel:yc+NCenterSubPixel + 1, xc-NCenterSubPixel:xc+NCenterSubPixel+1]
        x3 = np.linspace(0, 2 * NCenterSubPixel, 2* NCenterSubPixel + 1)
        y3 = np.linspace(0, 2 * NCenterSubPixel, 2*NCenterSubPixel+1)
        x3, y3 = np.meshgrid(x3, y3)
        
        





        centerImage = [[] for _ in xrange(len(xc))] # centerImage[Molecule#][Frame#]
    
        frameNTemp = int(self.rawImageFN.GetValue())
        print 'center image data for frame # ', frameNTemp
        for k in range(len(xpos)):            
            
            frameLargeTemp = frames[frameNTemp][yc[k]-NCenterSubPixel:yc[k]+NCenterSubPixel+1, xc[k]-NCenterSubPixel:xc[k]+NCenterSubPixel+1]
            #frameLargeTemp = np.transpose(np.rot90(frames[frameNTemp][xc[k]-NCenterSubPixel:xc[k]+NCenterSubPixel+1, yc[k]-NCenterSubPixel:yc[k]+NCenterSubPixel+1], 3)) # for trackpy 3.0
            
            
            #frameCenterTemp = frames[frameTemp][yc[k]-2:yc[k]+3, xc[k]-2:xc[k]+3]
            
            centerImage[k].append(frameLargeTemp) # it's not matching to the actual frame #, only first a few are saved



        self.centerImage = centerImage





        
        #print "centerImage \n", centerImage
        
        #print "centerImage[0] \n", centerImage[0]
        
        #print "centerImage[0][0] \n", centerImage[0][0]
        
        #print "centerImage[0][0][0]", centerImage[0][0][0]
        
        #print "centerImage[0][0][0][0]", centerImage[0][0][0][0]
        
        NFeatures = len(xpos)
  
       
        NCfig = int(  math.ceil((NFeatures/84.0)))
        self.fig_AllDataCenterImagesOnly = [[] for _ in xrange(NCfig)]
        
        
        
            



        x3 = np.linspace(0, 15, 16)
        y3 = np.linspace(0, 15, 16)
        x3, y3 = np.meshgrid(x3, y3)
        
        fsn = int(self.Nframe.GetValue()) # frame start number
        k = 0
        fn = 0
        for n in range(NCfig):
            #print 'n = ', n
            self.fig_AllDataCenterImagesOnly[n] = plt.figure('Raw_Data_Images'+ str(n), figsize = (18, 9))
            plt.clf()
            self.fig_AllDataCenterImagesOnly[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.36, hspace = 0.4)
            plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n'+ 'Frame#: '+ str(frameNTemp) + ',      ' + str(self.filenameOnly)
                ,ha='center', color='black', weight='normal', size='small')
            plt.figtext(0.05, 0.94, 'Total #: '+str(len(x_ImageJ)) , ha='left', color='black', weight='normal', size='small')
            
                


            for m in range(84):


                try:
                        
                    plt.subplot(6,14, m+1)
                    plt.title('M# ' + str(m+fn) + ',(' + str(int(x_ImageJ[m+fn])) + ',' + str(int(y_ImageJ[m+fn]))+')'  , fontsize=8, color = 'black')
                    plt.tick_params(labelsize=7)
                    #plt.imshow(centerImage[m+fn][0].reshape(15, 15), interpolation = 'None', cmap=plt.cm.jet, origin='bottom', extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    plt.imshow(centerImage[m+fn][0].reshape(15, 15), interpolation = 'None', cmap=plt.cm.jet, origin='upper', extent=(x3.min(), x3.max(), y3.min(), y3.max())) # trackpy 0.3
                    
                    
                except:
                    print '\nNo center image for ', 'M# = ', str(m + fn), '  F# ', str(fsn)


                
                k += 1
                if k%84 ==0:
                    fn += 84
                if k == NFeatures:
                    break



        plt.ion()
        plt.show()
        




    

    def ShowAnisotropyM(self, event):
                
        print '\n running ShowAnisotropyM'

        
        if self.IntensityFilter_Applied_yn == 'y':
            self.xposImageJFound_temp = self.xposImageJFound_temp_filtered
            self.yposImageJFound_temp = self.yposImageJFound_temp_filtered
        else: pass
            
        
        
        
        
        centerLength = int(self.centerPixelN.GetValue())
        if centerLength == 7 or centerLength == 9 or centerLength == 11 or centerLength == 13:
            pass
        else:
            MessageBox(self, 'Choose 7, 9, 11 or 13 for the center pixel #s', 'Error')
            return
        
        
        
        
        firstFrameN = int(self.rawImageFN1.GetValue())
        secondFrameN = int(self.rawImageFN2.GetValue())
        
        
        if secondFrameN <= self.TotalNframes-1:
            pass
        else:
            MessageBox(self, 'Change the 2nd image frame #', 'Error')
            return

        
        
        
        
        
        try: plt.close(self.fig_BG_ref)
        except: pass

        
        try: plt.close(self.fig_ExcitationM_Histogram)
        except: pass
    
        try:        
            for n in range(len(self.fig_AllDataCenterImages)):    
                plt.close(self.fig_AllDataCenterImages[n])
        except: pass

    
        
        self.MaxTfitRatio = 1 + float(self.MaxTFitGoodnessErr.GetValue())*0.01
        self.MinTfitRatio = 1 - float(self.MaxTFitGoodnessErr.GetValue())*0.01
        



        
        if self.cb_use_loaded_xy_yn.GetValue():
            if self.cb_xy_offset_yn.GetValue():
                #self.xposImageJFound = [x + int(float(self.OffSet_x.GetValue())) for x in self.Loaded_ImageJxpos]
                #self.yposImageJFound = [y - int(float(self.OffSet_y.GetValue())) for y in self.Loaded_ImageJypos]
                
                self.xposImageJFound = [x - int(float(self.OffSet_x.GetValue())) for x in self.Loaded_ImageJxpos] #trackpy 0.3
                self.yposImageJFound = [y + int(float(self.OffSet_y.GetValue())) for y in self.Loaded_ImageJypos] #trackpy 0.3
                
            else:
                self.xposImageJFound = self.Loaded_ImageJxpos
                self.yposImageJFound = self.Loaded_ImageJypos
           
        else:
            if self.cb_xy_offset_yn.GetValue():
                #self.xposImageJFound = [x + int(float(self.OffSet_x.GetValue())) for x in self.xposImageJFound_temp]
                #self.yposImageJFound = [y - int(float(self.OffSet_y.GetValue())) for y in self.yposImageJFound_temp]
                
                self.xposImageJFound = [x - int(float(self.OffSet_x.GetValue())) for x in self.xposImageJFound_temp] #trackpy 0.3
                self.yposImageJFound = [y + int(float(self.OffSet_y.GetValue())) for y in self.yposImageJFound_temp] #trackpy 0.3
                
            else:
                self.xposImageJFound = self.xposImageJFound_temp
                self.yposImageJFound = self.yposImageJFound_temp
                
            


                

        moreXtemp = []
        moreYtemp = []
        for n in range(len(self.MorePairs_x)):
            if int(self.MorePairs_x[n].GetValue()) != 0:
                #moreXtemp.append(int(self.MorePairs_x[n].GetValue()))
                #moreYtemp.append(self.N_xPixel - int(self.MorePairs_y[n].GetValue()) - 1)
            
                moreXtemp.append(int(self.MorePairs_x[n].GetValue()) - 0) #trackpy 0.3
                moreYtemp.append(int(self.MorePairs_y[n].GetValue()) + 0) #trackpy 0.3
                
               
        
        

        BG_ref_Xtemp = []
        BG_ref_Ytemp = []
        for n in range(len(self.BG_ref_x)):
            if int(self.BG_ref_x[n].GetValue()) != 0:
                BG_ref_Xtemp.append(int(self.BG_ref_x[n].GetValue()))
                BG_ref_Ytemp.append(self.N_xPixel - int(self.BG_ref_y[n].GetValue()))
            
        
        
        
        
        xposTemp = np.append(np.around(self.xposImageJFound, 0) , moreXtemp)
        yposTemp = np.append(np.around(self.yposImageJFound, 0) , moreYtemp)
        
        
        
        
        #######################
        # limiting x,y in the finding fieature window
        if self.cb_apply_xyLimits_yn.GetValue():
                
            x_limit_min = int(float(self.x_Limit_Min.GetValue()))
            x_limit_max = int(float(self.x_Limit_Max.GetValue()))
            y_limit_min = self.N_xPixel - int(float(self.y_Limit_Max.GetValue()))
            y_limit_max = self.N_xPixel - int(float(self.y_Limit_Min.GetValue()))
            
    
            print 'xposTemp ', xposTemp
            print 'x_limit_min ', x_limit_min
            print 'x_limit_max ', x_limit_max
            print 'y_limit_min ', y_limit_min
            print 'y_limit_max ', y_limit_max
    
            xpos = []
            ypos = []
            
            for n in range(len(xposTemp)):
                if x_limit_min <= xposTemp[n] <= x_limit_max and y_limit_min <= yposTemp[n] <= y_limit_max:
                    xpos.append(xposTemp[n])
                    ypos.append(yposTemp[n])

                
            Allxpos = np.array(xpos)                        
            Allypos = np.array(ypos)
                            
        else:
            Allxpos = xposTemp
            Allypos = yposTemp
        # limiting x,y in the finding fieature window
        #######################













        
        N_BG_ref = len(BG_ref_Xtemp)
        print 'N_BG_ref ', N_BG_ref

        
        Allxpos = np.append(Allxpos , BG_ref_Xtemp)
        Allypos = np.append(Allypos , BG_ref_Ytemp)
        
        
        
        print 'Allxpos ', Allxpos
        print 'Allypos ', Allypos
        
        

        radius = np.sqrt( (float(self.RefPoint_x1.GetValue()) - float(self.RefPoint_x2.GetValue()))**2 + (float(self.RefPoint_y1.GetValue()) - float(self.RefPoint_y2.GetValue()))**2 )/2.0
        refCenter_x0 = (float(self.RefPoint_x1.GetValue()) + float(self.RefPoint_x2.GetValue()))/2.0
        refCenter_y0 = (float(self.RefPoint_y1.GetValue()) + float(self.RefPoint_y2.GetValue()))/2.0
        
        theta0 = np.arccos((float(self.RefPoint_x1.GetValue())-refCenter_x0)/radius)
        theta0_by_y = np.arccos((float(self.RefPoint_y1.GetValue())-refCenter_y0)/radius)

        ref_dx =  refCenter_x0 - float(self.RefPoint_x1.GetValue()) # so frame[0] imageJ_x + ref_dx = rotation center x
        ref_dy = refCenter_y0 - float(self.RefPoint_y1.GetValue())
        dtheta = (3.14159)/( float(self.RefPoint_F2N.GetValue()) - 1 )







        dx1 = int(radius * np.cos(theta0 - dtheta * 0))
        dy1 = int(radius * np.sin(theta0 - dtheta * 0))
        
        secondFrameN = int(self.rawImageFN2.GetValue())
        dx2 = int(radius * np.cos(theta0 - dtheta * 90))
        dy2 = int(radius * np.sin(theta0 - dtheta * 90))
        
        #print 'dx1, dy1', dx1, dy1
        #print 'dx2, dy2', dx2, dy2
        
        #print 'refCenter_x0 + dx1 ', refCenter_x0 + dx1
        #print 'refCenter_y0 + dy1 ', refCenter_y0 + dy1

        distance_ref1 = np.sqrt( (float(self.RefPoint_x1.GetValue()) - (refCenter_x0 + dx1) )**2 + (float(self.RefPoint_y1.GetValue()) - (refCenter_y0 + dy1) )**2 )
        #print 'distance_ref1 ', distance_ref1
        
        if distance_ref1 > 2: # to avoid a degerate wrong theta0
            theta0 = -theta0
            dx1 = int(radius * np.cos(theta0 - dtheta * 0))
            dy1 = int(radius * np.sin(theta0 - dtheta * 0))
        
            distance_ref1 = np.sqrt( (float(self.RefPoint_x1.GetValue()) - (refCenter_x0 + dx1) )**2 + (float(self.RefPoint_y1.GetValue()) - (refCenter_y0 + dy1) )**2 )
            #print 'distance_ref1 2', distance_ref1
            


        '''
        print ' \n\n'
        print 'theta0 ', theta0
        print 'theta0_by_y ', theta0_by_y


        print 'ref_dx ', ref_dx
        print 'ref_dy ', ref_dy


        print 'self.RefPoint_x1.GetValue() ', self.RefPoint_x1.GetValue()
        print 'self.RefPoint_y1.GetValue() ', self.RefPoint_y1.GetValue()
        '''

        #exit()
        
        

        IntensityMethod = self.radio_IntensityMethod.GetSelection()

        
        
        if self.cb_TraceRotation_yn.GetValue() == True:            
            centerIntensityAveTrajectory, centerIntensityTrajectory, backgroundIntensityAveTrajectory, backgroundIntensityTrajectory = IntensityTrajectoryPlotsFuncRotForM(self.framesloaded, self.filepath, Allxpos, Allypos, radius, theta0, ref_dx, ref_dy, dtheta, IntensityMethod, centerLength)
        else:
            centerIntensityAveTrajectory, centerIntensityTrajectory, backgroundIntensityAveTrajectory, backgroundIntensityTrajectory = IntensityTrajectoryPlotsFuncForM(self.framesloaded, self.filepath, Allxpos, Allypos, IntensityMethod, centerLength)



        #print 'centerIntensityMaxTrajectory ', centerIntensityMaxTrajectory
        

        IntensityTrajectory_BG_ref = np.array(centerIntensityAveTrajectory[-N_BG_ref:])
        IntensityTrajectory_BG_ref_BS = np.array(centerIntensityAveTrajectory[-N_BG_ref:]) - np.array(backgroundIntensityAveTrajectory[-N_BG_ref:])

        
            
        #IntensityTrajectory_BG_ref_BS_Ave = np.average(IntensityTrajectory_BG_ref_BS, axis=0)
        
        

        
        self.fig_BG_ref = plt.figure('BG_ref')
        plt.title('Reference background data')        
        for n in range(len(IntensityTrajectory_BG_ref_BS)):
            if len(IntensityTrajectory_BG_ref[n]) ==1:
                plt.plot([IntensityTrajectory_BG_ref[n][0], IntensityTrajectory_BG_ref[n][0]], label='Ave')
                plt.plot([IntensityTrajectory_BG_ref_BS[n], IntensityTrajectory_BG_ref_BS[n][0]], label='Ave_BS')
                    
            else:
                plt.plot(IntensityTrajectory_BG_ref[n], label='Ave')
                plt.plot(IntensityTrajectory_BG_ref_BS[n], label='Ave_BS')
                
        plt.legend()
            
        
        
        
        
        IntensityTrajectory = np.array(centerIntensityTrajectory[:-N_BG_ref])
        self.backgroundIntensityTrajectory = np.array(backgroundIntensityTrajectory[:-N_BG_ref])
        
        
        
        
        

        self.backgroundIntensityTrajectory_Fit = []
        self.backgroundIntensityTrajectory_Fit_Min = []
        self.backgroundIntensityTrajectory_Fit_Max = []
        self.backgroundIntensityTrajectory_Fit_Ave = []
        T_in_frame = float(self.fitInitialGuess_T.GetValue())
        IntensityTrajectory_fit_x = np.arange(len(IntensityTrajectory[0]))
        for n in range(len(self.backgroundIntensityTrajectory)):
            A_guess = ( np.max(self.backgroundIntensityTrajectory[n]) - np.min(self.backgroundIntensityTrajectory[n]) )/2
            y0_guess = ( np.max(self.backgroundIntensityTrajectory[n]) + np.min(self.backgroundIntensityTrajectory[n]) )/2
            initial_guess = [A_guess, T_in_frame, 0, y0_guess]
            try:
                popt, pcov = curve_fit(Sine_fit_func, IntensityTrajectory_fit_x, self.backgroundIntensityTrajectory[n], p0=initial_guess )
                fit_temp = Sine_fit_func(IntensityTrajectory_fit_x, popt[0],popt[1],popt[2],popt[3])
            except:
                print 'BG fit failed for n: ', n
                fit_temp = Sine_fit_func(IntensityTrajectory_fit_x, 0.0, initial_guess[1], initial_guess[2], initial_guess[3] )
            self.backgroundIntensityTrajectory_Fit.append(fit_temp)
            self.backgroundIntensityTrajectory_Fit_Min.append(np.array([np.min(fit_temp)]))
            self.backgroundIntensityTrajectory_Fit_Max.append(np.array([np.max(fit_temp)]))
            self.backgroundIntensityTrajectory_Fit_Ave.append(np.array([np.mean(fit_temp)]))
        
        print '\n IntensityTrajectory\n', IntensityTrajectory
        print '\n np.array(self.backgroundIntensityTrajectory_Fit_Min) \n ', np.array(self.backgroundIntensityTrajectory_Fit_Min)

        
        
        if self.cb_BG_subtraction_yn.GetValue():
            #IntensityTrajectoryBS = np.array(centerIntensityMaxTrajectory[:-N_BG_ref]) - np.array(backgroundIntensityAveTrajectory[:-N_BG_ref]) - np.array(IntensityTrajectory_BG_ref_BS_Ave)
            BG_ref_correction_text = 'Background Subtraction: Yes, '
            
            if self.radio_BSMethod.GetSelection() == 0:
                IntensityTrajectoryBS = IntensityTrajectory - self.backgroundIntensityTrajectory
                BG_ref_correction_text += 'Frame by Frame, '
                self.BGS_type = 'FrameByFrame'
            elif self.radio_BSMethod.GetSelection() == 1:
                IntensityTrajectoryBS = IntensityTrajectory - np.array(self.backgroundIntensityTrajectory_Fit_Ave)
                BG_ref_correction_text += 'Fit Ave, '
                self.BGS_type = 'FitAve'
            
            
                                
        else:
            IntensityTrajectoryBS = IntensityTrajectory
            #BG_ref_correction_text = 'Background Substraction 9x9 edge: Total Ave'
            BG_ref_correction_text = 'Background Subtraction: No,'
            self.BGS_type = 'noBGS'
            
   
          
        if IntensityMethod == 0:
            BG_ref_correction_text += ' Max 5 Ave '
        elif IntensityMethod == 1:
            BG_ref_correction_text += ' Max 5 Med '
        elif IntensityMethod == 2:
            BG_ref_correction_text += ' Max 1 '
        elif IntensityMethod == 3:
            BG_ref_correction_text += ' Sum '
            
            
            
            
            
        
        
        Pol_Cor_factor = float(self.Excitation_Pol_correction_factor.GetValue())
        print '\nPol_Cor_factor = ', Pol_Cor_factor
        
        
        if self.cb_Excitation_Pol_correction_yn.GetValue():
            Excitation_Pol_correction_text = 'Excitation Polarization Correction: YES by ' + str(Pol_Cor_factor)
        else:
            Excitation_Pol_correction_text = 'Excitation Polarization Correction: NO'
            
            
            
        

        if self.cb_MaxTFitGoodnessErr_yn.GetValue():
            MaxT_Err_text_temp = str(self.MaxTFitGoodnessErr.GetValue())
        else:
            MaxT_Err_text_temp = 'na'

        if self.cb_MinAveIntensity_yn.GetValue():
            Min_Ave_I_text_temp = str(self.MinAveIntensity.GetValue())
        else:
            Min_Ave_I_text_temp = 'na'
                    
        
        
        
        self.AnalysisConditions = ( 'Chosen f # ' + str(self.Nframe.GetValue()) + '   Ave # f: ' + str(self.NFramesForAveraging.GetValue()) + 
        '   f size: '+ str(self.FeatureSize.GetValue()) +  '   Min I: ' + str(self.MinIntensity.GetValue()) 
        + ',   Total # Frames:' + str(self.TotalNframes) + ',   Threshold Ave I: ' + Min_Ave_I_text_temp
        + ',   MaxT%Err:' +  MaxT_Err_text_temp  )
        
        
        
        
        xposChosen = []
        self.xposChosen_GoodFits = []
        yposChosen = []
        self.yposChosen_GoodFits = []
        
        self.backgroundIntensityTrajectoryChosen = []
        self.backgroundIntensityTrajectory_FitChosen = []
        
        self.IntensityTrajectoryChosen = []
        self.IntensityTrajectoryChosenBS = []
        self.IntensityTrajectoryChosenBS_fit = []
        self.IntensityTrajectoryChosenBS_fit_M = []
        self.IntensityTrajectoryChosenBS_fit_M_GoodFits = []
        self.IntensityTrajectoryChosenBS_fit_RsquareFitGoodness = []
        self.IntensityTrajectoryChosenBS_fit_A = []
        self.IntensityTrajectoryChosenBS_fit_A_GoodFits = []
        self.IntensityTrajectoryChosenBS_fit_T = []
        self.IntensityTrajectoryChosenBS_fit_T_GoodFits = []
        self.IntensityTrajectoryChosenBS_fit_phi = []
        self.IntensityTrajectoryChosenBS_fit_phi_GoodFits = []
        self.IntensityTrajectoryChosenBS_fit_y0 = []
        self.IntensityTrajectoryChosenBS_fit_y0_GoodFits = []
        self.IntensityTrajectoryChosenBS_fit_MaxFrameN = []
        self.IntensityTrajectoryChosenBS_fit_MaxFrameN_GoodFit = []
        
        IntensityTrajectoryChosenBS_fit_x = np.arange(len(IntensityTrajectory[0]))
        
        
        T_in_frame = float(self.fitInitialGuess_T.GetValue())
        Min_AverageIntensity = float(self.MinAveIntensity.GetValue())
        
        
        print 'IntensityTrajectoryBS ', IntensityTrajectoryBS
        
        for n in range(len(IntensityTrajectoryBS)):
            AveIntensityTemp = np.mean(IntensityTrajectoryBS[n][:])
            if AveIntensityTemp >= Min_AverageIntensity or self.cb_MinAveIntensity_yn.GetValue() == False:
                                
                                
                                
                

                if self.radio_FitOptions.GetSelection() ==0:
                    
                    A_guess = ( np.max(IntensityTrajectoryBS[n]) - np.min(IntensityTrajectoryBS[n]) )/2
                    y0_guess = ( np.max(IntensityTrajectoryBS[n]) + np.min(IntensityTrajectoryBS[n]) )/2
                    
                    initial_guess = [A_guess, T_in_frame, 0, y0_guess]
                    
                    try:
                        popt, pcov = curve_fit(Sine_fit_func, IntensityTrajectoryChosenBS_fit_x, IntensityTrajectoryBS[n], p0=initial_guess )
                        fit_temp = Sine_fit_func(IntensityTrajectoryChosenBS_fit_x, popt[0],popt[1],popt[2],popt[3])
                        
                        A_temp = popt[0]
                        if A_temp < 0:
                            A_temp *= -1
                            
                        T_temp = popt[1]
                        phi_temp = (popt[2]/(3.141592*2))*180
                        y0_temp = popt[3]
                        MaxFrameN_temp = np.argmax(fit_temp)
                        
                        if popt[0] < 0:
                            phi_temp += 90
                        
                        while phi_temp >= 180:
                            phi_temp -= 180
                            
                        while phi_temp <= -180:
                            phi_temp += 180
                            
                        if phi_temp >= 90:
                            phi_temp = -(180 - phi_temp)
                            
                        if phi_temp <= -90:
                            phi_temp = 180 + phi_temp
                            
                        
                        print 'A, T, phi_temp, popt[2] = ', popt[0], T_temp, phi_temp, popt[2]
                        
                        #AmpErrTemp = np.sqrt(pcov[0,0])
                        #AmpErrPer = (AmpErrTemp/popt[0])*100.0
                        
                        RsquareFitGoodness = R_SqureFitGoodness(IntensityTrajectoryBS[n], fit_temp)
                        
                        
                    except:
                        print 'fit failed for n: ', n
                    
                        fit_temp = Sine_fit_func(IntensityTrajectoryChosenBS_fit_x, 0.0, initial_guess[1], initial_guess[2], initial_guess[3] )
                        
                        #AmpErrPer = 9999
                        A_temp = -999
                        T_temp = -999
                        phi_temp = -999
                        y0_temp = -999
                        
                        RsquareFitGoodness = -9999
                        MaxFrameN_temp = np.argmax(fit_temp)
                        
                    
                    

                elif self.radio_FitOptions.GetSelection() == 1:
                    print 'median y0 fit'
                    
                    A_temp = 0
                    T_temp = 0
                    phi_temp = 0
                    y0_temp = np.median(IntensityTrajectoryBS[n])
                    
                    fit_temp = IntensityTrajectoryBS[n]*0 + y0_temp
                    
                    
                    
                    RsquareFitGoodness = 0
                    MaxFrameN_temp = np.argmax(fit_temp)
                    
                
                
                
                    
                    
                    
                    
                    
                    
                    
                if self.cb_Excitation_Pol_correction_yn.GetValue():
                    Mtemp = (np.amax(fit_temp) - np.amin(fit_temp)       )/( np.amax(fit_temp) + np.amin(fit_temp) - Pol_Cor_factor*np.amax(fit_temp) )
                    Excitation_Pol_correction_text = 'Excitation Polarization Correction: YES by ' + str(Pol_Cor_factor)
                else:
                    Mtemp = (np.amax(fit_temp) - np.amin(fit_temp)       )/( np.amax(fit_temp) + np.amin(fit_temp))
                    Excitation_Pol_correction_text = 'Excitation Polarization Correction: NO'
                
                
                
                            
                if ( ( T_temp <= (float(self.fitInitialGuess_T.GetValue()))*self.MaxTfitRatio ) and ( T_temp >= (float(self.fitInitialGuess_T.GetValue()))*self.MinTfitRatio ) )\
                or ( self.cb_MaxTFitGoodnessErr_yn.GetValue() == False or self.radio_FitOptions.GetSelection() == 1):
                
                    xposChosen.append(Allxpos[n])
                    yposChosen.append(Allypos[n])
                    self.IntensityTrajectoryChosen.append(IntensityTrajectory[n])
                    self.IntensityTrajectoryChosenBS.append(IntensityTrajectoryBS[n])
                    self.backgroundIntensityTrajectoryChosen.append(self.backgroundIntensityTrajectory[n])
                    self.backgroundIntensityTrajectory_FitChosen.append(self.backgroundIntensityTrajectory_Fit[n])
    
                    
                    self.IntensityTrajectoryChosenBS_fit.append(fit_temp)
                    self.IntensityTrajectoryChosenBS_fit_M.append(Mtemp)
                    self.IntensityTrajectoryChosenBS_fit_RsquareFitGoodness.append(RsquareFitGoodness)
                    self.IntensityTrajectoryChosenBS_fit_A.append(A_temp)
                    self.IntensityTrajectoryChosenBS_fit_T.append(T_temp)
                    self.IntensityTrajectoryChosenBS_fit_phi.append(phi_temp)
                    self.IntensityTrajectoryChosenBS_fit_y0.append(y0_temp)
                    self.IntensityTrajectoryChosenBS_fit_MaxFrameN.append(MaxFrameN_temp)

                    
                    self.xposChosen_GoodFits.append(Allxpos[n])
                    self.yposChosen_GoodFits.append(Allypos[n])
                    self.IntensityTrajectoryChosenBS_fit_M_GoodFits.append(Mtemp)
                    self.IntensityTrajectoryChosenBS_fit_A_GoodFits.append(A_temp)
                    self.IntensityTrajectoryChosenBS_fit_T_GoodFits.append(T_temp)
                    self.IntensityTrajectoryChosenBS_fit_phi_GoodFits.append(phi_temp)
                    self.IntensityTrajectoryChosenBS_fit_y0_GoodFits.append(y0_temp)
                    self.IntensityTrajectoryChosenBS_fit_MaxFrameN_GoodFit.append(MaxFrameN_temp)
                    
                
         
        print 'IntensityTrajectoryChosenBS_fit_M ', self.IntensityTrajectoryChosenBS_fit_M
        print 'IntensityTrajectoryChosenBS_fit_M_GoodFits ', self.IntensityTrajectoryChosenBS_fit_M_GoodFits

        print 'Excitation_Pol_correction_text ', Excitation_Pol_correction_text
        
        Corrections = BG_ref_correction_text + ',       ' + Excitation_Pol_correction_text    
        Corrections2 = BG_ref_correction_text + '\n' + Excitation_Pol_correction_text + '\nF#s: ' + str(self.rawImageFN1.GetValue())+ ', ' + str(self.rawImageFN2.GetValue())
        self.Corrections = Corrections        
        self.Corrections2 = Corrections2

        
        NFeatures = len(self.IntensityTrajectoryChosen)
        #Nfig = int(  math.ceil((NFeatures/30.0)))
        
        TotalFeatureN_text = 'Total # ' + str(NFeatures)  + ' ,   Center: ' + str(centerLength) + 'x' +  str(centerLength) + ' ,   Background: ' + str(centerLength+2) + 'x' + str(centerLength+2)
        
        #self.fig_IntensityTrajectory = [[]] * Nfig
        #self.fig_ExcitationM_Trajectory = [[] for _ in xrange(Nfig)]
        









                    

        ############################################################################
        ##  center images  ##
        print '\n For RawImages'
        
        xpos = np.array(xposChosen)
        ypos = np.array(yposChosen)       
        
        self.xposChosen = xpos
        self.yposChosen = ypos
      
        ################################################
        if self.cb_TraceRotation_yn.GetValue() == False:
            x_ImageJ, y_ImageJ = xpos, ypos
            
        elif self.cb_TraceRotation_yn.GetValue() == True:
            x_ImageJ, y_ImageJ = xpos + ref_dx, ypos + ref_dy #making x,y as rotation centers
        
        #xc = y_ImageJ
        #yc = self.N_xPixel - x_ImageJ - 1
        
        xc = x_ImageJ - 0 # trackpy 0.3
        yc = y_ImageJ # trackpy 0.3
        ################################################
        
        
        
        filepath = self.filepath
        frameStart = int(self.Nframe.GetValue())
        frameEnd = frameStart 
        
     
        print 'frameStart ', frameStart
        print 'frameEnd ', frameEnd
        print 'data file = ', filepath
        
    
        
        frames = self.framesloaded
        
        TotalNframes = len(frames)
        print 'TotalNframes: ', TotalNframes

        
        NCenterSubPixel = 7  # number of pixels from the center
        
        
        
        #frameTemp0 = frame0[yc-NCenterSubPixel:yc+NCenterSubPixel + 1, xc-NCenterSubPixel:xc+NCenterSubPixel+1]
        x3 = np.linspace(0, 2 * NCenterSubPixel, 2* NCenterSubPixel + 1)
        y3 = np.linspace(0, 2 * NCenterSubPixel, 2* NCenterSubPixel + 1)
        x3, y3 = np.meshgrid(x3, y3)
        

        dx1 = int(radius * np.cos(theta0 - dtheta * firstFrameN))
        dy1 = int(radius * np.sin(theta0 - dtheta * firstFrameN))
        

        dx2 = int(radius * np.cos(theta0 - dtheta * secondFrameN))
        dy2 = int(radius * np.sin(theta0 - dtheta * secondFrameN))
        
        
        dx12 = int(radius * np.cos(theta0_by_y - dtheta * firstFrameN))
        dy12 = int(radius * np.sin(theta0_by_y - dtheta * firstFrameN))
        
        dx22 = int(radius * np.cos(theta0_by_y - dtheta * secondFrameN))
        dy22 = int(radius * np.sin(theta0_by_y - dtheta * secondFrameN))
        
        
        
        print 'dx1 ', dx1
        print 'dy1', dy1
        print 'dx2 ', dx2
        print 'dy2', dy2
        
        print 'dx12', dx12
        print 'dy12', dy12
        print 'dx22', dx22
        print 'dy22', dy22       
        
        print 'dtheta ', dtheta
        print 'radius ', radius
        print 'theta0 ', theta0

        self.centerImage = [[] for _ in xrange(len(xc))] # centerImage[Molecule#][Frame#]
        
        if self.cb_TraceRotation_yn.GetValue() == False:
            for frameTemp in [firstFrameN,secondFrameN]: # from the 0th frame to the end frame.
                print 'Not rotated center image data for frame # ', frameTemp
                for k in range(len(xpos)):            
                    frameLargeTemp = frames[frameTemp][yc[k]-NCenterSubPixel:yc[k]+NCenterSubPixel+1, xc[k]-NCenterSubPixel:xc[k]+NCenterSubPixel+1]
                    #frameCenterTemp = frames[frameTemp][yc[k]-2:yc[k]+3, xc[k]-2:xc[k]+3]
                    self.centerImage[k].append(frameLargeTemp) # it's not matching to the actual frame #, only first a few are saved
            
            
        elif self.cb_TraceRotation_yn.GetValue() == True:

            dx, dy = dx1, dy1
            print 'dx, dy 1', dx, dy

            for frameTemp in [firstFrameN, secondFrameN]: # from the 0th frame to the end frame.
                print 'Rotated center image data for frame # ', frameTemp

                for k in range(len(xpos)):
                    frameLargeTemp = frames[frameTemp][yc[k]-dx - NCenterSubPixel : yc[k]-dx + NCenterSubPixel + 1, xc[k]+dy - NCenterSubPixel : xc[k]+dy + NCenterSubPixel + 1]
                    #frameCenterTemp = frames[frameTemp][yc[k]-2:yc[k]+3, xc[k]-2:xc[k]+3]
                    self.centerImage[k].append(frameLargeTemp) # it's not matching to the actual frame #, only first a few are saved
                
                dx, dy = dx2, dy2
                print 'dx, dy 2', dx, dy
                    

        ##  center images  ##
        ############################################################################











        
        NFeatures = len(xpos)
        print '\n\nNFeatures: ', NFeatures
  
        NFeaturePerFig = 30.0
       
        NCfig = int(  math.ceil((NFeatures/float(NFeaturePerFig))))
        self.fig_AllDataCenterImages = [[] for _ in xrange(NCfig)]

        x3 = np.linspace(0, 15, 16)
        y3 = np.linspace(0, 15, 16)
        x3, y3 = np.meshgrid(x3, y3)
        
        #fsn = int(self.Nframe.GetValue()) # frame start number
        k = 0
        fn = 0
        for n in range(NCfig):
            #print 'n = ', n
            self.fig_AllDataCenterImages[n] = plt.figure('M_Data_Images'+ str(n), figsize = (18, 9))
            self.fig_AllDataCenterImages[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.36, hspace = 0.4)
            self.fig_AllDataCenterImages[n].text(0.02, 0.93, TotalFeatureN_text +'\n' +Corrections2, ha="left", va="bottom", size="small",color="black", weight='normal')
            
            plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly+ '\n' + self.AnalysisConditions ), ha='center', color='black', weight='normal', size='small')
                

            xlim_max = len(self.IntensityTrajectoryChosen[0])
            for m in range(int(NFeaturePerFig)):

                #try:
                #if self.IntensityTrajectoryChosenBS_fit_RsquareFitGoodness[m + fn] >= float(self.MinRsquareFitGoodness.GetValue()):
                if ( self.IntensityTrajectoryChosenBS_fit_T[m + fn] <= (float(self.RefPoint_F2N.GetValue()) - 1)*self.MaxTfitRatio ) and ( self.IntensityTrajectoryChosenBS_fit_T[m + fn] >= (float(self.RefPoint_F2N.GetValue()) - 1)*self.MinTfitRatio ):
                    titlecolorTemp = 'black'
                    weight='normal'
                else:
                    titlecolorTemp = 'red'
                    weight='normal'
                
                
                imageTemp = np.append(np.ravel(self.centerImage[m+fn][0]), np.ravel(self.centerImage[m+fn][1]))
                
                vminTemp = int(np.mean( np.sort(imageTemp, axis=None)[:5]) )
                vmaxTemp = int(np.mean( np.sort(imageTemp, axis=None)[-5:]) )
                
                
                plt.subplot(6,15,3*m+1)
                plt.title('# ' + str(m+fn) + ', (' + str(int(x_ImageJ[m+fn])) + ',' + str(int(y_ImageJ[m+fn]))+')'  , fontsize=9, color = titlecolorTemp, weight=weight)
                plt.tick_params(labelsize=7)
                #plt.imshow(self.centerImage[m+fn][0].reshape(15, 15), interpolation = 'None',vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom', extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                plt.imshow(self.centerImage[m+fn][0].reshape(15, 15), interpolation = 'None',vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='upper', extent=(x3.min(), x3.max(), y3.min(), y3.max())) # trackpy 0.3
                
                plt.subplot(6,15,3*m+2)
                plt.title('R2=' + str(np.around(self.IntensityTrajectoryChosenBS_fit_RsquareFitGoodness[m + fn],2)) + ',T=' + str(int(np.around(self.IntensityTrajectoryChosenBS_fit_T[m + fn],0)) )  , fontsize=9, color = titlecolorTemp, weight=weight)
                plt.tick_params(labelsize=7)
                #plt.imshow(self.centerImage[m+fn][1].reshape(15, 15), interpolation = 'None',vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                plt.imshow(self.centerImage[m+fn][1].reshape(15, 15), interpolation = 'None',vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='upper',extent=(x3.min(), x3.max(), y3.min(), y3.max())) # trackpy 0.3
                                      
                           
                plt.subplot(6,15, 3*m+3)
                plt.tick_params(labelsize=7)
                plt.xlim(0,xlim_max)
                maxFrameN = np.argmax(self.IntensityTrajectoryChosenBS_fit[m + fn])
                #print '\n:maxFrameN = ', maxFrameN 
                if maxFrameN >= 90:
                    maxFrameN = maxFrameN - 90
                #if maxFrameN >= 45:
                #    maxFrameN = maxFrameN - 45                    
                #print 'maxFrameN = ', maxFrameN 
                
                #plt.title( str(maxFrameN*2) + '$^\circ$'  + ', M=' + str(np.around(IntensityTrajectoryChosenBS_fit_M[m + fn],2) )  , fontsize=9, color = titlecolorTemp, weight=weight)
                plt.title( str(int(float(np.around(self.IntensityTrajectoryChosenBS_fit_phi[m + fn],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.IntensityTrajectoryChosenBS_fit_M[m + fn],2) )  , fontsize=9, color = titlecolorTemp, weight=weight)
                
                #print 'len(self.IntensityTrajectoryChosenBS[m + fn]) ', len(self.IntensityTrajectoryChosenBS[m + fn])
                #print 'self.IntensityTrajectoryChosenBS[m + fn]', self.IntensityTrajectoryChosenBS[m + fn]
                #print 'self.IntensityTrajectoryChosenBS[m + fn]*2', self.IntensityTrajectoryChosenBS[m + fn]*2
                
                if len(self.IntensityTrajectoryChosenBS[m + fn]) == 1:
                    plt.plot(([self.IntensityTrajectoryChosenBS[m + fn][0], self.IntensityTrajectoryChosenBS[m + fn][0]]), color = 'OrangeRed')
                    plt.plot(([self.backgroundIntensityTrajectoryChosen[m + fn][0], self.backgroundIntensityTrajectoryChosen[m + fn][0]]), color = 'Gray')
                    plt.plot(([self.backgroundIntensityTrajectory_FitChosen[m + fn][0], self.backgroundIntensityTrajectory_FitChosen[m + fn][0]]), color = 'Cyan')
                    plt.plot(([self.IntensityTrajectoryChosen[m + fn][0], self.IntensityTrajectoryChosen[m + fn][0]]), color = 'Black')
                else:
                    plt.plot((self.IntensityTrajectoryChosenBS[m + fn]), color = 'OrangeRed')
                    plt.plot((self.backgroundIntensityTrajectoryChosen[m + fn]), color = 'Gray')
                    plt.plot((self.backgroundIntensityTrajectory_FitChosen[m + fn]), color = 'Cyan')
                    plt.plot((self.IntensityTrajectoryChosen[m + fn]), color = 'Black')

                
                
                                
                
                
                
                if self.cb_includeFitGraphs_yn.GetValue() == True:
                    plt.plot((self.IntensityTrajectoryChosenBS_fit[m + fn]), color = 'b')
                    plt.locator_params(nbins=4)
                
                
                
                
                    
                #except:
                 #   print '\nNo center image for ', 'M# = ', str(m + fn), '  F# ', str(fsn)


                
                k += 1
                if k%int(NFeaturePerFig) ==0:
                    fn += int(NFeaturePerFig)
                if k == NFeatures:
                    break
      
      
     
        self.xpos_ImageJ = x_ImageJ
        self.ypos_ImageJ = y_ImageJ
        
        
        
        
        

        self.fig_ExcitationM_Histogram = plt.figure('histogram_M', figsize = (12,9))
        self.fig_ExcitationM_Histogram.subplots_adjust(top = 0.89, bottom = 0.07, left = 0.1, right = 0.9, wspace = 0.2, hspace = 0.3)
        self.fig_ExcitationM_Histogram.text(0.02, 0.97, TotalFeatureN_text +',       ' + Corrections, ha="left", va="bottom", size=9,color="black", weight='normal')
        plt.figtext(0.5, 0.92, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly+ '\n' + self.AnalysisConditions ), ha='center', color='black', weight='normal', size=9)

        #plt.title('Total # ' + str(len(IntensityTrajectoryChosenBS_fit_M)))
        plt.subplot(2,2,1)
        if len(self.IntensityTrajectoryChosenBS_fit_M) == 0:
            self.IntensityTrajectoryChosenBS_fit_M = [999,9999]            
        plt.hist(self.IntensityTrajectoryChosenBS_fit_M, bins=np.arange(0.0, 1.3, 0.1))
        plt.xlabel('Modulation Depth (M)')
        plt.xlim(0.0, 1.31)
        
        
        
        '''
        plt.subplot(2,2,2)
        if len(self.IntensityTrajectoryChosenBS_fit_M_GoodFits) == 0:
            self.IntensityTrajectoryChosenBS_fit_M_GoodFits = [999,9999]            

        plt.title('Total # Good T Fits: ' + str(len(self.IntensityTrajectoryChosenBS_fit_M_GoodFits)), size=10)
        plt.hist(self.IntensityTrajectoryChosenBS_fit_M_GoodFits, bins=np.arange(0.0, 1.3, 0.1))
        plt.xlabel('Modulation Depth (M) for Good T Fits')
        plt.xlim(0.0, 1.31)
        '''



        plt.subplot(2,2,2)
        if len(self.IntensityTrajectoryChosenBS_fit_phi) == 0 or len(self.IntensityTrajectoryChosenBS_fit_phi) == 1:
            self.IntensityTrajectoryChosenBS_fit_phi = [11, 111]
        plt.hist(self.IntensityTrajectoryChosenBS_fit_phi)
        plt.xlabel('phi (deg)')



        plt.subplot(2,2,3)
        plt.title('Intensity')
        if len(self.IntensityTrajectoryChosenBS_fit_y0) == 0 or len(self.IntensityTrajectoryChosenBS_fit_y0) == 1:
            self.IntensityTrajectoryChosenBS_fit_y0 = [11, 111]
        plt.hist(self.IntensityTrajectoryChosenBS_fit_y0)
        plt.xlabel('Intensity')
        plt.ylabel('Occurence')





        self.btn_Save_data_figures.Show()
        
        self.text_NFeatures.SetLabel(str(len(self.IntensityTrajectoryChosenBS)))
        
        
        self.btn_SelectData.Show()
        
        
        


        plt.ion()
        plt.show()












                
      
    def SaveFindFeatureData(self, event):
        prefix = self.FilePrefix.GetValue()
        figFileName = prefix + self.filenameOnly[:30]
        
        

        todayDate = time.strftime("%Y%m%d_%Hh%Mm")
        print "\nToday's Date: ", time.strftime("%Y-%m-%d_%Hh%Mm")
        
        
        self.fig_FindingFeatures.savefig(self.filedirectoryOnly + '\\' + figFileName + '_' + todayDate + '_0_FindingFeatureFigure.png')
        



        ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_0_FindingFeatureData(x,y,I).txt','w')
        for i in range(len(self.FindingFeature_xposImageJ)):
            ff.write(str(int(self.FindingFeature_xposImageJ[i])) + ' ' + str(int(self.FindingFeature_yposImageJ[i])) + ' ' + str(int(self.FindingFeature_IntensityData[i])) + '\n'   )
        ff.close()       


        



                
      
    def SaveFindFeatureData_filtered(self, event):
        prefix = self.FilePrefix.GetValue()
        figFileName = prefix + self.filenameOnly[:30]
        
        

        todayDate = time.strftime("%Y%m%d_%Hh%Mm")
        print "\nToday's Date: ", time.strftime("%Y-%m-%d_%Hh%Mm")
        
        
        self.fig_FindingFeatures_filtered.savefig(self.filedirectoryOnly + '\\' + figFileName + '_' + todayDate + '_0_FindingFeatureFigure_filtered.png')
        



        ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_0_FindingFeatureDataFiltered(x,y,I).txt','w')
        for i in range(len(self.FindingFeature_xposImageJ)):
            ff.write(str(int(self.xposImageJFound_temp_filtered[i])) + ' ' + str(int(self.yposImageJFound_temp_filtered[i])) + ' ' + str(int(self.intensity_temp_filtered[i])) + '\n'   )
        ff.close()       


        






        
                
      
    def SaveDataAndFigures(self, event):
        prefix = self.FilePrefix.GetValue()
        if prefix != '':
            prefix += '_'
        figFileName = prefix + self.filenameOnly[:40]
        
        

        todayDate = time.strftime("%Y%m%d_%Hh%Mm")
        print "\nToday's Date: ", time.strftime("%Y-%m-%d_%Hh%Mm")
        

        save_data_error = 0



        try:
            for n in range(len(self.fig_AllDataCenterImagesOnly)):    
                plt.close(self.fig_AllDataCenterImagesOnly[n])
        except:
            pass


        
        if self.cb_SaveFigs_yn.GetValue():          
        
            self.fig_ShowAfewFrames.savefig(self.filedirectoryOnly + '\\' + figFileName + '_' + todayDate + '_1_RawImageFrames.png')
            plt.close(self.fig_ShowAfewFrames)
            print 'Raw Frames figures are saved'
        
            try:
                self.fig_FindingFeatures.savefig(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_2_Molecules.png')
                plt.close(self.fig_FindingFeatures)
                print 'Feature finding figure is saved'
            except:
                pass


            try:
                self.fig_FindingFeatures_filtered.savefig(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_2_Molecules2.png')
                plt.close(self.fig_FindingFeatures_filtered)
                print 'Feature finding figure_filtered is saved'
            except:
                pass


            try:
                    
                self.fig_BG_ref.savefig(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_3_Background.png')
                plt.close(self.fig_BG_ref)
                print 'fig_BG_ref figure is saved'
            
            except:
                save_data_error += 1
                pass
            
    
    
    
            
            
    
    
    
    
            '''
            try:        
                for n in range(len(self.fig_ExcitationM_Trajectory)):    
                    self.fig_ExcitationM_Trajectory[n].savefig(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate + '_3_ExcitationM_' + str(n) + '.png')
                    #plt.close(self.fig_ExcitationM_Trajectory[n])
                #self.fig_LDAveHist.savefig(figFileName + '_' + todayDate + '_4_LD_Histogram.png')
                print 'fig_ExcitationM_Trajectory figures are saved'
                
            except:
                print 'fig_ExcitationM_Trajectory figures are NOT saved: error'  
                save_fig_error += 1
                
            ###'''
                
                
                
              
            try:        
                for n in range(len(self.fig_AllDataCenterImages)):    
                    self.fig_AllDataCenterImages[n].savefig(self.filedirectoryOnly + '\\' + figFileName + '_' + todayDate + '_4_DataImages_' + str(n) + '.png')
                    plt.close(self.fig_AllDataCenterImages[n])
                #self.fig_LDAveHist.savefig(figFileName + '_' + todayDate + '_4_LD_Histogram.png')
                print 'fig_AllDataCenterImages figures are saved'
                
            except:
                print 'fig_AllDataCenterImages figures are NOT saved: error'  
                save_data_error += 1
                
    
    
    
            self.fig_ExcitationM_Histogram.savefig(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_5_M_histogram.png')
            plt.close(self.fig_ExcitationM_Histogram)
            print 'fig_ExcitationM_Histogram figure is saved'
    
            
    
    
            plt.close(self.fig_BG_ref)
            
            

        
        if self.cb_SaveImageJ_xyData_yn.GetValue():
            print ''
            

            ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_6_xy_A_T(frames)_phi(deg)_y0_M.txt','w')
            
            for i in range(len(self.xpos_ImageJ)):
                
                ff.write(str(self.xposChosen[i]) + ' ' + str(self.yposChosen[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_A[i])
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_T[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_phi[i]) 
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_y0[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_M[i]) + '\n'   )    
                
            ff.close()       



        if self.cb_SaveImageJ_xyDataGoodFit_yn.GetValue():

            ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_7_xy_A_T(frames)_phi(deg)_y0_M_GoodTfits.txt','w')
            
            
            for i in range(len(self.xposChosen_GoodFits)):
                
                ff.write(str(self.xposChosen_GoodFits[i]) + ' ' + str(self.yposChosen_GoodFits[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_A_GoodFits[i])
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_T_GoodFits[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_phi_GoodFits[i]) 
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_y0_GoodFits[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_M_GoodFits[i]) + '\n'   )    
                
            ff.close()       
            
            
        if self.cb_SaveIntensityData_yn.GetValue():
            
            for n in range(len(self.IntensityTrajectoryChosenBS)):
                    
                ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_8_Intensity_' + str(n) + '.txt','w')
                
                
                for i in range(len(self.IntensityTrajectoryChosenBS[n])):
                    
                    ff.write(str(int(self.IntensityTrajectoryChosenBS[n][i])) + '\n')
                    
                ff.close()       
                
            
                       

        if self.cb_SaveIntensityFitData_yn.GetValue():
            
            for n in range(len(self.IntensityTrajectoryChosenBS)):
                    
                ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_9_Fit_' + str(n) + '.txt','w')
                

                for i in range(len(self.IntensityTrajectoryChosenBS_fit[n])):
                    
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_fit[n][i])) + '\n')
                    
                
                ff.close()       
                
            
   

                    

        if self.cb_Complete_Data_yn.GetValue():
            
            ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_10_CompleteDataAll' + '.txt','w')
            
            ff.write('%% x y A T phi y0 M Background_Subtraction_Type: '+ self.BGS_type + ' - rawIntensity , intensityBS , fit , BG , BGfit , image0 , image1' + '\n')
            
            for n in range(len(self.IntensityTrajectoryChosenBS)):
                

                ff.write('#' + str(n) + '\n')
                
                    

                    
                ff.write(str(self.xposChosen[n]) + ' ' + str(self.yposChosen[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_A[n])
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_T[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_phi[n]) 
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_y0[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_M[n]) + '\n'   )    
                



                for i in range(len(self.IntensityTrajectoryChosen[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosen[n][i])) + ' ')
                ff.write('\n')



                for i in range(len(self.IntensityTrajectoryChosenBS[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosenBS[n][i])) + ' ')
                ff.write('\n')
                
                for i in range(len(self.IntensityTrajectoryChosenBS_fit[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_fit[n][i])) + ' ')
                ff.write('\n')
                

                for i in range(len(self.backgroundIntensityTrajectoryChosen[n])):
                    ff.write(str(int(self.backgroundIntensityTrajectoryChosen[n][i])) + ' ')
                ff.write('\n')


                for i in range(len(self.backgroundIntensityTrajectory_FitChosen[n])):
                    ff.write(str(int(self.backgroundIntensityTrajectory_FitChosen[n][i])) + ' ')
                ff.write('\n')





                
                for i in range(len(self.centerImage[n][0])):
                    for k in range(len(self.centerImage[n][0][i])):
                        ff.write(str(int(self.centerImage[n][0][i][k])) + ' ')
                    ff.write(' , ')
                    
                ff.write('\n')
                
                
                for i in range(len(self.centerImage[n][1])):
                    for k in range(len(self.centerImage[n][1][i])):
                        ff.write(str(int(self.centerImage[n][1][i][k])) + ' ')
                    ff.write(' , ')
                    
                    
                
                    
                ff.write('\n\n')
                                
                                  
                    
            ff.close()       
                
            
  


                     



        
        
        if save_data_error == 0 :
            save_fig_text = 'All Data are saved'
            print '\n', save_fig_text, '\n'
        else:
            save_fig_text = 'NOT All Data are saved'
            print '\n NOT All Data are saved\n # error: ', save_data_error
            
        
        
                
        
        
        self.text_Save_figures.SetLabel(save_fig_text)
     
     
     
     
     
    def SelectData(self, event):
        #self.SelectData_window = wx.Frame(None, title='SelectData', size=(800, 1000), pos=(800,0))
        #self.SelectData_window.Show()
        

        self.SelectData_frame = wx.Frame(self, title='SelectData', size=(850, 1000), pos=(800,0))
        self.SelectData_window = wx.ScrolledWindow(self.SelectData_frame, -1)
        self.SelectData_window.SetScrollbars(1, 1, 800, int(len(self.xpos_ImageJ)/10)*25 + 200)
        self.SelectData_frame.Show()
        
        
        
        

        wx.StaticText(self.SelectData_window, -1, "Data File Prefix", pos=(160, 35))
        self.FilePrefixForSelected = wx.TextCtrl(self.SelectData_window, -1, str(FilePrefixInput), pos=(250, 30), size=(150, -1))
        
        
        btn_ShowFiguresSelectedData = wx.Button(self.SelectData_window, pos=(50,30), label="Show Figures")
        btn_ShowFiguresSelectedData.Bind(wx.EVT_BUTTON, self.SelectData_ShowFigures)
        
        self.cb_EachFigure_Selected_yn = wx.CheckBox(self.SelectData_window, -1, 'Show individual images', (50, 55))
        self.cb_EachFigure_Selected_yn.SetValue(False)
        
        
        
        
        
        self.btn_SaveSelectedData = wx.Button(self.SelectData_window, pos=(420,30), label="Save Selected Data")
        self.btn_SaveSelectedData.Bind(wx.EVT_BUTTON, self.SaveSelectedData)
        self.btn_SaveSelectedData.Hide()
        
        self.text_SaveSelectedData = wx.StaticText(self.SelectData_window, -1, "...", pos=(420, 60))
        

        self.cb_SaveFigs_Selected_yn = wx.CheckBox(self.SelectData_window, -1, 'Figures', (550, 10))
        self.cb_SaveFigs_Selected_yn.SetValue(True)
        self.cb_SaveImageJ_xyData_Selected_yn = wx.CheckBox(self.SelectData_window, -1, 'ImageJ x,y data & Fit: A, T, phase, y0, M', (550, 30))
        self.cb_SaveImageJ_xyData_Selected_yn.SetValue(False)
        self.cb_SaveIntensityData_Selected_yn = wx.CheckBox(self.SelectData_window, -1, 'Intensity Trajectory Data for all', (550, 50))
        self.cb_SaveIntensityData_Selected_yn.SetValue(False)
        self.cb_SaveIntensityFitData_Selected_yn = wx.CheckBox(self.SelectData_window, -1, 'Fit Trajectory Data for all', (550, 70))
        self.cb_SaveIntensityFitData_Selected_yn.SetValue(False)
        self.cb_Complete_Data_Selected_yn = wx.CheckBox(self.SelectData_window, -1, 'Complete Data File', (550, 90))
        self.cb_Complete_Data_Selected_yn.SetValue(True)


        
        
        wx.StaticText(self.SelectData_window, -1, 'Select Data', pos=(50, 110)) 
                
        #self.xpos_ImageJ = [[] for _ in xrange(150)]
        self.cb_Data_yn = [[] for _ in xrange(len(self.xpos_ImageJ))]
        rn = 0
        cn = 0
        countertemp = 0
        for n in range(len(self.xpos_ImageJ)):
            self.cb_Data_yn[n] = wx.CheckBox(self.SelectData_window, -1, '# '+str(n), (50 + 70*cn, 145 + 25*rn))
            self.cb_Data_yn[n].SetValue(True)
            countertemp += 1
            cn += 1
            if countertemp % 10 == 0:
                rn += 1
                cn = 0
            
        self.btn_SelectAll = wx.Button(self.SelectData_window, pos=(150,105), label="Select All")
        self.btn_SelectAll.Bind(wx.EVT_BUTTON, self.SelectDataAll)
        self.btn_SelectAll.Show()
        
        self.btn_SelectNone = wx.Button(self.SelectData_window, pos=(250,105), label="Deselect All")
        self.btn_SelectNone.Bind(wx.EVT_BUTTON, self.SelectDataNone)
        self.btn_SelectNone.Show()
        
        
    def SelectDataAll(self, event):
        print '\nSelect All '
        for n in range(len(self.cb_Data_yn)):
            self.cb_Data_yn[n].SetValue(True)
    def SelectDataNone(self, event):
        print '\nDeselect All '
        for n in range(len(self.cb_Data_yn)):
            self.cb_Data_yn[n].SetValue(False)
            
                        
            
    def SelectData_ShowFigures(self, event):
        print '\n SelectData_ShowFigures '
        
        try:
            for n in range(len(self.fig_AllDataCenterImages_Selected)):
                plt.close(self.fig_AllDataCenterImages_Selected[n])
                plt.close(self.fig_ExcitationM_Histogram_Selected[n])
                                
        except: pass
        
        
        
        
        self.NFeatures_Selected = 0
        self.x_ImageJ_Selected = []
        self.y_ImageJ_Selected = []
        self.centerImage_Selected = []
        
        self.IntensityTrajectoryChosen_Selected = []
        self.IntensityTrajectoryChosenBS_Selected = []
        self.IntensityTrajectoryChosenBS_fit_Selected = []
        self.backgroundIntensityTrajectoryChosen_Selected = []
        self.backgroundIntensityTrajectoryChosen_Fit_Selected = []
        
        self.IntensityTrajectoryChosenBS_fit_A_Selected = []
        self.IntensityTrajectoryChosenBS_fit_T_Selected = []
        self.IntensityTrajectoryChosenBS_fit_phi_Selected = []
        self.IntensityTrajectoryChosenBS_fit_y0_Selected = []
        
        
        self.IntensityTrajectoryChosenBS_fit_M_Selected = []
        
        self.IntensityTrajectoryChosenBS_fit_RsquareFitGoodness_Selected = []
        
        
        
        
        
        
        for n in range(len(self.cb_Data_yn)):
            if self.cb_Data_yn[n].GetValue():
                self.NFeatures_Selected += 1
                
                self.x_ImageJ_Selected.append(self.xposChosen[n])
                self.y_ImageJ_Selected.append(self.yposChosen[n])
                
                self.centerImage_Selected.append(self.centerImage[n])
                
                self.IntensityTrajectoryChosen_Selected.append(self.IntensityTrajectoryChosen[n])
                self.IntensityTrajectoryChosenBS_Selected.append(self.IntensityTrajectoryChosenBS[n])
                self.IntensityTrajectoryChosenBS_fit_Selected.append(self.IntensityTrajectoryChosenBS_fit[n])
                self.backgroundIntensityTrajectoryChosen_Selected.append(self.backgroundIntensityTrajectoryChosen[n])
                self.backgroundIntensityTrajectoryChosen_Fit_Selected.append(self.backgroundIntensityTrajectory_FitChosen[n])

                self.IntensityTrajectoryChosenBS_fit_A_Selected.append(self.IntensityTrajectoryChosenBS_fit_A[n])
                self.IntensityTrajectoryChosenBS_fit_T_Selected.append(self.IntensityTrajectoryChosenBS_fit_T[n])
                self.IntensityTrajectoryChosenBS_fit_phi_Selected.append(self.IntensityTrajectoryChosenBS_fit_phi[n])
                self.IntensityTrajectoryChosenBS_fit_y0_Selected.append(self.IntensityTrajectoryChosenBS_fit_y0[n])
                
                self.IntensityTrajectoryChosenBS_fit_M_Selected.append(self.IntensityTrajectoryChosenBS_fit_M[n])

                
                self.IntensityTrajectoryChosenBS_fit_RsquareFitGoodness_Selected.append(self.IntensityTrajectoryChosenBS_fit_RsquareFitGoodness[n])
                        
        
        
        self.IntensityTrajectoryChosenBS_fit_M_GoodFits_Selected = []
        self.IntensityTrajectoryChosenBS_fit_phi_GoodFits_Selected = []
        
        for n in range(len(self.x_ImageJ_Selected)):  #finding good fit data
            T_temp = self.IntensityTrajectoryChosenBS_fit_T_Selected[n]
            if ( T_temp <= (float(self.RefPoint_F2N.GetValue()) - 1)*self.MaxTfitRatio ) and ( T_temp >= (float(self.RefPoint_F2N.GetValue()) - 1)*self.MinTfitRatio ) :
                self.IntensityTrajectoryChosenBS_fit_M_GoodFits_Selected.append(self.IntensityTrajectoryChosenBS_fit_M_Selected[n])
                
        
                
        
        

        centerLength = int(self.centerPixelN.GetValue())

        TotalFeatureN_text = 'Total # ' + str(self.NFeatures_Selected) + ' ,   Center: ' + str(centerLength) + 'x' +  str(centerLength) + ' ,   Background: ' + str(centerLength+2) + 'x' + str(centerLength+2)
         



        
        if self.cb_EachFigure_Selected_yn.GetValue():
                   
      
            NFeaturePerFig = 30.0
           
            NCfig = int(  math.ceil((self.NFeatures_Selected/float(NFeaturePerFig))))
            self.fig_AllDataCenterImages_Selected = [[] for _ in xrange(NCfig)]
    
            x3 = np.linspace(0, 15, 16)
            y3 = np.linspace(0, 15, 16)
            x3, y3 = np.meshgrid(x3, y3)
            
            #fsn = int(self.Nframe.GetValue()) # frame start number
            k = 0
            fn = 0
            for n in range(NCfig):
                #print 'n = ', n
                self.fig_AllDataCenterImages_Selected[n] = plt.figure('Raw_Data_Images_M_s'+ str(n), figsize = (18, 9))
                self.fig_AllDataCenterImages_Selected[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.36, hspace = 0.4)
                self.fig_AllDataCenterImages_Selected[n].text(0.05, 0.93, TotalFeatureN_text +'\n' + self.Corrections2, ha="left", va="bottom", size="small",color="black", weight='normal')
                
                plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly+ '\n' + self.AnalysisConditions ), ha='center', color='black', weight='normal', size='small')
                    
    
    
                for m in range(int(NFeaturePerFig)):
    
                    #try:
                    #if self.IntensityTrajectoryChosenBS_fit_RsquareFitGoodness[m + fn] >= float(self.MinRsquareFitGoodness.GetValue()):
                    if ( self.IntensityTrajectoryChosenBS_fit_T_Selected[m + fn] <= (float(self.RefPoint_F2N.GetValue()) - 1)*self.MaxTfitRatio ) and ( self.IntensityTrajectoryChosenBS_fit_T_Selected[m + fn] >= (float(self.RefPoint_F2N.GetValue()) - 1)*self.MinTfitRatio ):
                        titlecolorTemp = 'black'
                        weight='normal'
                    else:
                        titlecolorTemp = 'red'
                        weight='normal'
                    
                    
                    imageTemp = np.append(np.ravel(self.centerImage_Selected[m+fn][0]), np.ravel(self.centerImage_Selected[m+fn][1]))
                    
                    vminTemp = int(np.mean( np.sort(imageTemp)[:5]) )
                    vmaxTemp = int(np.mean( np.sort(imageTemp)[-5:]) )
                    
                    
                    plt.subplot(6,15,3*m+1)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.x_ImageJ_Selected[m+fn])) + ',' + str(int(self.y_ImageJ_Selected[m+fn]))+')'  , fontsize=9, color = titlecolorTemp, weight=weight)
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.centerImage_Selected[m+fn][0].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    
                    plt.subplot(6,15,3*m+2)
                    plt.title('R2=' + str(np.around(self.IntensityTrajectoryChosenBS_fit_RsquareFitGoodness_Selected[m + fn],2)) + ',T=' + str(int(np.around(self.IntensityTrajectoryChosenBS_fit_T_Selected[m + fn],0)) )  , fontsize=9, color = titlecolorTemp, weight=weight)
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.centerImage_Selected[m+fn][1].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                                          
                               
                    plt.subplot(6,15, 3*m+3)
                    plt.tick_params(labelsize=7)
                    maxFrameN = np.argmax(self.IntensityTrajectoryChosenBS_fit_Selected[m + fn])
                    #print '\n:maxFrameN = ', maxFrameN 
                    if maxFrameN >= 90:
                        maxFrameN = maxFrameN - 90
                    #if maxFrameN >= 45:
                    #    maxFrameN = maxFrameN - 45                    
                    #print 'maxFrameN = ', maxFrameN 
                    
                    #plt.title( str(maxFrameN*2) + '$^\circ$'  + ', M=' + str(np.around(IntensityTrajectoryChosenBS_fit_M[m + fn],2) )  , fontsize=9, color = titlecolorTemp, weight=weight)
                    plt.title( str(int(float(np.around(self.IntensityTrajectoryChosenBS_fit_phi_Selected[m + fn],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.IntensityTrajectoryChosenBS_fit_M_Selected[m + fn],2) )  , fontsize=9, color = titlecolorTemp, weight=weight)
    
                    if len(self.IntensityTrajectoryChosenBS_Selected[m + fn]) == 1:
                        plt.plot(([self.IntensityTrajectoryChosenBS_Selected[m + fn][0], self.IntensityTrajectoryChosenBS_Selected[m + fn][0]]), color = 'OrangeRed')
                        plt.plot(([self.backgroundIntensityTrajectoryChosen_Selected[m + fn][0], self.backgroundIntensityTrajectoryChosen_Selected[m + fn][0]]), color = 'Gray')
                        plt.plot(([self.backgroundIntensityTrajectoryChosen_Fit_Selected[m + fn][0], self.backgroundIntensityTrajectoryChosen_Fit_Selected[m + fn][0]]), color = 'Cyan')
                        plt.plot(([self.IntensityTrajectoryChosen_Selected[m + fn][0], self.IntensityTrajectoryChosen_Selected[m + fn][0]]), color = 'Black')
                    else:
                        plt.plot((self.IntensityTrajectoryChosenBS_Selected[m + fn]), color = 'OrangeRed')
                        plt.plot((self.backgroundIntensityTrajectoryChosen_Selected[m + fn]), color = 'Gray')
                        plt.plot((self.backgroundIntensityTrajectoryChosen_Fit_Selected[m + fn]), color = 'Cyan')
                        plt.plot((self.IntensityTrajectoryChosen_Selected[m + fn]), color = 'Black')
        
                    
                    
    
                    
                    if self.cb_includeFitGraphs_yn.GetValue() == True:
                        plt.plot((self.IntensityTrajectoryChosenBS_fit_Selected[m + fn]), color = 'b')
                        plt.locator_params(nbins=4)
                    
                    
                    
                    
                        
                    #except:
                     #   print '\nNo center image for ', 'M# = ', str(m + fn), '  F# ', str(fsn)
    
    
                    
                    k += 1
                    if k%int(NFeaturePerFig) ==0:
                        fn += int(NFeaturePerFig)
                    if k == self.NFeatures_Selected:
                        break
             
                
        

        self.fig_ExcitationM_Histogram_Selected = plt.figure('histogram_M_S', figsize = (12,9))
        self.fig_ExcitationM_Histogram_Selected.subplots_adjust(top = 0.89, bottom = 0.07, left = 0.1, right = 0.9, wspace = 0.2, hspace = 0.3)
        self.fig_ExcitationM_Histogram_Selected.text(0.02, 0.97, TotalFeatureN_text +',       ' + self.Corrections, ha="left", va="bottom", size=9,color="black", weight='normal')
        plt.figtext(0.5, 0.92, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly+ '\n' + self.AnalysisConditions ), ha='center', color='black', weight='normal', size=9)

        #plt.title('Total # ' + str(len(IntensityTrajectoryChosenBS_fit_M)))
        plt.subplot(2,2,1)
        if len(self.IntensityTrajectoryChosenBS_fit_M_Selected) == 0:
            self.IntensityTrajectoryChosenBS_fit_M_Selected = [999,9999]            
        plt.hist(self.IntensityTrajectoryChosenBS_fit_M_Selected, bins=np.arange(0.0, 1.3, 0.1))
        plt.xlabel('Modulation Depth (M)')
        plt.xlim(0.0, 1.21)
        
        
        '''
        plt.subplot(2,2,2)
        if len(self.IntensityTrajectoryChosenBS_fit_M_GoodFits_Selected) == 0:
            self.IntensityTrajectoryChosenBS_fit_M_GoodFits_Selected = [999,9999]            

        plt.title('Total # Good T Fits: ' + str(len(self.IntensityTrajectoryChosenBS_fit_M_GoodFits_Selected)), size=10)
        plt.hist(self.IntensityTrajectoryChosenBS_fit_M_GoodFits_Selected, bins=np.arange(0.0, 1.3, 0.1))
        plt.xlabel('Modulation Depth (M) for Good T Fits')
        plt.xlim(0.0, 1.21)
        '''
        

        plt.subplot(2,2,2)
        if len(self.IntensityTrajectoryChosenBS_fit_phi_Selected) == 0 or len(self.IntensityTrajectoryChosenBS_fit_phi_Selected) == 1:
            self.IntensityTrajectoryChosenBS_fit_phi_Selected = [11, 111]
        plt.hist(self.IntensityTrajectoryChosenBS_fit_phi_Selected)
        plt.xlabel('phi (deg)')


 
        plt.subplot(2,2,3)
        plt.title('Intensity')
        if len(self.IntensityTrajectoryChosenBS_fit_y0_Selected) == 0 or len(self.IntensityTrajectoryChosenBS_fit_y0_Selected) == 1:
            self.IntensityTrajectoryChosenBS_fit_y0_Selected = [11, 111]
        plt.hist(self.IntensityTrajectoryChosenBS_fit_y0_Selected)
        plt.xlabel('Intensity')
        plt.ylabel('Occurence')





        self.btn_SaveSelectedData.Show()
        
        plt.ion()
        
        plt.show()
        
        print 'SelectData_ShowFigures done'
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    def SaveSelectedData(self, event):
        print ' '
        self.text_SaveSelectedData.SetLabel('Saving Data ....')
        
        prefix = self.FilePrefixForSelected.GetValue()
        if prefix != '':
            prefix += '_'
        figFileName = prefix + self.filenameOnly[:40]
        todayDate = time.strftime("%Y%m%d_%Hh%Mm")        



        if self.cb_SaveFigs_Selected_yn.GetValue():         
            
            
            
        
            self.fig_ShowAfewFrames.savefig(self.filedirectoryOnly + '\\' + figFileName + '_' + todayDate + '_1_RawImageFrames.png')
            plt.close(self.fig_ShowAfewFrames)
            print 'Raw Frames figures are saved'
        
            try:    
                self.fig_FindingFeatures.savefig(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_2_Molecules.png')
                plt.close(self.fig_FindingFeatures)
                print 'Feature finding figure is saved'
            except:
                pass


            try:
                self.fig_FindingFeatures_filtered.savefig(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_2_Molecules2.png')
                plt.close(self.fig_FindingFeatures_filtered)
                print 'Feature finding figure_filtered is saved'
            except:
                pass


            
            try:    
                self.fig_BG_ref.savefig(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_3_Background.png')
                plt.close(self.fig_BG_ref)
                print 'fig_BG_ref figure is saved'
            except:
                pass
            
    
            
            
            
            try:        
                for n in range(len(self.fig_AllDataCenterImages_Selected)):    
                    self.fig_AllDataCenterImages_Selected[n].savefig(self.filedirectoryOnly + '\\' + figFileName + '_' + todayDate + '_4_DataImages_' + str(n) + '.png')
                    plt.close(self.fig_AllDataCenterImages_Selected[n])
                #self.fig_LDAveHist.savefig(figFileName + '_' + todayDate + '_4_LD_Histogram.png')
                print 'fig_AllDataCenterImages figures are saved'
                
            except:
                print 'fig_AllDataCenterImages figures are NOT saved: error'  
                
                
    
    
    
            self.fig_ExcitationM_Histogram_Selected.savefig(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_5_M_histogram.png')
            plt.close(self.fig_ExcitationM_Histogram_Selected)
            print 'fig_ExcitationM_Histogram figure is saved'
    
            
    
    
            
        if self.cb_SaveImageJ_xyData_Selected_yn.GetValue():
            print ''
            

            ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_6_xy_A_T(frames)_phi(deg)_y0_M.txt','w')
            
            for i in range(len(self.x_ImageJ_Selected)):
                
                ff.write(str(self.x_ImageJ_Selected[i]) + ' ' + str(self.y_ImageJ_Selected[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_A_Selected[i])
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_T_Selected[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_phi_Selected[i]) 
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_y0_Selected[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_M_Selected[i]) + '\n'   )    
                
            ff.close()       

            
            
        if self.cb_SaveIntensityData_Selected_yn.GetValue():
            
            for n in range(len(self.IntensityTrajectoryChosenBS_Selected)):
                    
                ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_8_Intensity_' + str(n) + '.txt','w')
                
                
                for i in range(len(self.IntensityTrajectoryChosenBS_Selected[n])):
                    
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_Selected[n][i])) + '\n')
                    
                ff.close()       
                
            
                       

        if self.cb_SaveIntensityFitData_Selected_yn.GetValue():
            
            for n in range(len(self.IntensityTrajectoryChosenBS_fit_Selected)):
                    
                ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_9_Fit_' + str(n) + '.txt','w')
                
                
                for i in range(len(self.IntensityTrajectoryChosenBS_fit_Selected[n])):
                    
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_fit_Selected[n][i])) + '\n')
                    
                ff.close()       
                




        if self.cb_Complete_Data_Selected_yn.GetValue():
            
            
            ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_10_CompleteDataAll_Selected' + '.txt','w')
            
            ff.write('%% x y A T phi y0 M Background_Subtraction_Type: '+ self.BGS_type + ' - rawIntensity , intensityBS , fit , BG , BGfit , image0 , image1' + '\n')
            
            for n in range(len(self.IntensityTrajectoryChosenBS_Selected)):
                    

                
                ff.write('#' + str(n) + '\n')
                
                    

                    
                ff.write(str(self.x_ImageJ_Selected[n]) + ' ' + str(self.y_ImageJ_Selected[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_A_Selected[n])
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_T_Selected[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_phi_Selected[n]) 
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_y0_Selected[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_M_Selected[n]) + '\n'   )    
                



                for i in range(len(self.IntensityTrajectoryChosen_Selected[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosen_Selected[n][i])) + ' ')
                ff.write('\n')


                for i in range(len(self.IntensityTrajectoryChosenBS_Selected[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_Selected[n][i])) + ' ')
                ff.write('\n')
                

                for i in range(len(self.IntensityTrajectoryChosenBS_fit_Selected[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_fit_Selected[n][i])) + ' ')
                ff.write('\n')
                

                for i in range(len(self.backgroundIntensityTrajectoryChosen_Selected[n])):
                    ff.write(str(int(self.backgroundIntensityTrajectoryChosen_Selected[n][i])) + ' ')
                ff.write('\n')

                
                for i in range(len(self.backgroundIntensityTrajectoryChosen_Fit_Selected[n])):
                    ff.write(str(int(self.backgroundIntensityTrajectoryChosen_Fit_Selected[n][i])) + ' ')
                ff.write('\n')
  






              
                for i in range(len(self.centerImage_Selected[n][0])):
                    for k in range(len(self.centerImage_Selected[n][0][i])):
                        ff.write(str(int(self.centerImage_Selected[n][0][i][k])) + ' ')
                    ff.write(' , ')
                    
                ff.write('\n')
                
                
                for i in range(len(self.centerImage_Selected[n][1])):
                    for k in range(len(self.centerImage_Selected[n][1][i])):
                        ff.write(str(int(self.centerImage_Selected[n][1][i][k])) + ' ')
                    ff.write(' , ')
                    
                ff.write('\n\n')
                                
                                  
                    
            ff.close()       
            
        self.text_SaveSelectedData.SetLabel('Data Saved')
     
     
     
     

        
        
          
'''
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
'''





class OpenDataFile(wx.Frame):
    def __init__(self):
        print ' '
        """
        Create and show the Open FileDialog
        """
        dlg = wx.FileDialog(None, message="Choose a file", defaultFile="", wildcard=wildcardDataFile, style=wx.OPEN | wx.MULTIPLE) # | wx.CHANGE_DIR
            
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            print "You chose the following file(s):"
            for path in paths:
                print path
        print dlg
        print 'path = ', dlg.GetPath()
        print 'filename = ', dlg.GetFilename()
        
        
        
        
        self.fileName = dlg.GetFilename()
        self.fileDirectory = dlg.GetDirectory()
        
        
        fdata = open(paths[0], 'r')
        linestemp = fdata.readlines()
        fdata.close()
        
        
        print 'linestemp[0][0]:', linestemp[0][0]
        print '(linestemp[0].split())[0]:', (linestemp[0].split())[0]
        


        
        if (linestemp[0].split())[0] == '#0':
            
            self.Data_type_indicator = (linestemp[0].split())[0]
            
            print 'len(linestemp) ', len(linestemp)
            
            self.Nfeatures = int((len(linestemp))/7)
            print 'self.Nfeatures ', self.Nfeatures
            
            self.Loaded_ImageJxpos = []
            self.Loaded_ImageJypos = []
            self.Loaded_A = []
            self.Loaded_T = []
            self.Loaded_phi = []
            self.Loaded_y0 = []
            self.Loaded_M = []
            
            self.Loaded_CenterImages = []
            self.Loaded_Intensity = []
            self.Loaded_Intensity_Fit = []
            
            for n in range(self.Nfeatures):
                p1 = linestemp[7*n + 1].split()
                #print 'p1 ', p1
                
                self.Loaded_ImageJxpos.append(int(float(p1[0])) ) 
                self.Loaded_ImageJypos.append(int(float(p1[1])) )
                self.Loaded_A.append((float(p1[2])) )
                self.Loaded_T.append((float(p1[3])) )
                self.Loaded_phi.append((float(p1[4])) )
                self.Loaded_y0.append((float(p1[5])) )
                self.Loaded_M.append((float(p1[6])) )
                
                p2 = linestemp[7*n + 2].split()
                intensityTemp = []
                for k in p2: #Intensity data
                    intensityTemp.append(float(k)) 
                self.Loaded_Intensity.append(intensityTemp)
                    
                p3 = linestemp[7*n + 3].split()
                intensityFitTemp = []
                for k in p3: #Intensity Fit data
                    intensityFitTemp.append(float(k)) 
                self.Loaded_Intensity_Fit.append(intensityFitTemp)
                
 
                p5 = linestemp[7*n + 4].split()
                ImageTemp1 = []
                for k in range(15): #Intensity Image 1 data
                    temp = []
                    for p in range(15):
                        temp.append(float(p5[16*k + p])) 
                    ImageTemp1.append(temp)
                    
                    
                p6 = linestemp[7*n + 5].split()
                ImageTemp2 = []
                for k in range(15): #Intensity Image 1 data
                    temp = []
                    for p in range(15):
                        temp.append(float(p6[16*k + p])) 
                    ImageTemp2.append(temp)
                    
                self.Loaded_CenterImages.append([ImageTemp1, ImageTemp2])
                
            self.Loaded_Intensity_BG = np.array(self.Loaded_Intensity)*0
            self.Loaded_Intensity_BG_Fit = np.array(self.Loaded_Intensity_BG) * 0
            self.Loaded_Intensity_noBGS = np.array(self.Loaded_Intensity)*0
            self.BGS_type = 'n/a' 


        
        elif (linestemp[0].split())[0] == '##':    
            
            self.Data_type_indicator = (linestemp[0].split())[0]
            
            print 'len(linestemp) ', len(linestemp)
            
            self.Nfeatures = int((len(linestemp)-1)/8)
            print 'self.Nfeatures ', self.Nfeatures
            
            self.Loaded_ImageJxpos = []
            self.Loaded_ImageJypos = []
            self.Loaded_A = []
            self.Loaded_T = []
            self.Loaded_phi = []
            self.Loaded_y0 = []
            self.Loaded_M = []
            
            self.Loaded_CenterImages = []
            self.Loaded_Intensity = []
            self.Loaded_Intensity_Fit = []
            self.Loaded_Intensity_BG = []
            
            for n in range(self.Nfeatures):
                p1 = linestemp[8*n + 2].split()
                #print 'p1 ', p1
                
                self.Loaded_ImageJxpos.append(int(float(p1[0])) ) 
                self.Loaded_ImageJypos.append(int(float(p1[1])) )
                self.Loaded_A.append((float(p1[2])) )
                self.Loaded_T.append((float(p1[3])) )
                self.Loaded_phi.append((float(p1[4])) )
                self.Loaded_y0.append((float(p1[5])) )
                self.Loaded_M.append((float(p1[6])) )
                
                p2 = linestemp[8*n + 3].split()
                intensityTemp = []
                for k in p2: #Intensity data
                    intensityTemp.append(float(k)) 
                self.Loaded_Intensity.append(np.array(intensityTemp))
                    
                p3 = linestemp[8*n + 4].split()
                intensityFitTemp = []
                for k in p3: #Intensity Fit data
                    intensityFitTemp.append(float(k)) 
                self.Loaded_Intensity_Fit.append(np.array(intensityFitTemp))
                
 
                p4 = linestemp[8*n + 5].split()
                intensityFitTemp = []
                intensityBGTemp = []
                for k in p4: #Intensity BG data
                    intensityBGTemp.append(float(k)) 
                self.Loaded_Intensity_BG.append(np.array(intensityBGTemp))
       

               
                p5 = linestemp[8*n + 6].split()
                ImageTemp1 = []
                for k in range(15): #Intensity Image 1 data
                    temp = []
                    for p in range(15):
                        temp.append(float(p5[16*k + p])) 
                    ImageTemp1.append(temp)
                    
                p6 = linestemp[8*n + 7].split()
                ImageTemp2 = []
                for k in range(15): #Intensity Image 1 data
                    temp = []
                    for p in range(15):
                        temp.append(float(p6[16*k + p])) 
                    ImageTemp2.append(temp)
                    
                self.Loaded_CenterImages.append([ImageTemp1, ImageTemp2])
             
             
            self.Loaded_Intensity_BG_Fit = np.array(self.Loaded_Intensity_BG) * 0
            self.BGS_type = 'n/a'
            
            self.Loaded_Intensity_noBGS = np.array(self.Loaded_Intensity) + np.array(self.Loaded_Intensity_BG)
             


        elif (linestemp[0].split())[0] == '%%':
            
            self.Data_type_indicator = (linestemp[0].split())[0]
            
            self.BGS_type = (linestemp[0].split())[9]
            print 'self.BGS_type ', self.BGS_type
            
            print 'len(linestemp) ', len(linestemp)
            
            self.Nfeatures = int((len(linestemp)-1)/10)
            print 'self.Nfeatures ', self.Nfeatures
            
            self.Loaded_ImageJxpos = []
            self.Loaded_ImageJypos = []
            self.Loaded_A = []
            self.Loaded_T = []
            self.Loaded_phi = []
            self.Loaded_y0 = []
            self.Loaded_M = []
            
            
            self.Loaded_Intensity_noBGS = []
            self.Loaded_Intensity = []
            self.Loaded_Intensity_Fit = []
            self.Loaded_Intensity_BG = []
            self.Loaded_Intensity_BG_Fit = []
            self.Loaded_CenterImages = []
            
            
            
            for n in range(self.Nfeatures):
                p1 = linestemp[10*n + 2].split()
                #print 'p1 ', p1
                
                self.Loaded_ImageJxpos.append(int(float(p1[0])) ) 
                self.Loaded_ImageJypos.append(int(float(p1[1])) )
                self.Loaded_A.append((float(p1[2])) )
                self.Loaded_T.append((float(p1[3])) )
                self.Loaded_phi.append((float(p1[4])) )
                self.Loaded_y0.append((float(p1[5])) )
                self.Loaded_M.append((float(p1[6])) )

                
                lineElementsTemp = linestemp[10*n + 3].split()
                RawintensityTemp = []
                for k in lineElementsTemp: #Raw Intensity data
                    RawintensityTemp.append(float(k)) 
                self.Loaded_Intensity_noBGS.append(np.array(RawintensityTemp))


                lineElementsTemp = linestemp[10*n + 4].split()
                intensityTemp = []
                for k in lineElementsTemp: #Intensity data
                    intensityTemp.append(float(k)) 
                self.Loaded_Intensity.append(np.array(intensityTemp))


                lineElementsTemp = linestemp[10*n + 5].split()
                intensityFitTemp = []
                for k in lineElementsTemp: #Intensity Fit data
                    intensityFitTemp.append(float(k)) 
                self.Loaded_Intensity_Fit.append(np.array(intensityFitTemp))
                
 
                lineElementsTemp = linestemp[10*n + 6].split()
                intensityFitTemp = []
                intensityBGTemp = []
                for k in lineElementsTemp: #Intensity BG data
                    intensityBGTemp.append(float(k)) 
                self.Loaded_Intensity_BG.append(np.array(intensityBGTemp))
       
 
                lineElementsTemp = linestemp[10*n + 7].split()
                intensityFitTemp = []
                intensityBGfitTemp = []
                for k in lineElementsTemp: #Intensity BG data
                    intensityBGfitTemp.append(float(k)) 
                self.Loaded_Intensity_BG_Fit.append(np.array(intensityBGfitTemp))
       


               
               
               
                lineElementsTemp = linestemp[10*n + 8].split()
                ImageTemp1 = []
                for k in range(15): #Intensity Image 1 data
                    temp = []
                    for p in range(15):
                        temp.append(float(lineElementsTemp[16*k + p])) 
                    ImageTemp1.append(temp)
                    
                lineElementsTemp = linestemp[10*n + 9].split()
                ImageTemp2 = []
                for k in range(15): #Intensity Image 1 data
                    temp = []
                    for p in range(15):
                        temp.append(float(lineElementsTemp[16*k + p])) 
                    ImageTemp2.append(temp)
                    
                self.Loaded_CenterImages.append([ImageTemp1, ImageTemp2])
             
             


          
        else:
            MessageBox(self, 'Data type is not supported', 'error')

        






        self.Loaded_IntensityNoiseSTD = []
        self.Loaded_SBR = []

        for n in range(len(self.Loaded_Intensity)):
            STDtemp = np.std( np.array(self.Loaded_Intensity[n]) - np.array(self.Loaded_Intensity_Fit[n]) )
            #print n , ' STDtemp ', ' = ' ,STDtemp
            self.Loaded_IntensityNoiseSTD.append(STDtemp)
            
            SBRtemp = np.mean(self.Loaded_Intensity_noBGS[n])/np.min(self.Loaded_Intensity_BG_Fit[n])
            
            print 'SBRtemp ', SBRtemp
            self.Loaded_SBR.append(SBRtemp)

        self.SNR = np.array(self.Loaded_y0)/np.array(self.Loaded_IntensityNoiseSTD)
        


        print 'self.Loaded_SBR ', self.Loaded_SBR
        

  

          
        self.M_error = []
     
        IntensityTrajectory_fit_x = np.arange(len(self.Loaded_Intensity[0]))         
        for n in range(len(self.Loaded_Intensity)):
            A_guess = self.Loaded_y0[n]
            M_guess = self.Loaded_M[n]
            T_in_frame = self.Loaded_T[n]
            #phi_guess = self.Loaded_phi[n]
            
            initial_guess = [A_guess, M_guess, T_in_frame, 0]  # A, M, T, phi
            
            
            try:
                popt, pcov = curve_fit(ModulationDepth_M_func, IntensityTrajectory_fit_x, self.Loaded_Intensity[n], p0=initial_guess )
                #fit_temp = ModulationDepth_M_func(IntensityTrajectory_fit_x, popt[0],popt[1],popt[2],popt[3])
                
                err_temp = np.sqrt(np.diag(pcov))

                self.M_error.append(err_temp[1])
                #print '\nA_guess, M_guess, T_in_frame, phi_guess ', A_guess, M_guess, T_in_frame, phi_guess
                #print 'popt[0],popt[1],popt[2],popt[3] ', popt[0],popt[1],popt[2],popt[3]
                #print 'np.sqrt(np.diag(pcov)) ', np.sqrt(np.diag(pcov))

            except: 
                #print '\n fit error '
                self.M_error.append(0)
                pass
     

        #print 'self.M_error ', self.M_error

                      
                                        
        self.Loaded_CenterImages = np.array(self.Loaded_CenterImages)                    
                            
        
        









     
class DataAnalysisMainClass(wx.Frame):
    #ImageFilePath = ''
 
    #----------------------------------------------------------------------
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, "M Data Analysis", size=(1100, 1000), pos=(0,0))
        
        self.MainPanel = wx.ScrolledWindow(self, -1)
        self.MainPanel.SetScrollbars(1, 1, 800, 900)
        self.MainPanel.Show()

        
        
        
        self.btn_SwitchToOneChSM1 = wx.Button(self.MainPanel, pos=(250,5), label="Switch To Anisotropy M")
        self.btn_SwitchToOneChSM1.Bind(wx.EVT_BUTTON, self.SwitchToOneChSM1)

        
        self.btn_SwitchToDataCombining = wx.Button(self.MainPanel, pos=(450,5), label="Switch To Data Combining")
        self.btn_SwitchToDataCombining.Bind(wx.EVT_BUTTON, self.SwitchToDataCombining)
        
        
        btn_reset = wx.Button(self.MainPanel, pos=(800,5), label="Restarting Current GUI")
        btn_reset.Bind(wx.EVT_BUTTON, self.resetting)
    
    
    

        btn1 = wx.Button(self.MainPanel, pos=(30,40), label="Open #1 Data File")
        self.DataFilePath1 = btn1.Bind(wx.EVT_BUTTON, self.onOpenDataFile1)
        self.btn1text = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath1), pos=(150, 43))
        self.radio_file_1 = wx.RadioBox(self.MainPanel, -1, choices=['Ext', 'Emi', 'na' ], label='Type ', pos=(150, 63), size=wx.Size(130, 40), style=wx.RA_SPECIFY_COLS)
        self.radio_file_1.SetSelection(2)
        self.btn1text_2 = wx.StaticText(self.MainPanel, -1, '-', pos=(300, 83))


        self.btn1_plot = wx.Button(self.MainPanel, pos=(500,80), label="Plot M Histogram #1")
        self.btn1_plot.Bind(wx.EVT_BUTTON, self.Plot_Data1_Histogram)
        self.btn1_plot_label = wx.TextCtrl(self.MainPanel, -1, '', pos=(750, 80))

        self.btn1_plot_save_xyPlots = wx.Button(self.MainPanel, pos=(640,80), label="Save xy Plots")
        self.btn1_plot_save_xyPlots.Bind(wx.EVT_BUTTON, self.Plot_Data1_Histogram_save_xyPlots)
        



        
        
        btn2 = wx.Button(self.MainPanel, pos=(30,120), label="Open #2 Data File")
        self.DataFilePath2 = btn2.Bind(wx.EVT_BUTTON, self.onOpenDataFile2)
        self.btn2text = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath2), pos=(150, 123))
        self.radio_file_2 = wx.RadioBox(self.MainPanel, -1, choices=['Ext', 'Emi', 'na' ], label='Type ', pos=(150, 143), size=wx.Size(130, 40), style=wx.RA_SPECIFY_COLS)
        self.radio_file_2.SetSelection(2)
        self.btn2text_2 = wx.StaticText(self.MainPanel, -1, '-', pos=(300, 163))
        

        self.btn2_plot = wx.Button(self.MainPanel, pos=(500,160), label="Plot M Histogram #2")
        self.btn2_plot.Bind(wx.EVT_BUTTON, self.Plot_Data2_Histogram)
        self.btn2_plot_label = wx.TextCtrl(self.MainPanel, -1, '', pos=(750, 163))

        

        self.btn23_plot = wx.Button(self.MainPanel, pos=(740,120), label="Plot M Histograms Data 1,2")
        self.btn23_plot.Bind(wx.EVT_BUTTON, self.Plot_Data_2_3_Histogram)
        
        
        self.cb_publication_figures_yn = wx.CheckBox(self.MainPanel, -1, 'Pubication figures', (920, 125))
        self.cb_publication_figures_yn.SetValue(False)
        
        self.publication_featureN_text = wx.StaticText(self.MainPanel, -1, 'feature #:', pos=(920, 146))
        self.publication_featureN = wx.TextCtrl(self.MainPanel, -1, '118', pos=(980,143), size=(40,-1))

        
 

        
        
        btn3 = wx.Button(self.MainPanel, pos=(300,220), label="Open #3 Data File")
        self.DataFilePath3 = btn3.Bind(wx.EVT_BUTTON, self.onOpenDataFile3)
        self.btn3text = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath3), pos=(420, 223))
        self.cb_include_Data3_trajectories_yn = wx.CheckBox(self.MainPanel, -1, 'Include #3 Trajectories between #1 and #2', (50, 225))
        self.cb_include_Data3_trajectories_yn.SetValue(False)
        self.cb_include_Data3_fit_yn = wx.CheckBox(self.MainPanel, -1, 'Include #3 Data Fits', (90, 245))
        self.cb_include_Data3_fit_yn.SetValue(False)
        self.btn3text_2 = wx.StaticText(self.MainPanel, -1, '-', pos=(420, 243))
                
        
        
        btn4 = wx.Button(self.MainPanel, pos=(300,270), label="Open #4 Data File")
        self.DataFilePath4 = btn4.Bind(wx.EVT_BUTTON, self.onOpenDataFile4)
        self.btn4text = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath4), pos=(420, 273))
        self.cb_include_Data4_trajectories_yn = wx.CheckBox(self.MainPanel, -1, 'Include #4 Trajectories after #2', (50, 275))
        self.cb_include_Data4_trajectories_yn.SetValue(False)
        self.cb_include_Data4_fit_yn = wx.CheckBox(self.MainPanel, -1, 'Include #4 Data Fits', (90, 295))
        self.cb_include_Data4_fit_yn.SetValue(False)        
        self.btn4text_2 = wx.StaticText(self.MainPanel, -1, '-', pos=(420, 293))
                
        
        

        
        
        btn5 = wx.Button(self.MainPanel, pos=(300,340), label="Open #5 Data File")
        self.DataFilePath5 = btn5.Bind(wx.EVT_BUTTON, self.onOpenDataFile5)
        self.btn5text = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath5), pos=(420, 343))
        self.btn5text_2 = wx.StaticText(self.MainPanel, -1, '-', pos=(420, 363))
                
        
    
        btn6 = wx.Button(self.MainPanel, pos=(300,390), label="Open #6 Data File")
        self.DataFilePath6 = btn6.Bind(wx.EVT_BUTTON, self.onOpenDataFile6)
        self.btn6text = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath6), pos=(420, 393))
        self.btn6text_2 = wx.StaticText(self.MainPanel, -1, '-', pos=(420, 413))
                
        
        
                

        
        
        
        
        
        btn_plotDataTogether = wx.Button(self.MainPanel, pos=(20,450), label="Plot #1 #2 Together")
        self.DataFilePath2 = btn_plotDataTogether.Bind(wx.EVT_BUTTON, self.plotData)
        self.btn_plotData_text1 = wx.StaticText(self.MainPanel, -1, str('Total # '), pos=(15, 480))
        self.btn_plotData_text2 = wx.StaticText(self.MainPanel, -1, ' =>  # ', pos=(15, 500))




        



        self.cb_2nd_Plot_yn = wx.CheckBox(self.MainPanel, -1, '2nd Plots', (920, 455))
        self.cb_2nd_Plot_yn.SetValue(False)

        
        self.cb_90deg_offset_yn = wx.CheckBox(self.MainPanel, -1, '90 deg shift adjustment in phi', (150, 455))
        self.cb_90deg_offset_yn.SetValue(False)
        self.cb_fit_1_yn = wx.CheckBox(self.MainPanel, -1, 'Include Fit #1', (150, 475))
        self.cb_fit_1_yn.SetValue(True)
        self.cb_fit_2_yn = wx.CheckBox(self.MainPanel, -1, 'Include Fit #2', (150, 495))
        self.cb_fit_2_yn.SetValue(True)
        self.cb_individual_fig_yn = wx.CheckBox(self.MainPanel, -1, 'Individual feature\'s data', (150, 515))
        self.cb_individual_fig_yn.SetValue(False)
        


        
        
        
        
        self.cb_Data34_Trajectories_fig_yn = wx.CheckBox(self.MainPanel, -1, 'Separate #3,#4 Data Trajectories', (150, 565))
        self.cb_Data34_Trajectories_fig_yn.SetValue(False)
               
        
        
        
        self.cb_Data3_Trajectories_fig_yn = wx.CheckBox(self.MainPanel, -1, 'Separate #3 Data Trajectories', (150, 585))
        self.cb_Data3_Trajectories_fig_yn.SetValue(False)
        
        self.cb_Data4_Trajectories_fig_yn = wx.CheckBox(self.MainPanel, -1, 'Separate #4 Data Trajectories', (150, 605))
        self.cb_Data4_Trajectories_fig_yn.SetValue(False)
        
        
              
              
              
              

        self.cb_M1Range_yn = wx.CheckBox(self.MainPanel, -1, 'M1 Range', (370, 455))
        self.cb_M1Range_yn.SetValue(False)
        self.M1RangeValueMin = wx.TextCtrl(self.MainPanel, -1, '0.0', pos=(510,452), size=(40,-1))
        self.M1RangeValueMax = wx.TextCtrl(self.MainPanel, -1, '0.8', pos=(560,452), size=(40,-1))
        #self.M1Range = wx.SpinButton(self.MainPanel, -1, pos=(400, 232), size=(50,-1))
        #self.M1Range.SetRange(0, 1.1)  
        #self.M1Range.SetValue(0.7)   
        self.cb_M2Range_yn = wx.CheckBox(self.MainPanel, -1, 'M2 Range', (370, 475))
        self.cb_M2Range_yn.SetValue(False)
        self.M2RangeValueMin = wx.TextCtrl(self.MainPanel, -1, '0.0',  pos=(510,472), size=(40,-1))
        self.M2RangeValueMax = wx.TextCtrl(self.MainPanel, -1, '0.8',  pos=(560,472), size=(40,-1))
        
        
        
        self.cb_Intensity1Range_yn = wx.CheckBox(self.MainPanel, -1, 'Ave Intensity 1 Range', (370, 505))
        self.cb_Intensity1Range_yn.SetValue(False)
        self.Intensity1RangeValueMin = wx.TextCtrl(self.MainPanel, -1, '1000', pos=(510,502), size=(40,-1))
        self.Intensity1RangeValueMax = wx.TextCtrl(self.MainPanel, -1, '5000', pos=(560,502), size=(40,-1))
        self.cb_Intensity2Range_yn = wx.CheckBox(self.MainPanel, -1, 'Ave Intensity 2 Range', (370, 525))
        self.cb_Intensity2Range_yn.SetValue(False)
        self.Intensity2RangeValueMin = wx.TextCtrl(self.MainPanel, -1, '1000',  pos=(510,522), size=(40,-1))
        self.Intensity2RangeValueMax = wx.TextCtrl(self.MainPanel, -1, '5000',  pos=(560,522), size=(40,-1))
        
        
        
        
        self.cb_Fit_y0_1_Range_yn = wx.CheckBox(self.MainPanel, -1, 'Fit y0 1 Range', (640, 505))
        self.cb_Fit_y0_1_Range_yn.SetValue(False)
        self.Fit_y0_1RangeValueMin = wx.TextCtrl(self.MainPanel, -1, '1000', pos=(780,502), size=(40,-1))
        self.Fit_y0_1RangeValueMax = wx.TextCtrl(self.MainPanel, -1, '5000', pos=(830,502), size=(40,-1))
        
        self.cb_Fit_y0_2_Range_yn = wx.CheckBox(self.MainPanel, -1, 'Fit y0 2 Range', (640, 525))
        self.cb_Fit_y0_2_Range_yn.SetValue(False)
        self.Fit_y0_2RangeValueMin = wx.TextCtrl(self.MainPanel, -1, '1000',  pos=(780,522), size=(40,-1))
        self.Fit_y0_2RangeValueMax = wx.TextCtrl(self.MainPanel, -1, '5000',  pos=(830,522), size=(40,-1))
        
                
        
        
        
        self.cb_SNR_1Range_yn = wx.CheckBox(self.MainPanel, -1, 'SNR 1 Range', (370, 555))
        self.cb_SNR_1Range_yn.SetValue(False)
        self.SNR_1RangeValueMin = wx.TextCtrl(self.MainPanel, -1, '0', pos=(510,552), size=(40,-1))
        self.SNR_1RangeValueMax = wx.TextCtrl(self.MainPanel, -1, '9999', pos=(560,552), size=(40,-1))
        self.cb_SNR_2Range_yn = wx.CheckBox(self.MainPanel, -1, 'SNR 2 Range', (370, 575))
        self.cb_SNR_2Range_yn.SetValue(False)
        self.SNR_2RangeValueMin = wx.TextCtrl(self.MainPanel, -1, '0',  pos=(510,572), size=(40,-1))
        self.SNR_2RangeValueMax = wx.TextCtrl(self.MainPanel, -1, '9999',  pos=(560,572), size=(40,-1))
        
     
        
        self.cb_SBR_1Range_yn = wx.CheckBox(self.MainPanel, -1, 'SBR 1 Range', (370, 605))
        self.cb_SBR_1Range_yn.SetValue(False)
        self.SBR_1RangeValueMin = wx.TextCtrl(self.MainPanel, -1, '1.1', pos=(510,602), size=(40,-1))
        self.SBR_1RangeValueMax = wx.TextCtrl(self.MainPanel, -1, '9999', pos=(560,602), size=(40,-1))
        self.cb_SBR_2Range_yn = wx.CheckBox(self.MainPanel, -1, 'SBR 2 Range', (370, 625))
        self.cb_SBR_2Range_yn.SetValue(False)
        self.SBR_2RangeValueMin = wx.TextCtrl(self.MainPanel, -1, '1.1',  pos=(510,622), size=(40,-1))
        self.SBR_2RangeValueMax = wx.TextCtrl(self.MainPanel, -1, '9999',  pos=(560,622), size=(40,-1))
        
        






        
        btn_plotDataSeparate = wx.Button(self.MainPanel, pos=(20,650), label="Plot #1 #2 Separately")
        self.DataFilePath2 = btn_plotDataSeparate.Bind(wx.EVT_BUTTON, self.plotDataSeparate)
        self.btn_plotDataSeparate_text1 = wx.StaticText(self.MainPanel, -1, str('1. # '), pos=(15, 680))
        self.btn_plotDataSeparate_text2 = wx.StaticText(self.MainPanel, -1, str('2. # '), pos=(15, 700))


        
        
        
        
        wx.StaticText(self.MainPanel, -1, "File Name", pos=(50, 750))
        self.FileName = wx.TextCtrl(self.MainPanel, -1, ' ', pos=(50, 770), size=(500, -1))
        
        self.btn_SaveData = wx.Button(self.MainPanel, pos=(50,800), label="Save Figure")
        self.SaveData = self.btn_SaveData.Bind(wx.EVT_BUTTON, self.SaveData)
        self.btn_SaveData.Hide()
        
        self.SaveData_text = wx.StaticText(self.MainPanel, -1, str(' '), pos=(170, 803))





        self.btn_SelectData = wx.Button(self.MainPanel, pos=(50,900), label="Select Data")
        self.btn_SelectData.Bind(wx.EVT_BUTTON, self.SelectData)
        self.btn_SelectData.Hide()
        
        
        

        #MessageBox(self, 'Opening Data Analysis Window', 'Data Analysis')





















###############################################################################
#``````````````````````````````````````````````````````````````````````````````


    def SwitchToOneChSM1(self, event):
        
        Data1CHSM1Main = SingleChannelMainTracking()
        Data1CHSM1Main.Show()
        self.Close()


    def SwitchToDataCombining(self, event):
        
        DataCombiningMain = DataCombiningMainClass()
        DataCombiningMain.Show()
        self.Close()


     
     
    def resetting(self, event):
        
        print '\n\nResetting'
        try:
            plt.close()
            plt.close()
            plt.close()
            plt.close()
        except: pass
                
        try:
            for n in range(len(self.fig_DataPlots)):
                print 'fig_ExcitationM_Trajectory n = ', n
                plt.close()
        except: pass
        
        try:
            for n in range(len(self.fig_AllDataCenterImages_Loaded)):
                print 'fig_AllDataCenterImages n = ', n
                plt.close()
        except: pass
        
        
        try:
            self.SelectData_frame.Close()
        except: pass
            
            
        try:
            plt.close(self.fig_DataPlots_Selected)
        except: pass
 

      
        self.Close()                      
        frame2 = DataAnalysisMainClass()
        frame2.Show()
          
       
       
       


    def onOpenDataFile1(self, event):
        print  '\n clicked open file 1'
        
        self.data_1 = OpenDataFile()
        self.btn1text.SetLabel("#1: " + str(self.data_1.fileName))
        self.btn_SaveData.Hide()
        self.SaveData_text.SetLabel(' ')
        
        self.FileName.SetLabel(self.data_1.fileName[:40])
        print 'self.data_1 ', self.data_1
        
        self.btn1text_2.SetLabel('BGS type: ' + self.data_1.BGS_type + ' ,   Total: ' + str(len(self.data_1.Loaded_A)))
        
        self.btn_plotData_text1.SetLabel( 'Total # ' + str(len(self.data_1.Loaded_A)))



    def onOpenDataFile2(self, event):
        print  '\n clicked open file 2'
        
        self.data_2 = OpenDataFile()
        self.btn2text.SetLabel("#2: " + str(self.data_2.fileName))
        self.btn_SaveData.Hide()
        self.SaveData_text.SetLabel(' ')
        
        print 'self.data_2 ', self.data_2
        
        
        if len(self.data_1.Loaded_ImageJxpos) != len(self.data_2.Loaded_ImageJxpos):
            MessageBox(self, '#1, #2 data numbers are different' , 'Error')
            
        if self.data_1.BGS_type != self.data_2.BGS_type:
            MessageBox(self, '#1, #2 BGS (background subtraction) types are different' , 'Error')
            
        self.btn2text_2.SetLabel('BGS type: ' + self.data_2.BGS_type + ' ,   Total: ' + str(len(self.data_2.Loaded_A)))


    def onOpenDataFile3(self, event):
        print  '\n clicked open file 3'
        
        self.data_3 = OpenDataFile()
        self.btn3text.SetLabel("#3: " + str(self.data_3.fileName))
        self.btn_SaveData.Hide()
        self.SaveData_text.SetLabel(' ')
        
        print 'self.data_3 ', self.data_3
        
        if len(self.data_1.Loaded_ImageJxpos) != len(self.data_3.Loaded_ImageJxpos):
            MessageBox(self, '#1, #3 data numbers are different' , 'Error')
            
        if self.data_1.BGS_type != self.data_3.BGS_type:
            MessageBox(self, '#1, #3 BGS (background subtraction) types are different' , 'Error')
            
        self.btn3text_2.SetLabel('BGS type: ' + self.data_3.BGS_type + ' ,   Total: ' + str(len(self.data_3.Loaded_A)))


    def onOpenDataFile4(self, event):
        print  '\n clicked open file 4'
        self.data_4 = OpenDataFile()
        self.btn4text.SetLabel("#2: " + str(self.data_4.fileName))
        self.btn_SaveData.Hide()
        self.SaveData_text.SetLabel(' ')
        
        print 'self.data_4 ', self.data_4

        
        if len(self.data_1.Loaded_ImageJxpos) != len(self.data_4.Loaded_ImageJxpos):
            MessageBox(self, '#1, #4 data numbers are different' , 'Error')
            
        if self.data_1.BGS_type != self.data_4.BGS_type:
            MessageBox(self, '#1, #4 BGS (background subtraction) types are different' , 'Error')
            
        self.btn4text_2.SetLabel('BGS type: ' + self.data_4.BGS_type + ' ,   Total: ' + str(len(self.data_4.Loaded_A)))



    def onOpenDataFile5(self, event):
        print  '\n clicked open file 5'
        self.data_5 = OpenDataFile()
        self.btn5text.SetLabel("#2: " + str(self.data_5.fileName))
        self.btn_SaveData.Hide()
        self.SaveData_text.SetLabel(' ')
        
        print 'self.data_5 ', self.data_5

        
        if len(self.data_1.Loaded_ImageJxpos) != len(self.data_5.Loaded_ImageJxpos):
            MessageBox(self, '#1, #5 data numbers are different' , 'Error')
            
        if self.data_1.BGS_type != self.data_5.BGS_type:
            MessageBox(self, '#1, #5 BGS (background subtraction) types are different' , 'Error')
                        
        self.btn5text_2.SetLabel('BGS type: ' + self.data_5.BGS_type + ' ,   Total: ' + str(len(self.data_5.Loaded_A)))






    def onOpenDataFile6(self, event):
        print  '\n clicked open file 6'
        self.data_6 = OpenDataFile()
        self.btn6text.SetLabel("#2: " + str(self.data_6.fileName))
        self.btn_SaveData.Hide()
        self.SaveData_text.SetLabel(' ')
        
        print 'self.data_6 ', self.data_6

        
        if len(self.data_1.Loaded_ImageJxpos) != len(self.data_6.Loaded_ImageJxpos):
            MessageBox(self, '#1, #6 data numbers are different' , 'Error')
            
        if self.data_1.BGS_type != self.data_6.BGS_type:
            MessageBox(self, '#1, #6 BGS (background subtraction) types are different' , 'Error')
            
            
        self.btn6text_2.SetLabel('BGS type: ' + self.data_6.BGS_type + ' ,   Total: ' + str(len(self.data_6.Loaded_A)))


            
      
      
    def Plot_Data1_Histogram(self, event):
        print ''
        
        try:
            self.fig_Data1_Histogram.clf()
        except: pass
        

        self.fig_Data1_Histogram = plt.figure('Data1_M_Histogram', figsize = (27,12))
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(1110,30,600, 1000)
        #self.fig_Data1_Histogram.subplots_adjust(top = 0.89, bottom = 0.07, left = 0.05, right = 0.98, wspace = 0.4, hspace = 0.4)
        self.fig_Data1_Histogram.text(0.02, 0.93, '#1: ' + self.data_1.fileName  + '\nTotal #: ' + str(len(self.data_1.Loaded_M))\
        + ',    BGS_type: ' + self.data_1.BGS_type , size=12, color="black", weight='normal')
        
   
        plt.subplot(2,1,1)
        #plt.title('Median: ' + str(np.around(np.median(self.data_1.Loaded_M),2)), size = 16)
        plt.title('Median M: ' + str(np.around(np.median(self.data_1.Loaded_M),2)) + ' , intensity: ' + str(int(np.around(np.median(self.data_1.Loaded_y0),0))), size = 14)
        plt.hist(self.data_1.Loaded_M, bins=np.arange(0.0, 1.4, 0.1), label = '#1' + self.btn1_plot_label.GetValue(), alpha=1.0, color='red')
        plt.xlabel('Modulation Depth (M)', size=11)
        plt.xlim(0.0, 1.2)
        plt.legend(fontsize=10, loc='upper left')
        
        plt.subplot(2,1,2)
        #plt.title('Median: ' + str(np.around(np.median(self.data_1.Loaded_M),2)), size = 16)
        plt.title('Median intensity: ' + str(int(np.around(np.median(self.data_1.Loaded_y0),0))), size = 14)
        plt.hist(self.data_1.Loaded_y0, bins=np.arange(0, 12000, 500), label = '#1' + self.btn1_plot_label.GetValue(), alpha=1.0, color='red')
        plt.xlabel('Intensity (a.u.)', size=11)
        #plt.xlim(0.0, 1.2)
        plt.legend(fontsize=10, loc='upper left')
        
   
   
   
   
   
   
   
   
   

        if self.cb_publication_figures_yn.GetValue():
      
            figtemp = plt.figure('Data1_M_Histogram_publication', figsize = (27,12), facecolor = 'white')
      
            font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 14}
            matplotlib.rc('font', **font)
            majortick_size = 8
            majortick_width = 3
            plt.rc('axes', linewidth = 3)
            #matplotlib.rcParams['lines.linewidth'] = 4
            matplotlib.rcParams['xtick.major.size'] = majortick_size
            matplotlib.rcParams['xtick.major.width'] = majortick_width
            matplotlib.rcParams['xtick.minor.size'] = 5
            matplotlib.rcParams['xtick.minor.width'] = 4
            matplotlib.rcParams['ytick.major.size'] = majortick_size
            matplotlib.rcParams['ytick.major.width'] = majortick_width
            matplotlib.rcParams['ytick.minor.size'] = 5
            matplotlib.rcParams['ytick.minor.width'] = 4
            matplotlib.rcParams.update({'font.size': 22, 'family' : 'normal', 'weight' : 'bold',  'size': 8})
       
    
                
            
            
            
            ax_temp = figtemp.add_subplot(111)
            mngr = plt.get_current_fig_manager()
            mngr.window.setGeometry(1110,30,600, 500)
            #self.fig_Data1_Histogram.subplots_adjust(top = 0.89, bottom = 0.07, left = 0.05, right = 0.98, wspace = 0.4, hspace = 0.4)
            self.fig_Data1_Histogram.text(0.02, 0.93, '#1: ' + self.data_1.fileName  + '\nTotal #: ' + str(len(self.data_1.Loaded_M))\
            + ',    BGS_type: ' + self.data_1.BGS_type , size=12, color="black", weight='normal')
            
       
            plt.subplot(1,1,1)
            plt.title('Median M: ' + str(np.around(np.median(self.data_1.Loaded_M),2)) + ' , intensity: ' + str(np.around(np.median(self.data_1.Loaded_y0),0)), size = 14)
            plt.hist(self.data_1.Loaded_M, bins=np.arange(0.0, 1.4, 0.1), label = '#1' + self.btn1_plot_label.GetValue(), alpha=1.0, color='c', linewidth=1.5)
            plt.xlabel('Modulation Depth (M)', size=22, fontweight='bold')
            plt.ylabel('Occurence', size=22, fontweight='bold')
            plt.xlim(0.0, 1.2)
            
            plt.yticks(np.arange(0, 50, 10))
            
            #ax_temp.xaxis.set_tick_params(width=2)
            #ax_temp.yaxis.set_tick_params(width=2)
            
            #plt.legend(fontsize=10, loc='upper left')
            
            plt.tight_layout()
            
        plt.ion()
        plt.show()
            

    def Plot_Data1_Histogram_save_xyPlots(self, event):
        print ''
        

        todayDate = time.strftime("%Y%m%d_%Hh%Mm")


        
        self.fig_Data1_Histogram.savefig(self.data_1.fileDirectory + '\\' + str(self.data_1.fileName[:40]) + '_' + todayDate + '_hist' +'.png')
        plt.close(self.fig_Data1_Histogram)


        #################################################################################################            
        ### hist to xy plots in files.


        

        bin_width = 0.1
        hist_x, hist_y = histtoxy2(self.data_1.Loaded_M, 0.0, 1.3, bin_width, 'Y')
        
        #hist_x.append(1.0 + bin_width/2)
        #hist_y.append(hist_y[-1])
        
        x = np.array(hist_x) - bin_width/2
        y = np.array(hist_y)/float(np.max(hist_y))
        
        plt.figure('h1')
        plt.plot(x, y, '-ro')
        plt.bar(x,y, width=bin_width, alpha=0.7)
        plt.show()
        
        ff = open(self.data_1.fileDirectory + '\\' + str(self.data_1.fileName[:40]) + 'M_hist_xyPlots_' + todayDate +'.txt','w')
        for i in range(len(x)):
            ff.write(str(x[i]) + ' ' + str(y[i]) + '\n'   )    
        ff.close()            

        


        bin_width = 500
        hist_x, hist_y = histtoxy2(self.data_1.Loaded_y0, 0, 16000, bin_width, 'Y')
        
        #hist_x.append(1.0 + bin_width/2)
        #hist_y.append(hist_y[-1])
        
        x = np.array(hist_x) - bin_width/2
        y = np.array(hist_y)/float(np.max(hist_y))
        
        plt.figure('h2')
        plt.plot(x, y, '-ro')
        plt.bar(x,y, width=bin_width, alpha=0.7)
        plt.show()
        
        ff = open(self.data_1.fileDirectory + '\\' + str(self.data_1.fileName[:40]) + 'intensity_hist_xyPlots_' + todayDate +'.txt','w')
        for i in range(len(x)):
            ff.write(str(x[i]) + ' ' + str(y[i]) + '\n'   )    
        ff.close()            

        
        
        
        
        print 'hist to xy plots in files created.'
        
        ### hist to xy plots in files.            
        #################################################################################################
            
        
    
    
    
    


            
            
     
      
    def Plot_Data2_Histogram(self, event):
        print ''
        
        try:
            self.fig_Data2_Histogram.clf()
            
            self.fig_Data2_Histogram2.clf()
            
        except: pass
        
        '''

        self.fig_Data2_Histogram = plt.figure('Data2_M_Histogram', figsize = (27,12))
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(1110,530,600, 500)
        #self.fig_Data1_Histogram.subplots_adjust(top = 0.89, bottom = 0.07, left = 0.05, right = 0.98, wspace = 0.4, hspace = 0.4)
        self.fig_Data2_Histogram.text(0.02, 0.93, '#2: ' + self.data_2.fileName  + '\nTota2 #: ' + str(len(self.data_2.Loaded_M))\
        + ',    BGS_type: ' + self.data_1.BGS_type , size=12, color="black", weight='normal')
        
   
        plt.subplot(1,1,1)
        plt.title('Median: ' + str(np.around(np.median(self.data_2.Loaded_M),2)), size = 16)
        plt.hist(self.data_2.Loaded_M, bins=np.arange(0.0, 1.4, 0.1), label = '#2' + self.btn2_plot_label.GetValue(), alpha=1.0, color='blue')
        plt.xlabel('Modulation Depth (M)', size=11)
        plt.xlim(0.0, 1.2)
        plt.legend(fontsize=10, loc='upper left')
       
        '''
       
       
       
       
       

        self.fig_Data2_Histogram2 = plt.figure('Data2_M_Histogram2', figsize = (27,12))
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(1110,30,600, 1000)
        #self.fig_Data1_Histogram.subplots_adjust(top = 0.89, bottom = 0.07, left = 0.05, right = 0.98, wspace = 0.4, hspace = 0.4)
        self.fig_Data2_Histogram2.text(0.02, 0.93, '#1: ' + self.data_2.fileName  + '\nTotal #: ' + str(len(self.data_2.Loaded_M))\
        + ',    BGS_type: ' + self.data_2.BGS_type , size=12, color="black", weight='normal')
        
   
        plt.subplot(2,1,1)
        #plt.title('Median: ' + str(np.around(np.median(self.data_1.Loaded_M),2)), size = 16)
        plt.title('Median M: ' + str(np.around(np.median(self.data_2.Loaded_M),2)) + ' , intensity: ' + str(int(np.around(np.median(self.data_2.Loaded_y0),0))), size = 14)
        plt.hist(self.data_2.Loaded_M, bins=np.arange(0.0, 1.4, 0.1), label = '#2' + self.btn2_plot_label.GetValue(), alpha=1.0, color='blue')
        plt.xlabel('Modulation Depth (M)', size=11)
        plt.xlim(0.0, 1.2)
        plt.legend(fontsize=10, loc='upper left')
        
        plt.subplot(2,1,2)
        #plt.title('Median: ' + str(np.around(np.median(self.data_1.Loaded_M),2)), size = 16)
        plt.title('Median intensity: ' + str(int(np.around(np.median(self.data_2.Loaded_y0),0))), size = 14)
        plt.hist(self.data_2.Loaded_y0, bins=np.arange(0, 12000, 500), label = '#2' + self.btn2_plot_label.GetValue(), alpha=1.0, color='blue')
        plt.xlabel('Intensity (a.u.)', size=11)
        #plt.xlim(0.0, 1.2)
        plt.legend(fontsize=10, loc='upper left')
        
          
       
       
       
       
       
       
       
       
       
       
       
       
       
        
        
        plt.ion()
        plt.show()
        
        
        
      
    def Plot_Data_2_3_Histogram(self, event):
        print ''
        
   
        try:
            self.fig_Data23_Histogram.clf()
        except: pass
        
        

        self.fig_Data23_Histogram = plt.figure('Data_2,3_M_Histogram', figsize = (27,12), facecolor='white')
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(1110,230,600, 500)
        #self.fig_Data1_Histogram.subplots_adjust(top = 0.89, bottom = 0.07, left = 0.05, right = 0.98, wspace = 0.4, hspace = 0.4)
        
   
        plt.subplot(1,1,1)
        plt.hist(self.data_1.Loaded_M, bins=np.arange(0.0, 1.4, 0.1), label = '#1' + self.btn1_plot_label.GetValue(), alpha=.5, color='red')
        plt.hist(self.data_2.Loaded_M, bins=np.arange(0.0, 1.4, 0.1), label = '#2' + self.btn2_plot_label.GetValue(), alpha=.5, color='blue')
        plt.xlabel('Modulation Depth (M)', size=11)
        plt.xlim(0.0, 1.2)
        plt.yticks(np.arange(20, 100, 20))
        plt.legend(fontsize=22, loc='upper left',frameon=False)
        
        

        if self.cb_publication_figures_yn.GetValue():
            
        
            font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 14}
            matplotlib.rc('font', **font)
            majortick_size = 8
            majortick_width = 3
            plt.rc('axes', linewidth = 3)
            #matplotlib.rcParams['lines.linewidth'] = 4
            matplotlib.rcParams['xtick.major.size'] = majortick_size
            matplotlib.rcParams['xtick.major.width'] = majortick_width
            matplotlib.rcParams['xtick.minor.size'] = 5
            matplotlib.rcParams['xtick.minor.width'] = 4
            matplotlib.rcParams['ytick.major.size'] = majortick_size
            matplotlib.rcParams['ytick.major.width'] = majortick_width
            matplotlib.rcParams['ytick.minor.size'] = 5
            matplotlib.rcParams['ytick.minor.width'] = 4
            matplotlib.rcParams.update({'font.size': 22, 'family' : 'normal', 'weight' : 'bold',  'size': 8})
    
            
            fig_temp1 = plt.figure(facecolor='white')
            
            plt.rcParams.update({'mathtext.default': 'regular' })
            
            #plt.title('Excitation M')

            ax1 = fig_temp1.add_subplot(111)
            ax2 = ax1.twinx()
            
            
            
            #plt.tick_params(left='on',  top='off', right='off', bottom='on',labelleft='on', labeltop='off', labelright='on', labelbottom='on')
                                   
            ax1.hist(self.data_1.Loaded_M, bins=np.arange(0.0, 1.4, 0.1), label = '' + self.btn1_plot_label.GetValue(), alpha=.5, color='c', linewidth=1.5)
            ax2.hist(self.data_2.Loaded_M, bins=np.arange(0.0, 1.4, 0.1), label = '' + self.btn2_plot_label.GetValue(), alpha=.5, color='m', linewidth=1.5)
            
            #plt.xlabel('Modulation Depth (M)', size=11)
            
            plt.xlim(0.0, 1.2)
            
            ax1.legend(fontsize=20, loc='upper left', frameon=False)
            ax1.set_ylim(0, 82)
            
            
            ax1.tick_params(axis='x', labelsize=20, pad=10)
            ax1.tick_params(axis='y', labelsize=20, pad=5)
            ax1.yaxis.set_ticks([20, 40, 60, 80])
            
            
            ax2.set_ylim(0, 139)
            ax2.tick_params(axis='y', labelsize=20, right='on', labelright='on')
            ax2.yaxis.set_ticks([30, 60, 90, 120])
            ax2.legend(fontsize=20, loc='upper right', frameon=False)

              
            ax1.set_xlabel('Modulation Depth (M)', size=22, fontweight='bold')
            ax1.set_ylabel('Occurence', color='k', size=22, fontweight='bold')


            ax2.set_ylabel('', color='b', size = 22)
            

            
            '''
            plt.figure(facecolor='white')


            plt.tick_params(left='on',  top='off', right='off', bottom='on',
                                   labelleft='on', labeltop='off', labelright='off', labelbottom='on')
                                   
            plt.hist(self.data_1.Loaded_M, bins=np.arange(0.0, 1.4, 0.1), label = '#1' + self.btn1_plot_label.GetValue(), alpha=.5, color='red')
            plt.hist(self.data_2.Loaded_M, bins=np.arange(0.0, 1.4, 0.1), label = '#2' + self.btn2_plot_label.GetValue(), alpha=.5, color='blue')
            plt.xlabel('Modulation Depth (M)', size=11)
            plt.xlim(0.0, 1.2)
            plt.legend(fontsize=10, loc='upper left')
            '''
            
            
            
            
            
            plt.tight_layout()
            

        else:
            font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   : 12}
            matplotlib.rc('font', **font)
            majortick_size = 3
            majortick_width = 1
            plt.rc('axes', linewidth = 1)
            #matplotlib.rcParams['lines.linewidth'] = 4
            matplotlib.rcParams['xtick.major.size'] = majortick_size
            matplotlib.rcParams['xtick.major.width'] = majortick_width
            matplotlib.rcParams['xtick.minor.size'] = 2
            matplotlib.rcParams['xtick.minor.width'] = 1
            matplotlib.rcParams['ytick.major.size'] = majortick_size
            matplotlib.rcParams['ytick.major.width'] = majortick_width
            matplotlib.rcParams['ytick.minor.size'] = 2
            matplotlib.rcParams['ytick.minor.width'] = 1
            matplotlib.rcParams.update({'font.size': 12, 'family' : 'normal', 'weight' : 'normal',  'size': 8})
    
            
        plt.ion()
        plt.show()        
        
      
      

    def plotData(self, event):
        
        print ' '
   
   
        try:
            plt.close(self.fig_temp2)
        except: pass


        try:
            self.data_1.Loaded_ImageJxpos
        except:
            MessageBox(self, 'No #1 Data' , 'Error')
            return

        try:
            self.data_2.Loaded_ImageJxpos
        except:
            MessageBox(self, 'No #2 Data' , 'Error')
            return

       
            
            
            
            


        if len(self.data_1.Loaded_ImageJxpos) != len(self.data_2.Loaded_ImageJxpos):
            MessageBox(self, '# of Data is not matching\nCheck data', 'Error')
            return
            





        if self.cb_Data3_Trajectories_fig_yn.GetValue():
            try:
                self.data_3.Loaded_ImageJxpos
            except:
                MessageBox(self, 'No #3 Data' , 'Error')
                return
    
           


        if self.cb_Data4_Trajectories_fig_yn.GetValue():
            try:
                self.data_4.Loaded_ImageJxpos
            except:
                MessageBox(self, 'No #4 Data' , 'Error')
                return
    
           


























        self.btn_SaveData.Show()
        self.btn_SelectData.Show()        
        self.SaveData_text.SetLabel(' ')
        
        
        try:
            self.SelectData_frame.Destroy()
        except:pass
    

        try:
            plt.close(self.fig_DataPlots)
            plt.close(self.fig_DataPlots_2)
        except: pass
    
        
        try:
            for n in range(len(self.fig_AllDataCenterImages_Loaded)):
                plt.close(self.fig_AllDataCenterImages_Loaded[n])
        except: pass
        
        
        
        try:
            for n in range(len(self.fig_Data3_Trajectories_Loaded)):
                plt.close(self.fig_Data3_Trajectories_Loaded[n])
        except: pass     
        
        
        
        try:
            for n in range(len(self.fig_Data4_Trajectories_Loaded)):
                plt.close(self.fig_Data4_Trajectories_Loaded[n])
        except: pass     
        
               
        
        
        
        if self.radio_file_1.GetSelection() == 0:
            fileType_1 = 'Ext'
        elif self.radio_file_1.GetSelection() == 1:
            fileType_1 = 'Emi'
        elif self.radio_file_1.GetSelection() == 2:
            fileType_1 = ' '
            
        
        if self.radio_file_2.GetSelection() == 0:
            fileType_2 = 'Ext'
        elif self.radio_file_2.GetSelection() == 1:
            fileType_2 = 'Emi'
        elif self.radio_file_2.GetSelection() == 2:
            fileType_2 = ' '
            
            
        self.BGS_type_1 = self.data_1.BGS_type
        self.BGS_type_2 = self.data_2.BGS_type
            
            
            
        self.Loaded_ImageJxpos_1_Range = []    
        self.Loaded_ImageJypos_1_Range = []
        self.Loaded_ImageJxpos_2_Range = []
        self.Loaded_ImageJypos_2_Range = []
        self.Loaded_M_1_Range = []
        self.Loaded_M_2_Range = []
        self.Loaded_A_1_Range = []
        self.Loaded_A_2_Range = []
        self.Loaded_T_1_Range = []
        self.Loaded_T_2_Range = []
        self.Loaded_phi_1_Range = []
        self.Loaded_phi_2_Range = []
        self.Loaded_y0_1_Range = []
        self.Loaded_y0_2_Range = []
        
        self.Loaded_Intensity_1_Range = []
        self.Loaded_Intensity_2_Range = []
        self.Loaded_Intensity_Fit_1_Range = []
        self.Loaded_Intensity_Fit_2_Range = []

        self.Loaded_Intensity_noBGS_1_Range = []
        self.Loaded_Intensity_noBGS_2_Range = []
        self.Loaded_Intensity_BG_1_Range = []
        self.Loaded_Intensity_BG_2_Range = []
        self.Loaded_Intensity_BG_Fit_1_Range = []
        self.Loaded_Intensity_BG_Fit_2_Range = []
        
        self.Loaded_CenterImages_1_Range = []
        self.Loaded_CenterImages_2_Range = []


        self.M_error_1_Range = []
        self.M_error_2_Range = []
        self.SNR_1_Range = []
        self.SNR_2_Range = []
        self.SBR_1_Range = []
        self.SBR_2_Range = []
        
        
        
        self.M_error_3_Range = []
        self.M_error_4_Range = []
        self.SNR_3_Range = []
        self.SNR_4_Range = []
        self.SBR_3_Range = []
        self.SBR_4_Range = []
        



        
        
        self.Loaded_ImageJxpos_3_Range = []    
        self.Loaded_ImageJypos_3_Range = []
        self.Loaded_A_3_Range = []
        self.Loaded_T_3_Range = []
        self.Loaded_phi_3_Range = []
        self.Loaded_y0_3_Range = []
        self.Loaded_M_3_Range = []

        self.Loaded_Intensity_3_Range = []
        self.Loaded_Intensity_Fit_3_Range = []
        self.Loaded_Intensity_noBGS_3_Range = []
        self.Loaded_Intensity_BG_3_Range = []
        self.Loaded_Intensity_BG_Fit_3_Range = []
        self.Loaded_CenterImages_3_Range = []
        
        
        
        
        self.Loaded_ImageJxpos_4_Range = []    
        self.Loaded_ImageJypos_4_Range = []
        self.Loaded_A_4_Range = []
        self.Loaded_T_4_Range = []
        self.Loaded_phi_4_Range = []
        self.Loaded_y0_4_Range = []
        self.Loaded_M_4_Range = []

        self.Loaded_Intensity_4_Range = []
        self.Loaded_Intensity_Fit_4_Range = []
        self.Loaded_Intensity_noBGS_4_Range = []
        self.Loaded_Intensity_BG_4_Range = []
        self.Loaded_Intensity_BG_Fit_4_Range = []
        self.Loaded_CenterImages_4_Range = []
        
                
        
        self.Loaded_ImageJxpos_5_Range = []    
        self.Loaded_ImageJypos_5_Range = []
        self.Loaded_A_5_Range = []
        self.Loaded_T_5_Range = []
        self.Loaded_phi_5_Range = []
        self.Loaded_y0_5_Range = []
        self.Loaded_M_5_Range = []

        self.Loaded_Intensity_5_Range = []
        self.Loaded_Intensity_Fit_5_Range = []
        self.Loaded_Intensity_noBGS_5_Range = []
        self.Loaded_Intensity_BG_5_Range = []
        self.Loaded_Intensity_BG_Fit_5_Range = []
        self.Loaded_CenterImages_5_Range = []
        
                
        
        self.Loaded_ImageJxpos_6_Range = []    
        self.Loaded_ImageJypos_6_Range = []
        self.Loaded_A_6_Range = []
        self.Loaded_T_6_Range = []
        self.Loaded_phi_6_Range = []
        self.Loaded_y0_6_Range = []
        self.Loaded_M_6_Range = []

        self.Loaded_Intensity_6_Range = []
        self.Loaded_Intensity_Fit_6_Range = []
        self.Loaded_Intensity_noBGS_6_Range = []
        self.Loaded_Intensity_BG_6_Range = []
        self.Loaded_Intensity_BG_Fit_6_Range = []
        self.Loaded_CenterImages_6_Range = []
        
                
        


                   

        if self.cb_M1Range_yn.GetValue() == False and self.cb_M2Range_yn.GetValue() == False and self.cb_Intensity1Range_yn.GetValue() == False and self.cb_Intensity2Range_yn.GetValue() == False\
        and self.cb_Fit_y0_1_Range_yn.GetValue() == False and self.cb_Fit_y0_2_Range_yn.GetValue() == False and self.cb_SNR_1Range_yn.GetValue() == False and self.cb_SNR_2Range_yn.GetValue() == False \
        and self.cb_SBR_1Range_yn.GetValue() == False and self.cb_SBR_2Range_yn.GetValue() == False:
        
            self.Loaded_ImageJxpos_1_Range = self.data_1.Loaded_ImageJxpos
            self.Loaded_ImageJypos_1_Range = self.data_1.Loaded_ImageJypos
            self.Loaded_ImageJxpos_2_Range = self.data_2.Loaded_ImageJxpos
            self.Loaded_ImageJypos_2_Range = self.data_2.Loaded_ImageJypos
            self.Loaded_M_1_Range = self.data_1.Loaded_M
            self.Loaded_M_2_Range = self.data_2.Loaded_M
            self.Loaded_A_1_Range = self.data_1.Loaded_A
            self.Loaded_A_2_Range = self.data_2.Loaded_A
            self.Loaded_T_1_Range = self.data_1.Loaded_T
            self.Loaded_T_2_Range = self.data_2.Loaded_T
            self.Loaded_phi_1_Range = self.data_1.Loaded_phi
            self.Loaded_phi_2_Range = self.data_2.Loaded_phi
            self.Loaded_y0_1_Range = self.data_1.Loaded_y0
            self.Loaded_y0_2_Range = self.data_2.Loaded_y0
            self.Loaded_CenterImages_1_Range = self.data_1.Loaded_CenterImages
            self.Loaded_CenterImages_2_Range = self.data_2.Loaded_CenterImages
            self.Loaded_Intensity_1_Range = self.data_1.Loaded_Intensity
            self.Loaded_Intensity_2_Range = self.data_2.Loaded_Intensity
            self.Loaded_Intensity_noBGS_1_Range = self.data_1.Loaded_Intensity_noBGS
            self.Loaded_Intensity_noBGS_2_Range = self.data_2.Loaded_Intensity_noBGS
            
            
            self.Loaded_Intensity_BG_1_Range = self.data_1.Loaded_Intensity_BG
            self.Loaded_Intensity_BG_2_Range = self.data_2.Loaded_Intensity_BG
            self.Loaded_Intensity_BG_Fit_1_Range = self.data_1.Loaded_Intensity_BG_Fit
            self.Loaded_Intensity_BG_Fit_2_Range = self.data_2.Loaded_Intensity_BG_Fit
            
            self.Loaded_Intensity_Fit_1_Range = self.data_1.Loaded_Intensity_Fit
            self.Loaded_Intensity_Fit_2_Range = self.data_2.Loaded_Intensity_Fit
            self.M_error_1_Range = self.data_1.M_error
            self.M_error_2_Range = self.data_2.M_error
            self.SNR_1_Range = self.data_1.SNR
            self.SNR_2_Range = self.data_2.SNR
            self.SBR_1_Range = self.data_1.Loaded_SBR
            self.SBR_2_Range = self.data_2.Loaded_SBR


            
            try:
                self.Loaded_ImageJxpos_3_Range = self.data_3.Loaded_ImageJxpos
                self.Loaded_ImageJypos_3_Range = self.data_3.Loaded_ImageJypos
                self.Loaded_A_3_Range = self.data_3.Loaded_A
                self.Loaded_T_3_Range = self.data_3.Loaded_T
                self.Loaded_phi_3_Range = self.data_3.Loaded_phi
                self.Loaded_y0_3_Range = self.data_3.Loaded_y0
                self.Loaded_M_3_Range = self.data_3.Loaded_M
                
                self.Loaded_Intensity_noBGS_3_Range = self.data_3.Loaded_Intensity_noBGS
                self.Loaded_Intensity_3_Range = self.data_3.Loaded_Intensity
                self.Loaded_Intensity_Fit_3_Range = self.data_3.Loaded_Intensity_Fit
                self.Loaded_Intensity_BG_3_Range = self.data_3.Loaded_Intensity_BG
                self.Loaded_Intensity_BG_Fit_3_Range = self.data_3.Loaded_Intensity_BG_Fit
                self.Loaded_CenterImages_3_Range = self.data_3.Loaded_CenterImages
                
                self.M_error_3_Range = self.data_3.M_error
                self.SNR_3_Range = self.data_3.SNR
                self.SBR_3_Range = self.data_3.Loaded_SBR



            except: pass



            try:
                self.Loaded_ImageJxpos_4_Range = self.data_4.Loaded_ImageJxpos
                self.Loaded_ImageJypos_4_Range = self.data_4.Loaded_ImageJypos
                self.Loaded_A_4_Range = self.data_4.Loaded_A
                self.Loaded_T_4_Range = self.data_4.Loaded_T
                self.Loaded_phi_4_Range = self.data_4.Loaded_phi
                self.Loaded_y0_4_Range = self.data_4.Loaded_y0
                self.Loaded_M_4_Range = self.data_4.Loaded_M
                
                self.Loaded_Intensity_noBGS_4_Range = self.data_4.Loaded_Intensity_noBGS
                self.Loaded_Intensity_4_Range = self.data_4.Loaded_Intensity
                self.Loaded_Intensity_Fit_4_Range = self.data_4.Loaded_Intensity_Fit
                self.Loaded_CenterImages_4_Range = self.data_4.Loaded_CenterImages
                self.Loaded_Intensity_BG_4_Range = self.data_4.Loaded_Intensity_BG
                self.Loaded_Intensity_BG_Fit_4_Range = self.data_4.Loaded_Intensity_BG_Fit
                
                
                self.M_error_4_Range = self.data_4.M_error
                self.SNR_4_Range = self.data_4.SNR
                self.SBR_4_Range = self.data_4.Loaded_SBR




            except: pass
                



            try:
                self.Loaded_ImageJxpos_5_Range = self.data_5.Loaded_ImageJxpos
                self.Loaded_ImageJypos_5_Range = self.data_5.Loaded_ImageJypos
                self.Loaded_A_5_Range = self.data_5.Loaded_A
                self.Loaded_T_5_Range = self.data_5.Loaded_T
                self.Loaded_phi_5_Range = self.data_5.Loaded_phi
                self.Loaded_y0_5_Range = self.data_5.Loaded_y0
                self.Loaded_M_5_Range = self.data_5.Loaded_M
                
                self.Loaded_Intensity_noBGS_5_Range = self.data_5.Loaded_Intensity_noBGS
                self.Loaded_Intensity_5_Range = self.data_5.Loaded_Intensity
                self.Loaded_Intensity_Fit_5_Range = self.data_5.Loaded_Intensity_Fit
                self.Loaded_CenterImages_5_Range = self.data_5.Loaded_CenterImages
                self.Loaded_Intensity_BG_5_Range = self.data_5.Loaded_Intensity_BG
                self.Loaded_Intensity_BG_Fit_5_Range = self.data_5.Loaded_Intensity_BG_Fit
            except: pass
        
        
        

            try:
                self.Loaded_ImageJxpos_6_Range = self.data_6.Loaded_ImageJxpos
                self.Loaded_ImageJypos_6_Range = self.data_6.Loaded_ImageJypos
                self.Loaded_A_6_Range = self.data_6.Loaded_A
                self.Loaded_T_6_Range = self.data_6.Loaded_T
                self.Loaded_phi_6_Range = self.data_6.Loaded_phi
                self.Loaded_y0_6_Range = self.data_6.Loaded_y0
                self.Loaded_M_6_Range = self.data_6.Loaded_M
                
                self.Loaded_Intensity_noBGS_6_Range = self.data_6.Loaded_Intensity_noBGS
                self.Loaded_Intensity_6_Range = self.data_6.Loaded_Intensity
                self.Loaded_Intensity_Fit_6_Range = self.data_6.Loaded_Intensity_Fit
                self.Loaded_CenterImages_6_Range = self.data_6.Loaded_CenterImages
                self.Loaded_Intensity_BG_6_Range = self.data_6.Loaded_Intensity_BG
                self.Loaded_Intensity_BG_Fit_6_Range = self.data_6.Loaded_Intensity_BG_Fit
            except: pass
        
        
        
                    
                
                

                    
        else:
            for n in range(len(self.data_1.Loaded_M)):
                
                if self.cb_M1Range_yn.GetValue() == True:
                    if self.data_1.Loaded_M[n] >= float(self.M1RangeValueMin.GetValue()) and self.data_1.Loaded_M[n] <= float(self.M1RangeValueMax.GetValue()):
                        pass
                    else: continue
                    
                if self.cb_M2Range_yn.GetValue() == True:
                    if self.data_2.Loaded_M[n] >= float(self.M2RangeValueMin.GetValue()) and self.data_2.Loaded_M[n] <= float(self.M2RangeValueMax.GetValue()):
                        pass
                    else: continue
                


                if self.cb_Intensity1Range_yn.GetValue() == True:
                    if np.mean(self.data_1.Loaded_Intensity[n]) >= float(self.Intensity1RangeValueMin.GetValue()) and np.mean(self.data_1.Loaded_Intensity[n]) <= float(self.Intensity1RangeValueMax.GetValue()):
                        pass
                    else: continue
                    
                if self.cb_Intensity2Range_yn.GetValue() == True:
                    if np.mean(self.data_2.Loaded_Intensity[n]) >= float(self.Intensity2RangeValueMin.GetValue()) and np.mean(self.data_2.Loaded_Intensity[n]) <= float(self.Intensity2RangeValueMax.GetValue()):
                        pass
                    else: continue




                if self.cb_Fit_y0_1_Range_yn.GetValue() == True:
                    if self.data_1.Loaded_y0[n] >= float(self.Fit_y0_1RangeValueMin.GetValue()) and self.data_1.Loaded_y0[n] <= float(self.Fit_y0_1RangeValueMax.GetValue()):
                        pass
                    else: continue
                    



                if self.cb_Fit_y0_2_Range_yn.GetValue() == True:
                    if self.data_2.Loaded_y0[n] >= float(self.Fit_y0_2RangeValueMin.GetValue()) and self.data_2.Loaded_y0[n] <= float(self.Fit_y0_2RangeValueMax.GetValue()):
                        pass
                    else: continue
                    
  
        
        



                if self.cb_SNR_1Range_yn.GetValue() == True:
                    if (self.data_1.SNR[n]) >= float(self.SNR_1RangeValueMin.GetValue()) and (self.data_1.SNR[n]) <= float(self.SNR_1RangeValueMax.GetValue()):
                        pass
                    else: continue

                if self.cb_SNR_2Range_yn.GetValue() == True:
                    if (self.data_2.SNR[n]) >= float(self.SNR_2RangeValueMin.GetValue()) and (self.data_2.SNR[n]) <= float(self.SNR_2RangeValueMax.GetValue()):
                        pass
                    else: continue


                if self.cb_SBR_1Range_yn.GetValue() == True:
                    if (self.data_1.Loaded_SBR[n]) >= float(self.SBR_1RangeValueMin.GetValue()) and (self.data_1.Loaded_SBR[n]) <= float(self.SBR_1RangeValueMax.GetValue()):
                        pass
                    else: continue

                if self.cb_SBR_2Range_yn.GetValue() == True:
                    if (self.data_2.Loaded_SBR[n]) >= float(self.SBR_2RangeValueMin.GetValue()) and (self.data_2.Loaded_SBR[n]) <= float(self.SBR_2RangeValueMax.GetValue()):
                        pass
                    else: continue
               

               
                self.Loaded_ImageJxpos_1_Range.append(self.data_1.Loaded_ImageJxpos[n])
                self.Loaded_ImageJypos_1_Range.append(self.data_1.Loaded_ImageJypos[n])
                self.Loaded_ImageJxpos_2_Range.append(self.data_2.Loaded_ImageJxpos[n])
                self.Loaded_ImageJypos_2_Range.append(self.data_2.Loaded_ImageJypos[n])
                self.Loaded_M_1_Range.append(self.data_1.Loaded_M[n])
                self.Loaded_M_2_Range.append(self.data_2.Loaded_M[n])
                self.Loaded_A_1_Range.append(self.data_1.Loaded_A[n])
                self.Loaded_A_2_Range.append(self.data_2.Loaded_A[n])
                self.Loaded_T_1_Range.append(self.data_1.Loaded_T[n])
                self.Loaded_T_2_Range.append(self.data_2.Loaded_T[n])
                self.Loaded_phi_1_Range.append(self.data_1.Loaded_phi[n])
                self.Loaded_phi_2_Range.append(self.data_2.Loaded_phi[n])
                self.Loaded_y0_1_Range.append(self.data_1.Loaded_y0[n])
                self.Loaded_y0_2_Range.append(self.data_2.Loaded_y0[n])
                self.Loaded_CenterImages_1_Range.append(self.data_1.Loaded_CenterImages[n])
                self.Loaded_CenterImages_2_Range.append(self.data_2.Loaded_CenterImages[n])
                self.Loaded_Intensity_1_Range.append(self.data_1.Loaded_Intensity[n])
                self.Loaded_Intensity_2_Range.append(self.data_2.Loaded_Intensity[n])
                self.Loaded_Intensity_noBGS_1_Range.append(self.data_1.Loaded_Intensity_noBGS[n])
                self.Loaded_Intensity_noBGS_2_Range.append(self.data_2.Loaded_Intensity_noBGS[n])
                
                
                self.Loaded_Intensity_BG_1_Range.append(self.data_1.Loaded_Intensity_BG[n])
                self.Loaded_Intensity_BG_2_Range.append(self.data_2.Loaded_Intensity_BG[n])
                self.Loaded_Intensity_BG_Fit_1_Range.append(self.data_1.Loaded_Intensity_BG_Fit[n])
                self.Loaded_Intensity_BG_Fit_2_Range.append(self.data_2.Loaded_Intensity_BG_Fit[n])
                
                
                self.Loaded_Intensity_Fit_1_Range.append(self.data_1.Loaded_Intensity_Fit[n])
                self.Loaded_Intensity_Fit_2_Range.append(self.data_2.Loaded_Intensity_Fit[n])
                self.M_error_1_Range.append(self.data_1.M_error[n])
                self.M_error_2_Range.append(self.data_2.M_error[n])
                self.SNR_1_Range.append(self.data_1.SNR[n])
                self.SNR_2_Range.append(self.data_2.SNR[n])
                self.SBR_1_Range.append(self.data_1.Loaded_SBR[n])
                self.SBR_2_Range.append(self.data_2.Loaded_SBR[n])
                
                
                    
                
                try:
                    self.Loaded_ImageJxpos_3_Range.append(self.data_3.Loaded_ImageJxpos[n])
                    self.Loaded_ImageJypos_3_Range.append(self.data_3.Loaded_ImageJypos[n])
                    self.Loaded_A_3_Range.append(self.data_3.Loaded_A[n])
                    self.Loaded_T_3_Range.append(self.data_3.Loaded_T[n])
                    self.Loaded_phi_3_Range.append(self.data_3.Loaded_phi[n])
                    self.Loaded_y0_3_Range.append(self.data_3.Loaded_y0[n])
                    self.Loaded_M_3_Range.append(self.data_3.Loaded_M[n])
                    
                    self.Loaded_Intensity_noBGS_3_Range.append(self.data_3.Loaded_Intensity_noBGS[n])
                    self.Loaded_Intensity_3_Range.append(self.data_3.Loaded_Intensity[n])
                    self.Loaded_Intensity_Fit_3_Range.append(self.data_3.Loaded_Intensity_Fit[n])
                    self.Loaded_CenterImages_3_Range.append(self.data_3.Loaded_CenterImages[n])
                    self.Loaded_Intensity_BG_3_Range.append(self.data_3.Loaded_Intensity_BG[n])
                    self.Loaded_Intensity_BG_Fit_3_Range.append(self.data_3.Loaded_Intensity_BG_Fit[n])
                except: pass
    
    
    
                try:
                    self.Loaded_ImageJxpos_4_Range.append(self.data_4.Loaded_ImageJxpos[n])
                    self.Loaded_ImageJypos_4_Range.append(self.data_4.Loaded_ImageJypos[n])
                    self.Loaded_A_4_Range.append(self.data_4.Loaded_A[n])
                    self.Loaded_T_4_Range.append(self.data_4.Loaded_T[n])
                    self.Loaded_phi_4_Range.append(self.data_4.Loaded_phi[n])
                    self.Loaded_y0_4_Range.append(self.data_4.Loaded_y0[n])
                    self.Loaded_M_4_Range.append(self.data_4.Loaded_M[n])
                    
                    self.Loaded_Intensity_noBGS_4_Range.append(self.data_4.Loaded_Intensity_noBGS[n])
                    self.Loaded_Intensity_4_Range.append(self.data_4.Loaded_Intensity[n])
                    self.Loaded_Intensity_Fit_4_Range.append(self.data_4.Loaded_Intensity_Fit[n])
                    self.Loaded_CenterImages_4_Range.append(self.data_4.Loaded_CenterImages[n])
                    self.Loaded_Intensity_BG_4_Range.append(self.data_4.Loaded_Intensity_BG[n])
                    self.Loaded_Intensity_BG_Fit_4_Range.append(self.data_4.Loaded_Intensity_BG_Fit[n])
                except: pass
    
                    
                    
    
                try:
                    self.Loaded_ImageJxpos_5_Range.append(self.data_5.Loaded_ImageJxpos[n])
                    self.Loaded_ImageJypos_5_Range.append(self.data_5.Loaded_ImageJypos[n])
                    self.Loaded_A_5_Range.append(self.data_5.Loaded_A[n])
                    self.Loaded_T_5_Range.append(self.data_5.Loaded_T[n])
                    self.Loaded_phi_5_Range.append(self.data_5.Loaded_phi[n])
                    self.Loaded_y0_5_Range.append(self.data_5.Loaded_y0[n])
                    self.Loaded_M_5_Range.append(self.data_5.Loaded_M[n])
                    
                    self.Loaded_Intensity_noBGS_5_Range.append(self.data_5.Loaded_Intensity_noBGS[n])
                    self.Loaded_Intensity_5_Range.append(self.data_5.Loaded_Intensity[n])
                    self.Loaded_Intensity_Fit_5_Range.append(self.data_5.Loaded_Intensity_Fit[n])
                    self.Loaded_CenterImages_5_Range.append(self.data_5.Loaded_CenterImages[n])
                    self.Loaded_Intensity_BG_5_Range.append(self.data_5.Loaded_Intensity_BG[n])
                    self.Loaded_Intensity_BG_Fit_5_Range.append(self.data_5.Loaded_Intensity_BG_Fit[n])
                except: pass
            
            
            
    
                try:
                    self.Loaded_ImageJxpos_6_Range.append(self.data_6.Loaded_ImageJxpos[n])
                    self.Loaded_ImageJypos_6_Range.append(self.data_6.Loaded_ImageJypos[n])
                    self.Loaded_A_6_Range.append(self.data_6.Loaded_A[n])
                    self.Loaded_T_6_Range.append(self.data_6.Loaded_T[n])
                    self.Loaded_phi_6_Range.append(self.data_6.Loaded_phi[n])
                    self.Loaded_y0_6_Range.append(self.data_6.Loaded_y0[n])
                    self.Loaded_M_6_Range.append(self.data_6.Loaded_M[n])
                    
                    self.Loaded_Intensity_noBGS_6_Range.append(self.data_6.Loaded_Intensity_noBGS[n])
                    self.Loaded_Intensity_6_Range.append(self.data_6.Loaded_Intensity[n])
                    self.Loaded_Intensity_Fit_6_Range.append(self.data_6.Loaded_Intensity_Fit[n])
                    self.Loaded_CenterImages_6_Range.append(self.data_6.Loaded_CenterImages[n])
                    self.Loaded_Intensity_BG_6_Range.append(self.data_6.Loaded_Intensity_BG[n])
                    self.Loaded_Intensity_BG_Fit_6_Range.append(self.data_6.Loaded_Intensity_BG_Fit[n])
                except: pass
            
            
            
        
                    
                    
      
      



               
                 
                 
                 
                 

        if self.cb_SBR_1Range_yn.GetValue() == True or self.cb_SBR_2Range_yn.GetValue() == True:
            if self.data_1.Data_type_indicator != '%%':
                MessageBox(self, 'SBR is not supported with current data set\nUncheck SBR', 'Error')
                return

 

        if len(self.Loaded_ImageJxpos_1_Range) == 0:
            MessageBox(self, '# of data is zero', 'Error')
            return

        if len(self.Loaded_ImageJxpos_1_Range) == 1:
            MessageBox(self, '# of data is one', 'Error')
            return







        print '\n\n ## intensity ratio #2/#1'
        for n in range(len(self.Loaded_y0_1_Range)):
            print '# ' + str(n) + ' : ' , self.Loaded_y0_2_Range[n]/self.Loaded_y0_1_Range[n]















        data_conditions = ': Conditions : '

        if self.cb_M1Range_yn.GetValue() == True:
            data_conditions += '\n' + self.M1RangeValueMin.GetValue() +' <= M1 <= '+ self.M1RangeValueMax.GetValue() + ' ,   '
        
        if self.cb_M2Range_yn.GetValue() == True:
            data_conditions += self.M2RangeValueMin.GetValue() +' <= M2 <= '+ self.M2RangeValueMax.GetValue() + ' ,   '
        
        
        if self.cb_Intensity1Range_yn.GetValue() == True:
            data_conditions += '\n' + self.Intensity1RangeValueMin.GetValue() +' <= Intensity1 <= '+ self.Intensity1RangeValueMax.GetValue() + ' ,   '
        
        if self.cb_Intensity2Range_yn.GetValue() == True:
            data_conditions += self.Intensity2RangeValueMin.GetValue() +' <= Intensity2 <= '+ self.Intensity2RangeValueMax.GetValue() + ' ,   '
        
        
        
        if self.cb_Fit_y0_1_Range_yn.GetValue() == True:
            data_conditions += '\n' + self.Fit_y0_1RangeValueMin.GetValue() +' <= Fit:y0_#1 <= '+ self.Fit_y0_1RangeValueMax.GetValue() + ' ,   '
        

        if self.cb_Fit_y0_2_Range_yn.GetValue() == True:
            data_conditions += self.Fit_y0_2RangeValueMin.GetValue() +' <= Fit:y0_#2 <= '+ self.Fit_y0_2RangeValueMax.GetValue() + ' ,   '
        


        
        

        if self.cb_SNR_1Range_yn.GetValue() == True:
            data_conditions +='\n' +  self.SNR_1RangeValueMin.GetValue() +' <= SNR1 <= '+ self.SNR_1RangeValueMax.GetValue() + ' ,   '
            
        if self.cb_SNR_2Range_yn.GetValue() == True:
            data_conditions += self.SNR_2RangeValueMin.GetValue() +' <= SNR1 <= '+ self.SNR_2RangeValueMax.GetValue() + ' ,   '




        if self.cb_SBR_1Range_yn.GetValue() == True:
            data_conditions += '\n' +  self.SBR_1RangeValueMin.GetValue() +' <= SBR1 <= '+ self.SBR_1RangeValueMax.GetValue() + ' ,   '

        if self.cb_SBR_2Range_yn.GetValue() == True:
            data_conditions += self.SBR_2RangeValueMin.GetValue() +' <= SBR2 <= '+ self.SBR_2RangeValueMax.GetValue()            
            










        self.btn_plotData_text2.SetLabel(' =>  # ' + str(len(self.Loaded_M_1_Range)))
        
        
        
        
        
        


   

        if self.cb_publication_figures_yn.GetValue():
            title_text1 = '1. ' + self.data_1.fileName + '\n2. ' + self.data_2.fileName
            title_text2 = data_conditions
            M_before_bleaching_to_plotfunc = self.Loaded_M_1_Range
            M_after_bleaching_to_plotfunc = self.Loaded_M_2_Range
            Intensity_before_threshold_to_plotfunc = self.Loaded_y0_1_Range
            Intensity_after_threshold_to_plotfunc = self.Loaded_y0_2_Range
            dataplot(title_text1,title_text2, M_before_bleaching_to_plotfunc, M_after_bleaching_to_plotfunc, Intensity_before_threshold_to_plotfunc, Intensity_after_threshold_to_plotfunc)
            






   
        """       
        if self.cb_publication_figures_yn.GetValue():




        
            font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 14}
            matplotlib.rc('font', **font)
            majortick_size = 8
            majortick_width = 3
            plt.rc('axes', linewidth = 3)
            #matplotlib.rcParams['lines.linewidth'] = 4
            matplotlib.rcParams['xtick.major.size'] = majortick_size
            matplotlib.rcParams['xtick.major.width'] = majortick_width
            matplotlib.rcParams['xtick.minor.size'] = 5
            matplotlib.rcParams['xtick.minor.width'] = 4
            matplotlib.rcParams['ytick.major.size'] = majortick_size
            matplotlib.rcParams['ytick.major.width'] = majortick_width
            matplotlib.rcParams['ytick.minor.size'] = 5
            matplotlib.rcParams['ytick.minor.width'] = 4
            matplotlib.rcParams.update({'font.size': 22, 'family' : 'normal', 'weight' : 'bold',  'size': 8})
    
            
            fig_temp1 = plt.figure('temp1', facecolor='white')
            

            
            #plt.title('Excitation M')

            ax1 = fig_temp1.add_subplot(111)
            ax2 = ax1.twinx()
            
            
            
            plt.tick_params(left='on',  top='off', right='off', bottom='off', labelleft='on', labeltop='off', labelright='off', labelbottom='off')
                                   
            ax1.hist(self.Loaded_M_1_Range, bins=np.arange(0.0, 1.4, 0.1), label = '' + self.btn1_plot_label.GetValue(), alpha=.6, color='m', linewidth=1.5)
            #ax2.hist(self.Loaded_M_2_Range, bins=np.arange(0.0, 1.4, 0.1), label = '' + self.btn2_plot_label.GetValue(), alpha=.6, color='c', linewidth=1.5)
            
            #plt.xlabel('Modulation Depth (M)', size=11)

            plt.xlim(0.0, 1.2)
            
            ax1.legend(fontsize=22, loc='upper left', frameon=False)
            ax2.legend(fontsize=22, loc='upper right', frameon=False)
            
            ax1.set_ylim(0, 44)
            ax2.set_ylim(0, 44)
            
            ax1.set_yticks([0,20,40])
            
            
            
            ax1.tick_params(axis='x', labelsize=20)
            ax1.tick_params(axis='y', labelsize=20)
            #ax2.tick_params(axis='y', labelsize=20)
              
            ax1.set_xlabel('Modulation Depth (M)', size=22, fontweight='bold')
            ax1.set_ylabel('Occurence', color='k', size=22, fontweight='bold')
            ax2.set_ylabel('', color='b', size = 22)
            
            #plt.legend()
            
            plt.tight_layout()
            plt.show()
            
            
            
            
            '''
            #################################################################################################            
            ### hist to xy plots in files.

            todayDate = time.strftime("%Y%m%d_%Hh%Mm")

            
    
            bin_width = 0.1
            hist_x, hist_y = histtoxy2(self.Loaded_M_1_Range, 0.0, 1.0, bin_width, 'Y')
            
            #hist_x.append(1.0 + bin_width/2)
            #hist_y.append(hist_y[-1])
            
            x = np.array(hist_x) - bin_width/2
            y = np.array(hist_y)/float(np.max(hist_y))
            
            plt.figure('h1')
            plt.plot(x, y, '-ro')
            plt.bar(x,y, width=bin_width, alpha=0.7)
            plt.show()
            
            ff = open(self.data_1.fileDirectory + '\\' + 'CF_ext1_before_bleaching_bin_xyPlots_' + str(bin_width) + '_' + todayDate +'.txt','w')
            for i in range(len(x)):
                ff.write(str(x[i]) + ' ' + str(y[i]) + '\n'   )    
            ff.close()            

            
            
            bin_width = 0.1
            hist_x, hist_y = histtoxy2(self.Loaded_M_2_Range, 0.0, 1.0, bin_width, 'Y')
            
            #hist_x.append(1.0 + bin_width/2)
            #hist_y.append(hist_y[-1])
            
            x = np.array(hist_x) - bin_width/2
            y = np.array(hist_y)/float(np.max(hist_y))
            
            plt.figure('h2')
            plt.plot(x, y, '-ro')
            plt.bar(x,y, width=bin_width, alpha=0.7)
            plt.show()
            
            ff = open(self.data_1.fileDirectory + '\\' + 'CF_ext2_after_bleaching_bin_xyPlots_' + str(bin_width) + '_' + todayDate +'.txt','w')
            for i in range(len(x)):
                ff.write(str(x[i]) + ' ' + str(y[i]) + '\n'   )    
            ff.close()            
            
            
            
            print 'hist to xy plots in files created.'
            
            ### hist to xy plots in files.            
            #################################################################################################
            '''    
    
    



            
            fN_temp = int(self.publication_featureN.GetValue())
            
            
            figtemp = plt.figure('pf', figsize = (27,12), facecolor = 'white')
            figtemp.subplots_adjust(top = 0.9, bottom = 0.1, left = 0.05, right = 0.98, wspace = 0.5, hspace = 0.5)
    
            figtemp.text(0.02, 0.96, '#: ' + str(fN_temp) + ' ,   ' + self.data_1.fileName  , size=10, color="black", weight='normal')    
    

            mngr = plt.get_current_fig_manager()
            mngr.window.setGeometry(400,300,1500, 600)
            
            

          
            
            
            
            font = {'family' : 'normal',
            'weight' : 'bold',
            'size'   : 24}
            matplotlib.rc('font', **font)
            majortick_size = 8
            majortick_width = 3
            plt.rc('axes', linewidth = 3)
            #matplotlib.rcParams['lines.linewidth'] = 4
            matplotlib.rcParams['xtick.major.size'] = majortick_size
            matplotlib.rcParams['xtick.major.width'] = majortick_width
            matplotlib.rcParams['xtick.minor.size'] = 5
            matplotlib.rcParams['xtick.minor.width'] = 4
            matplotlib.rcParams['ytick.major.size'] = majortick_size
            matplotlib.rcParams['ytick.major.width'] = majortick_width
            matplotlib.rcParams['ytick.minor.size'] = 5
            matplotlib.rcParams['ytick.minor.width'] = 4
            matplotlib.rcParams.update({'font.size': 20, 'family' : 'normal', 'weight' : 'normal',  'size': 20})
       
    
                
            scalefontsize = 16
            #scalefontweight = 'bold'
            
               
            
            
            xaxis_deg = np.arange(0, 360, 2)
            
            imageTemp = np.append(np.ravel(self.Loaded_CenterImages_1_Range[fN_temp][0]), np.ravel(self.Loaded_CenterImages_1_Range[fN_temp][1]))
            vminTemp = int(np.mean( np.sort(imageTemp[1], axis=None)[:5]) * 0.5 )
            vmaxTemp = int(np.mean( np.sort(imageTemp[1], axis=None)[-5:]) * 2 )
  
  
            x3 = np.linspace(0, 15, 16)
            y3 = np.linspace(0, 15, 16)
            x3, y3 = np.meshgrid(x3, y3)
            

            plt.subplot(2,5, 1)
            plt.title('T = ' + str(int(np.around(self.Loaded_T_1_Range[fN_temp],0)) )  , fontsize=16, color = 'red', weight= 'normal')
            plt.tick_params(labelsize=scalefontsize)
            plt.imshow(self.Loaded_CenterImages_1_Range[fN_temp][1].reshape(15, 15), interpolation = 'None', vmin = 1000, vmax=10000, cmap=plt.cm.jet, origin='bottom',
                       extent=(x3.min(), x3.max(), y3.min(), y3.max()))
            plt.xlabel('pixel')                       
            plt.xticks(np.arange(0, 18, 3))
            plt.yticks(np.arange(0, 18, 3))
                     


                     
            #print '[self.Loaded_Intensity_1_Range[fN_temp][-1]]  ', [self.Loaded_Intensity_1_Range[fN_temp][-1]]
            #print 'len [self.Loaded_Intensity_1_Range[fN_temp]]  ', len(self.Loaded_Intensity_1_Range[fN_temp] + [self.Loaded_Intensity_1_Range[fN_temp][-1]])
            
            
                        
            
            
            print ' len(xaxis_deg) ', len(xaxis_deg)
            
            
          
            plt.subplot(2,5,2)
            plt.tick_params(labelsize=scalefontsize)
            plt.title( 'M = ' + str(np.around(self.Loaded_M_1_Range[fN_temp],2) )  , fontsize=16, color = 'red', weight='normal')
            plt.plot(xaxis_deg, self.Loaded_Intensity_1_Range[fN_temp], color = 'OrangeRed', linewidth = 2)
            if self.cb_fit_1_yn.GetValue():
                plt.plot(xaxis_deg, self.Loaded_Intensity_Fit_1_Range[fN_temp], color = 'b')
            plt.xlim(0,360)
            plt.locator_params(nbins=4)
            plt.xlabel('Rotation angle (deg)')

            plt.subplot(2,5,3)
            plt.tick_params(labelsize=scalefontsize)
            plt.plot(self.Loaded_Intensity_3_Range[fN_temp], color = 'green')
            if self.cb_include_Data3_fit_yn.GetValue():
                plt.plot(self.Loaded_Intensity_Fit_3_Range[fN_temp], color = 'b')
                plt.title( str(int(float(np.around(self.Loaded_phi_3_Range[fN_temp],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.Loaded_M_3_Range[fN_temp],2) )  , fontsize=16, color = 'blue', weight='normal')
            plt.locator_params(nbins=4)
            plt.xlim(0,1000)
            plt.xticks(np.arange(0, 1000, 500))
            plt.xlabel('Frame #')
          
          
            plt.subplot(2,5,4)
            plt.tick_params(labelsize=scalefontsize)
            plt.title('M = ' + str(np.around(self.Loaded_M_2_Range[fN_temp],2) )  , fontsize=16, color = 'blue', weight='normal')
            plt.plot(xaxis_deg, self.Loaded_Intensity_2_Range[fN_temp], color = 'OrangeRed', linewidth = 2)
            plt.plot(xaxis_deg, self.Loaded_Intensity_Fit_2_Range[fN_temp], color = 'b')
            plt.xlim(0,360)
            plt.locator_params(nbins=4)
            plt.xlabel('Rotation angle (deg)')            
          
            plt.subplot(2,5,5)
            plt.tick_params(labelsize=scalefontsize)
            plt.plot((self.Loaded_Intensity_4_Range[fN_temp]), color = 'limegreen')
            if self.cb_include_Data4_fit_yn.GetValue():
                plt.plot(self.Loaded_Intensity_Fit_4_Range[fN_temp], color = 'b')
                plt.title( str(int(float(np.around(self.Loaded_phi_4_Range[fN_temp],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.Loaded_M_4_Range[fN_temp],2) )  , fontsize=16, color = 'blue', weight='normal')
            plt.locator_params(nbins=4)
            plt.xlim(0,2000)
            plt.xticks(np.arange(0, 2000, 1000))
            plt.xlabel('Frame #')            
        






            plt.subplot(2,5, 6)
            #plt.title('T = ' + str(int(np.around(self.Loaded_T_1_Range[fN_temp],0)) )  , fontsize=16, color = 'red', weight= 'normal')
            plt.tick_params(labelsize=scalefontsize)
            plt.imshow(self.Loaded_CenterImages_1_Range[fN_temp][1].reshape(15, 15), interpolation = 'None', vmin = 1000, vmax=10000, cmap=plt.cm.jet, origin='bottom',
                       extent=(x3.min(), x3.max(), y3.min(), y3.max()))
            plt.xlabel('pixel', fontsize=20)                       
            plt.ylabel('pixel', fontsize=20)        
            plt.xticks(np.arange(0, 18, 3))
            plt.yticks(np.arange(0, 18, 3))
          
            plt.subplot(2,5,7)
            plt.tick_params(labelsize=scalefontsize)
            plt.title( 'M = ' + str(np.around(self.Loaded_M_1_Range[fN_temp],2) )  , fontsize=16, color = 'red', weight='normal')
            plt.plot(xaxis_deg, self.Loaded_Intensity_1_Range[fN_temp], color = 'OrangeRed', linewidth = 2)
            if self.cb_fit_1_yn.GetValue():
                plt.plot(xaxis_deg, self.Loaded_Intensity_Fit_1_Range[fN_temp], color = 'b')
            plt.xlim(0,360)
            plt.ylim(0,2000)
            plt.locator_params(nbins=4)
            plt.xlabel('Rotation angle '+ r"$\phi \ (^\circ)$", fontsize=18)     
            plt.ylabel('Intesnsity (a.u.)', fontsize=18)
            plt.xticks([0, 90, 180, 270, 360])
            plt.yticks([0, 1000, 2000])
            



             
            bleachingIntensityTrajectoryBS_1st = np.array(self.Loaded_Intensity_3_Range[fN_temp]) - np.array(self.Loaded_Intensity_BG_3_Range[fN_temp])
            bleachingIntensityTrajectoryBS_1st_per = 100*bleachingIntensityTrajectoryBS_1st / bleachingIntensityTrajectoryBS_1st[0]
            
            plt.subplot(2,5,8)
            #plt.title('1st photobleaching', fontsize=16, color = 'k', weight='normal')
            plt.tick_params(labelsize=scalefontsize)
            #plt.plot(self.Loaded_Intensity_3_Range[fN_temp], color = 'green')
            #plt.plot(bleachingIntensityTrajectoryBS_1st, color = 'b')
            
            plt.plot(np.arange(len(bleachingIntensityTrajectoryBS_1st_per))*0.2, bleachingIntensityTrajectoryBS_1st_per, color = 'limegreen')
            plt.locator_params(nbins=3)
            plt.xlim(0,200)
            plt.xticks([0, 100, 200])
            plt.xlabel(r"$^\circ$" + ' Time (s) ' + r"$^\circ$", fontsize=18)
            plt.ylim(0,120)
            plt.ylabel('Intensity (%)', fontsize=18)
            
            
          
            plt.subplot(2,5,9)
            plt.tick_params(labelsize=scalefontsize)
            plt.title('M = ' + str(np.around(self.Loaded_M_2_Range[fN_temp],2) )  , fontsize=16, color = 'red', weight='normal')
            plt.plot(xaxis_deg, self.Loaded_Intensity_2_Range[fN_temp], color = 'OrangeRed', linewidth = 2)
            plt.plot(xaxis_deg, self.Loaded_Intensity_Fit_2_Range[fN_temp], color = 'b')
            plt.xlim(0,360)
            plt.ylim(0, 7200)
            plt.locator_params(nbins=4)
            plt.xlabel('Rotation angle '+ r"$\phi \ (^\circ)$", fontsize=18)     
            plt.ylabel('Intesnsity (a.u.)', fontsize=18)
            plt.xticks([0, 90, 180, 270, 360])
            plt.yticks([0, 2000, 4000, 6000])
            
          
          
 
             
            bleachingIntensityTrajectoryBS_2nd = np.array(self.Loaded_Intensity_4_Range[fN_temp]) - np.array(self.Loaded_Intensity_BG_4_Range[fN_temp])
            bleachingIntensityTrajectoryBS_2nd_per = 100*bleachingIntensityTrajectoryBS_2nd / bleachingIntensityTrajectoryBS_2nd[0]
                     
          
            plt.subplot(2,5,10)
            #plt.title('2nd photobleaching', fontsize=16, color = 'k', weight='normal')
            plt.tick_params(labelsize=scalefontsize)
            plt.plot(np.arange(len(bleachingIntensityTrajectoryBS_2nd_per))*0.2, bleachingIntensityTrajectoryBS_2nd_per, color = 'limegreen')
            plt.locator_params(nbins=4)
            plt.xlim(0,400)
            plt.xticks([0, 200, 400])
            plt.xlabel(r"$^\circ$" + ' Time (s) ' + r"$^\circ$", fontsize=18)
            plt.ylim(0,120)
            plt.yticks([0,50,100])
            plt.ylabel('Intensity (%)', fontsize=18)
        









            #plt.tight_layout()

            

            plt.show()
            
        





            
            return




        else:
            pass
            

        """    
            
            


                   

        self.Loaded_Intensity_Fit_1_Range_CB = np.array(self.Loaded_y0_1_Range) - np.abs(self.Loaded_A_1_Range)  # CB == constant background
        self.Loaded_Intensity_Fit_2_Range_CB = np.array(self.Loaded_y0_2_Range) - np.abs(self.Loaded_A_2_Range)
        
        self.Loaded_Intensity_Fit_1_Range_Max = np.array(self.Loaded_y0_1_Range) + np.abs(self.Loaded_A_1_Range)  
        self.Loaded_Intensity_Fit_2_Range_Max = np.array(self.Loaded_y0_2_Range) + np.abs(self.Loaded_A_2_Range)
        
        self.Loaded_Intensity_Fit_1_Range_Min = np.array(self.Loaded_y0_1_Range) - np.abs(self.Loaded_A_1_Range)  
        self.Loaded_Intensity_Fit_2_Range_Min = np.array(self.Loaded_y0_2_Range) - np.abs(self.Loaded_A_2_Range)
        
    
 
    
    
    

    
        self.Loaded_Intensity_Fit_y0_Range_Diff = np.array(self.Loaded_y0_2_Range) -  np.array(self.Loaded_y0_1_Range)
        self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage = ( np.array(self.Loaded_y0_2_Range) - np.array(self.Loaded_y0_1_Range) )*100/np.array(self.Loaded_y0_1_Range)
        
        
        self.Loaded_Intensity_Fit_Range_CB_Diff_Percentage = (self.Loaded_Intensity_Fit_2_Range_CB - self.Loaded_Intensity_Fit_1_Range_CB)*100/self.Loaded_Intensity_Fit_1_Range_CB
        self.Loaded_Intensity_Fit_Range_A_Diff_Percentage = ( np.array(self.Loaded_A_2_Range) - np.array(self.Loaded_A_1_Range) )*100/np.array(self.Loaded_A_1_Range)
        









        self.Loaded_IntensityNoiseSTD_1_Range = []
        self.Loaded_IntensityNoiseSTD_2_Range = []
        for n in range(len(self.Loaded_Intensity_1_Range)):
            STDtemp1 = np.std( np.array(self.Loaded_Intensity_1_Range[n]) - np.array(self.Loaded_Intensity_Fit_1_Range[n]) )
            #print n , ' STDtemp1 ', ' = ' ,STDtemp1
            STDtemp2 = np.std( np.array(self.Loaded_Intensity_2_Range[n]) - np.array(self.Loaded_Intensity_Fit_2_Range[n]) )
            #print n , ' STDtemp2 ', ' = ' ,STDtemp2
            
            self.Loaded_IntensityNoiseSTD_1_Range.append(STDtemp1)
            self.Loaded_IntensityNoiseSTD_2_Range.append(STDtemp2)




        self.SNR_1_Range = np.array(self.Loaded_y0_1_Range)/np.array(self.Loaded_IntensityNoiseSTD_1_Range)
        self.SNR_2_Range = np.array(self.Loaded_y0_2_Range)/np.array(self.Loaded_IntensityNoiseSTD_2_Range)
        self.SNR_diff_Range = self.SNR_2_Range - self.SNR_1_Range


            

        self.Loaded_phi_diff_Range = np.array(self.Loaded_phi_2_Range) - np.array(self.Loaded_phi_1_Range)
        self.Loaded_phi_diff_abs_Range = np.abs(self.Loaded_phi_diff_Range)
        
        for n in range(len(self.Loaded_phi_diff_abs_Range)):
            if self.Loaded_phi_diff_abs_Range[n] >= 90:
                self.Loaded_phi_diff_abs_Range[n] = 180 - self.Loaded_phi_diff_abs_Range[n]
                
        
        if self.cb_90deg_offset_yn.GetValue():
            self.Loaded_phi_diff_abs_Range = 90 - self.Loaded_phi_diff_abs_Range
            self.Loaded_phi_diff_abs_text = ' 90 deg shift applied'
        else:
            self.Loaded_phi_diff_abs_text = ' '
            




        
        self.Loaded_M_diff_Range = np.array(self.Loaded_M_2_Range) - np.array(self.Loaded_M_1_Range)
        


        self.Loaded_Ndata_Range = len(self.Loaded_M_diff_Range)
        
        
        
        
        
        
        
        
        
        
        
        

        self.fig_DataPlots = plt.figure('Data_Plots', figsize = (27,12), facecolor='white')
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(110,30,1800, 1000)
        self.fig_DataPlots.subplots_adjust(top = 0.89, bottom = 0.07, left = 0.05, right = 0.98, wspace = 0.4, hspace = 0.4)
        self.fig_DataPlots.text(0.02, 0.93, '#1: ' + self.data_1.fileName + '\n' + '#2: ' + self.data_2.fileName + '\nTotal #: ' + str(self.Loaded_Ndata_Range)\
        + ',    BGS_type: ' + self.data_1.BGS_type + ', ' + self.data_2.BGS_type, size=12, color="black", weight='normal')
        
        self.fig_DataPlots.text(0.62, 0.93, data_conditions, size=12, color="red", weight='normal')

        
        
        #self.Data_Plots.text(0.02, 0.97, TotalFeatureN_text +',       ' + self.Corrections, ha="left", va="bottom", size=9,color="black", weight='normal')
        #plt.figtext(0.5, 0.92, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly+ '\n' + self.AnalysisConditions ), ha='center', color='black', weight='normal', size=9)

        #plt.title('Total # ' + str(len(IntensityTrajectoryChosenBS_fit_M)))
        plt.subplot(4,6,1)
        plt.title('Median: #1: ' + str(np.around(np.median(self.Loaded_M_1_Range),2)) + ' ,  #2: ' + str(np.around(np.median(self.Loaded_M_2_Range),2)), size = 10)
        plt.hist(self.Loaded_M_1_Range, bins=np.arange(0.0, 1.4, 0.1), label = '#1.'+fileType_1, alpha=1.0, color='red')
        plt.hist(self.Loaded_M_2_Range, bins=np.arange(0.0, 1.4, 0.1), label = '#2.'+fileType_2, alpha=0.3, color='blue')
        plt.xlabel('Modulation Depth (M)', size=11)
        plt.xlim(0.0, 1.31)
        plt.legend(fontsize=10, loc='upper left')
        
        
        
        plt.subplot(4,6,2)
        plt.title('t ', size=10)
        plt.scatter(self.Loaded_M_1_Range, self.Loaded_phi_diff_abs_Range, color='red', label = '#1.'+fileType_1)
        plt.scatter(self.Loaded_M_2_Range, self.Loaded_phi_diff_abs_Range, color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('phi diff ABS (M1, M2)', size=11)
        plt.legend(fontsize=10, loc='upper left')
        plt.xlim(-0.1, 1.21)
        
        
        
        plt.subplot(4,6,3, aspect='equal')
        plt.title('M2 vs M1 ', size=10)
        #plt.scatter(self.Loaded_M_1_Range, self.Loaded_M_2_Range, color='black', s = 4 )
        plt.errorbar(self.Loaded_M_1_Range, self.Loaded_M_2_Range, xerr = self.M_error_1_Range, yerr = self.M_error_2_Range, color='black', linestyle="None" )
        plt.plot([0,2],[0,2], 'r')
        plt.locator_params(nbins=5)
        plt.xlabel('M1', size=11)
        plt.ylabel('M2', size=11)
        plt.xlim(0.0, 1.21)
        plt.ylim(0.0, 1.21)
        
        
        

        plt.subplot(4,6,4)
        plt.title('M change vs M', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.Loaded_M_diff_Range, s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(self.Loaded_M_2_Range, self.Loaded_M_diff_Range, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.ylabel('M change  (M1, M2)', size=11)
        plt.xlabel('M', size=11)
        plt.legend(fontsize=9, loc='upper left')
        
        

        plt.subplot(4,6,5)
        plt.title('Fit: average intensity y0 vs M Change', size=11)
        plt.scatter(self.Loaded_M_diff_Range, self.Loaded_y0_1_Range, s = 6,  color='red', label = '#1.'+fileType_1)
        plt.xlabel('M Change', size=11)
        plt.ylabel('Fit y0', size=11)
        #plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')
        plt.locator_params(nbins=5)



        
        plt.subplot(4,6,6)
        plt.title('Fit: average intensity change % VS M Change', size=11)
        plt.scatter(self.Loaded_M_diff_Range, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage, s = 6,  color='red', label = '#1.'+fileType_1)
        plt.xlabel('M Change', size=11)
        plt.ylabel('Fit y0 diff ABS %', size=11)
        #plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')
        plt.locator_params(nbins=5)
        
        
        
        
        plt.subplot(4,6,7)
        plt.title('phi abs difference,' + self.Loaded_phi_diff_abs_text, size=11)
        plt.hist(self.Loaded_phi_diff_abs_Range, bins=np.arange(0,100,10))
        plt.xlabel('|phi2 - phi1| (deg)', size=11)
        plt.xlim(0,100)
        

        plt.subplot(4,6,8)
        plt.title('M difference', size=11)
        plt.hist(self.Loaded_M_diff_Range, bins=np.arange(-1,1.1,0.1))
        plt.xlabel('M2 - M1', size=11)
        
        
        
        
        
        plt.subplot(4,6,9)
        plt.title('Fit: average intensity y0', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.Loaded_y0_1_Range, s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(self.Loaded_M_2_Range, self.Loaded_y0_2_Range, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('Fit y0', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')

        
        
        plt.subplot(4,6,10)
        plt.title('Fit: average intensity change', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.Loaded_Intensity_Fit_y0_Range_Diff, s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(self.Loaded_M_2_Range, self.Loaded_Intensity_Fit_y0_Range_Diff, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('Fit y0 diff ABS (M1, M2)', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=9, loc='upper left')



        
        
        
        A = np.vstack([self.Loaded_M_1_Range, np.ones(len(self.Loaded_M_1_Range))]).T        
        m, c = np.linalg.lstsq(A, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage)[0]        
        x_temp = np.array([-0.2, 1.2])
        
        plt.subplot(4,6,11)
        plt.title('Fit: average intensity change %', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage, s = 6,  color='red', label = '#1.'+fileType_1)
        plt.plot(x_temp, m*x_temp + c, 'c', label='Fit: slope='+str(np.around(m,1)), linewidth = 2)
        plt.xlabel('Modulation Depth M1', size=11)
        plt.ylabel('Fit y0 diff ABS % (M1)', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')



        #for n in range(len(self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage)):
            #print '#,  %', n, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage[n]



        
        A = np.vstack([self.Loaded_M_2_Range, np.ones(len(self.Loaded_M_2_Range))]).T        
        m, c = np.linalg.lstsq(A, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage)[0]        
        x_temp = np.array([-0.2, 1.2])
        
        plt.subplot(4,6,12)
        plt.title('Fit: average intensity change %', size=11)
        plt.scatter(self.Loaded_M_2_Range, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.plot(x_temp, m*x_temp + c, 'c', label='Fit: slope='+str(np.around(m,1)), linewidth = 2)
        plt.xlabel('Modulation Depth M2', size=11)
        plt.ylabel('Fit y0 diff ABS % (M2)', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')

        
        
        
        
        plt.subplot(4,6,13)
        plt.title('Fit: average intensity y0', size=11)
        plt.scatter(self.Loaded_y0_1_Range, self.Loaded_y0_2_Range, s = 6,  color='red', label = ' ')
        
        plt.xlabel('#1 Fit y0', size=11)
        plt.ylabel('#2 Fit y0', size=11)
        #plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')
        plt.locator_params(nbins=5)

        





        plt.subplot(4,6,14)
        plt.title('Fit: modulation A change', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.Loaded_Intensity_Fit_Range_A_Diff_Percentage, s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(self.Loaded_M_2_Range, self.Loaded_Intensity_Fit_Range_A_Diff_Percentage, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('Fit A diff % (M1, M2)', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')
        plt.locator_params(nbins=5)

        
        

        plt.subplot(4,6,15)
        plt.title('Fit: constant intensity change', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.Loaded_Intensity_Fit_Range_CB_Diff_Percentage, s = 6, color='red', label = '#1.'+fileType_1)
        plt.scatter(self.Loaded_M_2_Range, self.Loaded_Intensity_Fit_Range_CB_Diff_Percentage, s = 6, color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('Fit CB diff % (M1, M2)', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')
        plt.locator_params(nbins=5)

        
        
                
        

        
        
        

        
        
        plt.subplot(4,6,16)
        plt.title('NoiseSTD vs M', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.Loaded_IntensityNoiseSTD_1_Range,  s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(self.Loaded_M_2_Range, self.Loaded_IntensityNoiseSTD_2_Range, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('M', size=11)
        plt.legend(fontsize=9, loc='upper left')
        plt.locator_params(nbins=5)
        
        
        
        
        
        A = np.vstack([self.Loaded_y0_1_Range, np.ones(len(self.Loaded_y0_1_Range))]).T        
        m, c = np.linalg.lstsq(A, self.Loaded_IntensityNoiseSTD_1_Range)[0]        
        x_temp = np.arange(-1000, 20000, 10)
    
        A_guess = m*10
        B_guess = np.min(self.Loaded_IntensityNoiseSTD_1_Range)
        
        initial_guess = [A_guess, B_guess]
        
        try:
            popt, pcov = curve_fit(Noise_func, self.Loaded_y0_1_Range, self.Loaded_IntensityNoiseSTD_1_Range, p0=initial_guess )
            fit_temp = Noise_func(x_temp, popt[0],popt[1])
            #err_temp = np.sqrt(np.diag(pcov))
        except: 
            print '\n fit error '
            fit_temp = x_temp *0.0
            popt = [0, 0]
            pass
   


        plt.subplot(4,6,17)
        plt.title('NoiseSTD vs Ave Intensity (y0)', size=11)
        plt.scatter(self.Loaded_y0_1_Range, self.Loaded_IntensityNoiseSTD_1_Range,  s = 5,  color='red', label = '#1.'+fileType_1)
        plt.plot(x_temp,  m*x_temp + c, 'c', label='Linear Fit: slope='+str(np.around(m,1)), linewidth = 2)
        plt.plot(x_temp, fit_temp, 'm', label='Noise Fit: A,B= '+str(np.around(popt[0],1))+','+str(np.around(popt[1],1)), linewidth = 2)
        #plt.scatter(self.Loaded_y0_2_Range, self.Loaded_IntensityNoiseSTD_2, s = 2,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Ave Intensity', size=11)
        plt.legend(fontsize=9, loc='upper left')
        plt.xlim(0, np.max(self.Loaded_y0_1_Range))
        plt.ylim(0, np.max(self.Loaded_IntensityNoiseSTD_1_Range))
        plt.locator_params(nbins=5)
        #plt.xticks(np.arange(0, np.max(self.Loaded_y0_1_Range), )
        









   
        NetNoiseSTD = []
        for n in range(len(self.Loaded_M_1_Range)):
            dI_temp = Noise_func(self.Loaded_y0_1_Range[n], popt[0],popt[1])
            NetNoiseTemp = self.Loaded_IntensityNoiseSTD_1_Range[n] - dI_temp
            NetNoiseSTD.append(NetNoiseTemp)
            
            
        A = np.vstack([self.Loaded_M_1_Range, np.ones(len(self.Loaded_M_1_Range))]).T        
        m, c = np.linalg.lstsq(A, NetNoiseSTD)[0]        
        x_temp = np.arange(-0, 1.2, 0.1)
        
        
        plt.subplot(4,6,18)
        plt.title('Renormalized NetNoiseSTD vs M', size=11)
        plt.scatter(self.Loaded_M_1_Range, NetNoiseSTD,  s = 6,  color='red', label = '#1.'+fileType_1)
        plt.plot(x_temp,  m*x_temp + c, 'c', label='Linear Fit: slope='+str(np.around(m,1)), linewidth = 2)
        #plt.scatter(self.Loaded_M_2_Range, self.Loaded_IntensityNoiseSTD_2, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('M', size=11)
        plt.legend(fontsize=9, loc='upper left')
        plt.locator_params(nbins=5)





        
        plt.subplot(4,6,19)
        plt.title('SBR vs M', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.SBR_1_Range, s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(self.Loaded_M_2_Range, self.SBR_2_Range, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('M', size=11)
        plt.ylabel('SBR', size=11)
        plt.legend(fontsize=9, loc='upper right')
        plt.locator_params(nbins=5)
        
        





        plt.subplot(4,6,20)
        plt.hist(self.SNR_1_Range, label = '#1.'+fileType_1, alpha=1.0, color='red')
        plt.hist(self.SNR_2_Range, label = '#2.'+fileType_2, alpha=0.3, color='blue')
        plt.xlabel('SNR', size=11)
        #plt.xlim(0.0, 1.21)
        plt.legend(fontsize=10, loc='upper right')
        
        
        plt.subplot(4,6,21)
        plt.title('SNR vs M', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.SNR_1_Range,  s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(self.Loaded_M_2_Range, self.SNR_2_Range, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('M', size=11)
        plt.ylabel('SNR', size=11)
        plt.legend(fontsize=9, loc='upper right')
        
        




        
        plt.subplot(4,6,22)
        plt.title('SNR change vs M change', size=11)
        plt.scatter(self.Loaded_M_diff_Range, self.SNR_diff_Range,  s = 6,  color='red', label = '#1.'+fileType_1)
        plt.xlabel('M change', size=11)
        plt.ylabel('SNR change', size=11)
        plt.legend(fontsize=9, loc='upper left')
        plt.locator_params(nbins=5)
        
        



        
        plt.subplot(4,6,23)
        plt.title('M_error vs SNR', size=11)
        plt.scatter(self.SNR_1_Range, self.M_error_1_Range,  s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(self.SNR_2_Range, self.M_error_2_Range, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('SNR', size=11)
        plt.ylabel('M_error', size=11)
        plt.legend(fontsize=9, loc='upper left')
        
        


        
        plt.subplot(4,6,24)
        plt.title('M_error vs 1/SNR', size=11)
        plt.scatter(1/self.SNR_1_Range, self.M_error_1_Range,  s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(1/self.SNR_2_Range, self.M_error_2_Range, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('1/SNR', size=11)
        plt.ylabel('M_error', size=11)
        plt.legend(fontsize=9, loc='upper left')
        plt.locator_params(nbins=6)
        
        






































        if self.cb_2nd_Plot_yn.GetValue() or self.cb_publication_figures_yn.GetValue():
    
            self.fig_DataPlots_2 = plt.figure('Data_Plots_2', figsize = (27,12), facecolor = 'white')
            mngr = plt.get_current_fig_manager()
            mngr.window.setGeometry(210,30,1800, 1000)
            self.fig_DataPlots_2.subplots_adjust(top = 0.89, bottom = 0.07, left = 0.05, right = 0.98, wspace = 0.4, hspace = 0.4)
            self.fig_DataPlots_2.text(0.02, 0.93, '#1.' + self.data_1.fileName + '\n' + '#2.' + self.data_2.fileName + '\nTotal #: ' + str(self.Loaded_Ndata_Range)\
            + ',    BGS_type: ' + self.data_1.BGS_type + ', ' + self.data_2.BGS_type, size=12, color="black", weight='normal')
            
            self.fig_DataPlots_2.text(0.62, 0.93, data_conditions, size=12, color="red", weight='normal')
    
            
            
            
    
            plt.subplot(4,6,1)
            plt.hist(self.Loaded_M_1_Range, bins=np.arange(0.0, 1.4, 0.1), label = '#1.'+fileType_1, alpha=1.0, color='red')
            #plt.hist(self.Loaded_M_2_Range, bins=np.arange(0.0, 1.4, 0.1), label = '#2.'+fileType_2, alpha=0.3, color='blue')
            plt.xlabel('Modulation Depth (M)', size=11)
            plt.xlim(0.0, 1.31)
            plt.legend(fontsize=10, loc='upper left')
            
    
            plt.subplot(4,6,2)
            #plt.hist(self.Loaded_M_1_Range, bins=np.arange(0.0, 1.4, 0.1), label = '#1.'+fileType_1, alpha=1.0, color='red')
            plt.hist(self.Loaded_M_2_Range, bins=np.arange(0.0, 1.4, 0.1), label = '#2.'+fileType_2, alpha=0.3, color='blue')
            plt.xlabel('Modulation Depth (M)', size=11)
            plt.xlim(0.0, 1.31)
            plt.legend(fontsize=10, loc='upper left')
            
            
    
    
    
            x_temp = np.array([-0.2, 1.3])
            initial_guess = [-0.5, 0.5]
    
            popt, pcov = opt.curve_fit(Linear_func, self.Loaded_M_1_Range, self.Loaded_M_diff_Range, p0 = initial_guess)
            #print '\nFIONA popt: ', popt
            #print '\nFIONA pcov: \n', pcov
            Afit = popt[0]
            Bfit = popt[1]
            AfitErr = np.sqrt(pcov[0,0])                    
            #BfitErr = np.sqrt(pcov[1,1])
            
            
            
            
            
            plt.subplot(4,6,3)
            plt.title('M change vs M', size=11)
            plt.scatter(self.Loaded_M_1_Range, self.Loaded_M_diff_Range, s = 6,  color='red', label = '#1.'+fileType_1)
            plt.plot(x_temp, Afit*x_temp + Bfit, 'c', label='Fit: slope='+str(np.around(Afit,2)) + ' Err: ' + str(np.around(AfitErr,2)), linewidth = 2)
            
            
            plt.xlim(0, 1.2)        
            plt.ylabel('M change  (M1, M2)', size=11)
            plt.xlabel('M', size=11)
            plt.legend(fontsize=9, loc='upper left')
            
            
            
            
            
            
            self.Loaded_M_1_bins_Range = np.arange(0.0, 1.2, 0.1)        
            self.Loaded_M_1_bins_elements_Range = [[] for _ in xrange(len(self.Loaded_M_1_bins_Range))]
            self.Loaded_M_1_bins_elements_increased_Range = [[] for _ in xrange(len(self.Loaded_M_1_bins_Range))]
            self.Loaded_M_1_bins_elements_decreased_Range = [[] for _ in xrange(len(self.Loaded_M_1_bins_Range))]
            self.Loaded_M_1_bins_numbers_of_elements_Range = [0 for _ in xrange(len(self.Loaded_M_1_bins_Range))]
            self.Loaded_M_1_bins_numbers_of_elements_increased_Range = [0 for _ in xrange(len(self.Loaded_M_1_bins_Range))]
            self.Loaded_M_1_bins_numbers_of_elements_decreased_Range = [0 for _ in xrange(len(self.Loaded_M_1_bins_Range))]
            self.Loaded_M_1_bins_Ave_M_diff_Range = [0 for _ in xrange(len(self.Loaded_M_1_bins_Range))]
            
            
            for n in range(len(self.Loaded_M_1_Range)):
    
                M_temp = self.Loaded_M_1_Range[n]
                Mdiff_temp = self.Loaded_M_diff_Range[n]
                
                for m in range(len(self.Loaded_M_1_bins_Range)):
                    if m*0.1 <= M_temp <= m*0.1 + 0.1:
                        self.Loaded_M_1_bins_elements_Range[m].append(Mdiff_temp)
                        if Mdiff_temp > 0:
                            self.Loaded_M_1_bins_elements_increased_Range[m].append(Mdiff_temp)
                        if Mdiff_temp < 0:
                            self.Loaded_M_1_bins_elements_decreased_Range[m].append(Mdiff_temp)
                        
                        
            for n in range(len(self.Loaded_M_1_bins_Range)):
                self.Loaded_M_1_bins_numbers_of_elements_Range[n] = len(self.Loaded_M_1_bins_elements_Range[n])
                self.Loaded_M_1_bins_numbers_of_elements_increased_Range[n] = len(self.Loaded_M_1_bins_elements_increased_Range[n])
                self.Loaded_M_1_bins_numbers_of_elements_decreased_Range[n] = len(self.Loaded_M_1_bins_elements_decreased_Range[n])
                self.Loaded_M_1_bins_Ave_M_diff_Range[n] = np.mean(self.Loaded_M_1_bins_elements_Range[n])
                
                
                
                
                
            #print 'self.Loaded_M_1_bins_Range ', self.Loaded_M_1_bins_Range
            #print 'self.Loaded_M_1_bins_elements_Range ', self.Loaded_M_1_bins_elements_Range
            #print 'self.Loaded_M_1_bins_Ave_M_diff_Range ', self.Loaded_M_1_bins_Ave_M_diff_Range
            
    
    
    
                
            plt.subplot(4,6,4)
            plt.title('Ave M change vs M', size=11)
            #plt.scatter(self.Loaded_M_1_bins_Range, self.Loaded_M_1_bins_Ave_M_diff_Range, s = 6,  color='red', label = '#1.'+fileType_1)
            plt.bar(self.Loaded_M_1_bins_Range, self.Loaded_M_1_bins_Ave_M_diff_Range, 0.1)
            
            plt.xlim(0, 1.2)
            plt.ylabel('Ave M change  (M1, M2)', size=11)
            plt.xlabel('M1', size=11)
            plt.legend(fontsize=9, loc='upper right')
            
            
            
            
            plt.subplot(4,6,5)
            plt.title('# vs M', size=11)
            plt.bar(self.Loaded_M_1_bins_Range, self.Loaded_M_1_bins_numbers_of_elements_Range, 0.1)
            
            plt.xlim(0, 1.2)
            plt.ylabel('#', size=11)
            plt.xlabel('M1', size=11)
            plt.legend(fontsize=9, loc='upper right')
            
            
                    
                           
            
            
            plt.subplot(4,6,6)
            plt.title('# vs M', size=11)
            plt.bar(self.Loaded_M_1_bins_Range, self.Loaded_M_1_bins_numbers_of_elements_increased_Range, 0.1, color='m', alpha=0.4, label='increased')
            plt.bar(self.Loaded_M_1_bins_Range, self.Loaded_M_1_bins_numbers_of_elements_decreased_Range, 0.1, color='c', alpha=0.4, label='decreased')
            
            plt.xlim(0, 1.2)
            plt.ylabel('#', size=11)
            plt.xlabel('M1', size=11)
            plt.legend(fontsize=9, loc='upper right')
            
            
                    
                           
                    
            
            
            
            
            
            plt.subplot(4,6,7)
            plt.title('Ave Intensity (fit: y0)', size=11)
            plt.hist(self.Loaded_y0_1_Range,  label = '#1.'+fileType_1, alpha=0.5, color='red')
            #plt.hist(self.Loaded_y0_2_Range,  label = '#2.'+fileType_2, alpha=0.5, color='blue')
            
            plt.xlabel('Fit: Ave Intensity (y0)', size=11)
            plt.legend(fontsize=10, loc='upper left')
            plt.locator_params(nbins=5)
            
            
    
            
            
            
            plt.subplot(4,6,8)
            plt.title('Ave Intensity (fit: y0)', size=11)
            #plt.hist(self.Loaded_y0_1_Range,  label = '#1.'+fileType_1, alpha=0.5, color='red')
            plt.hist(self.Loaded_y0_2_Range,  label = '#2.'+fileType_2, alpha=0.5, color='blue')
            
            plt.xlabel('Fit: Ave Intensity (y0)', size=11)
            plt.legend(fontsize=10, loc='upper left')
            plt.locator_params(nbins=5)
            
            
    
            
            
            
            plt.subplot(4,6,9)
            plt.title('Fit: average intensity y0', size=11)
            plt.scatter(self.Loaded_M_1_Range, self.Loaded_y0_1_Range, s = 6,  color='red', label = '#1.'+fileType_1)
            #plt.scatter(self.Loaded_M_2_Range, self.Loaded_y0_2_Range, s = 6,  color='blue', label = '#2.'+fileType_2)
            plt.xlabel('Modulation Depth M1', size=11)
            plt.ylabel('Fit y0', size=11)
            plt.xlim(-0.2, 1.21)
            plt.legend(fontsize=10, loc='upper left')
    
            
            
            
                    
            
            plt.subplot(4,6,10)
            plt.title('Fit: average intensity y0', size=11)
            #plt.scatter(self.Loaded_M_1_Range, self.Loaded_y0_1_Range, s = 6,  color='red', label = '#1.'+fileType_1)
            plt.scatter(self.Loaded_M_2_Range, self.Loaded_y0_2_Range, s = 6,  color='blue', label = '#2.'+fileType_2)
            plt.xlabel('Modulation Depth M2', size=11)
            plt.ylabel('Fit y0', size=11)
            plt.xlim(-0.2, 1.21)
            plt.legend(fontsize=10, loc='upper left')
    
            

            self.Loaded_Intensity_BG_1_Range_median = []
            self.Loaded_Intensity_BG_2_Range_median = []
            for n in range(len(self.Loaded_Intensity_BG_1_Range)):
                self.Loaded_Intensity_BG_1_Range_median.append(np.median(self.Loaded_Intensity_BG_1_Range[n]))
                self.Loaded_Intensity_BG_2_Range_median.append(np.median(self.Loaded_Intensity_BG_2_Range[n]))
                
                
                        
    
            plt.subplot(4,6,11)
            plt.title('Median: #1: ' + str(np.around(np.median(self.Loaded_Intensity_BG_1_Range),0)) + '\nMedian: #2: ' + str(np.around(np.median(self.Loaded_Intensity_BG_2_Range),0)), size = 10)
            plt.hist(self.Loaded_Intensity_BG_1_Range_median, label = '#1', alpha=.4, color='c')
            plt.hist(self.Loaded_Intensity_BG_2_Range_median, label = '#2', alpha=.4, color='m')
            plt.xlabel('BG intensity', size=11)
            #plt.xlim(0.0, 1.31)
            plt.legend(fontsize=9, loc='upper right',frameon=False)
            plt.locator_params(nbins=4)
        
                    
                    
                    
                    
            
            plt.subplot(4,6,12)
            plt.title('Phi diff abs vs M1', size=11)
            #plt.hist(self.Loaded_y0_1_Range,  label = '#1.'+fileType_1, alpha=0.5, color='red')
            plt.scatter(self.Loaded_M_1_Range, self.Loaded_phi_diff_abs_Range, s = 6,  color='red', label = '')
            
            plt.xlabel('M1', size=11)
            plt.legend(fontsize=10, loc='upper left')
            plt.locator_params(nbins=5)
            
 












           


    
            plt.subplot(4,6,13)
            plt.title('Median: #1: ' + str(np.around(np.median(self.Loaded_phi_1_Range),2)), size = 10)
            plt.hist(self.Loaded_phi_1_Range, label = '#1', alpha=1.0, color='red')
            plt.xlabel('phi (deg)', size=11)
            #plt.xlim(0.0, 1.31)
            plt.legend(fontsize=9, loc='upper right',frameon=False)
        

            
        
            
            plt.subplot(4,6,14)
            plt.title('Median: #2: ' + str(np.around(np.median(self.Loaded_phi_2_Range),2)), size = 10)
            plt.hist(self.Loaded_phi_2_Range, label = '#2', alpha=1.0, color='blue')
            plt.xlabel('phi (deg)', size=11)
            #plt.xlim(0.0, 1.31)
            plt.legend(fontsize=9, loc='upper right',frameon=False)
            

            
    
    
            
            
            plt.subplot(4,6,15)
            plt.title('Fit: average intensity change', size=11)
            plt.scatter(self.Loaded_M_1_Range, self.Loaded_Intensity_Fit_y0_Range_Diff, s = 6,  color='red', label = '#1.'+fileType_1)
            #plt.scatter(self.Loaded_M_2_Range, self.Loaded_Intensity_Fit_y0_Range_Diff, s = 6,  color='blue', label = '#2.'+fileType_2)
            plt.xlabel('Modulation Depth M1', size=11)
            plt.ylabel('Fit y0 diff ABS (M1)', size=11)
            plt.xlim(-0.2, 1.21)
            plt.legend(fontsize=9, loc='upper left')
    
    
    
    
            
            
            plt.subplot(4,6,16)
            plt.title('Fit: average intensity change', size=11)
            #plt.scatter(self.Loaded_M_1_Range, self.Loaded_Intensity_Fit_y0_Range_Diff, s = 6,  color='red', label = '#1.'+fileType_1)
            plt.scatter(self.Loaded_M_2_Range, self.Loaded_Intensity_Fit_y0_Range_Diff, s = 6,  color='blue', label = '#2.'+fileType_2)
            plt.xlabel('Modulation Depth M2', size=11)
            plt.ylabel('Fit y0 diff ABS (M2)', size=11)
            plt.xlim(-0.2, 1.21)
            plt.legend(fontsize=9, loc='upper left')
    
    
    
    
            #A = np.vstack([self.Loaded_M_1_Range, np.ones(len(self.Loaded_M_1_Range))]).T        
            #m, c = np.linalg.lstsq(A, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage)[0]        
            x_temp = np.array([-0.2, 1.3])
            initial_guess = [0.0, -80.0]
    
            popt, pcov = opt.curve_fit(Linear_func, self.Loaded_M_1_Range, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage, p0 = initial_guess)
            #print '\nFIONA popt: ', popt
            #print '\nFIONA pcov: \n', pcov
            Afit = popt[0]
            Bfit = popt[1]
            
            
            AfitErr = np.sqrt(pcov[0,0])                    
            #BfitErr = np.sqrt(pcov[1,1])
            
            
            
            plt.subplot(4,6,17)
            plt.title('Fit: average intensity change %', size=11)
            plt.scatter(self.Loaded_M_1_Range, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage, s = 6,  color='red', label = '#1.'+fileType_1)
            plt.plot(x_temp, Afit*x_temp + Bfit, 'c', label='Fit: slope='+str(np.around(Afit,1)) + ' Err: ' + str(np.around(AfitErr,1)), linewidth = 2)
            plt.xlabel('Modulation Depth M1', size=11)
            plt.ylabel('Fit y0 diff ABS % (M1)', size=11)
            plt.xlim(-0.2, 1.21)
            plt.legend(fontsize=10, loc='upper left')
    
    
    
    
            intesnsity_ratio_y0 = np.array(self.Loaded_y0_2_Range)/(np.array(self.Loaded_y0_1_Range).astype(float))
    
            FitMinIntesnityRatio = (self.Loaded_Intensity_Fit_2_Range_Min.astype(float))/(self.Loaded_Intensity_Fit_1_Range_Min.astype(float))
            FitMaxIntesnityRatio = (self.Loaded_Intensity_Fit_2_Range_Max.astype(float))/(self.Loaded_Intensity_Fit_1_Range_Max.astype(float))
                                        
            
            plt.subplot(4,6,19)
            plt.title('Fit: average intensity y0 ratio', size=11)
            #plt.scatter(self.Loaded_M_1_Range, self.Loaded_y0_1_Range, s = 6,  color='red', label = '#1.'+fileType_1)
            plt.scatter(self.Loaded_M_1_Range, intesnsity_ratio_y0, s = 6,  color='blue', label = '')
            plt.xlabel('Modulation Depth M1', size=11)
            plt.ylabel('Fit y0 ratio:  #2/#1', size=11)
            plt.xlim(-0.2, 1.21)
            #plt.legend(fontsize=10, loc='upper right',frameon=False)
    
          

            
            plt.subplot(4,6,20)
            plt.title('Fit: Imin, Imax intensity ratio', size=11)
            plt.scatter(self.Loaded_M_1_Range, FitMinIntesnityRatio, s = 12,  color='c', label = 'Imin ratio')
            plt.scatter(self.Loaded_M_1_Range, FitMaxIntesnityRatio, s = 12,  color='m', label = 'Imax ratio')
            plt.xlabel('Modulation Depth M1', size=11)
            plt.ylabel('Fit intensity ratio #2/#1', size=11)
            plt.xlim(-0.2, 1.21)
            plt.legend(fontsize=10, loc='upper left',frameon=False)
    
              
    

            plt.subplot(4,6,21)
            plt.title('Fit: IminRatio/ImaxRatio', size=11)
            plt.scatter(self.Loaded_M_1_Range, FitMinIntesnityRatio/FitMaxIntesnityRatio, s = 12,  color='c', label = '')
            plt.xlabel('Modulation Depth M1', size=11)
            plt.ylabel('Ratio', size=11)
            plt.xlim(-0.2, 1.21)
            #plt.legend(fontsize=10, loc='upper left',frameon=False)
    
              
    



                    
            
            plt.subplot(4,6,22)
            plt.title('Median: ' + str(np.around(np.median(intesnsity_ratio_y0),2)) , size=11)
            bins = np.arange(0, 10, 0.2)
            plt.hist(intesnsity_ratio_y0, bins=bins)
            plt.xlabel('intensity ratio #2/#1', size=11)
            plt.ylabel('Occurence', size=11)
            #plt.legend(fontsize=10, loc='upper right',frameon=False)
    
          


    
    
    
    
    
    
            plt.subplot(4,6,23, aspect='equal')
            plt.title('M2 vs M1 ', size=10)
            
            colors = np.array(self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage)
            plt.scatter(self.Loaded_M_1_Range, self.Loaded_M_2_Range, c=colors, s=50, marker='o', linewidth='0', alpha=0.5 )
            plt.colorbar()
    
    
            
            
            plt.plot([0,2],[0,2], 'r')
            plt.xlabel('M1', size=11)
            plt.ylabel('M2', size=11)
            plt.xlim(0.0, 1.21)
            plt.ylim(0.0, 1.21)
            
            
    
    
            plt.ion()
            plt.show()
            
            
    
















   

        if self.cb_publication_figures_yn.GetValue() == True:


            #'''
 
            plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='on',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='on') # labels along the bottom edge are off
            plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left='on',      # ticks along the bottom edge are off
            right='off',         # ticks along the top edge are off
            labelleft='on') # labels along the bottom edge are off
            
            font = {'family' : 'normal','weight' : 'normal', 'size'   : 14}
            matplotlib.rc('font', **font)
            majortick_size = 5
            majortick_width = 2
            plt.rc('axes', linewidth = 2)
            #matplotlib.rcParams['lines.linewidth'] = 4
            matplotlib.rcParams['xtick.major.size'] = majortick_size
            matplotlib.rcParams['xtick.major.width'] = majortick_width
            matplotlib.rcParams['xtick.minor.size'] = 5
            matplotlib.rcParams['xtick.minor.width'] = 4
            matplotlib.rcParams['ytick.major.size'] = majortick_size
            matplotlib.rcParams['ytick.major.width'] = majortick_width
            matplotlib.rcParams['ytick.minor.size'] = 5
            matplotlib.rcParams['ytick.minor.width'] = 4
            matplotlib.rcParams.update({'font.size': 20, 'family' : 'normal', 'weight' : 'normal',  'size': 8})
            #'''
                   
            
            
            
            

            
            
            
            
            
            self.fig_temp2 = plt.figure('fig_temp2',figsize = (16,8), facecolor='white')
            self.fig_temp2.subplots_adjust(top = 0.89, bottom = 0.12, left = 0.1, right = 0.98, wspace = 0.4, hspace = 0.4)
        
            self.fig_temp2.text(0.02, 0.93, '#1: ' + self.data_1.fileName + '\n' + '#2: ' + self.data_2.fileName + '\nTotal #: ' + str(self.Loaded_Ndata_Range)\
            + ',    BGS_type: ' + self.data_1.BGS_type + ', ' + self.data_2.BGS_type, size=12, color="black", weight='normal')
            
            self.fig_temp2.text(0.62, 0.93, data_conditions, size=12, color="red", weight='normal')
    



             
            plt.subplot(2,3,1)
            plt.rcParams.update({'mathtext.default': 'regular' })
            
            plt.title('median: '+str(np.around(np.median(self.Loaded_y0_1_Range),1)))

            #plt.title('Average intensity ', size=18)
            plt.plot([0,2],[0,0], 'k--', linewidth=1.4)
            plt.scatter(self.Loaded_M_1_Range, self.Loaded_y0_1_Range, s = 50,  color='red', label = '#1.'+fileType_1, alpha=0.5)
            #plt.scatter(self.Loaded_M_2_Range, self.Loaded_y0_2_Range, s = 6,  color='blue', label = '#2.'+fileType_2)
            plt.xlabel('$M$', size=22)
            plt.ylabel('$I_{ave}$ (a.u.)', size=20)
            
            plt.xticks([0, 0.5, 1.0])
            plt.xlim(-0.01, 1.2)
            
            
            plt.yticks([0, 2000, 4000, 6000, 8000])
            plt.ylim(-1000, 9000)
            
            #plt.legend(fontsize=10, loc='upper left')
    
            
        
        
        
         
            
            plt.subplot(2,3,2)
            #plt.title('Average intensity change', size=18)
            plt.scatter(self.Loaded_M_1_Range, self.Loaded_Intensity_Fit_y0_Range_Diff, s = 6,  color='red', label = '#1.'+fileType_1)
            #plt.scatter(self.Loaded_M_2_Range, self.Loaded_Intensity_Fit_y0_Range_Diff, s = 6,  color='blue', label = '#2.'+fileType_2)
            plt.xlabel('$M$', size=22)
            plt.ylabel(r'$\Delta I_{ave} \, (a.u.)$', size=20)
            
            
            plt.xticks([0, 0.5, 1.0])
            plt.xlim(-0.01, 1.2)
            
            
            plt.yticks([-6000, -4000, -2000, 0])
            #plt.legend(fontsize=9, loc='upper left')
    
    


            
            
        
            
            plt.subplot(2,3,3)
            plt.title('Fit: slope='+str(np.around(Afit,1)) + ' ,   Err: ' + str(np.around(AfitErr,1)), size=18)
            plt.plot([0,2],[0,0], 'k--', linewidth=1.4)
            plt.scatter(self.Loaded_M_1_Range, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage, s = 50,  color='r', label = '', alpha=0.5)
            plt.plot(x_temp, Afit*x_temp + Bfit, 'c', label='', linewidth = 2)
            #plt.xlabel('$M_{ex}$', size=22)
            plt.xlabel('$M$', size=22)
            plt.ylabel(r'$\Delta I_{ave} \, (\%)$', size=20)
            
            plt.xticks([0, 0.5, 1.0])
            plt.xlim(-0.01, 1.2)
            
            plt.ylim(-110, 60)
            plt.yticks([-100, -50, 0, 50])
            
            #plt.legend(fontsize=12, loc='upper left',frameon=False)
            





    
            
        
            
            plt.subplot(2,3,4)
            #plt.title('Relative average intensity change', size=18)




    
            plt.subplot(2,3,5, aspect='equal')
            colors = np.array(self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage[::-1])
            plt.scatter(self.Loaded_M_1_Range[::-1], self.Loaded_M_2_Range[::-1], c=colors, s=60, marker='o', linewidth='0', alpha=0.6, cmap=plt.cm.spectral )
            #plt.colorbar(ticks=[-100, -50, 0, 20])
            #plt.colorbar()
            cbar = plt.colorbar(ticks=[-100, -75, -50, -25, 0, 25, 50, 75, 100])
            cbar.set_clim(-100, 100)
            plt.plot([0,2],[0,2], 'r', linewidth=2.5)
            plt.xlabel('$M_{before}$', size=22)
            plt.ylabel('$M_{after}$', size=22)     
            plt.xlim(0.0, 1.21)
            plt.ylim(0.0, 1.21)
            plt.xticks([0, 0.5, 1.0])
            plt.yticks([0, 0.5, 1.0])
            
            
            

          
            
            plt.subplot(2,3,6, aspect='equal')
            #plt.title('M2 vs M1 ', size=10)
            #plt.scatter(self.Loaded_M_1_Range, self.Loaded_M_2_Range, color='black', s = 4 )
            plt.errorbar(self.Loaded_M_1_Range, self.Loaded_M_2_Range, xerr = self.M_error_1_Range, yerr = self.M_error_2_Range, color='black', linestyle="None" )
            plt.plot([0,2],[0,2], 'r')
            plt.xlabel('$M_{before}$', size=22)
            plt.ylabel('$M_{after}$', size=22)
            plt.xlim(0.0, 1.21)
            plt.ylim(0.0, 1.21)
            
            plt.xticks([0, 0.5, 1.0])
            plt.yticks([0, 0.5, 1.0])
            
            
        
        
        
        
            plt.ion()
            plt.show()
        
               
            todayDate = time.strftime("%Y%m%d_%Hh%Mm%Ss")
            
            
            self.fig_temp2.savefig('fig_temp2_' + todayDate + '.png', dpi=500)
            
                 
        
        
        
        
        
        
        
        
        
        
        

            xl = 5
            yl = 4

            plt.figure(figsize=(xl,yl), facecolor='white')
            plt.hist(self.Loaded_M_1_Range, bins=np.arange(0.0, 1.4, 0.1), label = 'Before', alpha=.6, color='c', linewidth=2)
            plt.hist(self.Loaded_M_2_Range, bins=np.arange(0.0, 1.4, 0.1), label = 'After', alpha=.6, color='m', linewidth=2)
            plt.xlabel('$M$', size=22)
            plt.ylabel('Occurence', size=20)
            plt.xticks([0, 0.5, 1.0])
            plt.xlim(0, 1.21)

            #plt.ylim(0, 52) # for CF
            #plt.yticks([0, 10, 20, 30, 40]) # for CF
            plt.ylim(0, 89)  # for toluene
            plt.yticks([0, 20, 40, 60, 80])  # for toluene

            plt.legend(fontsize=12, loc='upper right',frameon=False)
            plt.tight_layout()
            

            ra = 1.259
   
            plt.figure(figsize=(xl*ra, yl*ra), facecolor='white')
            plt.subplot(1,1,1, aspect='equal')
            colors = np.array(self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage[::-1])
            plt.scatter(self.Loaded_M_1_Range[::-1], self.Loaded_M_2_Range[::-1], c=colors, s=60, marker='o', linewidth='0', alpha=0.6, cmap=plt.cm.spectral )
            #plt.colorbar(ticks=[-100, -50, 0, 20])
            #plt.colorbar()
            cbar = plt.colorbar(ticks=[-100, -75, -50, -25, 0, 25, 50, 75, 100])
            cbar.set_clim(-100, 100)
            plt.plot([0,2],[0,2], 'r', linewidth=2.5)
            plt.xlabel('$M_{before}$', size=22)
            plt.ylabel('$M_{after}$', size=22)     
            plt.xlim(0.0, 1.21)
            plt.ylim(0.0, 1.21)
            plt.xticks([0, 0.5, 1.0])
            plt.yticks([0, 0.5, 1.0])
            plt.tight_layout()
            

    
    





        
 #           try:
  
            '''  
            fig_temp3 = plt.figure(figsize = (16,8), facecolor='white')
            fig_temp3.subplots_adjust(top = 0.89, bottom = 0.12, left = 0.1, right = 0.98, wspace = 0.4, hspace = 0.4)
        
            fig_temp3.text(0.02, 0.93, '#3: ' + self.data_3.fileName + '\n' + '#4: ' + self.data_4.fileName + '\nTotal #: ' + str(self.Loaded_Ndata_Range)\
            + ',    BGS_type: ' + self.data_3.BGS_type + ', ' + self.data_4.BGS_type, size=12, color="black", weight='normal')
            
            fig_temp3.text(0.62, 0.93, data_conditions, size=12, color="red", weight='normal')
    



             
            plt.subplot(2,3,1)
            plt.rcParams.update({'mathtext.default': 'regular' })

            #plt.title('Average intensity ', size=18)
            plt.scatter(self.Loaded_M_3_Range, self.Loaded_y0_3_Range, s = 6,  color='red', label = '#3')
            #plt.scatter(self.Loaded_M_2_Range, self.Loaded_y0_2_Range, s = 6,  color='blue', label = '#2.'+fileType_2)
            plt.xlabel('$M_{ex}$', size=22)
            plt.ylabel('$I_{ave}$ (a.u.)', size=20)
            
            plt.xticks([0, 0.5, 1.0])
            plt.xlim(-0.01, 1.21)
            
            
            plt.yticks([0, 2000, 4000, 6000, 8000])
            
            #plt.legend(fontsize=10, loc='upper left')
    
            
        
        
            self.Loaded_Intensity_Fit_y0_Range_Diff_34 = np.array(self.Loaded_y0_4_Range) -  np.array(self.Loaded_y0_3_Range)
            self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage_34 = ( np.array(self.Loaded_y0_4_Range) - np.array(self.Loaded_y0_3_Range) )*100/np.array(self.Loaded_y0_3_Range)
         
            
            plt.subplot(2,3,2)
            #plt.title('Average intensity change', size=18)
            plt.scatter(self.Loaded_M_3_Range, self.Loaded_Intensity_Fit_y0_Range_Diff_34, s = 6,  color='red', label = '#3')
            #plt.scatter(self.Loaded_M_2_Range, self.Loaded_Intensity_Fit_y0_Range_Diff, s = 6,  color='blue', label = '#2.'+fileType_2)
            plt.xlabel('$M_{ex}$', size=22)
            plt.ylabel(r'$\Delta I_{ave} \, (a.u.)$', size=20)
            
            plt.xticks([0, 0.5, 1.0])
            plt.xlim(-0.01, 1.21)
            
            plt.yticks([-6000, -4000, -2000, 0])
            #plt.legend(fontsize=9, loc='upper left')
    
    
    
    

            x_temp = np.array([-0.2, 1.3])
            initial_guess = [0.0, -80.0]
    
            popt, pcov = opt.curve_fit(Linear_func, self.Loaded_M_3_Range, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage_34, p0 = initial_guess)
            #print '\nFIONA popt: ', popt
            #print '\nFIONA pcov: \n', pcov
            Afit_34 = popt[0]
            Bfit_34 = popt[1]
            
            
            AfitErr_34 = np.sqrt(pcov[0,0])               
    
        
            
            plt.subplot(2,3,3)
            #plt.title('Relative average intensity change', size=18)
            plt.scatter(self.Loaded_M_3_Range, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage_34, s = 4,  color='r', label = '')
            plt.plot(x_temp, Afit_34*x_temp + Bfit_34, 'c', label='Fit: slope='+str(np.around(Afit_34,1)) + ' ,   Err: ' + str(np.around(AfitErr_34,1)), linewidth = 2)
            plt.xlabel('$M_{ex}$', size=22)
            plt.ylabel(r'$\Delta I_{ave} \, (\%)$', size=20)
            
            plt.xticks([0, 0.5, 1.0])
            plt.xlim(-0.01, 1.21)
            
            plt.ylim(-110, 60)
            plt.yticks([-100, -50, 0, 50])
            
            plt.legend(fontsize=12, loc='upper left',frameon=False)
            



    
            
            
            plt.subplot(2,3,4, aspect='equal')
            #plt.title('M2 vs M1 ', size=10)
            #plt.scatter(self.Loaded_M_1_Range, self.Loaded_M_2_Range, color='black', s = 4 )
            plt.errorbar(self.Loaded_M_3_Range, self.Loaded_M_4_Range, xerr = self.M_error_3_Range, yerr = self.M_error_4_Range, color='black', linestyle="None" )
            plt.plot([0,2],[0,2], 'r')
            plt.xlabel('$M_{before}$', size=22)
            plt.ylabel('$M_{after}$', size=22)
            plt.xlim(0.0, 1.21)
            plt.ylim(0.0, 1.21)
            
            plt.xticks([0, 0.5, 1.0])
            plt.yticks([0, 0.5, 1.0])
            
            
        



    
            plt.subplot(2,3,5, aspect='equal')
            
            
            colors = np.array(self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage_34)
            plt.scatter(self.Loaded_M_3_Range, self.Loaded_M_4_Range, c=colors, s=50, marker='o', linewidth='0', alpha=0.5 )
            plt.colorbar()
            #cbar = plt.colorbar(ticks=[-100, -50, 0, 50])
            #cbar.set_clim(-100, 50)
            
            
            
            plt.plot([0,2],[0,2], 'r')
            
            plt.xlabel('$M_{before}$', size=22)
            plt.ylabel('$M_{after}$', size=22)     
            
            plt.xlim(0.0, 1.21)
            plt.ylim(0.0, 1.21)
            
            plt.xticks([0, 0.5, 1.0])
            plt.yticks([0, 0.5, 1.0])
            
                
                
       
    
                #plt.tight_layout()
            '''
            
            





































        if self.cb_Data34_Trajectories_fig_yn.GetValue():






            NFeatures = len(self.Loaded_M_1_Range)
      
            NFeaturePerFig = 60.0
           
            NCfig = int(  math.ceil((NFeatures/float(NFeaturePerFig))))
            self.fig_Data34_Trajectories_Loaded = [[] for _ in xrange(NCfig)]
    

            k = 0
            fn = 0
            for n in range(NCfig):
                #print 'n = ', n
                self.fig_Data34_Trajectories_Loaded[n] = plt.figure('Data_3_4_Trajectories_'+ str(n), figsize = (20, 11))
                self.fig_Data34_Trajectories_Loaded[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.36, hspace = 0.4)
                
                self.fig_Data34_Trajectories_Loaded[n].text(0.03, 0.97, 'Total #: ' + str(NFeatures), ha="left", va="bottom", size="small",color="black", weight='normal')
                
                plt.figtext(0.2, 0.94, '3.' + str(self.data_3.fileName) + '\n' + '4.' + str(self.data_4.fileName), ha='left', color='black', weight='normal', size='small')
                
                plt.figtext(0.72, 0.93, data_conditions, size=10, color="red", weight='normal')
                    
    
    
                ColCounter = 0
                RawCounter = 0
                for m in range(int(NFeaturePerFig)):
    
  
                    
                    plt.subplot(8,15, m+1 + RawCounter*15)
                    plt.tick_params(labelsize=7)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_3_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_3_Range[m+fn]))+')'  , fontsize=9, color = 'red', weight='normal')
                    plt.plot((self.Loaded_Intensity_3_Range[m + fn]), color = 'OrangeRed')
                    plt.locator_params(nbins=4)
                
                    


        
                    plt.subplot(8,15, m+1 + RawCounter*15 + 15)
                    plt.tick_params(labelsize=7)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_4_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_4_Range[m+fn]))+')'  , fontsize=9, color = 'blue', weight='normal')
                    plt.plot((self.Loaded_Intensity_4_Range[m + fn]), color = 'darkcyan')
                    plt.locator_params(nbins=4)
                
                    
                    
                    
           
           
           
                  
                    ColCounter += 1
                    if ColCounter%15 == 0:
                        RawCounter += 1
                        
                        
    
                    
                    k += 1
                    if k%int(NFeaturePerFig) ==0:
                        fn += int(NFeaturePerFig)
                    if k == NFeatures:
                        break
          
          
          



































        if self.cb_Data3_Trajectories_fig_yn.GetValue():

            NFeatures = len(self.Loaded_M_1_Range)
            NFeaturePerFig = 90.0
            NCfig = int(  math.ceil((NFeatures/float(NFeaturePerFig))))
            self.fig_Data3_Trajectories_Loaded = [[] for _ in xrange(NCfig)]
    
            x3 = np.linspace(0, 15, 16)
            y3 = np.linspace(0, 15, 16)
            x3, y3 = np.meshgrid(x3, y3)

            k = 0
            fn = 0
            for n in range(NCfig):
                self.fig_Data3_Trajectories_Loaded[n] = plt.figure('Data_3_Trajectories_'+ str(n), figsize = (21, 11))
                self.fig_Data3_Trajectories_Loaded[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.36, hspace = 0.6)
                plt.figtext(0.02, 0.94, '#3: ' + str(self.data_3.fileName) + '\n', ha='left', color='black', weight='normal', size='small')
                plt.figtext(0.62, 0.935, data_conditions, size=10, color="red", weight='normal')
                    
                for m in range(int(NFeaturePerFig)):
                    plt.subplot(6,15, m+1)
                    plt.tick_params(labelsize=7)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_3_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_3_Range[m+fn]))+\
                    ')\nM1:' + str(np.around(self.Loaded_M_1_Range[m+fn],2)) + ', M2:' + str(np.around(self.Loaded_M_2_Range[m+fn],2)), fontsize=8, color = 'green', weight='normal')
                    plt.plot((self.Loaded_Intensity_3_Range[m + fn]), color = 'OrangeRed')
                    plt.locator_params(nbins=4)                
                    
                    k += 1
                    if k%int(NFeaturePerFig) ==0:
                        fn += int(NFeaturePerFig)
                    if k == NFeatures:
                        break
          
                






        if self.cb_Data4_Trajectories_fig_yn.GetValue():

            NFeatures = len(self.Loaded_M_1_Range)
            NFeaturePerFig = 90.0
            NCfig = int(  math.ceil((NFeatures/float(NFeaturePerFig))))
            self.fig_Data4_Trajectories_Loaded = [[] for _ in xrange(NCfig)]
    
            x3 = np.linspace(0, 15, 16)
            y3 = np.linspace(0, 15, 16)
            x3, y3 = np.meshgrid(x3, y3)

            k = 0
            fn = 0
            for n in range(NCfig):
                self.fig_Data4_Trajectories_Loaded[n] = plt.figure('Data_4_Trajectories_'+ str(n), figsize = (21, 11))
                self.fig_Data4_Trajectories_Loaded[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.36, hspace = 0.6)
                plt.figtext(0.02, 0.94, '#4: ' + str(self.data_4.fileName) + '\n', ha='left', color='black', weight='normal', size='small')
                plt.figtext(0.62, 0.935, data_conditions, size=10, color="red", weight='normal')
                    
                for m in range(int(NFeaturePerFig)):
                    plt.subplot(6,15, m+1)
                    plt.tick_params(labelsize=7)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_4_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_4_Range[m+fn]))+\
                    ')\nM1:' + str(np.around(self.Loaded_M_1_Range[m+fn],2)) + ', M2:' + str(np.around(self.Loaded_M_2_Range[m+fn],2)), fontsize=8, color = 'green', weight='normal')
                    plt.plot((self.Loaded_Intensity_4_Range[m + fn]), color = 'OrangeRed')
                    plt.locator_params(nbins=4)                
                    
                    k += 1
                    if k%int(NFeaturePerFig) ==0:
                        fn += int(NFeaturePerFig)
                    if k == NFeatures:
                        break
          
                



































        if self.cb_individual_fig_yn.GetValue():
            pass
        else: return





        if self.cb_include_Data3_trajectories_yn.GetValue() == False and self.cb_include_Data4_trajectories_yn.GetValue() == True:
            MessageBox(self, 'Error: include #3 bleaching trajectories', 'Error')





        if self.cb_include_Data3_trajectories_yn.GetValue() == False and self.cb_include_Data4_trajectories_yn.GetValue() == False:
    
            
            NFeatures = len(self.Loaded_M_1_Range)
      
            NFeaturePerFig = 15.0
           
            NCfig = int(  math.ceil((NFeatures/float(NFeaturePerFig))))
            self.fig_AllDataCenterImages_Loaded = [[] for _ in xrange(NCfig)]
    
            x3 = np.linspace(0, 15, 16)
            y3 = np.linspace(0, 15, 16)
            x3, y3 = np.meshgrid(x3, y3)
            
            #fsn = int(self.Nframe.GetValue()) # frame start number
            k = 0
            fn = 0
            for n in range(NCfig):
                #print 'n = ', n
                self.fig_AllDataCenterImages_Loaded[n] = plt.figure('Raw_Data_Images_M_Loaded'+ str(n), figsize = (18, 9))
                self.fig_AllDataCenterImages_Loaded[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.36, hspace = 0.4)
                self.fig_AllDataCenterImages_Loaded[n].text(0.03, 0.97, 'Total #: ' + str(NFeatures) + ',    BGS_type: ' + self.data_1.BGS_type + ', ' + self.data_2.BGS_type \
                , ha="left", va="bottom", size="small",color="black", weight='normal')
                
                plt.figtext(0.2, 0.94, '1.' + str(self.data_1.fileName) + '\n' + '2.' + str(self.data_2.fileName), ha='left', color='black', weight='normal', size='small')
                
                plt.figtext(0.72, 0.93, data_conditions, size=10, color="red", weight='normal')
                    
    
    
                ColCounter = 0
                RawCounter = 0
                for m in range(int(NFeaturePerFig)):
    
                    
                    imageTemp = np.append(np.ravel(self.Loaded_CenterImages_1_Range[m+fn][0]), np.ravel(self.Loaded_CenterImages_1_Range[m+fn][1]))
                    vminTemp = int(np.mean( np.sort(imageTemp, axis=None)[:5]) )
                    vmaxTemp = int(np.mean( np.sort(imageTemp, axis=None)[-5:]) )
                    
                    plt.subplot(6,15, 3*m+1 + RawCounter*15)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_1_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_1_Range[m+fn]))+')'  , fontsize=9, color = 'red', weight='normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_1_Range[m+fn][0].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    
                    plt.subplot(6,15, 3*m+2 + RawCounter*15)
                    plt.title('T = ' + str(int(np.around(self.Loaded_T_1_Range[m + fn],0)) )  , fontsize=9, color = 'red', weight= 'normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_1_Range[m+fn][1].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                  
                    plt.subplot(6,15, 3*m+3 + RawCounter*15)
                    plt.tick_params(labelsize=7)
                    plt.title( str(int(float(np.around(self.Loaded_phi_1_Range[m + fn],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.Loaded_M_1_Range[m + fn],2) )  , fontsize=9, color = 'red', weight='normal')
                    plt.plot((self.Loaded_Intensity_1_Range[m + fn]), color = 'OrangeRed')
                    plt.plot((self.Loaded_Intensity_noBGS_1_Range[m + fn]), color = 'black')
                    plt.plot((self.Loaded_Intensity_BG_1_Range[m + fn]), color = 'gray')
 
                    if self.cb_fit_1_yn.GetValue():
                        plt.plot((self.Loaded_Intensity_Fit_1_Range[m + fn]), color = 'b')
                        plt.plot((self.Loaded_Intensity_BG_Fit_1_Range[m + fn]), color = 'cyan')
                    plt.locator_params(nbins=4)
                
                    
                    
                    
                    imageTemp = np.append(np.ravel(self.Loaded_CenterImages_2_Range[m+fn][0]), np.ravel(self.Loaded_CenterImages_2_Range[m+fn][1]))
                    vminTemp = int(np.mean( np.sort(imageTemp, axis=None)[:5]) )
                    vmaxTemp = int(np.mean( np.sort(imageTemp, axis=None)[-5:]) )
                    
                    #print '\nm+fn ', m+fn
                    #print 'vminTemp', vminTemp
                    #print 'vmaxTemp ', vmaxTemp
                    
                    plt.subplot(6,15, 3*m+1 + RawCounter*15 + 15)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_2_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_2_Range[m+fn]))+')'  , fontsize=9, color = 'blue', weight='bold')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_2_Range[m+fn][0].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    
                    plt.subplot(6,15, 3*m+2 + RawCounter*15 + 15)
                    plt.title('T = ' + str(int(np.around(self.Loaded_T_2_Range[m + fn],0)) )  , fontsize=9, color = 'blue', weight= 'bold')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_2_Range[m+fn][1].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                  
                    plt.subplot(6,15, 3*m+3 + RawCounter*15 + 15)
                    plt.tick_params(labelsize=7)
                    plt.title( str(int(float(np.around(self.Loaded_phi_2_Range[m + fn],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.Loaded_M_2_Range[m + fn],2) )  , fontsize=9, color = 'blue', weight='bold')
                    plt.plot((self.Loaded_Intensity_2_Range[m + fn]), color = 'OrangeRed')
                    plt.plot((self.Loaded_Intensity_noBGS_2_Range[m + fn]), color = 'black')
                    plt.plot((self.Loaded_Intensity_BG_2_Range[m + fn]), color = 'gray')
                    if self.cb_fit_2_yn.GetValue():
                        plt.plot((self.Loaded_Intensity_Fit_2_Range[m + fn]), color = 'b')
                        plt.plot((self.Loaded_Intensity_BG_Fit_2_Range[m + fn]), color = 'cyan')
                    plt.locator_params(nbins=4)
                
                    
                    
                    
                    
                  
                    ColCounter += 1
                    if ColCounter%5 == 0:
                        RawCounter += 1
                        
                        
    
                    
                    k += 1
                    if k%int(NFeaturePerFig) ==0:
                        fn += int(NFeaturePerFig)
                    if k == NFeatures:
                        break
          
          
          
          
      



        elif self.cb_include_Data3_trajectories_yn.GetValue() == True and self.cb_include_Data4_trajectories_yn.GetValue() == False:
    
            
            NFeatures = len(self.Loaded_M_1_Range)
      
            NFeaturePerFig = 10.0
           
            NCfig = int(  math.ceil((NFeatures/float(NFeaturePerFig))))
            self.fig_AllDataCenterImages_Loaded = [[] for _ in xrange(NCfig)]
    
            x3 = np.linspace(0, 15, 16)
            y3 = np.linspace(0, 15, 16)
            x3, y3 = np.meshgrid(x3, y3)
            
            #fsn = int(self.Nframe.GetValue()) # frame start number
            k = 0
            fn = 0
            for n in range(NCfig):
                #print 'n = ', n
                self.fig_AllDataCenterImages_Loaded[n] = plt.figure('Raw_Data_Images_M_Loaded'+ str(n), figsize = (18, 9))
                self.fig_AllDataCenterImages_Loaded[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.36, hspace = 0.4)
                #self.fig_AllDataCenterImages_Loaded[n].text(0.03, 0.97, 'Total #: ' + str(NFeatures) , ha="left", va="bottom", size="small",color="black", weight='normal')
                self.fig_AllDataCenterImages_Loaded[n].text(0.03, 0.97, 'Total #: ' + str(NFeatures) + ',    BGS_type: ' + self.data_1.BGS_type + ', ' + self.data_2.BGS_type \
                , ha="left", va="bottom", size="small",color="black", weight='normal')
                
                plt.figtext(0.2, 0.94, '1.' + str(self.data_1.fileName) + '\n' + '2.' + str(self.data_2.fileName), ha='left', color='black', weight='normal', size='small')
                plt.figtext(0.72, 0.93, data_conditions, size=10, color="red", weight='normal')
                    
        
    
    
                ColCounter = 0
                RawCounter = 0
                for m in range(int(NFeaturePerFig)):
    
                    
                    imageTemp = np.append(np.ravel(self.Loaded_CenterImages_1_Range[m+fn][0]), np.ravel(self.Loaded_CenterImages_1_Range[m+fn][1]))
                    vminTemp = int(np.mean( np.sort(imageTemp, axis=None)[:5]) )
                    vmaxTemp = int(np.mean( np.sort(imageTemp, axis=None)[-5:]) )
                    
                    plt.subplot(6,15, 3*m+1 + RawCounter*30)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_1_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_1_Range[m+fn]))+')'  , fontsize=9, color = 'red', weight='normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_1_Range[m+fn][0].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    
                    plt.subplot(6,15, 3*m+2 + RawCounter*30)
                    plt.title('T = ' + str(int(np.around(self.Loaded_T_1_Range[m + fn],0)) )  , fontsize=9, color = 'red', weight= 'normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_1_Range[m+fn][1].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                  
                    plt.subplot(6,15, 3*m+3 + RawCounter*30)
                    plt.tick_params(labelsize=7)
                    plt.title( str(int(float(np.around(self.Loaded_phi_1_Range[m + fn],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.Loaded_M_1_Range[m + fn],2) )  , fontsize=9, color = 'red', weight='normal')
                    plt.plot((self.Loaded_Intensity_1_Range[m + fn]), color = 'OrangeRed')
                    if self.cb_fit_1_yn.GetValue():
                        plt.plot((self.Loaded_Intensity_Fit_1_Range[m + fn]), color = 'b')
                    plt.locator_params(nbins=4)
                
                    
                    




                    imageTemp = np.append(np.ravel(self.Loaded_CenterImages_3_Range[m+fn][0]), np.ravel(self.Loaded_CenterImages_3_Range[m+fn][1]))
                    vminTemp = int(np.mean( np.sort(imageTemp, axis=None)[:5]) )
                    vmaxTemp = int(np.mean( np.sort(imageTemp, axis=None)[-5:]) )
                    
                    plt.subplot(6,15, 3*m+1 + RawCounter*30 + 15)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_3_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_3_Range[m+fn]))+')'  , fontsize=9, color = 'green', weight='bold')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_3_Range[m+fn][0].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    
                    plt.subplot(6,15, 3*m+2 + RawCounter*30 + 15)
                    if self.cb_include_Data3_fit_yn.GetValue():
                        plt.title('T = ' + str(int(np.around(self.Loaded_T_3_Range[m + fn],0)) )  , fontsize=9, color = 'blue', weight= 'normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_3_Range[m+fn][1].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                  
                    plt.subplot(6,15, 3*m+3 + RawCounter*30 + 15)
                    plt.tick_params(labelsize=7)
                    plt.plot((self.Loaded_Intensity_3_Range[m + fn]), color = 'green')
                    if self.cb_include_Data3_fit_yn.GetValue():
                        plt.plot((self.Loaded_Intensity_Fit_3_Range[m + fn]), color = 'b')
                        plt.title( str(int(float(np.around(self.Loaded_phi_3_Range[m + fn],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.Loaded_M_3_Range[m + fn],2) )  , fontsize=9, color = 'blue', weight='normal')
                    plt.locator_params(nbins=4)
                
                    
                    
                    




                    
                    imageTemp = np.append(np.ravel(self.Loaded_CenterImages_2_Range[m+fn][0]), np.ravel(self.Loaded_CenterImages_2_Range[m+fn][1]))
                    vminTemp = int(np.mean( np.sort(imageTemp, axis=None)[:5]) )
                    vmaxTemp = int(np.mean( np.sort(imageTemp, axis=None)[-5:]) )
                    
                    plt.subplot(6,15, 3*m+1 + RawCounter*30 + 30)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_2_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_2_Range[m+fn]))+')'  , fontsize=9, color = 'blue', weight='bold')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_2_Range[m+fn][0].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    
                    plt.subplot(6,15, 3*m+2 + RawCounter*30 + 30)
                    plt.title('T = ' + str(int(np.around(self.Loaded_T_2_Range[m + fn],0)) )  , fontsize=9, color = 'blue', weight= 'bold')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_2_Range[m+fn][1].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                  
                    plt.subplot(6,15, 3*m+3 + RawCounter*30 + 30)
                    plt.tick_params(labelsize=7)
                    plt.title( str(int(float(np.around(self.Loaded_phi_2_Range[m + fn],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.Loaded_M_2_Range[m + fn],2) )  , fontsize=9, color = 'blue', weight='bold')
                    plt.plot((self.Loaded_Intensity_2_Range[m + fn]), color = 'OrangeRed')
                    if self.cb_fit_2_yn.GetValue():
                        plt.plot((self.Loaded_Intensity_Fit_2_Range[m + fn]), color = 'b')
                    plt.locator_params(nbins=4)
                
                    
                    
                    
                  
                    ColCounter += 1
                    if ColCounter%5 == 0:
                        RawCounter += 1
                        
                        
    
                    
                    k += 1
                    if k%int(NFeaturePerFig) ==0:
                        fn += int(NFeaturePerFig)
                    if k == NFeatures:
                        break
          
                











        elif self.cb_include_Data3_trajectories_yn.GetValue() == True and self.cb_include_Data4_trajectories_yn.GetValue() == True:
    
            
            NFeatures = len(self.Loaded_M_1_Range)
      
            NFeaturePerFig = 10.0
           
            NCfig = int(  math.ceil((NFeatures/float(NFeaturePerFig))))
            self.fig_AllDataCenterImages_Loaded = [[] for _ in xrange(NCfig)]
    
            x3 = np.linspace(0, 15, 16)
            y3 = np.linspace(0, 15, 16)
            x3, y3 = np.meshgrid(x3, y3)
            
            #fsn = int(self.Nframe.GetValue()) # frame start number
            k = 0
            fn = 0
            for n in range(NCfig):
                #print 'n = ', n
                self.fig_AllDataCenterImages_Loaded[n] = plt.figure('Raw_Data_Images_M_Loaded'+ str(n), figsize = (18, 9.2))
                self.fig_AllDataCenterImages_Loaded[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.40, hspace = 0.46)
                #self.fig_AllDataCenterImages_Loaded[n].text(0.03, 0.97, 'Total #: ' + str(NFeatures) , ha="left", va="bottom", size="small",color="black", weight='normal')
                self.fig_AllDataCenterImages_Loaded[n].text(0.03, 0.97, 'Total #: ' + str(NFeatures) + ',    BGS_type: ' + self.data_1.BGS_type + ', ' + self.data_2.BGS_type \
                , ha="left", va="bottom", size="small",color="black", weight='normal')
                
                
                plt.figtext(0.2, 0.94, '1.' + str(self.data_1.fileName) + '\n' + '2.' + str(self.data_2.fileName), ha='left', color='black', weight='normal', size='small')
                plt.figtext(0.72, 0.93, data_conditions, size=10, color="red", weight='normal')
    
    
                ColCounter = 0
                RawCounter = 0
                for m in range(int(NFeaturePerFig)):
    
                    
                    imageTemp = np.append(np.ravel(self.Loaded_CenterImages_1_Range[m+fn][0]), np.ravel(self.Loaded_CenterImages_1_Range[m+fn][1]))
                    vminTemp = int(np.mean( np.sort(imageTemp, axis=None)[:5]) )
                    vmaxTemp = int(np.mean( np.sort(imageTemp, axis=None)[-5:]) )
                    
                    plt.subplot(8,15, 3*m+1 + RawCounter*45)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_1_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_1_Range[m+fn]))+')'  , fontsize=9, color = 'red', weight='normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_1_Range[m+fn][0].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    
                    plt.subplot(8,15, 3*m+2 + RawCounter*45)
                    plt.title('T = ' + str(int(np.around(self.Loaded_T_1_Range[m + fn],0)) )  , fontsize=9, color = 'red', weight= 'normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_1_Range[m+fn][1].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                  
                    plt.subplot(8,15, 3*m+3 + RawCounter*45)
                    plt.tick_params(labelsize=7)
                    plt.title( str(int(float(np.around(self.Loaded_phi_1_Range[m + fn],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.Loaded_M_1_Range[m + fn],2) )  , fontsize=9, color = 'red', weight='normal')
                    plt.plot((self.Loaded_Intensity_1_Range[m + fn]), color = 'OrangeRed')
                    if self.cb_fit_1_yn.GetValue():
                        plt.plot((self.Loaded_Intensity_Fit_1_Range[m + fn]), color = 'b')
                    plt.locator_params(nbins=4)
                
                    
                    




                    imageTemp = np.append(np.ravel(self.Loaded_CenterImages_3_Range[m+fn][0]), np.ravel(self.Loaded_CenterImages_3_Range[m+fn][1]))
                    vminTemp = int(np.mean( np.sort(imageTemp, axis=None)[:5]) )
                    vmaxTemp = int(np.mean( np.sort(imageTemp, axis=None)[-5:]) )
                    
                    plt.subplot(8,15, 3*m+1 + RawCounter*45 + 15)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_3_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_3_Range[m+fn]))+')'  , fontsize=9, color = 'green', weight='normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_3_Range[m+fn][0].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    
                    plt.subplot(8,15, 3*m+2 + RawCounter*45 + 15)
                    if self.cb_include_Data3_fit_yn.GetValue():
                        plt.title('T = ' + str(int(np.around(self.Loaded_T_3_Range[m + fn],0)) )  , fontsize=9, color = 'blue', weight= 'normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_3_Range[m+fn][1].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                  
                    plt.subplot(8,15, 3*m+3 + RawCounter*45 + 15)
                    plt.tick_params(labelsize=7)
                    plt.plot((self.Loaded_Intensity_3_Range[m + fn]), color = 'green')
                    if self.cb_include_Data3_fit_yn.GetValue():
                        plt.plot((self.Loaded_Intensity_Fit_3_Range[m + fn]), color = 'b')
                        plt.title( str(int(float(np.around(self.Loaded_phi_3_Range[m + fn],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.Loaded_M_3_Range[m + fn],2) )  , fontsize=9, color = 'blue', weight='normal')
                    plt.locator_params(nbins=4)
                
                    




                    
                    imageTemp = np.append(np.ravel(self.Loaded_CenterImages_2_Range[m+fn][0]), np.ravel(self.Loaded_CenterImages_2_Range[m+fn][1]))
                    vminTemp = int(np.mean( np.sort(imageTemp, axis=None)[:5]) )
                    vmaxTemp = int(np.mean( np.sort(imageTemp, axis=None)[-5:]) )
                    
                    plt.subplot(8,15, 3*m+1 + RawCounter*45 + 30)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_2_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_2_Range[m+fn]))+')'  , fontsize=9, color = 'blue', weight='normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_2_Range[m+fn][0].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    
                    plt.subplot(8,15, 3*m+2 + RawCounter*45 + 30)
                    plt.title('T = ' + str(int(np.around(self.Loaded_T_2_Range[m + fn],0)) )  , fontsize=9, color = 'blue', weight= 'normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_2_Range[m+fn][1].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                  
                    plt.subplot(8,15, 3*m+3 + RawCounter*45 + 30)
                    plt.tick_params(labelsize=7)
                    plt.title( str(int(float(np.around(self.Loaded_phi_2_Range[m + fn],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.Loaded_M_2_Range[m + fn],2) )  , fontsize=9, color = 'blue', weight='normal')
                    plt.plot((self.Loaded_Intensity_2_Range[m + fn]), color = 'OrangeRed')
                    if self.cb_fit_2_yn.GetValue():
                        plt.plot((self.Loaded_Intensity_Fit_2_Range[m + fn]), color = 'b')
                    plt.locator_params(nbins=4)
                
                    
                    
                    
                    
                    
                    


                    imageTemp = np.append(np.ravel(self.Loaded_CenterImages_4_Range[m+fn][0]), np.ravel(self.Loaded_CenterImages_4_Range[m+fn][1]))
                    vminTemp = int(np.mean( np.sort(imageTemp, axis=None)[:5]) )
                    vmaxTemp = int(np.mean( np.sort(imageTemp, axis=None)[-5:]) )
                    
                    plt.subplot(8,15, 3*m+1 + RawCounter*45 + 45)
                    plt.title('# ' + str(m+fn) + ', (' + str(int(self.Loaded_ImageJxpos_4_Range[m+fn])) + ',' + str(int(self.Loaded_ImageJypos_4_Range[m+fn]))+')'  , fontsize=9, color = 'green', weight='normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_4_Range[m+fn][0].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    
                    plt.subplot(8,15, 3*m+2 + RawCounter*45 + 45)
                    if self.cb_include_Data4_fit_yn.GetValue():
                        plt.title('T = ' + str(int(np.around(self.Loaded_T_4_Range[m + fn],0)) )  , fontsize=9, color = 'blue', weight= 'normal')
                    plt.tick_params(labelsize=7)
                    plt.imshow(self.Loaded_CenterImages_4_Range[m+fn][1].reshape(15, 15), interpolation = 'None', vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='bottom',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                  
                    plt.subplot(8,15, 3*m+3 + RawCounter*45 + 45)
                    plt.tick_params(labelsize=7)
                    plt.plot((self.Loaded_Intensity_4_Range[m + fn]), color = 'limegreen')
                    if self.cb_include_Data4_fit_yn.GetValue():
                        plt.plot((self.Loaded_Intensity_Fit_4_Range[m + fn]), color = 'b')
                        plt.title( str(int(float(np.around(self.Loaded_phi_4_Range[m + fn],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.Loaded_M_4_Range[m + fn],2) )  , fontsize=9, color = 'blue', weight='normal')
                    plt.locator_params(nbins=4)
                
                    
                    
                    
                    
                    
                    
                  
                    ColCounter += 1
                    if ColCounter%5 == 0:
                        RawCounter += 1
                        
                        
    
                    
                    k += 1
                    if k%int(NFeaturePerFig) ==0:
                        fn += int(NFeaturePerFig)
                    if k == NFeatures:
                        break
          
                
                
        

        plt.ion()
        plt.show()
        
        

                
                
                




                
                
                
                
                
        




    def SaveData(self, event):
        print ' '
        prefix = self.FileName.GetValue()


        if prefix != '':
            prefix += '_'




        todayDate = time.strftime("%Y%m%d_%Hh%Mm")
        
        figFileName = prefix
        
        
        print "\nToday's Date: ", time.strftime("%Y-%m-%d_%Hh%Mm")
        
        
    
        self.fig_DataPlots.savefig(self.data_1.fileDirectory + '\\' + figFileName + '_' + todayDate + '_Data_1' +'.png')
        plt.close(self.fig_DataPlots)


        
        try:
            self.fig_DataPlots_2.savefig(self.data_1.fileDirectory + '\\' + figFileName + '_' + todayDate + '_Data_2' +'.png')
            plt.close(self.fig_DataPlots_2)
        except: pass



        
        try:
            for n in range(len(self.fig_AllDataCenterImages_Loaded)):    
                self.fig_AllDataCenterImages_Loaded[n].savefig(self.data_1.fileDirectory + '\\' + figFileName + '_' + todayDate + '_DataImages_' + str(n) + '.png')
                plt.close(self.fig_AllDataCenterImages_Loaded[n])
                print 'Raw Frames figures are saved'
        except:
            print 'no individual figures '
            pass




        

  
  
        if self.cb_Data34_Trajectories_fig_yn.GetValue():
            try:
                for n in range(len(self.fig_Data34_Trajectories_Loaded)):    
                    self.fig_Data34_Trajectories_Loaded[n].savefig(self.data_1.fileDirectory + '\\' + figFileName + '_' + todayDate + '_Data34_Trajectories_' + str(n) + '.png')
                    plt.close(self.fig_Data34_Trajectories_Loaded[n])
                    print 'fig_Data34_Trajectories_Loaded are saved n: ' , n
            except:
                print 'no fig_Data34_Trajectories_Loaded '
                pass

      

  
  
        if self.cb_Data3_Trajectories_fig_yn.GetValue():
            try:
                for n in range(len(self.fig_Data3_Trajectories_Loaded)):    
                    self.fig_Data3_Trajectories_Loaded[n].savefig(self.data_1.fileDirectory + '\\' + figFileName + '_' + todayDate + '_Data3_Trajectories_' + str(n) + '.png')
                    plt.close(self.fig_Data3_Trajectories_Loaded[n])
                    print 'fig_Data3_Trajectories_Loaded are saved'
            except:
                print 'no fig_Data3_Trajectories_Loaded '
                pass

      



  
        if self.cb_Data4_Trajectories_fig_yn.GetValue():
            try:
                for n in range(len(self.fig_Data4_Trajectories_Loaded)):    
                    self.fig_Data4_Trajectories_Loaded[n].savefig(self.data_1.fileDirectory + '\\' + figFileName + '_' + todayDate + '_Data4_Trajectories_' + str(n) + '.png')
                    plt.close(self.fig_Data4_Trajectories_Loaded[n])
                    print 'fig_Data4_Trajectories_Loaded are saved'
            except:
                print 'no fig_Data4_Trajectories_Loaded '
                pass

      
        
        
        
        self.SaveData_text.SetLabel('Saved ')


























     
    def SelectData(self, event):
        
        
        self.SelectData_frame = wx.Frame(self, title='SelectData', size=(950, 1000), pos=(800,0))
        #self.SelectData_frame.Show()
        
        self.SelectData_window = wx.ScrolledWindow(self.SelectData_frame, -1)
        self.SelectData_window.SetScrollbars(1, 1, 800, int(len(self.Loaded_ImageJxpos_1_Range)/10)*25 + 400)
        
        self.SelectData_frame.Show()



        wx.StaticText(self.SelectData_window, -1, "# of Letters to be removed in each file name ", pos=(170, 5))
        wx.StaticText(self.SelectData_window, -1, "At Beginning: ", pos=(170, 25))
        self.NofLettersToRemoveFinleNameForSelected_Beginning = wx.TextCtrl(self.SelectData_window, -1, '10', pos=(250, 23), size=(30, -1))
        wx.StaticText(self.SelectData_window, -1, "At End: ", pos=(320, 25))
        self.NofLettersToRemoveFinleNameForSelected_End = wx.TextCtrl(self.SelectData_window, -1, '35', pos=(360, 23), size=(30, -1))
        
        

         
        wx.StaticText(self.SelectData_window, -1, "Data File Prefix: ", pos=(170, 55))
        self.FilePrefixForSelected = wx.TextCtrl(self.SelectData_window, -1, str(FilePrefixInput), pos=(260, 50), size=(150, -1))
        
        
        btn_ShowFiguresSelectedData = wx.Button(self.SelectData_window, pos=(50,30), label="Plot Selected Data")
        btn_ShowFiguresSelectedData.Bind(wx.EVT_BUTTON, self.SelectData_ShowFigures)
        

        
        self.btn_SaveSelectedData = wx.Button(self.SelectData_window, pos=(470,30), label="Save Selected Data #1 #2")
        self.btn_SaveSelectedData.Bind(wx.EVT_BUTTON, self.SaveSelectedData12)
        self.btn_SaveSelectedData.Hide()
        self.text_SaveSelectedData = wx.StaticText(self.SelectData_window, -1, "...", pos=(470, 60))
        

        self.cb_SaveFigs_Selected_yn = wx.CheckBox(self.SelectData_window, -1, 'Figures', (650, 10))
        self.cb_SaveFigs_Selected_yn.SetValue(True)
        self.cb_SaveImageJ_xyData_Selected_yn = wx.CheckBox(self.SelectData_window, -1, 'ImageJ x,y data & Fit: A, T, phase, y0, M', (650, 30))
        self.cb_SaveImageJ_xyData_Selected_yn.SetValue(False)
        self.cb_SaveIntensityData_Selected_yn = wx.CheckBox(self.SelectData_window, -1, 'Intensity Trajectory Data for all', (650, 50))
        self.cb_SaveIntensityData_Selected_yn.SetValue(False)
        self.cb_SaveIntensityFitData_Selected_yn = wx.CheckBox(self.SelectData_window, -1, 'Fit Trajectory Data for all', (650, 70))
        self.cb_SaveIntensityFitData_Selected_yn.SetValue(False)
        self.cb_Complete_Data_Selected_yn = wx.CheckBox(self.SelectData_window, -1, 'Complete Data File', (650, 90))
        self.cb_Complete_Data_Selected_yn.SetValue(True)






        self.cb_Save_Data1_yn = wx.CheckBox(self.SelectData_window, -1, '#1', (150, 123))
        self.cb_Save_Data1_yn.SetValue(True)
        self.cb_Save_Data2_yn = wx.CheckBox(self.SelectData_window, -1, '#2', (200, 123))
        self.cb_Save_Data2_yn.SetValue(True)
        self.cb_Save_Data3_yn = wx.CheckBox(self.SelectData_window, -1, '#3', (250, 123))
        self.cb_Save_Data3_yn.SetValue(True)
        self.cb_Save_Data4_yn = wx.CheckBox(self.SelectData_window, -1, '#4', (300, 123))
        self.cb_Save_Data4_yn.SetValue(True)
        self.cb_Save_Data5_yn = wx.CheckBox(self.SelectData_window, -1, '#5', (350, 123))
        self.cb_Save_Data5_yn.SetValue(True)
        self.cb_Save_Data6_yn = wx.CheckBox(self.SelectData_window, -1, '#6', (400, 123))
        self.cb_Save_Data6_yn.SetValue(True)





        self.btn_SaveSelectedDataForAll = wx.Button(self.SelectData_window, pos=(470,120), label="Save Selected Data For All - Complete Data Files Only")
        self.btn_SaveSelectedDataForAll.Bind(wx.EVT_BUTTON, self.SaveSelectedDataForAll)
        self.btn_SaveSelectedDataForAll.Hide()
        self.text_SaveSelectedDataForAll = wx.StaticText(self.SelectData_window, -1, "...", pos=(470, 150))
        





        
        
        wx.StaticText(self.SelectData_window, -1, 'Select Data', pos=(50, 203)) 
                
        #self.xpos_ImageJ = [[] for _ in xrange(150)]
        self.cb_Data_yn = [[] for _ in xrange(len(self.Loaded_ImageJxpos_1_Range))]
        rn = 0
        cn = 0
        countertemp = 0
        for n in range(len(self.Loaded_ImageJxpos_1_Range)):
            self.cb_Data_yn[n] = wx.CheckBox(self.SelectData_window, -1, '# '+str(n), (50 + 70*cn, 235 + 25*rn))
            self.cb_Data_yn[n].SetValue(True)
            countertemp += 1
            cn += 1
            if countertemp % 10 == 0:
                rn += 1
                cn = 0
            
        self.btn_SelectAll = wx.Button(self.SelectData_window, pos=(150,180), label="Select All")
        self.btn_SelectAll.Bind(wx.EVT_BUTTON, self.SelectDataAll)
        self.btn_SelectAll.Show()
        
        self.btn_SelectNone = wx.Button(self.SelectData_window, pos=(250,180), label="Deselect All")
        self.btn_SelectNone.Bind(wx.EVT_BUTTON, self.SelectDataNone)
        self.btn_SelectNone.Show()
        
        
        self.btn_SaveSelection = wx.Button(self.SelectData_window, pos=(440,180), label="Save Selections")
        self.btn_SaveSelection.Bind(wx.EVT_BUTTON, self.SaveSelection)
        self.btn_SaveSelection.Show()
        
        self.btn_LoadSelection = wx.Button(self.SelectData_window, pos=(550,180), label="Load Selections")
        self.btn_LoadSelection.Bind(wx.EVT_BUTTON, self.LoadSelection)
        self.btn_LoadSelection.Show()
        


        
        
        
    def SelectDataAll(self, event):
        print '\nSelect All '
        for n in range(len(self.cb_Data_yn)):
            self.cb_Data_yn[n].SetValue(True)
    def SelectDataNone(self, event):
        print '\nDeselect All '
        for n in range(len(self.cb_Data_yn)):
            self.cb_Data_yn[n].SetValue(False)
            
            
            
        
        

    def SelectData_ShowFigures(self, event):
        print '\n SelectData_ShowFigures '
        
        try:
            plt.close(self.fig_DataPlots_Selected)
        except: pass
        
        



        
        
        
        
        self.x_ImageJ_Selected_1 = []
        self.y_ImageJ_Selected_1 = []
        self.centerImage_Selected_1 = []

        self.IntensityTrajectoryChosen_Selected_1 = []
        self.IntensityTrajectoryChosenBS_Selected_1 = []
        #self.IntensityTrajectoryBGChosen_Selected_1 = []
        self.IntensityTrajectoryChosenBS_fit_Selected_1 = []
        self.backgroundIntensityTrajectoryChosen_Selected_1 = []
        self.backgroundIntensityTrajectoryChosen_Fit_Selected_1 = []
        
        self.IntensityTrajectoryChosenBS_fit_A_Selected_1 = []
        self.IntensityTrajectoryChosenBS_fit_T_Selected_1 = []
        self.IntensityTrajectoryChosenBS_fit_phi_Selected_1 = []
        self.IntensityTrajectoryChosenBS_fit_y0_Selected_1 = []
        self.IntensityTrajectoryChosenBS_fit_M_Selected_1 = []
        
        
        self.x_ImageJ_Selected_2 = []
        self.y_ImageJ_Selected_2 = []
        self.centerImage_Selected_2 = []

        self.IntensityTrajectoryChosen_Selected_2 = []
        self.IntensityTrajectoryChosenBS_Selected_2 = []
        #self.IntensityTrajectoryBGChosen_Selected_2 = []
        self.IntensityTrajectoryChosenBS_fit_Selected_2 = []
        self.backgroundIntensityTrajectoryChosen_Selected_2 = []
        self.backgroundIntensityTrajectoryChosen_Fit_Selected_2 = []

        
        self.IntensityTrajectoryChosenBS_fit_Selected_2 = []
        self.IntensityTrajectoryChosenBS_fit_A_Selected_2 = []
        self.IntensityTrajectoryChosenBS_fit_T_Selected_2 = []
        self.IntensityTrajectoryChosenBS_fit_phi_Selected_2 = []
        self.IntensityTrajectoryChosenBS_fit_y0_Selected_2 = []
        self.IntensityTrajectoryChosenBS_fit_M_Selected_2 = []
        



        self.NFeatures_Selected = 0
        
        for n in range(len(self.cb_Data_yn)):
            if self.cb_Data_yn[n].GetValue():
                self.NFeatures_Selected += 1
                
                self.x_ImageJ_Selected_1.append(self.Loaded_ImageJxpos_1_Range[n])
                self.y_ImageJ_Selected_1.append(self.Loaded_ImageJypos_1_Range[n])
                self.centerImage_Selected_1.append(self.Loaded_CenterImages_1_Range[n])
                
                self.IntensityTrajectoryChosen_Selected_1.append(self.Loaded_Intensity_noBGS_1_Range[n])
                self.IntensityTrajectoryChosenBS_Selected_1.append(self.Loaded_Intensity_1_Range[n])
                #self.IntensityTrajectoryBGChosen_Selected_1.append(self.Loaded_Intensity_BG_1[n])
                self.IntensityTrajectoryChosenBS_fit_Selected_1.append(self.Loaded_Intensity_Fit_1_Range[n])
                self.backgroundIntensityTrajectoryChosen_Selected_1.append(self.Loaded_Intensity_BG_1_Range[n])
                self.backgroundIntensityTrajectoryChosen_Fit_Selected_1.append(self.Loaded_Intensity_BG_Fit_1_Range[n])
                     
                self.IntensityTrajectoryChosenBS_fit_A_Selected_1.append(self.Loaded_A_1_Range[n])
                self.IntensityTrajectoryChosenBS_fit_T_Selected_1.append(self.Loaded_T_1_Range[n])
                self.IntensityTrajectoryChosenBS_fit_phi_Selected_1.append(self.Loaded_phi_1_Range[n])
                self.IntensityTrajectoryChosenBS_fit_y0_Selected_1.append(self.Loaded_y0_1_Range[n])
                self.IntensityTrajectoryChosenBS_fit_M_Selected_1.append(self.Loaded_M_1_Range[n])
                
                        


        
                self.x_ImageJ_Selected_2.append(self.Loaded_ImageJxpos_2_Range[n])
                self.y_ImageJ_Selected_2.append(self.Loaded_ImageJypos_2_Range[n])
              
                self.IntensityTrajectoryChosenBS_fit_A_Selected_2.append(self.Loaded_A_2_Range[n])
                self.IntensityTrajectoryChosenBS_fit_T_Selected_2.append(self.Loaded_T_2_Range[n])
                self.IntensityTrajectoryChosenBS_fit_phi_Selected_2.append(self.Loaded_phi_2_Range[n])
                self.IntensityTrajectoryChosenBS_fit_y0_Selected_2.append(self.Loaded_y0_2_Range[n])
                self.IntensityTrajectoryChosenBS_fit_M_Selected_2.append(self.Loaded_M_2_Range[n])
 
                self.IntensityTrajectoryChosen_Selected_2.append(self.Loaded_Intensity_noBGS_2_Range[n])
                self.IntensityTrajectoryChosenBS_Selected_2.append(self.Loaded_Intensity_2_Range[n])
                #self.IntensityTrajectoryBGChosen_Selected_2.append(self.Loaded_Intensity_BG_2[n])
                self.IntensityTrajectoryChosenBS_fit_Selected_2.append(self.Loaded_Intensity_Fit_2_Range[n])
                self.backgroundIntensityTrajectoryChosen_Selected_2.append(self.Loaded_Intensity_BG_2_Range[n])
                self.backgroundIntensityTrajectoryChosen_Fit_Selected_2.append(self.Loaded_Intensity_BG_Fit_2_Range[n])
                self.centerImage_Selected_2.append(self.Loaded_CenterImages_2_Range[n])
                          
                
   
                
                
                                


        if self.radio_file_1.GetSelection() == 0:
            fileType_1 = 'Ext'
        elif self.radio_file_1.GetSelection() == 1:
            fileType_1 = 'Emi'
        elif self.radio_file_1.GetSelection() == 2:
            fileType_1 = ' '
            
        
        if self.radio_file_2.GetSelection() == 0:
            fileType_2 = 'Ext'
        elif self.radio_file_2.GetSelection() == 1:
            fileType_2 = 'Emi'
        elif self.radio_file_2.GetSelection() == 2:
            fileType_2 = ' '
            


        self.Loaded_Intensity_Fit_1_Range_CB_Selected = np.array(self.IntensityTrajectoryChosenBS_fit_y0_Selected_1) - np.abs(self.IntensityTrajectoryChosenBS_fit_A_Selected_1)  # CB == constant background
        self.Loaded_Intensity_Fit_2_Range_CB_Selected = np.array(self.IntensityTrajectoryChosenBS_fit_y0_Selected_2) - np.abs(self.IntensityTrajectoryChosenBS_fit_A_Selected_2)
        
        self.Loaded_Intensity_Fit_1_Range_Max_Selected = np.array(self.IntensityTrajectoryChosenBS_fit_y0_Selected_1) + np.abs(self.IntensityTrajectoryChosenBS_fit_A_Selected_1)  # CB == constant background
        self.Loaded_Intensity_Fit_2_Range_Max_Selected = np.array(self.IntensityTrajectoryChosenBS_fit_y0_Selected_2) + np.abs(self.IntensityTrajectoryChosenBS_fit_A_Selected_2)
        
        
        self.Loaded_Intensity_Fit_y0_Range_Diff_Selected = np.array(self.IntensityTrajectoryChosenBS_fit_y0_Selected_2) -  np.array(self.IntensityTrajectoryChosenBS_fit_y0_Selected_1)
        self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage_Selected = ( np.array(self.IntensityTrajectoryChosenBS_fit_y0_Selected_2) - np.array(self.IntensityTrajectoryChosenBS_fit_y0_Selected_1) )*100/np.array(self.IntensityTrajectoryChosenBS_fit_y0_Selected_1)
        
        
        self.Loaded_Intensity_Fit_Range_CB_Diff_Percentage_Selected = (self.Loaded_Intensity_Fit_2_Range_CB_Selected - self.Loaded_Intensity_Fit_1_Range_CB_Selected)*100/self.Loaded_Intensity_Fit_1_Range_CB_Selected
        self.Loaded_Intensity_Fit_Range_A_Diff_Percentage_Selected = ( np.array(self.IntensityTrajectoryChosenBS_fit_A_Selected_2) - np.array(self.IntensityTrajectoryChosenBS_fit_A_Selected_1) )*100/np.array(self.IntensityTrajectoryChosenBS_fit_A_Selected_1)
        
            

        phi_diff_Selected = np.array(self.IntensityTrajectoryChosenBS_fit_phi_Selected_2) - np.array(self.IntensityTrajectoryChosenBS_fit_phi_Selected_1)
        phi_diff_abs_Selected = np.abs(phi_diff_Selected)
        
        for n in range(len(phi_diff_abs_Selected)):
            if phi_diff_abs_Selected[n] >= 90:
                phi_diff_abs_Selected[n] = 180 - phi_diff_abs_Selected[n]
                
        
        if self.cb_90deg_offset_yn.GetValue():
            phi_diff_abs_Selected = 90 - phi_diff_abs_Selected
            phi_diff_abs_text = ' 90 deg shift applied'
        else:
            phi_diff_abs_text = ' '
            
        
        M_diff_Selected = np.array(self.IntensityTrajectoryChosenBS_fit_M_Selected_2) - np.array(self.IntensityTrajectoryChosenBS_fit_M_Selected_1)
        
        Ndata_Selected = len(M_diff_Selected)




        self.fig_DataPlots_Selected = plt.figure('Data_Plots_Selected', figsize = (24,12))
        #plt.clf()
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(710,30,1500, 800)
        self.fig_DataPlots_Selected.subplots_adjust(top = 0.89, bottom = 0.07, left = 0.05, right = 0.98, wspace = 0.3, hspace = 0.5)
        self.fig_DataPlots_Selected.text(0.02, 0.93, '#1.' + self.data_1.fileName + '\n' + '#2.' + self.data_2.fileName + '\nTotal #: ' + str(Ndata_Selected), size=12, color="black", weight='normal')
        #self.Data_Plots.text(0.02, 0.97, TotalFeatureN_text +',       ' + self.Corrections, ha="left", va="bottom", size=9,color="black", weight='normal')
        #plt.figtext(0.5, 0.92, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly+ '\n' + self.AnalysisConditions ), ha='center', color='black', weight='normal', size=9)

        #plt.title('Total # ' + str(len(IntensityTrajectoryChosenBS_fit_M)))
        plt.subplot(3,4,1)
                   
        plt.hist(self.IntensityTrajectoryChosenBS_fit_M_Selected_1, bins=np.arange(0.0, 1.3, 0.1), label = '#1.'+fileType_1, alpha=1.0, color='red')
        plt.hist(self.IntensityTrajectoryChosenBS_fit_M_Selected_2, bins=np.arange(0.0, 1.3, 0.1), label = '#2.'+fileType_2, alpha=0.3, color='blue')
        plt.xlabel('Modulation Depth (M)', size=11)
        plt.xlim(0.0, 1.21)
        plt.legend(fontsize=10, loc='upper left')
        
        
        
        plt.subplot(3,4,2)
        plt.title('t ', size=10)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_1, phi_diff_abs_Selected, color='red', label = '#1.'+fileType_1)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_2, phi_diff_abs_Selected, color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('phi diff ABS (M1, M2)', size=11)
        plt.legend(fontsize=10, loc='upper left')
        plt.xlim(-0.1, 1.21)
        
        
        
        plt.subplot(3,4,3, aspect='equal')
        plt.title('M2 vs M1 ', size=10)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_1, self.IntensityTrajectoryChosenBS_fit_M_Selected_2, color='black' )
        plt.plot([0,2],[0,2], 'r')
        plt.xlabel('M1', size=11)
        plt.ylabel('M2', size=11)
        plt.xlim(0.0, 1.21)
        plt.ylim(0.0, 1.21)
        
        

        plt.subplot(3,4,4)
        plt.title('Fit: average intensity change', size=11)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_1, self.Loaded_Intensity_Fit_y0_Range_Diff_Selected, s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_2, self.Loaded_Intensity_Fit_y0_Range_Diff_Selected, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('Fit y0 diff ABS (M1, M2)', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')

        
        
        

        plt.subplot(3,4,5)
        plt.title('phi abs difference,' + phi_diff_abs_text, size=11)
        plt.hist(phi_diff_abs_Selected, bins=np.arange(0,100,10))
        plt.xlabel('|phi2 - phi1| (deg)', size=11)
        plt.xlim(0,100)



        plt.subplot(3,4,6)
        plt.title('M difference', size=11)
        plt.hist(M_diff_Selected, bins=np.arange(-1,1.1,0.1))
        plt.xlabel('M2 - M1', size=11)
                
        
        
        
        plt.subplot(3,4,7)
        plt.title('Fit: average intensity y0', size=11)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_1, self.IntensityTrajectoryChosenBS_fit_y0_Selected_1, s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_2, self.IntensityTrajectoryChosenBS_fit_y0_Selected_2, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('Fit y0', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')

        
        
        plt.subplot(3,4,8)
        plt.title('Fit: average intensity change %', size=11)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_1, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage_Selected, s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_2, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage_Selected, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('Fit y0 diff ABS % (M1, M2)', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')

        
        
        
        plt.subplot(3,4,9)
        plt.title('Fit: average intensity change %', size=11)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_1, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage_Selected, s = 6,  color='red', label = '#1.'+fileType_1)
        #plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_2, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M1', size=11)
        plt.ylabel('Fit y0 diff ABS % (M1)', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')

        
        
        
        plt.subplot(3,4,10)
        plt.title('Fit: average intensity change %', size=11)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_2, self.Loaded_Intensity_Fit_y0_Range_Diff_Percentage_Selected, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M2', size=11)
        plt.ylabel('Fit y0 diff ABS % (M2)', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')


        


        plt.subplot(3,4,11)
        plt.title('Fit: constant intensity change', size=11)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_1, self.Loaded_Intensity_Fit_Range_CB_Diff_Percentage_Selected, s = 6, color='red', label = '#1.'+fileType_1)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_2, self.Loaded_Intensity_Fit_Range_CB_Diff_Percentage_Selected, s = 6, color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('Fit CB diff % (M1, M2)', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')

        

        plt.subplot(3,4,12)
        plt.title('Fit: modulation A change', size=11)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_1, self.Loaded_Intensity_Fit_Range_A_Diff_Percentage_Selected, s = 6,  color='red', label = '#1.'+fileType_1)
        plt.scatter(self.IntensityTrajectoryChosenBS_fit_M_Selected_2, self.Loaded_Intensity_Fit_Range_A_Diff_Percentage_Selected, s = 6,  color='blue', label = '#2.'+fileType_2)
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('Fit A diff % (M1, M2)', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=10, loc='upper left')

        
        
        
        self.btn_SaveSelectedData.Show()
        self.btn_SaveSelectedDataForAll.Show()
        
        plt.ion()
        plt.show()
        
        print 'SelectData_ShowFigures done'
                
                
                
                
                
                
    def SaveSelection(self, event):
        print ''
        prefix = self.FilePrefixForSelected.GetValue()
        
        if prefix != '':
            prefix += '_'

        todayDate = time.strftime("%Y%m%d_%Hh%Mm")  
          
        ff = open(self.data_1.fileDirectory + '\\' + 'DataSelection_' + prefix + self.data_1.fileName[:30] + todayDate +  '.txt','w')
        
        for n in range(len(self.cb_Data_yn)):
            
            if self.cb_Data_yn[n].GetValue():
                ff.write(str(n) + ' ' + 'y\n')    
            else:
                ff.write(str(n) + ' ' + 'n\n')    
                
        ff.close()           
    
    
    def LoadSelection(self, event):
        print ''
        dlg = wx.FileDialog(
            self, message="Choose a file",
            defaultFile="",
            wildcard=wildcardDataFile,
            style=wx.OPEN | wx.MULTIPLE) # | wx.CHANGE_DIR
            
        if dlg.ShowModal() == wx.ID_OK:
            paths = dlg.GetPaths()
            print "You chose the following file(s):"
            for path in paths:
                print path
        print dlg
        print 'path = ', dlg.GetPath()
        print 'filename = ', dlg.GetFilename()
        self.Show()
        
        #filenameTemp = dlg.GetFilename()
        self.filenameLoadedSelection = dlg.GetFilename()
        self.filedirectoryOnly_LoadedSelection = dlg.GetDirectory()
        
        
        fdata = open(paths[0], 'r')
        linestemp = fdata.readlines()
        fdata.close()
        
        self.Loaded_SelectionNunmer = []
        self.Loaded_SelectionYN = []
        
        for line in linestemp:
            p = line.split()
            self.Loaded_SelectionNunmer.append(int(float(p[0])) ) 
            self.Loaded_SelectionYN.append(p[1]) 
  
  
      
        if len(self.Loaded_SelectionNunmer) != len(self.cb_Data_yn):
            print 'data numbers are different'
            MessageBox(self, 'data numbers are different', 'Error')
            
        else:
            for n in range(len(self.cb_Data_yn)):
                if self.Loaded_SelectionYN[n] == 'y':
                    self.cb_Data_yn[n].SetValue(True)
                else:
                    self.cb_Data_yn[n].SetValue(False)
                    
      
      
      
      
      
    def SaveSelectedData12(self, event):
        print ' '
        self.text_SaveSelectedData.SetLabel('Saving Data ....')
        
        prefix = self.FilePrefixForSelected.GetValue()
        
        if prefix != '':
            prefix += '_'


        N_lettersToOmitFileName_begining = int(float(self.NofLettersToRemoveFinleNameForSelected_Beginning.GetValue()))
        N_lettersToOmitFileName_end = int(float(self.NofLettersToRemoveFinleNameForSelected_End.GetValue()))


        
        figFileName_1 = prefix + self.data_1.fileName[N_lettersToOmitFileName_begining:-(4+N_lettersToOmitFileName_end)]
        figFileName_2 = prefix + self.data_2.fileName[N_lettersToOmitFileName_begining:-(4+N_lettersToOmitFileName_end)]
        todayDate = time.strftime("%Y%m%d_%Hh%Mm")        



        if self.cb_SaveFigs_Selected_yn.GetValue():         
  
            self.fig_DataPlots_Selected.savefig(self.data_1.fileDirectory + '\\' + figFileName_1 +'_' + todayDate +  'DataPlots.png')
            plt.close(self.fig_DataPlots_Selected)
            print 'fig_DataPlots_Selected figure is saved'
    
            
    
    
            
        if self.cb_SaveImageJ_xyData_Selected_yn.GetValue():
            print ''
            ff = open(self.data_1.fileDirectory + '\\' + figFileName_1 +'_' + todayDate +  '_6_xy_A_T(frames)_phi(deg)_y0_M_Selected.txt','w')
            for i in range(len(self.x_ImageJ_Selected_1)):
                ff.write(str(self.x_ImageJ_Selected_1[i]) + ' ' + str(self.y_ImageJ_Selected_1[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_A_Selected_1[i])
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_T_Selected_1[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_phi_Selected_1[i]) 
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_y0_Selected_1[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_M_Selected_1[i]) + '\n'   )    
            ff.close()       
        
        if self.cb_SaveIntensityData_Selected_yn.GetValue():
            for n in range(len(self.IntensityTrajectoryChosenBS_Selected_1)):
                ff = open(self.data_1.fileDirectory + '\\' + figFileName_1 +'_' + todayDate +  '_8_Intensity_' + str(n) + 'Selected.txt','w')
                for i in range(len(self.IntensityTrajectoryChosenBS_Selected_1[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_Selected_1[n][i])) + '\n')
                ff.close()       
    
        if self.cb_SaveIntensityFitData_Selected_yn.GetValue():
            for n in range(len(self.IntensityTrajectoryChosenBS_fit_Selected)):
                ff = open(self.data_1.fileDirectory + '\\' + figFileName_1 +'_' + todayDate +  '_9_Fit_' + str(n) + 'Selected.txt','w')
                for i in range(len(self.IntensityTrajectoryChosenBS_fit_Selected_1[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_fit_Selected_1[n][i])) + '\n')
                ff.close()       
    
    
    
    
        if self.cb_Complete_Data_Selected_yn.GetValue():
            
            ff = open(self.data_1.fileDirectory + '\\' + figFileName_1 +'_' + todayDate +  '_10_CompleteDataAll_Selected' + '.txt','w')
            ff.write('%% x y A T phi y0 M Background_Subtraction_Type: '+ self.BGS_type_1 + ' - rawIntensity , intensityBS , fit , BG , BGfit , image0 , image1' + '\n')
            for n in range(len(self.IntensityTrajectoryChosenBS_Selected_1)):
                ff.write('#' + str(n) + '\n')
                    
                ff.write(str(self.x_ImageJ_Selected_1[n]) + ' ' + str(self.y_ImageJ_Selected_1[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_A_Selected_1[n])
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_T_Selected_1[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_phi_Selected_1[n]) 
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_y0_Selected_1[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_M_Selected_1[n]) + '\n'   )    
                
                for i in range(len(self.IntensityTrajectoryChosen_Selected_1[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosen_Selected_1[n][i])) + ' ')
                ff.write('\n')
                for i in range(len(self.IntensityTrajectoryChosenBS_Selected_1[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_Selected_1[n][i])) + ' ')
                ff.write('\n')  
                for i in range(len(self.IntensityTrajectoryChosenBS_fit_Selected_1[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_fit_Selected_1[n][i])) + ' ')
                ff.write('\n')  
                for i in range(len(self.backgroundIntensityTrajectoryChosen_Selected_1[n])):
                    ff.write(str(int(self.backgroundIntensityTrajectoryChosen_Selected_1[n][i])) + ' ')
                ff.write('\n')
                for i in range(len(self.backgroundIntensityTrajectoryChosen_Fit_Selected_1[n])):
                    ff.write(str(int(self.backgroundIntensityTrajectoryChosen_Fit_Selected_1[n][i])) + ' ')
                ff.write('\n')
                for i in range(len(self.centerImage_Selected_1[n][0])):
                    for k in range(len(self.centerImage_Selected_1[n][0][i])):
                        ff.write(str(int(self.centerImage_Selected_1[n][0][i][k])) + ' ')
                    ff.write(' , ')
                ff.write('\n')
                for i in range(len(self.centerImage_Selected_1[n][1])):
                    for k in range(len(self.centerImage_Selected_1[n][1][i])):
                        ff.write(str(int(self.centerImage_Selected_1[n][1][i][k])) + ' ')
                    ff.write(' , ')
                ff.write('\n\n')
                          
            ff.close()       
     
     
     




     
     

            
        if self.cb_SaveImageJ_xyData_Selected_yn.GetValue():
            print ''
            ff = open(self.filedirectoryOnly_2 + '\\' + figFileName_2 +'_' + todayDate +  '_6_xy_A_T(frames)_phi(deg)_y0_M_Selected.txt','w')
            for i in range(len(self.x_ImageJ_Selected_2)):
                ff.write(str(self.x_ImageJ_Selected_2[i]) + ' ' + str(self.y_ImageJ_Selected_2[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_A_Selected_2[i])
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_T_Selected_2[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_phi_Selected_2[i]) 
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_y0_Selected_2[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_M_Selected_2[i]) + '\n'   )    
            ff.close()       
        
        if self.cb_SaveIntensityData_Selected_yn.GetValue():
            for n in range(len(self.IntensityTrajectoryChosenBS_Selected_2)):
                ff = open(self.filedirectoryOnly_2 + '\\' + figFileName_2 +'_' + todayDate +  '_8_Intensity_' + str(n) + 'Selected.txt','w')
                for i in range(len(self.IntensityTrajectoryChosenBS_Selected_2[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_Selected_2[n][i])) + '\n')
                ff.close()       
    
        if self.cb_SaveIntensityFitData_Selected_yn.GetValue():
            for n in range(len(self.IntensityTrajectoryChosenBS_fit_Selected)):
                ff = open(self.filedirectoryOnly_2 + '\\' + figFileName_2 +'_' + todayDate +  '_9_Fit_' + str(n) + 'Selected.txt','w')
                for i in range(len(self.IntensityTrajectoryChosenBS_fit_Selected_2[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_fit_Selected_2[n][i])) + '\n')
                ff.close()       
    
    

    
        if self.cb_Complete_Data_Selected_yn.GetValue():
            
            ff = open(self.data_2.fileDirectory + '\\' + figFileName_2 +'_' + todayDate +  '_10_CompleteDataAll_Selected' + '.txt','w')
            ff.write('%% x y A T phi y0 M Background_Subtraction_Type: '+ self.BGS_type_2 + ' - rawIntensity , intensityBS , fit , BG , BGfit , image0 , image1' + '\n')
            for n in range(len(self.IntensityTrajectoryChosenBS_Selected_2)):
                ff.write('#' + str(n) + '\n')
                    
                ff.write(str(self.x_ImageJ_Selected_2[n]) + ' ' + str(self.y_ImageJ_Selected_2[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_A_Selected_2[n])
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_T_Selected_2[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_phi_Selected_2[n]) 
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_y0_Selected_2[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_M_Selected_2[n]) + '\n'   )    
                
                for i in range(len(self.IntensityTrajectoryChosen_Selected_2[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosen_Selected_2[n][i])) + ' ')
                ff.write('\n')
                for i in range(len(self.IntensityTrajectoryChosenBS_Selected_2[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_Selected_2[n][i])) + ' ')
                ff.write('\n')  
                for i in range(len(self.IntensityTrajectoryChosenBS_fit_Selected_2[n])):
                    ff.write(str(int(self.IntensityTrajectoryChosenBS_fit_Selected_2[n][i])) + ' ')
                ff.write('\n')  
                for i in range(len(self.backgroundIntensityTrajectoryChosen_Selected_2[n])):
                    ff.write(str(int(self.backgroundIntensityTrajectoryChosen_Selected_2[n][i])) + ' ')
                ff.write('\n')
                for i in range(len(self.backgroundIntensityTrajectoryChosen_Fit_Selected_2[n])):
                    ff.write(str(int(self.backgroundIntensityTrajectoryChosen_Fit_Selected_2[n][i])) + ' ')
                ff.write('\n')
                for i in range(len(self.centerImage_Selected_2[n][0])):
                    for k in range(len(self.centerImage_Selected_2[n][0][i])):
                        ff.write(str(int(self.centerImage_Selected_2[n][0][i][k])) + ' ')
                    ff.write(' , ')
                ff.write('\n')
                for i in range(len(self.centerImage_Selected_2[n][1])):
                    for k in range(len(self.centerImage_Selected_2[n][1][i])):
                        ff.write(str(int(self.centerImage_Selected_2[n][1][i][k])) + ' ')
                    ff.write(' , ')
                ff.write('\n\n')
                          
            ff.close()       
            
        self.text_SaveSelectedData.SetLabel('Data saved')
     


















      
    def SaveSelectedDataForAll(self, event):
        print '\n\n SaveSelectedDataForAll'
        
        
 
        if self.cb_Save_Data1_yn.GetValue():
            try:
                self.data_1.BGS_type
            except:
                MessageBox(self, 'No #1 data' , 'Error')
                return
                

        if self.cb_Save_Data2_yn.GetValue():
            try:
                self.data_2.BGS_type
            except:
                MessageBox(self, 'No #2 data' , 'Error')
                return

            if len(self.data_1.Loaded_ImageJxpos) != len(self.data_2.Loaded_ImageJxpos):
                MessageBox(self, '#1, #2 data numbers are different' , 'Error')
                return
                
            




        if self.cb_Save_Data3_yn.GetValue():
            try:
                self.data_3.BGS_type
            except:
                MessageBox(self, 'No #3 data' , 'Error')
                return
                
                
            if len(self.data_1.Loaded_ImageJxpos) != len(self.data_3.Loaded_ImageJxpos):
                MessageBox(self, '#1, # data numbers are different' , 'Error')
                return
                
            

        if self.cb_Save_Data4_yn.GetValue():
            try:
                self.data_4.BGS_type
            except:
                MessageBox(self, 'No #4 data' , 'Error')
                return
                
                
            if len(self.data_1.Loaded_ImageJxpos) != len(self.data_4.Loaded_ImageJxpos):
                MessageBox(self, '#1, #4 data numbers are different' , 'Error')
                return
                

        if self.cb_Save_Data5_yn.GetValue():
            try:
                self.data_5.BGS_type
            except:
                MessageBox(self, 'No #5 data' , 'Error')
                return
                
                
            if len(self.data_1.Loaded_ImageJxpos) != len(self.data_5.Loaded_ImageJxpos):
                MessageBox(self, '#1, #5 data numbers are different' , 'Error')
                return
                
            

        if self.cb_Save_Data6_yn.GetValue():
            try:
                self.data_6.BGS_type
            except:
                MessageBox(self, 'No #6 data' , 'Error')
                return
                
            if len(self.data_1.Loaded_ImageJxpos) != len(self.data_6.Loaded_ImageJxpos):
                MessageBox(self, '#1, #6 data numbers are different' , 'Error')
                return
                
            


       
        
        self.text_SaveSelectedDataForAll.SetLabel('Saving Data ....')
        

        

        for n in range(6):
            
            print '\n\n Saving Data #', n+1
            
            if n == 0:
                if self.cb_Save_Data1_yn.GetValue():
                        
                    FileName_temp = self.data_1.fileName
                    filedirectoryOnly_temp = self.data_1.fileDirectory
                    prefix = self.FilePrefixForSelected.GetValue()
                    
                    if prefix != '':
                        prefix += '_'

                    N_lettersToOmitFileName_begining = int(float(self.NofLettersToRemoveFinleNameForSelected_Beginning.GetValue()))
                    N_lettersToOmitFileName_end = int(float(self.NofLettersToRemoveFinleNameForSelected_End.GetValue()))
                    filename_temp = prefix + FileName_temp[N_lettersToOmitFileName_begining:-(4+N_lettersToOmitFileName_end)]
                    BGS_type_temp = self.data_1.BGS_type
                
                    self.x_ImageJ_Range_temp = self.Loaded_ImageJxpos_1_Range
                    self.y_ImageJ_Range_temp = self.Loaded_ImageJypos_1_Range
                    self.IntensityTrajectoryChosenBGS_fit_A_Range_temp = self.Loaded_A_1_Range
                    self.IntensityTrajectoryChosenBGS_fit_T_Range_temp = self.Loaded_T_1_Range
                    self.IntensityTrajectoryChosenBGS_fit_phi_Range_temp = self.Loaded_phi_1_Range
                    self.IntensityTrajectoryChosenBGS_fit_y0_Range_temp = self.Loaded_y0_1_Range
                    self.IntensityTrajectoryChosenBGS_fit_M_Range_temp = self.Loaded_M_1_Range
                    
                    self.IntensityTrajectory_noBGS_Range_temp = self.Loaded_Intensity_noBGS_1_Range
                    self.IntensityTrajectory_BGS_Range_temp = self.Loaded_Intensity_1_Range
                    self.IntensityTrajectory_BGS_fit_Range_temp = self.Loaded_Intensity_Fit_1_Range
                    self.backgroundIntensityTrajectory_BG_Range_temp = self.Loaded_Intensity_BG_1_Range
                    self.backgroundIntensityTrajectory_BG_fit_Range_temp = self.Loaded_Intensity_BG_Fit_1_Range
                    self.centerImages_Range_temp = self.Loaded_CenterImages_1_Range
                    
                else: continue
                    
            if n == 1:
                if self.cb_Save_Data2_yn.GetValue():
    
                    FileName_temp = self.data_2.fileName
                    filedirectoryOnly_temp = self.data_2.fileDirectory
                    prefix = self.FilePrefixForSelected.GetValue()
                                        
                    if prefix != '':
                        prefix += '_'

                    N_lettersToOmitFileName_begining = int(float(self.NofLettersToRemoveFinleNameForSelected_Beginning.GetValue()))
                    N_lettersToOmitFileName_end = int(float(self.NofLettersToRemoveFinleNameForSelected_End.GetValue()))
                    filename_temp = prefix + FileName_temp[N_lettersToOmitFileName_begining:-(4+N_lettersToOmitFileName_end)]
                    BGS_type_temp = self.data_2.BGS_type
                
                    self.x_ImageJ_Range_temp = self.Loaded_ImageJxpos_2_Range
                    self.y_ImageJ_Range_temp = self.Loaded_ImageJypos_2_Range
                    self.IntensityTrajectoryChosenBGS_fit_A_Range_temp = self.Loaded_A_2_Range
                    self.IntensityTrajectoryChosenBGS_fit_T_Range_temp = self.Loaded_T_2_Range
                    self.IntensityTrajectoryChosenBGS_fit_phi_Range_temp = self.Loaded_phi_2_Range
                    self.IntensityTrajectoryChosenBGS_fit_y0_Range_temp = self.Loaded_y0_2_Range
                    self.IntensityTrajectoryChosenBGS_fit_M_Range_temp = self.Loaded_M_2_Range
                    
                    self.IntensityTrajectory_noBGS_Range_temp = self.Loaded_Intensity_noBGS_2_Range
                    self.IntensityTrajectory_BGS_Range_temp = self.Loaded_Intensity_2_Range
                    self.IntensityTrajectory_BGS_fit_Range_temp = self.Loaded_Intensity_Fit_2_Range
                    self.backgroundIntensityTrajectory_BG_Range_temp = self.Loaded_Intensity_BG_2_Range
                    self.backgroundIntensityTrajectory_BG_fit_Range_temp = self.Loaded_Intensity_BG_Fit_2_Range
                    self.centerImages_Range_temp = self.Loaded_CenterImages_2_Range
                
                else: continue

                    
            if n == 2:

                if self.cb_Save_Data3_yn.GetValue():
    
                    FileName_temp = self.data_3.fileName
                    filedirectoryOnly_temp = self.data_3.fileDirectory
                    prefix = self.FilePrefixForSelected.GetValue()                    
                    if prefix != '':
                        prefix += '_'

                    N_lettersToOmitFileName_begining = int(float(self.NofLettersToRemoveFinleNameForSelected_Beginning.GetValue()))
                    N_lettersToOmitFileName_end = int(float(self.NofLettersToRemoveFinleNameForSelected_End.GetValue()))
                    filename_temp = prefix + FileName_temp[N_lettersToOmitFileName_begining:-(4+N_lettersToOmitFileName_end)]
                    BGS_type_temp = self.data_3.BGS_type
                
                    self.x_ImageJ_Range_temp = self.Loaded_ImageJxpos_3_Range
                    self.y_ImageJ_Range_temp = self.Loaded_ImageJypos_3_Range
                    self.IntensityTrajectoryChosenBGS_fit_A_Range_temp = self.Loaded_A_3_Range
                    self.IntensityTrajectoryChosenBGS_fit_T_Range_temp = self.Loaded_T_3_Range
                    self.IntensityTrajectoryChosenBGS_fit_phi_Range_temp = self.Loaded_phi_3_Range
                    self.IntensityTrajectoryChosenBGS_fit_y0_Range_temp = self.Loaded_y0_3_Range
                    self.IntensityTrajectoryChosenBGS_fit_M_Range_temp = self.Loaded_M_3_Range
                    
                    self.IntensityTrajectory_noBGS_Range_temp = self.Loaded_Intensity_noBGS_3_Range
                    self.IntensityTrajectory_BGS_Range_temp = self.Loaded_Intensity_3_Range
                    self.IntensityTrajectory_BGS_fit_Range_temp = self.Loaded_Intensity_Fit_3_Range
                    self.backgroundIntensityTrajectory_BG_Range_temp = self.Loaded_Intensity_BG_3_Range
                    self.backgroundIntensityTrajectory_BG_fit_Range_temp = self.Loaded_Intensity_BG_Fit_3_Range
                    self.centerImages_Range_temp = self.Loaded_CenterImages_3_Range
                else:
                    continue
                
            if n == 3:

                if self.cb_Save_Data4_yn.GetValue():
    
                    FileName_temp = self.data_4.fileName
                    filedirectoryOnly_temp = self.data_4.fileDirectory
                    prefix = self.FilePrefixForSelected.GetValue()                    
                    if prefix != '':
                        prefix += '_'

                    N_lettersToOmitFileName_begining = int(float(self.NofLettersToRemoveFinleNameForSelected_Beginning.GetValue()))
                    N_lettersToOmitFileName_end = int(float(self.NofLettersToRemoveFinleNameForSelected_End.GetValue()))
                    filename_temp = prefix + FileName_temp[N_lettersToOmitFileName_begining:-(4+N_lettersToOmitFileName_end)]
                    BGS_type_temp = self.data_4.BGS_type
                
                    self.x_ImageJ_Range_temp = self.Loaded_ImageJxpos_4_Range
                    self.y_ImageJ_Range_temp = self.Loaded_ImageJypos_4_Range
                    self.IntensityTrajectoryChosenBGS_fit_A_Range_temp = self.Loaded_A_4_Range
                    self.IntensityTrajectoryChosenBGS_fit_T_Range_temp = self.Loaded_T_4_Range
                    self.IntensityTrajectoryChosenBGS_fit_phi_Range_temp = self.Loaded_phi_4_Range
                    self.IntensityTrajectoryChosenBGS_fit_y0_Range_temp = self.Loaded_y0_4_Range
                    self.IntensityTrajectoryChosenBGS_fit_M_Range_temp = self.Loaded_M_4_Range
                    
                    self.IntensityTrajectory_noBGS_Range_temp = self.Loaded_Intensity_noBGS_4_Range
                    self.IntensityTrajectory_BGS_Range_temp = self.Loaded_Intensity_4_Range
                    self.IntensityTrajectory_BGS_fit_Range_temp = self.Loaded_Intensity_Fit_4_Range
                    self.backgroundIntensityTrajectory_BG_Range_temp = self.Loaded_Intensity_BG_4_Range
                    self.backgroundIntensityTrajectory_BG_fit_Range_temp = self.Loaded_Intensity_BG_Fit_4_Range
                    self.centerImages_Range_temp = self.Loaded_CenterImages_4_Range
                    
                else: continue
        
            if n == 4:
                
                if self.cb_Save_Data5_yn.GetValue():
                    
                    FileName_temp = self.data_5.fileName
                    filedirectoryOnly_temp = self.data_5.fileDirectory
                    prefix = self.FilePrefixForSelected.GetValue()                    
                    if prefix != '':
                        prefix += '_'

                    N_lettersToOmitFileName_begining = int(float(self.NofLettersToRemoveFinleNameForSelected_Beginning.GetValue()))
                    N_lettersToOmitFileName_end = int(float(self.NofLettersToRemoveFinleNameForSelected_End.GetValue()))
                    filename_temp = prefix + FileName_temp[N_lettersToOmitFileName_begining:-(4+N_lettersToOmitFileName_end)]
                    BGS_type_temp = self.data_5.BGS_type
                
                    self.x_ImageJ_Range_temp = self.Loaded_ImageJxpos_5_Range
                    self.y_ImageJ_Range_temp = self.Loaded_ImageJypos_5_Range
                    self.IntensityTrajectoryChosenBGS_fit_A_Range_temp = self.Loaded_A_5_Range
                    self.IntensityTrajectoryChosenBGS_fit_T_Range_temp = self.Loaded_T_5_Range
                    self.IntensityTrajectoryChosenBGS_fit_phi_Range_temp = self.Loaded_phi_5_Range
                    self.IntensityTrajectoryChosenBGS_fit_y0_Range_temp = self.Loaded_y0_5_Range
                    self.IntensityTrajectoryChosenBGS_fit_M_Range_temp = self.Loaded_M_5_Range
                    
                    self.IntensityTrajectory_noBGS_Range_temp = self.Loaded_Intensity_noBGS_5_Range
                    self.IntensityTrajectory_BGS_Range_temp = self.Loaded_Intensity_5_Range
                    self.IntensityTrajectory_BGS_fit_Range_temp = self.Loaded_Intensity_Fit_5_Range
                    self.backgroundIntensityTrajectory_BG_Range_temp = self.Loaded_Intensity_BG_5_Range
                    self.backgroundIntensityTrajectory_BG_fit_Range_temp = self.Loaded_Intensity_BG_Fit_5_Range
                    self.centerImages_Range_temp = self.Loaded_CenterImages_5_Range
                    
                else: continue
                        
                        
            if n == 5:
                
                if self.cb_Save_Data6_yn.GetValue():
    
                    FileName_temp = self.data_6.fileName
                    filedirectoryOnly_temp = self.data_6.fileDirectory
                    prefix = self.FilePrefixForSelected.GetValue()                    
                    if prefix != '':
                        prefix += '_'

                    N_lettersToOmitFileName_begining = int(float(self.NofLettersToRemoveFinleNameForSelected_Beginning.GetValue()))
                    N_lettersToOmitFileName_end = int(float(self.NofLettersToRemoveFinleNameForSelected_End.GetValue()))
                    filename_temp = prefix + FileName_temp[N_lettersToOmitFileName_begining:-(4+N_lettersToOmitFileName_end)]
                    BGS_type_temp = self.data_6.BGS_type
                
                    self.x_ImageJ_Range_temp = self.Loaded_ImageJxpos_6_Range
                    self.y_ImageJ_Range_temp = self.Loaded_ImageJypos_6_Range
                    self.IntensityTrajectoryChosenBGS_fit_A_Range_temp = self.Loaded_A_6_Range
                    self.IntensityTrajectoryChosenBGS_fit_T_Range_temp = self.Loaded_T_6_Range
                    self.IntensityTrajectoryChosenBGS_fit_phi_Range_temp = self.Loaded_phi_6_Range
                    self.IntensityTrajectoryChosenBGS_fit_y0_Range_temp = self.Loaded_y0_6_Range
                    self.IntensityTrajectoryChosenBGS_fit_M_Range_temp = self.Loaded_M_6_Range
                    
                    self.IntensityTrajectory_noBGS_Range_temp = self.Loaded_Intensity_noBGS_6_Range
                    self.IntensityTrajectory_BGS_Range_temp = self.Loaded_Intensity_6_Range
                    self.IntensityTrajectory_BGS_fit_Range_temp = self.Loaded_Intensity_Fit_6_Range
                    self.backgroundIntensityTrajectory_BG_Range_temp = self.Loaded_Intensity_BG_6_Range
                    self.backgroundIntensityTrajectory_BG_fit_Range_temp = self.Loaded_Intensity_BG_Fit_6_Range
                    self.centerImages_Range_temp = self.Loaded_CenterImages_6_Range
                    
                else: continue
                    
            



            self.x_ImageJ_temp = []
            self.y_ImageJ_temp = []
            self.IntensityTrajectoryChosenBGS_fit_A_temp = []
            self.IntensityTrajectoryChosenBGS_fit_T_temp = []
            self.IntensityTrajectoryChosenBGS_fit_phi_temp = []
            self.IntensityTrajectoryChosenBGS_fit_y0_temp = []
            self.IntensityTrajectoryChosenBGS_fit_M_temp = []
            
            self.IntensityTrajectory_noBGS_temp = []
            self.IntensityTrajectory_BGS_temp = []
            self.IntensityTrajectory_BGS_fit_temp = []
            self.backgroundIntensityTrajectory_BG_temp = []
            self.backgroundIntensityTrajectory_BG_fit_temp = []
            self.centerImages_temp = []
            
    
    
            for n in range(len(self.cb_Data_yn)):
                if self.cb_Data_yn[n].GetValue():
                    
                    self.x_ImageJ_temp.append(self.x_ImageJ_Range_temp[n])
                    self.y_ImageJ_temp.append(self.y_ImageJ_Range_temp[n])
                    self.IntensityTrajectoryChosenBGS_fit_A_temp.append(self.IntensityTrajectoryChosenBGS_fit_A_Range_temp[n])
                    self.IntensityTrajectoryChosenBGS_fit_T_temp.append(self.IntensityTrajectoryChosenBGS_fit_T_Range_temp[n])
                    self.IntensityTrajectoryChosenBGS_fit_phi_temp.append(self.IntensityTrajectoryChosenBGS_fit_phi_Range_temp[n])
                    self.IntensityTrajectoryChosenBGS_fit_y0_temp.append(self.IntensityTrajectoryChosenBGS_fit_y0_Range_temp[n])
                    self.IntensityTrajectoryChosenBGS_fit_M_temp.append(self.IntensityTrajectoryChosenBGS_fit_M_Range_temp[n])
                    
                    self.IntensityTrajectory_noBGS_temp.append(self.IntensityTrajectory_noBGS_Range_temp[n])
                    self.IntensityTrajectory_BGS_temp.append(self.IntensityTrajectory_BGS_Range_temp[n])
                    self.IntensityTrajectory_BGS_fit_temp.append(self.IntensityTrajectory_BGS_fit_Range_temp[n])
                    self.backgroundIntensityTrajectory_BG_temp.append(self.backgroundIntensityTrajectory_BG_Range_temp[n])
                    self.backgroundIntensityTrajectory_BG_fit_temp.append(self.backgroundIntensityTrajectory_BG_fit_Range_temp[n])
                    self.centerImages_temp.append(self.centerImages_Range_temp[n])
                    
         
    
    
    
            filedirectoryOnly = filedirectoryOnly_temp
            filename = filename_temp
    
            
            xpos_ImageJ = self.x_ImageJ_temp
            ypos_ImageJ = self.y_ImageJ_temp
            A = self.IntensityTrajectoryChosenBGS_fit_A_temp
            T = self.IntensityTrajectoryChosenBGS_fit_T_temp
            phi = self.IntensityTrajectoryChosenBGS_fit_phi_temp
            y0 = self.IntensityTrajectoryChosenBGS_fit_y0_temp
            M = self.IntensityTrajectoryChosenBGS_fit_M_temp
            
            I_noBGS = self.IntensityTrajectory_noBGS_temp
            I_BGS = self.IntensityTrajectory_BGS_temp
            I_BGS_fit = self.IntensityTrajectory_BGS_fit_temp
            I_BG = self.backgroundIntensityTrajectory_BG_temp
            I_BG_fit = self.backgroundIntensityTrajectory_BG_fit_temp
            CenterImages = self.centerImages_temp
    
            
            cb_SaveImageJ_xyData_yn = False
            cb_Complete_Data_yn = True
            
            BGS_type = BGS_type_temp
    
    
            CreateDataFiles(filedirectoryOnly, filename, xpos_ImageJ, ypos_ImageJ, A, T, phi, y0, M, I_noBGS, I_BGS, I_BGS_fit, I_BG, I_BG_fit, CenterImages, cb_SaveImageJ_xyData_yn, cb_Complete_Data_yn, BGS_type)
    
    


            
            
            
        self.text_SaveSelectedDataForAll.SetLabel('Data saved')
     





















      

    def plotDataSeparate(self, event):
        
        print '\n\nplotDataSeparate '
   
   

        try:
            self.data_1.Loaded_ImageJxpos
        except:
            MessageBox(self, 'No #1 Data' , 'Error')
            return
              
            
            
        
        try:
            plt.close(self.fig_DataPlots_Separate)
            for n in range(len(self.fig_AllDataCenterImages_Loaded)):
                plt.close(self.fig_AllDataCenterImages_Loaded[n])
        except:
            print ' except: no previous figures'
            pass
        
        
        
        

        '''
        self.BGS_type_1 = self.data_1.BGS_type
        self.Loaded_ImageJxpos_1_Range = self.data_1.Loaded_ImageJxpos
        self.Loaded_ImageJypos_1_Range = self.data_1.Loaded_ImageJypos
        self.Loaded_M_1_Range = self.data_1.Loaded_M
        self.Loaded_A_1_Range = self.data_1.Loaded_A
        self.Loaded_T_1_Range = self.data_1.Loaded_T
        self.Loaded_phi_1_Range = self.data_1.Loaded_phi
        self.Loaded_y0_1_Range = self.data_1.Loaded_y0
        self.Loaded_CenterImages_1_Range = self.data_1.Loaded_CenterImages
        self.Loaded_Intensity_1_Range = self.data_1.Loaded_Intensity
        self.Loaded_Intensity_noBGS_1_Range = self.data_1.Loaded_Intensity_noBGS
        self.Loaded_Intensity_BG_1_Range = self.data_1.Loaded_Intensity_BG
        self.Loaded_Intensity_BG_Fit_1_Range = self.data_1.Loaded_Intensity_BG_Fit
        self.Loaded_Intensity_Fit_1_Range = self.data_1.Loaded_Intensity_Fit
        self.M_error_1_Range = self.data_1.M_error
        self.SNR_1_Range = self.data_1.SNR
        self.SBR_1_Range = self.data_1.Loaded_SBR
        '''

        
        
        

            
            
        self.BGS_type_1 = self.data_1.BGS_type
            
            
        self.Loaded_ImageJxpos_1_Range = []    
        self.Loaded_ImageJypos_1_Range = []
        self.Loaded_M_1_Range = []
        self.Loaded_A_1_Range = []
        self.Loaded_T_1_Range = []
        self.Loaded_phi_1_Range = []
        self.Loaded_y0_1_Range = []
        
        self.Loaded_Intensity_1_Range = []
        self.Loaded_Intensity_Fit_1_Range = []
        self.Loaded_Intensity_noBGS_1_Range = []
        self.Loaded_Intensity_BG_1_Range = []
        self.Loaded_Intensity_BG_Fit_1_Range = []
        self.Loaded_CenterImages_1_Range = []
        self.M_error_1_Range = []
        self.SNR_1_Range = []
        self.SBR_1_Range = []
        
            
        


                   

        if self.cb_M1Range_yn.GetValue() == False and self.cb_Intensity1Range_yn.GetValue() == False and self.cb_Fit_y0_1_Range_yn.GetValue() == False\
        and self.cb_SNR_1Range_yn.GetValue() == False and self.cb_SBR_1Range_yn.GetValue() == False :
        
            self.Loaded_ImageJxpos_1_Range = self.data_1.Loaded_ImageJxpos
            self.Loaded_ImageJypos_1_Range = self.data_1.Loaded_ImageJypos
            self.Loaded_M_1_Range = self.data_1.Loaded_M
            self.Loaded_A_1_Range = self.data_1.Loaded_A
            self.Loaded_T_1_Range = self.data_1.Loaded_T
            self.Loaded_phi_1_Range = self.data_1.Loaded_phi
            self.Loaded_y0_1_Range = self.data_1.Loaded_y0
            self.Loaded_CenterImages_1_Range = self.data_1.Loaded_CenterImages
            self.Loaded_Intensity_1_Range = self.data_1.Loaded_Intensity
            self.Loaded_Intensity_noBGS_1_Range = self.data_1.Loaded_Intensity_noBGS
            
            self.Loaded_Intensity_BG_1_Range = self.data_1.Loaded_Intensity_BG
            self.Loaded_Intensity_BG_Fit_1_Range = self.data_1.Loaded_Intensity_BG_Fit
            self.Loaded_Intensity_Fit_1_Range = self.data_1.Loaded_Intensity_Fit
            self.M_error_1_Range = self.data_1.M_error
            self.SNR_1_Range = self.data_1.SNR
            self.SBR_1_Range = self.data_1.Loaded_SBR





                
                

                    
        else:
            for n in range(len(self.data_1.Loaded_M)):
                
                if self.cb_M1Range_yn.GetValue() == True:
                    if self.data_1.Loaded_M[n] >= float(self.M1RangeValueMin.GetValue()) and self.data_1.Loaded_M[n] <= float(self.M1RangeValueMax.GetValue()):
                        pass
                    else: continue
                
                if self.cb_Intensity1Range_yn.GetValue() == True:
                    if np.mean(self.data_1.Loaded_Intensity[n]) >= float(self.Intensity1RangeValueMin.GetValue()) and np.mean(self.data_1.Loaded_Intensity[n]) <= float(self.Intensity1RangeValueMax.GetValue()):
                        pass
                    else: continue
                
                if self.cb_Fit_y0_1_Range_yn.GetValue() == True:
                    if self.data_1.Loaded_y0[n] >= float(self.Fit_y0_1RangeValueMin.GetValue()) and self.data_1.Loaded_y0[n] <= float(self.Fit_y0_1RangeValueMax.GetValue()):
                        pass
                    else: continue
                
                if self.cb_SNR_1Range_yn.GetValue() == True:
                    if (self.data_1.SNR[n]) >= float(self.SNR_1RangeValueMin.GetValue()) and (self.data_1.SNR[n]) <= float(self.SNR_1RangeValueMax.GetValue()):
                        pass
                    else: continue

                if self.cb_SBR_1Range_yn.GetValue() == True:
                    if (self.data_1.Loaded_SBR[n]) >= float(self.SBR_1RangeValueMin.GetValue()) and (self.data_1.Loaded_SBR[n]) <= float(self.SBR_1RangeValueMax.GetValue()):
                        pass
                    else: continue

               
                self.Loaded_ImageJxpos_1_Range.append(self.data_1.Loaded_ImageJxpos[n])
                self.Loaded_ImageJypos_1_Range.append(self.data_1.Loaded_ImageJypos[n])
                self.Loaded_M_1_Range.append(self.data_1.Loaded_M[n])
                self.Loaded_A_1_Range.append(self.data_1.Loaded_A[n])
                self.Loaded_T_1_Range.append(self.data_1.Loaded_T[n])
                self.Loaded_phi_1_Range.append(self.data_1.Loaded_phi[n])
                self.Loaded_y0_1_Range.append(self.data_1.Loaded_y0[n])
                self.Loaded_CenterImages_1_Range.append(self.data_1.Loaded_CenterImages[n])
                self.Loaded_Intensity_1_Range.append(self.data_1.Loaded_Intensity[n])
                self.Loaded_Intensity_noBGS_1_Range.append(self.data_1.Loaded_Intensity_noBGS[n])
                
                self.Loaded_Intensity_BG_1_Range.append(self.data_1.Loaded_Intensity_BG[n])
                self.Loaded_Intensity_BG_Fit_1_Range.append(self.data_1.Loaded_Intensity_BG_Fit[n])
                                
                self.Loaded_Intensity_Fit_1_Range.append(self.data_1.Loaded_Intensity_Fit[n])
                self.M_error_1_Range.append(self.data_1.M_error[n])
                self.SNR_1_Range.append(self.data_1.SNR[n])
                self.SBR_1_Range.append(self.data_1.Loaded_SBR[n])








        if len(self.Loaded_ImageJxpos_1_Range) == 0:
            MessageBox(self, '# of data is zero', 'Error')
            return




        data_conditions = ': Conditions : '

        if self.cb_M1Range_yn.GetValue() == True:
            data_conditions += '\n' + self.M1RangeValueMin.GetValue() +' <= M1 <= '+ self.M1RangeValueMax.GetValue() + ' ,   '
        
        if self.cb_M2Range_yn.GetValue() == True:
            data_conditions += self.M2RangeValueMin.GetValue() +' <= M2 <= '+ self.M2RangeValueMax.GetValue() + ' ,   '
        
        
        if self.cb_Intensity1Range_yn.GetValue() == True:
            data_conditions += '\n' + self.Intensity1RangeValueMin.GetValue() +' <= Intensity1 <= '+ self.Intensity1RangeValueMax.GetValue() + ' ,   '
        
        if self.cb_Intensity2Range_yn.GetValue() == True:
            data_conditions += self.Intensity2RangeValueMin.GetValue() +' <= Intensity2 <= '+ self.Intensity2RangeValueMax.GetValue() + ' ,   '
        
        
        
        if self.cb_Fit_y0_1_Range_yn.GetValue() == True:
            data_conditions += '\n' + self.Fit_y0_1RangeValueMin.GetValue() +' <= Fit:y0_#1 <= '+ self.Fit_y0_1RangeValueMax.GetValue() + ' ,   '
        

        if self.cb_Fit_y0_2_Range_yn.GetValue() == True:
            data_conditions += self.Fit_y0_2RangeValueMin.GetValue() +' <= Fit:y0_#2 <= '+ self.Fit_y0_2RangeValueMax.GetValue() + ' ,   '
        


        
        

        if self.cb_SNR_1Range_yn.GetValue() == True:
            data_conditions +='\n' +  self.SNR_1RangeValueMin.GetValue() +' <= SNR1 <= '+ self.SNR_1RangeValueMax.GetValue() + ' ,   '
            
        if self.cb_SNR_2Range_yn.GetValue() == True:
            data_conditions += self.SNR_2RangeValueMin.GetValue() +' <= SNR1 <= '+ self.SNR_2RangeValueMax.GetValue() + ' ,   '




        if self.cb_SBR_1Range_yn.GetValue() == True:
            data_conditions += '\n' +  self.SBR_1RangeValueMin.GetValue() +' <= SBR1 <= '+ self.SBR_1RangeValueMax.GetValue() + ' ,   '

        if self.cb_SBR_2Range_yn.GetValue() == True:
            data_conditions += self.SBR_2RangeValueMin.GetValue() +' <= SBR2 <= '+ self.SBR_2RangeValueMax.GetValue()            
            










        self.btn_plotDataSeparate_text1.SetLabel('1. #: ' + str(len(self.Loaded_M_1_Range)))
        
        
        
        
      



               
                 


           



                   

        self.Loaded_Intensity_Fit_1_Range_CB = np.array(self.Loaded_y0_1_Range) - np.abs(self.Loaded_A_1_Range)  # CB == constant background
        self.Loaded_Intensity_Fit_1_Range_Max = np.array(self.Loaded_y0_1_Range) + np.abs(self.Loaded_A_1_Range)  # CB == constant background

        self.Loaded_IntensityNoiseSTD_1_Range = []
        for n in range(len(self.Loaded_Intensity_1_Range)):
            STDtemp1 = np.std( np.array(self.Loaded_Intensity_1_Range[n]) - np.array(self.Loaded_Intensity_Fit_1_Range[n]) )
            print n , ' STDtemp1 ', ' = ' ,STDtemp1
            self.Loaded_IntensityNoiseSTD_1_Range.append(STDtemp1)
        self.SNR_1_Range = np.array(self.Loaded_y0_1_Range)/np.array(self.Loaded_IntensityNoiseSTD_1_Range)





        

        try:
            self.data_2
            plotText = '#1: ' + self.data_1.fileName + '\n' + '#2: ' + self.data_2.fileName + '\nBGS_type: ' + self.data_1.BGS_type + ', ' + self.data_2.BGS_type
            print 'self.data_2 exists'
            pass
        except:
            print '\n######################\n no data#2 \n####################\n'
            plotText = '#1: ' + self.data_1.fileName + '\n' + '\nBGS_type: ' + self.data_1.BGS_type
            
        

        self.fig_DataPlots_Separate = plt.figure('Data_Plots_Separate', figsize = (27,12), facecolor='white')
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(110,30,1800, 1000)
        self.fig_DataPlots_Separate.subplots_adjust(top = 0.89, bottom = 0.07, left = 0.05, right = 0.98, wspace = 0.4, hspace = 0.4)
        self.fig_DataPlots_Separate.text(0.02, 0.95, plotText , size=12, color="black", weight='normal')
        self.fig_DataPlots_Separate.text(0.62, 0.93, data_conditions, size=12, color="red", weight='normal')

        
        
        plt.subplot(4,6,1)
        plt.title('Median: #1: ' + str(np.around(np.median(self.Loaded_M_1_Range),2)), size = 10)
        plt.hist(self.Loaded_M_1_Range, bins=np.arange(0.0, 1.4, 0.1), label = '#1', alpha=1.0, color='red')
        plt.xlabel('Modulation Depth (M)', size=11)
        plt.xlim(0.0, 1.31)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        
        
        

        
        plt.subplot(4,6,2)
        plt.title('Median: #1: ' + str(np.around(np.median(self.Loaded_phi_1_Range),2)), size = 10)
        plt.hist(self.Loaded_phi_1_Range, label = '#1', alpha=1.0, color='red')
        plt.xlabel('phi (deg)', size=11)
        #plt.xlim(0.0, 1.31)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        


        
        plt.subplot(4,6,3)
        plt.title('Fit: average intensity y0', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.Loaded_y0_1_Range, s = 6,  color='red', label = '#1')
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('Fit y0', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=9, loc='upper right',frameon=False)

        

        
        plt.subplot(4,6,5)
        plt.title('NoiseSTD vs M', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.Loaded_IntensityNoiseSTD_1_Range,  s = 6,  color='red', label = '#1')
        plt.xlabel('M', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        plt.locator_params(nbins=5)
        
        
        
        A = np.vstack([self.Loaded_y0_1_Range, np.ones(len(self.Loaded_y0_1_Range))]).T        
        m, c = np.linalg.lstsq(A, self.Loaded_IntensityNoiseSTD_1_Range)[0]        
        x_temp = np.arange(-1000, 20000, 10)
    
        A_guess = m*10
        B_guess = np.min(self.Loaded_IntensityNoiseSTD_1_Range)
        
        initial_guess = [A_guess, B_guess]
        
        try:
            popt, pcov = curve_fit(Noise_func, self.Loaded_y0_1_Range, self.Loaded_IntensityNoiseSTD_1_Range, p0=initial_guess )
            fit_temp = Noise_func(x_temp, popt[0],popt[1])
            #err_temp = np.sqrt(np.diag(pcov))
        except: 
            print '\n fit error '
            fit_temp = x_temp *0.0
            popt = [0, 0]
            pass
   


        plt.subplot(4,6,6)
        plt.title('NoiseSTD vs Ave Intensity (y0)', size=11)
        plt.scatter(self.Loaded_y0_1_Range, self.Loaded_IntensityNoiseSTD_1_Range,  s = 5,  color='red', label = '#1')
        plt.plot(x_temp,  m*x_temp + c, 'c', label='Linear Fit: slope='+str(np.around(m,1)), linewidth = 2)
        plt.plot(x_temp, fit_temp, 'm', label='Noise Fit: A,B= '+str(np.around(popt[0],1))+','+str(np.around(popt[1],1)), linewidth = 2)
        plt.xlabel('Ave Intensity', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        plt.xlim(0, np.max(self.Loaded_y0_1_Range))
        plt.ylim(0, np.max(self.Loaded_IntensityNoiseSTD_1_Range))
        plt.locator_params(nbins=5)
        

   
        NetNoiseSTD = []
        for n in range(len(self.Loaded_M_1_Range)):
            dI_temp = Noise_func(self.Loaded_y0_1_Range[n], popt[0],popt[1])
            NetNoiseTemp = self.Loaded_IntensityNoiseSTD_1_Range[n] - dI_temp
            NetNoiseSTD.append(NetNoiseTemp)
            
            
        A = np.vstack([self.Loaded_M_1_Range, np.ones(len(self.Loaded_M_1_Range))]).T        
        m, c = np.linalg.lstsq(A, NetNoiseSTD)[0]        
        x_temp = np.arange(-0, 1.2, 0.1)
        
        
        plt.subplot(4,6,7)
        plt.title('Renormalized NetNoiseSTD vs M', size=11)
        plt.scatter(self.Loaded_M_1_Range, NetNoiseSTD,  s = 6,  color='red', label = '#1')
        plt.plot(x_temp,  m*x_temp + c, 'c', label='Linear Fit: slope='+str(np.around(m,1)), linewidth = 2)
        plt.xlabel('M', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        plt.locator_params(nbins=5)


        
        plt.subplot(4,6,8)
        plt.title('SBR vs M', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.SBR_1_Range, s = 6,  color='red', label = '#1')
        plt.xlabel('M', size=11)
        plt.ylabel('SBR', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        plt.locator_params(nbins=5)
        
        

        plt.subplot(4,6,9)
        plt.hist(self.SNR_1_Range, label = '#1', alpha=1.0, color='red')
        plt.xlabel('SNR', size=11)
        #plt.xlim(0.0, 1.21)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        
        
        plt.subplot(4,6,10)
        plt.title('SNR vs M', size=11)
        plt.scatter(self.Loaded_M_1_Range, self.SNR_1_Range,  s = 6,  color='red', label = '#1')
        plt.xlabel('M', size=11)
        plt.ylabel('SNR', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        
        
        plt.subplot(4,6,11)
        plt.title('M_error vs SNR', size=11)
        plt.scatter(self.SNR_1_Range, self.M_error_1_Range,  s = 6,  color='red', label = '#1')
        plt.xlabel('SNR', size=11)
        plt.ylabel('M_error', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        
        


        
        plt.subplot(4,6,12)
        plt.title('M_error vs 1/SNR', size=11)
        plt.scatter(1/self.SNR_1_Range, self.M_error_1_Range,  s = 6,  color='red', label = '#1')
        plt.xlabel('1/SNR', size=11)
        plt.ylabel('M_error', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        plt.locator_params(nbins=6)
        
        

        plt.ion()
        plt.show()
        
        










        try:
            self.data_2
            print 'self.data_2 exists'
            pass
        except:
            print '\n######################\n no data#2 \n####################\n'
            return
            
            


                            


            
        self.BGS_type_2 = self.data_2.BGS_type
            
            
        self.Loaded_ImageJxpos_2_Range = []    
        self.Loaded_ImageJypos_2_Range = []
        self.Loaded_M_2_Range = []
        self.Loaded_A_2_Range = []
        self.Loaded_T_2_Range = []
        self.Loaded_phi_2_Range = []
        self.Loaded_y0_2_Range = []
        
        self.Loaded_Intensity_2_Range = []
        self.Loaded_Intensity_Fit_2_Range = []
        self.Loaded_Intensity_noBGS_2_Range = []
        self.Loaded_Intensity_BG_2_Range = []
        self.Loaded_Intensity_BG_Fit_2_Range = []
        self.Loaded_CenterImages_2_Range = []
        self.M_error_2_Range = []
        self.SNR_2_Range = []
        self.SBR_2_Range = []
        
            
        


                   

        if self.cb_M2Range_yn.GetValue() == False and self.cb_Intensity2Range_yn.GetValue() == False and self.cb_Fit_y0_2_Range_yn.GetValue() == False\
        and self.cb_SNR_2Range_yn.GetValue() == False and self.cb_SBR_2Range_yn.GetValue() == False :
        
            self.Loaded_ImageJxpos_2_Range = self.data_2.Loaded_ImageJxpos
            self.Loaded_ImageJypos_2_Range = self.data_2.Loaded_ImageJypos
            self.Loaded_M_2_Range = self.data_2.Loaded_M
            self.Loaded_A_2_Range = self.data_2.Loaded_A
            self.Loaded_T_2_Range = self.data_2.Loaded_T
            self.Loaded_phi_2_Range = self.data_2.Loaded_phi
            self.Loaded_y0_2_Range = self.data_2.Loaded_y0
            self.Loaded_CenterImages_2_Range = self.data_2.Loaded_CenterImages
            self.Loaded_Intensity_2_Range = self.data_2.Loaded_Intensity
            self.Loaded_Intensity_noBGS_2_Range = self.data_2.Loaded_Intensity_noBGS
            
            self.Loaded_Intensity_BG_2_Range = self.data_2.Loaded_Intensity_BG
            self.Loaded_Intensity_BG_Fit_2_Range = self.data_2.Loaded_Intensity_BG_Fit
            self.Loaded_Intensity_Fit_2_Range = self.data_2.Loaded_Intensity_Fit
            self.M_error_2_Range = self.data_2.M_error
            self.SNR_2_Range = self.data_2.SNR
            self.SBR_2_Range = self.data_2.Loaded_SBR





                
                

                    
        else:
            for n in range(len(self.data_2.Loaded_M)):
                
                if self.cb_M2Range_yn.GetValue() == True:
                    if self.data_2.Loaded_M[n] >= float(self.M2RangeValueMin.GetValue()) and self.data_2.Loaded_M[n] <= float(self.M2RangeValueMax.GetValue()):
                        pass
                    else: continue
                
                if self.cb_Intensity2Range_yn.GetValue() == True:
                    if np.mean(self.data_2.Loaded_Intensity[n]) >= float(self.Intensity2RangeValueMin.GetValue()) and np.mean(self.data_2.Loaded_Intensity[n]) <= float(self.Intensity2RangeValueMax.GetValue()):
                        pass
                    else: continue
                
                if self.cb_Fit_y0_2_Range_yn.GetValue() == True:
                    if self.data_2.Loaded_y0[n] >= float(self.Fit_y0_2RangeValueMin.GetValue()) and self.data_2.Loaded_y0[n] <= float(self.Fit_y0_2RangeValueMax.GetValue()):
                        pass
                    else: continue
                
                if self.cb_SNR_2Range_yn.GetValue() == True:
                    if (self.data_2.SNR[n]) >= float(self.SNR_2RangeValueMin.GetValue()) and (self.data_2.SNR[n]) <= float(self.SNR_2RangeValueMax.GetValue()):
                        pass
                    else: continue

                if self.cb_SBR_2Range_yn.GetValue() == True:
                    if (self.data_2.Loaded_SBR[n]) >= float(self.SBR_2RangeValueMin.GetValue()) and (self.data_2.Loaded_SBR[n]) <= float(self.SBR_2RangeValueMax.GetValue()):
                        pass
                    else: continue

               
                self.Loaded_ImageJxpos_2_Range.append(self.data_2.Loaded_ImageJxpos[n])
                self.Loaded_ImageJypos_2_Range.append(self.data_2.Loaded_ImageJypos[n])
                self.Loaded_M_2_Range.append(self.data_2.Loaded_M[n])
                self.Loaded_A_2_Range.append(self.data_2.Loaded_A[n])
                self.Loaded_T_2_Range.append(self.data_2.Loaded_T[n])
                self.Loaded_phi_2_Range.append(self.data_2.Loaded_phi[n])
                self.Loaded_y0_2_Range.append(self.data_2.Loaded_y0[n])
                self.Loaded_CenterImages_2_Range.append(self.data_2.Loaded_CenterImages[n])
                self.Loaded_Intensity_2_Range.append(self.data_2.Loaded_Intensity[n])
                self.Loaded_Intensity_noBGS_2_Range.append(self.data_2.Loaded_Intensity_noBGS[n])
                
                self.Loaded_Intensity_BG_2_Range.append(self.data_2.Loaded_Intensity_BG[n])
                self.Loaded_Intensity_BG_Fit_2_Range.append(self.data_2.Loaded_Intensity_BG_Fit[n])
                                
                self.Loaded_Intensity_Fit_2_Range.append(self.data_2.Loaded_Intensity_Fit[n])
                self.M_error_2_Range.append(self.data_2.M_error[n])
                self.SNR_2_Range.append(self.data_2.SNR[n])
                self.SBR_2_Range.append(self.data_2.Loaded_SBR[n])
















        self.btn_plotDataSeparate_text2.SetLabel('2. #: ' + str(len(self.Loaded_M_2_Range)))




        if len(self.Loaded_ImageJxpos_2_Range) == 0:
            MessageBox(self, '# of data2 is zero', 'Error')
            return


                    
        
        
        self.Loaded_Intensity_Fit_2_Range_CB = np.array(self.Loaded_y0_2_Range) - np.abs(self.Loaded_A_2_Range)
        self.Loaded_Intensity_Fit_2_Range_Max = np.array(self.Loaded_y0_2_Range) + np.abs(self.Loaded_A_2_Range)

        self.Loaded_IntensityNoiseSTD_2_Range = []
        for n in range(len(self.Loaded_Intensity_2_Range)):
            STDtemp2 = np.std( np.array(self.Loaded_Intensity_2_Range[n]) - np.array(self.Loaded_Intensity_Fit_2_Range[n]) )
            print n , ' STDtemp2 ', ' = ' ,STDtemp2
            self.Loaded_IntensityNoiseSTD_2_Range.append(STDtemp2)
        self.SNR_2_Range = np.array(self.Loaded_y0_2_Range)/np.array(self.Loaded_IntensityNoiseSTD_2_Range)

      



        
        plt.subplot(4,6,13)
        plt.title('Median: #2: ' + str(np.around(np.median(self.Loaded_M_2_Range),2)), size = 10)
        plt.hist(self.Loaded_M_2_Range, bins=np.arange(0.0, 1.4, 0.1), label = '#2', alpha=1.0, color='blue')
        plt.xlabel('Modulation Depth (M)', size=11)
        plt.xlim(0.0, 1.31)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        


        
        plt.subplot(4,6,14)
        plt.title('Median: #2: ' + str(np.around(np.median(self.Loaded_phi_2_Range),2)), size = 10)
        plt.hist(self.Loaded_phi_2_Range, label = '#2', alpha=1.0, color='blue')
        plt.xlabel('phi (deg)', size=11)
        #plt.xlim(0.0, 1.31)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        


        
        
        plt.subplot(4,6,15)
        plt.title('Fit: average intensity y0', size=11)
        plt.scatter(self.Loaded_M_2_Range, self.Loaded_y0_2_Range, s = 6,  color='blue', label = '#2')
        plt.xlabel('Modulation Depth M', size=11)
        plt.ylabel('Fit y0', size=11)
        plt.xlim(-0.2, 1.21)
        plt.legend(fontsize=9, loc='upper right',frameon=False)

        
        
        plt.subplot(4,6,17)
        plt.title('NoiseSTD vs M', size=11)
        plt.scatter(self.Loaded_M_2_Range, self.Loaded_IntensityNoiseSTD_2_Range,  s = 6,  color='blue', label = '#2')
        plt.xlabel('M', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        plt.locator_params(nbins=5)
        
        
        
        A = np.vstack([self.Loaded_y0_2_Range, np.ones(len(self.Loaded_y0_2_Range))]).T        
        m, c = np.linalg.lstsq(A, self.Loaded_IntensityNoiseSTD_2_Range)[0]        
        x_temp = np.arange(-1000, 20000, 10)
    
        A_guess = m*10
        B_guess = np.min(self.Loaded_IntensityNoiseSTD_2_Range)
        
        initial_guess = [A_guess, B_guess]
        
        try:
            popt, pcov = curve_fit(Noise_func, self.Loaded_y0_2_Range, self.Loaded_IntensityNoiseSTD_2_Range, p0=initial_guess )
            fit_temp = Noise_func(x_temp, popt[0],popt[1])
            #err_temp = np.sqrt(np.diag(pcov))
        except: 
            print '\n fit error '
            fit_temp = x_temp *0.0
            popt = [0, 0]
            pass
   


        plt.subplot(4,6,18)
        plt.title('NoiseSTD vs Ave Intensity (y0)', size=11)
        plt.scatter(self.Loaded_y0_2_Range, self.Loaded_IntensityNoiseSTD_2_Range,  s = 5,  color='blue', label = '#2')
        plt.plot(x_temp,  m*x_temp + c, 'c', label='Linear Fit: slope='+str(np.around(m,1)), linewidth = 2)
        plt.plot(x_temp, fit_temp, 'm', label='Noise Fit: A,B= '+str(np.around(popt[0],1))+','+str(np.around(popt[1],1)), linewidth = 2)
        plt.xlabel('Ave Intensity', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        plt.xlim(0, np.max(self.Loaded_y0_2_Range))
        plt.ylim(0, np.max(self.Loaded_IntensityNoiseSTD_2_Range))
        plt.locator_params(nbins=5)
        

   
        NetNoiseSTD = []
        for n in range(len(self.Loaded_M_2_Range)):
            dI_temp = Noise_func(self.Loaded_y0_2_Range[n], popt[0],popt[1])
            NetNoiseTemp = self.Loaded_IntensityNoiseSTD_2_Range[n] - dI_temp
            NetNoiseSTD.append(NetNoiseTemp)
            
            
        A = np.vstack([self.Loaded_M_2_Range, np.ones(len(self.Loaded_M_2_Range))]).T        
        m, c = np.linalg.lstsq(A, NetNoiseSTD)[0]        
        x_temp = np.arange(-0, 1.2, 0.1)
        
        
        plt.subplot(4,6,19)
        plt.title('Renormalized NetNoiseSTD vs M', size=11)
        plt.scatter(self.Loaded_M_2_Range, NetNoiseSTD,  s = 6,  color='blue', label = '#2')
        plt.plot(x_temp,  m*x_temp + c, 'c', label='Linear Fit: slope='+str(np.around(m,1)), linewidth = 2)
        plt.xlabel('M', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        plt.locator_params(nbins=5)


        
        plt.subplot(4,6,20)
        plt.title('SBR vs M', size=11)
        plt.scatter(self.Loaded_M_2_Range, self.SBR_2_Range, s = 6,  color='blue', label = '#2')
        plt.xlabel('M', size=11)
        plt.ylabel('SBR', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        plt.locator_params(nbins=5)
        
        

        plt.subplot(4,6,21)
        plt.hist(self.SNR_2_Range, label = '#2', alpha=1.0, color='blue')
        plt.xlabel('SNR', size=11)
        #plt.xlim(0.0, 1.21)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        
        
        plt.subplot(4,6,22)
        plt.title('SNR vs M', size=11)
        plt.scatter(self.Loaded_M_2_Range, self.SNR_2_Range,  s = 6,  color='blue', label = '#2')
        plt.xlabel('M', size=11)
        plt.ylabel('SNR', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        
        
        plt.subplot(4,6,23)
        plt.title('M_error vs SNR', size=11)
        plt.scatter(self.SNR_2_Range, self.M_error_2_Range,  s = 6,  color='blue', label = '#2')
        plt.xlabel('SNR', size=11)
        plt.ylabel('M_error', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        
        


        
        plt.subplot(4,6,24)
        plt.title('M_error vs 1/SNR', size=11)
        plt.scatter(1/self.SNR_2_Range, self.M_error_2_Range,  s = 6,  color='blue', label = '#2')
        plt.xlabel('1/SNR', size=11)
        plt.ylabel('M_error', size=11)
        plt.legend(fontsize=9, loc='upper right',frameon=False)
        plt.locator_params(nbins=6)
        
        
        plt.ion()
        plt.show()
        
        
        


















        
          
'''
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
###########################################################################################################################     
Combining Data Files
'''


               





class DataCombiningMainClass(wx.Frame):
    #ImageFilePath = ''
 
    #----------------------------------------------------------------------
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "Data Combining", size=(1200, 1000), pos=(0,0))
                          
        self.MainPanel = wx.Panel(self, wx.ID_ANY)
        
        



        self.btn_SwitchToDataAnalysis = wx.Button(self.MainPanel, pos=(250,5), label="Switch To Data Analysis")
        self.btn_SwitchToDataAnalysis.Bind(wx.EVT_BUTTON, self.SwitchToDataAnalysis)
        


        btn_reset = wx.Button(self.MainPanel, pos=(800,5), label="Restarting Current GUI")
        btn_reset.Bind(wx.EVT_BUTTON, self.resetting)
            
        
    

        btn_data1 = wx.Button(self.MainPanel, pos=(50,30), label="Open #1 Data File")
        self.DataFilePath_data1 = btn_data1.Bind(wx.EVT_BUTTON, self.onOpenDataFile_data1)
        self.btn_text_data1 = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath_data1), pos=(170, 33))
        self.btn_text_2_data1 = wx.StaticText(self.MainPanel, -1, ' -', pos=(50, 63))
     
     
        
        btn_data2 = wx.Button(self.MainPanel, pos=(50,100), label="Open #2 Data File")
        self.cb_data2_yn = wx.CheckBox(self.MainPanel, -1, '', (30, 105))
        self.cb_data2_yn.SetValue(False)   
        self.DataFilePath_data2 = btn_data2.Bind(wx.EVT_BUTTON, self.onOpenDataFile_data2)
        self.btn_text_data2 = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath_data2), pos=(170, 103))
        self.btn_text_2_data2 = wx.StaticText(self.MainPanel, -1, ' -', pos=(50, 133))
        
        
        btn_data3 = wx.Button(self.MainPanel, pos=(50,170), label="Open #3 Data File")
        self.cb_data3_yn = wx.CheckBox(self.MainPanel, -1, '', (30, 175))
        self.cb_data3_yn.SetValue(False)   
        self.DataFilePath_data3 = btn_data3.Bind(wx.EVT_BUTTON, self.onOpenDataFile_data3)
        self.btn_text_data3 = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath_data3), pos=(170, 173))
        self.btn_text_2_data3 = wx.StaticText(self.MainPanel, -1, ' -', pos=(50, 203))
        
                
        
        btn_data4 = wx.Button(self.MainPanel, pos=(50,240), label="Open #4 Data File")
        self.cb_data4_yn = wx.CheckBox(self.MainPanel, -1, '', (30, 245))
        self.cb_data4_yn.SetValue(False)   
        self.DataFilePath_data4 = btn_data4.Bind(wx.EVT_BUTTON, self.onOpenDataFile_data4)
        self.btn_text_data4 = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath_data4), pos=(170, 243))
        self.btn_text_2_data4 = wx.StaticText(self.MainPanel, -1, ' -', pos=(50, 273))


        
        btn_data5 = wx.Button(self.MainPanel, pos=(50,310), label="Open #5 Data File")
        self.cb_data5_yn = wx.CheckBox(self.MainPanel, -1, '', (30, 315))
        self.cb_data5_yn.SetValue(False)   
        self.DataFilePath_data5 = btn_data5.Bind(wx.EVT_BUTTON, self.onOpenDataFile_data5)
        self.btn_text_data5 = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath_data5), pos=(170, 313))
        self.btn_text_2_data5 = wx.StaticText(self.MainPanel, -1, ' -', pos=(50, 343))
        
        
                
              
        btn_data6 = wx.Button(self.MainPanel, pos=(50,380), label="Open #6 Data File")
        self.cb_data6_yn = wx.CheckBox(self.MainPanel, -1, '', (30, 385))
        self.cb_data6_yn.SetValue(False)   
        self.DataFilePath_data6 = btn_data6.Bind(wx.EVT_BUTTON, self.onOpenDataFile_data6)
        self.btn_text_data6 = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath_data6), pos=(170, 383))
        self.btn_text_2_data6 = wx.StaticText(self.MainPanel, -1, ' -', pos=(50, 413))
        
        
                
              
        btn_data7 = wx.Button(self.MainPanel, pos=(50,450), label="Open #7 Data File")
        self.cb_data7_yn = wx.CheckBox(self.MainPanel, -1, '', (30, 455))
        self.cb_data7_yn.SetValue(False)   
        self.DataFilePath_data7 = btn_data7.Bind(wx.EVT_BUTTON, self.onOpenDataFile_data7)
        self.btn_text_data7 = wx.StaticText(self.MainPanel, -1, str(self.DataFilePath_data7), pos=(170, 453))
        self.btn_text_2_data7 = wx.StaticText(self.MainPanel, -1, ' -', pos=(50, 483))
        
        
                

        
                
      
                                
        
        
        
        
        
        
         

        wx.StaticText(self.MainPanel, -1, "File Name", pos=(50, 650))
        self.FileName = wx.TextCtrl(self.MainPanel, -1, ' ', pos=(50, 670), size=(500, -1))
        
               
        btn_combineData = wx.Button(self.MainPanel, pos=(50,700), label="Combine and Save Data Files")
        btn_combineData.Bind(wx.EVT_BUTTON, self.CombineData)
        self.SaveData_text = wx.StaticText(self.MainPanel, -1, str('_'), pos=(50, 733))
        
      
        self.cb_SaveImageJ_xyData_yn = wx.CheckBox(self.MainPanel, -1, 'ImageJ (x,y) ,   Fit: A, T, phase, y0, M', (300, 705))
        self.cb_SaveImageJ_xyData_yn.SetValue(False)
        
        self.cb_Complete_Data_yn = wx.CheckBox(self.MainPanel, -1, 'Complete Data File', (300, 725))
        self.cb_Complete_Data_yn.SetValue(True)


                        

        
        
        

        #MessageBox(self, 'Opening Data Analysis Window', 'Data Analysis')



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         

    def SwitchToDataAnalysis(self, enven):
        DataAnalysisMain = DataAnalysisMainClass()
        DataAnalysisMain.Show()
        self.Close()
        

     
     
    def resetting(self, event):
        
        print '\n\nResetting'
 

      
        self.Close()                      
        frame2 = DataCombiningMainClass()
        frame2.Show()
          
       
       

    def onOpenDataFile_data1(self, event):
        print  '\n clicked open file 1'
        
        self.SaveData_text.SetLabel(' ')
        
        self.data_1 = OpenDataFile()
        self.btn_text_data1.SetLabel("#1: " + str(self.data_1.fileName))
        
        self.btn_text_2_data1.SetLabel(' # of data: ' + str(len(self.data_1.Loaded_ImageJxpos)) + ' ,   BGS type: ' + str(self.data_1.BGS_type))
        
        
        self.FileName.SetLabel('Combined_' + self.data_1.fileName[:-4][:40])
        print 'self.data_1 ', self.data_1



    def onOpenDataFile_data2(self, event):
        print  '\n clicked open file 2'
        
        self.SaveData_text.SetLabel(' ')
        
        self.data_2 = OpenDataFile()
        self.btn_text_data2.SetLabel("#2: " + str(self.data_2.fileName))
        
        self.btn_text_2_data2.SetLabel(' # of data: ' + str(len(self.data_2.Loaded_ImageJxpos)) + ' ,   BGS type: ' + str(self.data_2.BGS_type))
        

        
        self.cb_data2_yn.SetValue(True)
        
        print 'self.data_2 ', self.data_2




    def onOpenDataFile_data3(self, event):
        print  '\n clicked open file 3'
        
        self.SaveData_text.SetLabel(' ')
        
        self.data_3 = OpenDataFile()
        self.btn_text_data3.SetLabel("#3: " + str(self.data_3.fileName))
        
        self.btn_text_2_data3.SetLabel(' # of data: ' + str(len(self.data_3.Loaded_ImageJxpos)) + ' ,   BGS type: ' + str(self.data_3.BGS_type))
        

        
        self.cb_data3_yn.SetValue(True)
        
        print 'self.data_3 ', self.data_3




    def onOpenDataFile_data4(self, event):
        print  '\n clicked open file 4'
        
        self.SaveData_text.SetLabel(' ')
        
        self.data_4 = OpenDataFile()
        self.btn_text_data4.SetLabel("#4: " + str(self.data_4.fileName))
        
        self.btn_text_2_data4.SetLabel(' # of data: ' + str(len(self.data_4.Loaded_ImageJxpos)) + ' ,   BGS type: ' + str(self.data_4.BGS_type))
        

        
        self.cb_data4_yn.SetValue(True)
        
        print 'self.data_4 ', self.data_4





    def onOpenDataFile_data5(self, event):
        print  '\n clicked open file 5'
        
        self.SaveData_text.SetLabel(' ')
        
        self.data_5 = OpenDataFile()
        self.btn_text_data5.SetLabel("#5: " + str(self.data_5.fileName))
        
        self.btn_text_2_data5.SetLabel(' # of data: ' + str(len(self.data_5.Loaded_ImageJxpos)) + ' ,   BGS type: ' + str(self.data_5.BGS_type))
        

        
        self.cb_data5_yn.SetValue(True)
        
        print 'self.data_5 ', self.data_5





    def onOpenDataFile_data6(self, event):
        print  '\n clicked open file 6'
        
        self.SaveData_text.SetLabel(' ')
        
        self.data_6 = OpenDataFile()
        self.btn_text_data6.SetLabel("#6: " + str(self.data_6.fileName))
        
        self.btn_text_2_data6.SetLabel(' # of data: ' + str(len(self.data_6.Loaded_ImageJxpos)) + ' ,   BGS type: ' + str(self.data_6.BGS_type))
        

        
        self.cb_data6_yn.SetValue(True)
        
        print 'self.data_6 ', self.data_6





    def onOpenDataFile_data7(self, event):
        print  '\n clicked open file 7'
        
        self.SaveData_text.SetLabel(' ')
        
        self.data_7 = OpenDataFile()
        self.btn_text_data7.SetLabel("#7: " + str(self.data_7.fileName))
        
        self.btn_text_2_data7.SetLabel(' # of data: ' + str(len(self.data_7.Loaded_ImageJxpos)) + ' ,   BGS type: ' + str(self.data_7.BGS_type))
        

        
        self.cb_data7_yn.SetValue(True)
        
        print 'self.data_7 ', self.data_7












    def CombineData(self, event):
        print ''
        
        try:
            self.data_1.BGS_type
        except:
            MessageBox(self, 'No #1 Data', 'Error')
            return
            


    

        self.Combined_ImageJxpos = []
        self.Combined_ImageJypos = []
        self.Combined_A = []
        self.Combined_T = []
        self.Combined_phi = []
        self.Combined_y0 = []
        self.Combined_M = []
        
        self.Combined_Intensity_noBGS = []
        self.Combined_Intensity_BGS = []
        self.Combined_Intensity_BGS_Fit = []
        self.Combined_Intensity_BG = []
        self.Combined_Intensity_BG_Fit = []
        self.Combined_CenterImages = []
        
            

        self.Combined_ImageJxpos.extend(self.data_1.Loaded_ImageJxpos)
        self.Combined_ImageJypos.extend(self.data_1.Loaded_ImageJypos)
        self.Combined_A.extend(self.data_1.Loaded_A)
        self.Combined_T.extend(self.data_1.Loaded_T)
        self.Combined_phi.extend(self.data_1.Loaded_phi)
        self.Combined_y0.extend(self.data_1.Loaded_y0)
        self.Combined_M.extend(self.data_1.Loaded_M)
        
        self.Combined_Intensity_noBGS.extend(self.data_1.Loaded_Intensity_noBGS)
        self.Combined_Intensity_BGS.extend(self.data_1.Loaded_Intensity)
        self.Combined_Intensity_BGS_Fit.extend(self.data_1.Loaded_Intensity_Fit)
        self.Combined_Intensity_BG.extend(self.data_1.Loaded_Intensity_BG)
        self.Combined_Intensity_BG_Fit.extend(self.data_1.Loaded_Intensity_BG_Fit)
        self.Combined_CenterImages.extend(self.data_1.Loaded_CenterImages)
        
        
        
        SaveData_text = 'Data #1'
        
        if self.cb_data2_yn.GetValue():
            if self.data_1.BGS_type != self.data_2.BGS_type:
                MessageBox(self, 'BGS types are not matching', 'Error' )
                print '\n\n '
                print 'self.data_1.BGS_type ', self.data_1.BGS_type
                print 'self.data_2.BGS_type ', self.data_2.BGS_type
                return
            
            self.Combined_ImageJxpos.extend(self.data_2.Loaded_ImageJxpos)
            self.Combined_ImageJypos.extend(self.data_2.Loaded_ImageJypos)
            self.Combined_A.extend(self.data_2.Loaded_A)
            self.Combined_T.extend(self.data_2.Loaded_T)
            self.Combined_phi.extend(self.data_2.Loaded_phi)
            self.Combined_y0.extend(self.data_2.Loaded_y0)
            self.Combined_M.extend(self.data_2.Loaded_M)
            
            self.Combined_Intensity_noBGS.extend(self.data_2.Loaded_Intensity_noBGS)
            self.Combined_Intensity_BGS.extend(self.data_2.Loaded_Intensity)
            self.Combined_Intensity_BGS_Fit.extend(self.data_2.Loaded_Intensity_Fit)
            self.Combined_Intensity_BG.extend(self.data_2.Loaded_Intensity_BG)
            self.Combined_Intensity_BG_Fit.extend(self.data_2.Loaded_Intensity_BG_Fit)
            self.Combined_CenterImages.extend(self.data_2.Loaded_CenterImages)
            
            SaveData_text += '  +  Data #2'
            
       
        
        if self.cb_data3_yn.GetValue():
            if self.data_1.BGS_type != self.data_3.BGS_type:
                MessageBox(self, 'BGS types are not matching', 'Error' )
                print '\n\n '
                print 'self.data_1.BGS_type ', self.data_1.BGS_type
                print 'self.data_3.BGS_type ', self.data_3.BGS_type
                return
                
            self.Combined_ImageJxpos.extend(self.data_3.Loaded_ImageJxpos)
            self.Combined_ImageJypos.extend(self.data_3.Loaded_ImageJypos)
            self.Combined_A.extend(self.data_3.Loaded_A)
            self.Combined_T.extend(self.data_3.Loaded_T)
            self.Combined_phi.extend(self.data_3.Loaded_phi)
            self.Combined_y0.extend(self.data_3.Loaded_y0)
            self.Combined_M.extend(self.data_3.Loaded_M)
            
            self.Combined_Intensity_noBGS.extend(self.data_3.Loaded_Intensity_noBGS)
            self.Combined_Intensity_BGS.extend(self.data_3.Loaded_Intensity)
            self.Combined_Intensity_BGS_Fit.extend(self.data_3.Loaded_Intensity_Fit)
            self.Combined_Intensity_BG.extend(self.data_3.Loaded_Intensity_BG)
            self.Combined_Intensity_BG_Fit.extend(self.data_3.Loaded_Intensity_BG_Fit)
            self.Combined_CenterImages.extend(self.data_3.Loaded_CenterImages)
            
            SaveData_text += '  +  Data #3'
            
            
            
        
        if self.cb_data4_yn.GetValue():
            if self.data_1.BGS_type != self.data_4.BGS_type:
                MessageBox(self, 'BGS types are not matching', 'Error' )
                print '\n\n '
                print 'self.data_1.BGS_type ', self.data_1.BGS_type
                print 'self.data_4.BGS_type ', self.data_4.BGS_type
                return
                
            self.Combined_ImageJxpos.extend(self.data_4.Loaded_ImageJxpos)
            self.Combined_ImageJypos.extend(self.data_4.Loaded_ImageJypos)
            self.Combined_A.extend(self.data_4.Loaded_A)
            self.Combined_T.extend(self.data_4.Loaded_T)
            self.Combined_phi.extend(self.data_4.Loaded_phi)
            self.Combined_y0.extend(self.data_4.Loaded_y0)
            self.Combined_M.extend(self.data_4.Loaded_M)
            
            self.Combined_Intensity_noBGS.extend(self.data_4.Loaded_Intensity_noBGS)
            self.Combined_Intensity_BGS.extend(self.data_4.Loaded_Intensity)
            self.Combined_Intensity_BGS_Fit.extend(self.data_4.Loaded_Intensity_Fit)
            self.Combined_Intensity_BG.extend(self.data_4.Loaded_Intensity_BG)
            self.Combined_Intensity_BG_Fit.extend(self.data_4.Loaded_Intensity_BG_Fit)
            self.Combined_CenterImages.extend(self.data_4.Loaded_CenterImages)
            
            SaveData_text += '  +  Data #4'
            
            
            
            
            
        
        if self.cb_data5_yn.GetValue():
            if self.data_1.BGS_type != self.data_5.BGS_type:
                MessageBox(self, 'BGS types are not matching', 'Error' )
                print '\n\n '
                print 'self.data_1.BGS_type ', self.data_1.BGS_type
                print 'self.data_5.BGS_type ', self.data_5.BGS_type
                return
                
            self.Combined_ImageJxpos.extend(self.data_5.Loaded_ImageJxpos)
            self.Combined_ImageJypos.extend(self.data_5.Loaded_ImageJypos)
            self.Combined_A.extend(self.data_5.Loaded_A)
            self.Combined_T.extend(self.data_5.Loaded_T)
            self.Combined_phi.extend(self.data_5.Loaded_phi)
            self.Combined_y0.extend(self.data_5.Loaded_y0)
            self.Combined_M.extend(self.data_5.Loaded_M)
            
            self.Combined_Intensity_noBGS.extend(self.data_5.Loaded_Intensity_noBGS)
            self.Combined_Intensity_BGS.extend(self.data_5.Loaded_Intensity)
            self.Combined_Intensity_BGS_Fit.extend(self.data_5.Loaded_Intensity_Fit)
            self.Combined_Intensity_BG.extend(self.data_5.Loaded_Intensity_BG)
            self.Combined_Intensity_BG_Fit.extend(self.data_5.Loaded_Intensity_BG_Fit)
            self.Combined_CenterImages.extend(self.data_5.Loaded_CenterImages)
            
            SaveData_text += '  +  Data #5'
            
            
            
            
        
        if self.cb_data6_yn.GetValue():
            if self.data_1.BGS_type != self.data_6.BGS_type:
                MessageBox(self, 'BGS types are not matching', 'Error' )
                print '\n\n '
                print 'self.data_1.BGS_type ', self.data_1.BGS_type
                print 'self.data_6.BGS_type ', self.data_6.BGS_type
                return
                
            self.Combined_ImageJxpos.extend(self.data_6.Loaded_ImageJxpos)
            self.Combined_ImageJypos.extend(self.data_6.Loaded_ImageJypos)
            self.Combined_A.extend(self.data_6.Loaded_A)
            self.Combined_T.extend(self.data_6.Loaded_T)
            self.Combined_phi.extend(self.data_6.Loaded_phi)
            self.Combined_y0.extend(self.data_6.Loaded_y0)
            self.Combined_M.extend(self.data_6.Loaded_M)
            
            self.Combined_Intensity_noBGS.extend(self.data_6.Loaded_Intensity_noBGS)
            self.Combined_Intensity_BGS.extend(self.data_6.Loaded_Intensity)
            self.Combined_Intensity_BGS_Fit.extend(self.data_6.Loaded_Intensity_Fit)
            self.Combined_Intensity_BG.extend(self.data_6.Loaded_Intensity_BG)
            self.Combined_Intensity_BG_Fit.extend(self.data_6.Loaded_Intensity_BG_Fit)
            self.Combined_CenterImages.extend(self.data_6.Loaded_CenterImages)
            
            SaveData_text += '  +  Data #6'
            
            
            
            
        
        if self.cb_data7_yn.GetValue():
            if self.data_1.BGS_type != self.data_7.BGS_type:
                MessageBox(self, 'BGS types are not matching', 'Error' )
                print '\n\n '
                print 'self.data_1.BGS_type ', self.data_1.BGS_type
                print 'self.data_7.BGS_type ', self.data_7.BGS_type
                return
                
            self.Combined_ImageJxpos.extend(self.data_7.Loaded_ImageJxpos)
            self.Combined_ImageJypos.extend(self.data_7.Loaded_ImageJypos)
            self.Combined_A.extend(self.data_7.Loaded_A)
            self.Combined_T.extend(self.data_7.Loaded_T)
            self.Combined_phi.extend(self.data_7.Loaded_phi)
            self.Combined_y0.extend(self.data_7.Loaded_y0)
            self.Combined_M.extend(self.data_7.Loaded_M)
            
            self.Combined_Intensity_noBGS.extend(self.data_7.Loaded_Intensity_noBGS)
            self.Combined_Intensity_BGS.extend(self.data_7.Loaded_Intensity)
            self.Combined_Intensity_BGS_Fit.extend(self.data_7.Loaded_Intensity_Fit)
            self.Combined_Intensity_BG.extend(self.data_7.Loaded_Intensity_BG)
            self.Combined_Intensity_BG_Fit.extend(self.data_7.Loaded_Intensity_BG_Fit)
            self.Combined_CenterImages.extend(self.data_7.Loaded_CenterImages)
            
            SaveData_text += '  +  Data #7'
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        filedirectoryOnly = self.data_1.fileDirectory
        filename = self.FileName.GetLabel()
        
        xpos_ImageJ = self.Combined_ImageJxpos
        ypos_ImageJ = self.Combined_ImageJypos
        A = self.Combined_A
        T = self.Combined_T
        phi = self.Combined_phi
        y0 = self.Combined_y0
        M = self.Combined_M
        
        I_noBGS = self.Combined_Intensity_noBGS
        I_BGS = self.Combined_Intensity_BGS
        I_BGS_fit = self.Combined_Intensity_BGS_Fit
        I_BG = self.Combined_Intensity_BG
        I_BG_fit = self.Combined_Intensity_BG_Fit
        CenterImages = self.Combined_CenterImages

        
        cb_SaveImageJ_xyData_yn = self.cb_SaveImageJ_xyData_yn.GetValue()
        cb_Complete_Data_yn = self.cb_Complete_Data_yn.GetValue()
        
        
        
        BGS_type = self.data_1.BGS_type
        
        
        
        if cb_SaveImageJ_xyData_yn or cb_Complete_Data_yn:
            pass
        else:
            MessageBox(self, 'Check at least a box to save', 'Error')
            return


        CreateDataFiles(filedirectoryOnly, filename, xpos_ImageJ, ypos_ImageJ, A, T, phi, y0, M, I_noBGS, I_BGS, I_BGS_fit, I_BG, I_BG_fit, CenterImages, cb_SaveImageJ_xyData_yn, cb_Complete_Data_yn, BGS_type)

            
        
        self.SaveData_text.SetLabel('Combined Data Saved \n\n' + SaveData_text + '\n\nTotal # of data: ' + str(len(xpos_ImageJ)))
        
        
        
        
        
        
        
        
        
        
        
        
  













#----------------------------------------------------------------------
# Run the program
   
if __name__ == "__main__":
    app = wx.App(False)
    
    if StartingPanel == 2:
        frame = TwoChannelMainTracking()
        frame.Show()
    elif StartingPanel == 1:
        frame = SingleChannelMainTracking()
        frame.Show()

    elif StartingPanel == 3:
        frame = DataAnalysisMainClass()
        frame.Show()

    else:
        print 'Nothing'
        pass
        
    app.MainLoop()
    
    
    
        
        

#MessageBox(self, 'Opening Data Analysis Window', 'Message')



       
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
