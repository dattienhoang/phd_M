# -*- coding: utf-8 -*-
"""
Created on Tue May 19 21:29:08 2015

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

import pandas as pd
import trackpy as tp
import pims


import os
import time
import math

import wx
from sys import exit

from pandas import DataFrame, Series  # for convenience


wildcard = "Tiff file (*.tif; *.tiff)|*.tif;*.tiff|" \
         "All files (*.*)|*.*"
 
wildcardDataFile = "txt file (*.txt)|*.txt|" \
         "All files (*.*)|*.*"




#######################################################
# feature find default values
ThrethholdIntensity = 5000000         
PixelSize = 7 # typically 7


ColorScaleMin = 1100
ColorScaleMax = 2200

FilePrefixInput = ''

#######################################################




#######################################################
# M analysis default values  ##########################
FitPeriodGuess = 90   # in frames

ThresholdAveIntensity = 0

MinRsquareFitGoodness = 0.7 # in %




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





manual_openFile = 'not yet'



manual_showM = 'not yet'

manual_saveData = 'not yet'








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
    
    




    


  
def IntensityTrajectoryPlotsFuncForM(framesLoaded, filename, xpos, ypos, IntensityMethod, centerPixelN, edgePixelN):
    print '\n Starting IntensityTrajectoryPlotsFunc for M' 
    ################################################
    x_ImageJ, y_ImageJ = np.array(xpos), np.array(ypos)
    


    print 'x_ImageJ', x_ImageJ
    print 'y_ImageJ', y_ImageJ
    print 'intensity method', IntensityMethod
    
    print 'filename = ', filename
    
    frames = framesLoaded
    
    Nframes = len(frames)
    print 'Nframes = ', Nframes
    
    N_xPixel = frames.frame_shape[0]
    N_yPixel = frames.frame_shape[1]
    
    
    print N_xPixel
    print 'type(N_xPixel)', type(N_xPixel)    
    print N_yPixel
    print 'type(N_yPixel)', type(N_yPixel)
    
    #xc = y_ImageJ
    #yc = N_xPixel - x_ImageJ -1

    xc = x_ImageJ # for trackpy 3.0
    yc = y_ImageJ # for trackpy 3.0
            
    ################################################




    
    
    frameStart = 0
    frameEnd = len(frames) - 1
    
    
    
    print '# of frames ', len(frames)    
    


    centerHalfPixelN = (centerPixelN-1)/2 
    edgeHalfPixelN = (edgePixelN-1)/2
    
    NofEdgePixels = edgePixelN**2 - (edgePixelN-2)**2
    
    
    # calculating each frame for all features => 
    #centerIntensityAveTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))] 
    centerIntensityTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    edgeIntensityTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    #backgroundIntensityMedianTrajectory = [np.zeros(Nframes) for _ in xrange(len(xc))]
    
    
    try:
        edgeArray = frames[0][yc[0]-edgeHalfPixelN:yc[0]+edgeHalfPixelN+1, xc[0]-edgeHalfPixelN:xc[0]+edgeHalfPixelN+1] * 0
        print 'edgeArray zeros try \n', edgeArray
        for n in range(edgePixelN):
            for k in range(edgePixelN):
                if n == 0 or n == edgePixelN-1 or k == 0 or k==edgePixelN-1:
                    edgeArray[n][k] = 1
        
    except:
        edgeArray = frames[0][yc[1]-edgeHalfPixelN:yc[1]+edgeHalfPixelN+1, xc[1]-edgeHalfPixelN:xc[1]+edgeHalfPixelN+1] * 0
        print 'edgeArray zeros try \n', edgeArray
        for n in range(edgePixelN):
            for k in range(edgePixelN):
                if n == 0 or n == edgePixelN-1 or k == 0 or k==edgePixelN-1:
                    edgeArray[n][k] = 1
                    
                    
                    
    print 'edgeArray \n ' , edgeArray
    
    

    
    
    for k in range(frameStart,frameEnd+1):
        print 'frame # ', k
        frametemp = frames[k]
        for n in range(len(xc)):
            
            print 'yc[n]', yc[n]
            print 'xc[n]', xc[n]
            
            
            frameLargeTemp = frametemp[yc[n]-edgeHalfPixelN:yc[n]+edgeHalfPixelN+1, xc[n]-edgeHalfPixelN:xc[n]+edgeHalfPixelN+1]
            frameCenterTemp = frametemp[yc[n]-centerHalfPixelN:yc[n]+centerHalfPixelN+1, xc[n]-centerHalfPixelN:xc[n]+centerHalfPixelN+1]


            
            print '\n\nframeLargeTemp \n', frameLargeTemp
            print '\nframeCenterTemp \n', frameCenterTemp


            

            #frameCenterTemp2 = frametemp[yc[n]-1:yc[n]+2, xc[n]-1:xc[n]+2]
            
            
            #centerIntensityAveTrajectory[n][k] = (np.sum(frameCenterTemp))/(centerPixelN**2)
            
            try:
                frameBackgroundTemp = frameLargeTemp * edgeArray
            except:
                frameBackgroundTemp = frameLargeTemp * 0
                
            
            #print 'frameBackgroundTemp \n', frameBackgroundTemp
            
            #print 'B_Max5Ave: ', int(np.mean( np.sort(frameBackgroundTemp.ravel())[-5:]) )
            
            
            
            
            if IntensityMethod == 0:
                centerIntensityTrajectory[n][k] = int(np.mean(frameCenterTemp))
                edgeIntensityTrajectory[n][k] = int(np.mean( np.sort(frameBackgroundTemp.ravel())[-NofEdgePixels:]) ) #  
                
            
            elif IntensityMethod == 1:
                centerIntensityTrajectory[n][k] = int(np.median(frameCenterTemp))
                edgeIntensityTrajectory[n][k] = int(np.median( np.sort(frameBackgroundTemp.ravel())[-NofEdgePixels:]) ) # 
                
            
            elif IntensityMethod == 2:
                centerIntensityTrajectory[n][k] = int(np.amax(frameCenterTemp))
                edgeIntensityTrajectory[n][k] = int(np.amax( np.sort(frameBackgroundTemp.ravel())[-NofEdgePixels:]) ) #  
                
                
                
    print 'End of Intensity Trajectory func \n'
    
    
    return centerIntensityTrajectory, edgeIntensityTrajectory
    

      
      
      
      
      
      
      
      
      
      
      
      
  
  


     
class Aggregate_M_analysis(wx.Frame):
    #ImageFilePath = ''
 
    #----------------------------------------------------------------------
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY,
                          "Aggregate_M_analysis", size=(810, 1000), pos=(0,0))
                          
        self.MainPanel = wx.Panel(self, wx.ID_ANY)
        
        self.MainPanel.Bind(wx.EVT_PAINT, self.on_paint)
        
        
        
        
        self.btn_Manual = wx.Button(self.MainPanel, pos=(700,5), label="Manual")
        self.btn_Manual.Bind(wx.EVT_BUTTON, self.ShowManual)
        self.btn_Manual.Hide()
        
        




                

                
        btn2 = wx.Button(self.MainPanel, pos=(50,30), label="Open Image File")
        self.ImageFilePath = btn2.Bind(wx.EVT_BUTTON, self.onOpenImageFile)
        self.btn2text = wx.StaticText(self.MainPanel, -1, str(self.ImageFilePath), pos=(170, 33))
        
        

        
        self.btn_AveragingFrames = wx.Button(self.MainPanel, pos=(50,70), label="Averaging Frames")
        self.btn_AveragingFrames.Bind(wx.EVT_BUTTON, self.AveragingFrames)
        self.btn_AveragingFrames.Hide()
        
        
        
        
        
        wx.StaticText(self.MainPanel, -1, "From", pos=(180, 75))
        self.averagingFrameN1_begin = wx.SpinCtrl(self.MainPanel, -1, pos=(220, 70), size=(50,-1))
        self.averagingFrameN1_begin.SetRange(0, 9999999)  
        self.averagingFrameN1_begin.SetValue(0)
        wx.StaticText(self.MainPanel, -1, "---     To", pos=(300, 75))
        self.averagingFrameN1_end = wx.SpinCtrl(self.MainPanel, -1, pos=(360, 70), size=(50,-1))
        self.averagingFrameN1_end.SetRange(0, 9999999)  
        self.averagingFrameN1_end.SetValue(0)
        
        
        
        
        self.btn_PixelHistogram = wx.Button(self.MainPanel, pos=(450,40), label="Pixel Histogram")
        self.btn_PixelHistogram.Bind(wx.EVT_BUTTON, self.PixelHistogram)
        self.btn_PixelHistogram.Hide()
        
        wx.StaticText(self.MainPanel, -1, "Frame #:", pos=(550, 43))
        self.PixelHistogram_FrameN = wx.SpinCtrl(self.MainPanel, -1, pos=(600, 40), size=(50,-1))
        self.PixelHistogram_FrameN.SetRange(0, 9999999)  
        self.PixelHistogram_FrameN.SetValue(0)
                
        
        self.btn_PixelHistogramSaveData = wx.Button(self.MainPanel, pos=(680,40), label="Save Data")
        self.btn_PixelHistogramSaveData.Bind(wx.EVT_BUTTON, self.SaveData_PixelHistogram)
        self.btn_PixelHistogramSaveData.Hide()
        
        
        self.text_PixelHistogramSaveData = wx.StaticText(self.MainPanel, -1, "--  ", pos=(765, 43))
        
        
                
        
        
        
        
        self.btn_CalculateMedianAve = wx.Button(self.MainPanel, pos=(450,70), label="Calculate Median, Average Trajectories")
        self.btn_CalculateMedianAve.Bind(wx.EVT_BUTTON, self.CalculateMedianAve)
        self.btn_CalculateMedianAve.Hide()
        
        
        
        self.btn_CalculateMedianAveSaveData = wx.Button(self.MainPanel, pos=(680,70), label="Save a data file")
        self.btn_CalculateMedianAveSaveData.Bind(wx.EVT_BUTTON, self.SaveData_CalculateMedianAve)
        self.btn_CalculateMedianAveSaveData.Hide()
        
        self.text_CalculateMedianAveSaveFile = wx.StaticText(self.MainPanel, -1, "--", pos=(700, 93))
        
        
        





        self.btn_ShowBackgroundImage = wx.Button(self.MainPanel, pos=(50,110), label="Show Background Image")
        self.btn_ShowBackgroundImage.Bind(wx.EVT_BUTTON, self.showBackgroundImage)
        self.btn_ShowBackgroundImage.Hide()
        wx.StaticText(self.MainPanel, -1, 'Enter center (x, y) for background\n from the Averaging Frames window', pos=(200, 110))
        self.BG_x = wx.SpinCtrl(self.MainPanel, -1, "234", pos=(390, 110), size=(50,-1))
        self.BG_y = wx.SpinCtrl(self.MainPanel, -1, "169", pos=(450, 110), size=(50,-1))
        wx.StaticText(self.MainPanel, -1, ',        # of pixels (odd number): ', pos=(510, 113))
        self.BG_pixel_size = wx.SpinCtrl(self.MainPanel, -1, pos=(670, 110), size=(50,-1))
        self.BG_pixel_size.SetRange(0, 9999999)  
        self.BG_pixel_size.SetValue(45)





        self.btn_ShowBackgroundIntensity = wx.Button(self.MainPanel, pos=(50,150), label="Calculate Global Background Intensity")
        self.btn_ShowBackgroundIntensity.Bind(wx.EVT_BUTTON, self.ShowBackgroundIntensity)
        self.btn_ShowBackgroundIntensity.Hide()
        wx.StaticText(self.MainPanel, -1, 'Dark Counts : ', pos=(300, 153))
        self.DartCount = wx.SpinCtrl(self.MainPanel, -1, pos=(380, 150), size=(60,-1))
        self.DartCount.SetRange(0, 9999999)  
        self.DartCount.SetValue(1100)
        
        
        
        
        
        
        
        

        wx.StaticText(self.MainPanel, -1, '             (x, y) for M analysis\n from the Averaging Frames window', pos=(10, 200))
        self.MorePairs_x = [[] for _ in xrange(15)] 
        self.MorePairs_y = [[] for _ in xrange(15)]
        for n in range(len(self.MorePairs_x)):
            ColN = n%5
            RowN = int(n/5)
            wx.StaticText(self.MainPanel, -1, str(n) + '.', pos=(210 + 100*ColN, 203 + 30*RowN))
            self.MorePairs_x[n] = wx.TextCtrl(self.MainPanel, -1, "-", pos=(230 + 100*ColN, 200 + 30*RowN), size=(30,-1))
            self.MorePairs_y[n] = wx.TextCtrl(self.MainPanel, -1, "-", pos=(263 + 100*ColN, 200 + 30*RowN), size=(30,-1))
            
            
        
        


        
        



        self.btn_ShowRawImages = wx.Button(self.MainPanel, pos=(50,300), label="Show Raw Images")
        self.btn_ShowRawImages.Bind(wx.EVT_BUTTON, self.ShowRawImages)
        self.btn_ShowRawImages.Hide()
        wx.StaticText(self.MainPanel, -1, "Raw Data Images F#:", pos=(180, 304))
        self.rawImageFN = wx.SpinCtrl(self.MainPanel, -1, pos=(300, 300), size=(50,-1))
        self.rawImageFN.SetRange(0, 100000)  
        self.rawImageFN.SetValue(0)
        




        self.radio_IntensityMethod = wx.RadioBox(self.MainPanel, -1, choices=['Average', 'Median', 'Max 1'], label='Intensity Method', pos=(50, 400), size=wx.Size(120, 80), style=wx.RA_SPECIFY_ROWS)
        self.radio_IntensityMethod.SetSelection(1)
        
        self.radio_BSMethod = wx.RadioBox(self.MainPanel, -1, choices=['Global frame by frame', 'Individual frame by frame'], label='Background Subtraction Method', pos=(200, 400), size=wx.Size(180, 98), style=wx.RA_SPECIFY_ROWS)
        self.radio_BSMethod.SetSelection(0)
        
        
                
        
        wx.StaticText(self.MainPanel, -1, 'Center Pixel #', pos=(50, 513))  
        self.centerPixelN = wx.SpinCtrl(self.MainPanel, -1, pos=(130, 510), size=(50,-1))
        self.centerPixelN.SetRange(1, 7)  
        self.centerPixelN.SetValue(3)


        
        wx.StaticText(self.MainPanel, -1, 'Edge Pixel #', pos=(250, 513))  
        self.edgePixelN = wx.SpinCtrl(self.MainPanel, -1, pos=(330, 510), size=(50,-1))
        self.edgePixelN.SetRange(9, 49)  
        self.edgePixelN.SetValue(21)


                
        


        
        self.cb_Excitation_Pol_correction_yn = wx.CheckBox(self.MainPanel, -1, 'Excitation Polarization Correction by', (50, 450))
        self.cb_Excitation_Pol_correction_yn.SetValue(False)
        self.cb_Excitation_Pol_correction_yn.Hide()
        self.Excitation_Pol_correction_factor = wx.TextCtrl(self.MainPanel, -1, "0.1", pos=(280, 476), size=(30, -1))
        self.Excitation_Pol_correction_factor.Hide()


        
        
        
        
                
        wx.StaticText(self.MainPanel, -1, 'Fit initial guess for period in frames:', pos=(466, 453))  
        self.fitInitialGuess_T = wx.TextCtrl(self.MainPanel, -1, str(FitPeriodGuess), pos=(650, 450), size=(30, -1))

        

        self.btn_ShowAnisotropyM = wx.Button(self.MainPanel, pos=(50,550), label="Show Anisotropy M")
        self.btn_ShowAnisotropyM.Bind(wx.EVT_BUTTON, self.ShowAnisotropyM)
        self.btn_ShowAnisotropyM.Hide()
        #wx.StaticText(self.MainPanel, -1, "Anisotropy M & Intensity", pos=(180, 554))
        
        
        wx.StaticText(self.MainPanel, -1, 'F#s: ', pos=(170, 553))  
        self.rawImageFN1 = wx.SpinCtrl(self.MainPanel, -1, pos=(200, 550), size=(50,-1))
        self.rawImageFN1.SetRange(0, 100000)  
        self.rawImageFN1.SetValue(0)
        
        self.rawImageFN2 = wx.SpinCtrl(self.MainPanel, -1, pos=(260, 550), size=(50,-1))
        self.rawImageFN2.SetRange(0, 100000)  
        self.rawImageFN2.SetValue(49)
        
        
        
        
        self.cb_includeFit_yn = wx.CheckBox(self.MainPanel, -1, 'Include Fit Graphs', (50, 580))
        self.cb_includeFit_yn.SetValue(True)

        self.cb_BG_subtraction_yn = wx.CheckBox(self.MainPanel, -1, 'Background Subtraction', (50, 600))
        self.cb_BG_subtraction_yn.SetValue(True)

        
        self.text_NFeatures = wx.StaticText(self.MainPanel, -1, ' ', pos=(50, 650)) 
        
        
        

        
                
        wx.StaticText(self.MainPanel, -1, "Data Files Prefix", pos=(50, 702))
        self.FilePrefix = wx.TextCtrl(self.MainPanel, -1, str(FilePrefixInput), pos=(150, 700), size=(180, -1))
        
        self.btn_Save_data_figures = wx.Button(self.MainPanel, pos=(50,740), label="Save Data")  
        self.btn_Save_data_figures.Hide()
        self.btn_Save_data_figures.Bind(wx.EVT_BUTTON, self.SaveDataAndFigures)
        self.text_Save_figures = wx.StaticText(self.MainPanel, -1, ' ', pos=(170, 790)) 
        
        
        self.cb_SaveFigs_yn = wx.CheckBox(self.MainPanel, -1, 'Figures', (350,722))
        self.cb_SaveFigs_yn.SetValue(True)
        self.cb_SaveImageJ_xyData_yn = wx.CheckBox(self.MainPanel, -1, 'ImageJ x,y data & Fit: A, M, T, phase', (350, 742))
        self.cb_SaveImageJ_xyData_yn.SetValue(False)
        self.cb_SaveIntensityData_yn = wx.CheckBox(self.MainPanel, -1, 'Intensity Trajectory Data for all', (350, 762))
        self.cb_SaveIntensityData_yn.SetValue(False)
        self.cb_SaveIntensityFitData_yn = wx.CheckBox(self.MainPanel, -1, 'Fit Trajectory Data for all', (350, 782))
        self.cb_SaveIntensityFitData_yn.SetValue(False)



        self.cb_Complete_Data_yn = wx.CheckBox(self.MainPanel, -1, 'Complete Data File', (350, 842))
        self.cb_Complete_Data_yn.SetValue(True)













     #----------------------------------------------------------------------


        

    def on_paint(self, event):
        dc = wx.PaintDC(event.GetEventObject())
        dc.Clear()
        dc.SetPen(wx.Pen("BLACK", 4))
        dc.DrawLine(0, 340, 800, 340)  
        dc.DrawLine(0, 690, 800, 690)  
        
        
        
        
    

    

    def ShowManual(self, event):
        self.manual_window = wx.Frame(None, title='Manual ', size=(800, 1000), pos=(800,0))
        self.manual_window.Show()
        
        wx.StaticText(self.manual_window, -1, manual_openFile, pos=(50, 30)) 
        wx.StaticText(self.manual_window, -1, manual_showM, pos=(50, 350)) 
        wx.StaticText(self.manual_window, -1, manual_saveData, pos=(50, 700)) 
        
        




        
        
    def onOpenImageFile(self, event):
        
        print '\n\n starting open a image file'


        
        self.btn_AveragingFrames.Hide()
        self.btn_ShowBackgroundImage.Hide()
        self.btn_ShowBackgroundIntensity.Hide()
        self.btn_ShowRawImages.Hide()
        self.btn_ShowAnisotropyM.Hide()
        self.btn_Save_data_figures.Hide()
        self.btn_PixelHistogram.Hide()
        self.btn_PixelHistogramSaveData.Hide()
        self.btn_CalculateMedianAve.Hide()
        self.btn_CalculateMedianAveSaveData.Hide()
        
        
        self.text_Save_figures.SetLabel(' ')
        self.text_NFeatures.SetLabel(' ')
        self.text_CalculateMedianAveSaveFile.SetLabel('--')
        self.text_PixelHistogramSaveData.SetLabel('--  ')
        



        tp_ver = tp.__version__[:3]
        if float(tp_ver[:3]) != 0.3:
            MessageBox(self, 'The trackpy version is not 0.3.#\nUpdate trackpy', 'Error')
            return





        
        try:
            plt.close(self.fig_ShowAfewFrames)
        except: pass
        
        
        try:
            plt.close(self.fig_AveragingFrames)
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
    
        try: plt.close(self.fig_MedianAverageIntensity)
        except: pass
            


  
                
        """
        Create and show the Open FileDialog
        """
          
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
        
        else:
            print 'no file is selected'
            return
            
        
        
        
        
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
        
        self.ShowAFewFrames()
        
        
        
        
        
        
        
        
    
    def ShowAFewFrames(self):
                
        frames = pims.TiffStack(self.filepath)
                
        self.framesLoaded = frames

        
        #self.TotalNframes = frames._count
        
        self.TotalNframes = len(frames)
        
        
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
            #plt.hold(True)
            #plt.imshow(np.rot90(frames[n],3), origin='upper')
            plt.imshow(frames[n], origin='upper')  # for trackpy 3.0
            #plt.imshow(frames[n], vmin=1000, vmax=16000)
        
        
        self.N_xPixel = frames.frame_shape[0]
        self.N_yPixel = frames.frame_shape[1]
        
        

        
        
        
        self.averagingFrameN1_begin.SetRange(0, self.TotalNframes - 1)  
        self.averagingFrameN1_end.SetRange(0, self.TotalNframes - 1)  
        
        self.BG_x.SetRange(0, self.N_xPixel - 1)  
        self.BG_y.SetRange(0, self.N_yPixel - 1)  
        
        

        plt.ion()
        plt.show()
        
        
        self.btn_AveragingFrames.Show()
        self.btn_CalculateMedianAve.Show()
        self.btn_PixelHistogram.Show()
  
        
        
        
        
     
     
     
     
        
        

    def AveragingFrames(self, event):
        


        frames = self.framesLoaded
        
        self.aveFrameN1_begin = int(self.averagingFrameN1_begin.GetValue())
        self.aveFrameN1_end = int(self.averagingFrameN1_end.GetValue())
        
        
        if self.aveFrameN1_begin > self.aveFrameN1_end:
            MessageBox(self, 'Check the range of the frames', 'error')
            return


        
        #framesSummed = np.rot90(frames[0].T, 2) * 0
        framesSummed = frames[0] * 0 # for trackpy 3.0



        for n in range(self.aveFrameN1_begin, self.aveFrameN1_end+1):
            #framesSummed += np.rot90(frames[n].T, 2)
            framesSummed += frames[n] # for trackpy 3.0
            
        framesSummed /= (self.aveFrameN1_end - self.aveFrameN1_begin + 1)
        
        self.framesSummed = framesSummed    


        self.fig_AveragingFrames = plt.figure('Averaging Frames')
        plt.clf()
        plt.title(str(os.path.basename(self.filepath)) + '\nAverage from frames:  ' + str(self.aveFrameN1_begin) +'  -  ' + str(self.aveFrameN1_end))
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(810,30,1000, 800) 


        #plt.imshow(framesSummed, origin='bottom', interpolation = 'None', cmap=plt.cm.jet)  
        
        plt.imshow(framesSummed, origin='upper', interpolation = 'None', cmap=plt.cm.jet)  # trackpy 0.3
        
        
        plt.colorbar()
        
        


      









   
        
        self.N_xPixel = frames.frame_shape[0]
        self.N_yPixel = frames.frame_shape[1]
        print 'self.N_xPixel: ', self.N_xPixel
        print 'self.N_yPixel: ', self.N_yPixel        
        
         
         
        
        #self.xposImageJFound_temp = np.array(self.xpos_temp) - 1 # the -1 is for correction between ImageJ and Python (x,y)
        #self.yposImageJFound_temp = self.N_xPixel - np.array(self.ypos_temp) - 1 # the -1 is for correction between ImageJ and Python (x,y)
        
        
     
        
        plt.ion()
        plt.show()

        
        
        self.btn_ShowBackgroundImage.Show()
        self.btn_ShowBackgroundIntensity.Hide()
                
        self.btn_ShowAnisotropyM.Hide()
        self.btn_ShowRawImages.Hide()
        self.btn_Save_data_figures.Hide()
        
        
        self.text_NFeatures.SetLabel(' ')
        
        



    def PixelHistogram(self, event):
        print '\n PixelHistogram'
        
        frames = self.framesLoaded
        
        self.frameN = int(self.PixelHistogram_FrameN.GetValue())
        
        self.frame_chosen = frames[self.frameN]
        print 'frame_chosen \n', self.frame_chosen
        
        
        try:
            plt.close(self.fig_PixelIntensityHistogram)
        except: pass
            
            
        self.fig_PixelIntensityHistogram = plt.figure('Pixel Intensity Histogram')
        plt.clf()
        plt.title(str(os.path.basename(self.filepath)))
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(1010,30,800, 600) 
        plt.hist(self.frame_chosen.ravel(), label='F#: ' + str(self.frameN))
        
        plt.ylabel('Occurence')
        plt.xlabel('Intensity')
        plt.legend(fontsize=16, loc='upper right',frameon=False)
        
        
        plt.ion()
        plt.show()


        
        self.btn_PixelHistogramSaveData.Show()
        

        
        

            

    def SaveData_PixelHistogram(self, event):
        print '\nCalculate Median Ave Save File'
        
        print ' '
        
        
        figFileName = self.filenameOnly[:40] + '_F#' + str(self.frameN)
        todayDate = time.strftime("%Y%m%d_%Hh%Mm")        



        ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_Pixel_Intesntiy.txt','w')
                
        for i in range(len(self.frame_chosen.ravel())):
            ff.write(str(self.frame_chosen.ravel()[i]) + '\n')
        ff.close()
        
        
        self.text_PixelHistogramSaveData.SetLabel('saved')
    
    




    def CalculateMedianAve(self, event):
        print '\nCalculate Median Ave'
        
        frames = self.framesLoaded
        self.AllPixelMedianTrajectory = []
        self.AllPixelAverageTrajectory = []
        
        for n in range(len(frames)):
            print 'frame n = ', n
            self.AllPixelMedianTrajectory.append(int(np.median(frames[n])))
            self.AllPixelAverageTrajectory.append(int(np.around(np.mean(frames[n]),0)))
            
            
        self.fig_MedianAverageIntensity = plt.figure('Median_Average_Intesntity_Trajectories')
        plt.clf()
        plt.title(str(os.path.basename(self.filepath)))
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(1010,30,800, 600) 
        plt.plot(self.AllPixelMedianTrajectory, label='Median')
        plt.plot(self.AllPixelAverageTrajectory, label='Average')
        plt.ylabel('Intesnsity')
        plt.xlabel('Frame #')
        plt.legend(fontsize=16, loc='upper right',frameon=False)
        
        
        plt.ion()
        plt.show()


        
        self.btn_CalculateMedianAveSaveData.Show()
        
              


            

    def SaveData_CalculateMedianAve(self, event):
        print '\nCalculate Median Ave Save File'
        
        print ' '
        
        
        prefix = ''
        figFileName = prefix + '_' + self.filenameOnly[:40]
        todayDate = time.strftime("%Y%m%d_%Hh%Mm")        



        ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_Median_Average.txt','w')
        ff.write('Frame# Median Average\n')
        
        for i in range(len(self.AllPixelMedianTrajectory)):
            ff.write(str(i+1) + ' ' + str(self.AllPixelMedianTrajectory[i]) + ' ' + str(self.AllPixelAverageTrajectory[i]) + '\n')
        ff.close()
        
        
        self.text_CalculateMedianAveSaveFile.SetLabel('saved')
    
    
    
    


    def showBackgroundImage(self, event):
        print '\n starting showBackgroundImage'



        framesSummed = self.framesSummed    
        self.NBGpixel = int(float(self.BG_pixel_size.GetValue()))
        
        if self.NBGpixel%2 == 0:
            MessageBox(self, 'The background pixel number must be an odd number', 'Error')
            return
        
        NBGpixel_oneSide = (self.NBGpixel - 1)/2


        ################################################
        #BG_x_ImageJ, BG_y_ImageJ = xpos, ypos
        
        #xc = y_ImageJ
        #yc = self.N_xPixel - x_ImageJ - 1
        ################################################
        
        
        
        BG_x = int(self.BG_x.GetValue())
        BG_y = int(self.BG_y.GetValue())
        
        #self.BG_x_ImageJ = BG_x
        #self.BG_y_ImageJ = self.N_xPixel - BG_y - 1
        
        self.BG_x_ImageJ = BG_x   # trackpy 0.3
        self.BG_y_ImageJ = BG_y   # trackpy 0.3
        
        
        
        print 'self.BG_x_ImageJ ', self.BG_x_ImageJ
        print 'self.BG_y_ImageJ ', self.BG_y_ImageJ
        
        
        
        

        
        self.BG_frameSummed = framesSummed[BG_y - NBGpixel_oneSide : BG_y + NBGpixel_oneSide + 1, BG_x - NBGpixel_oneSide : BG_x + NBGpixel_oneSide + 1]

        print 'self.BG_frameSummed \n', self.BG_frameSummed


        
        self.fig_GlobalBackground = plt.figure('Global Background Image and Intensity', figsize = (18, 9))
        mngr = plt.get_current_fig_manager()
        mngr.window.setGeometry(820,30,800, 900) 
        
          
        self.gs = gridspec.GridSpec(20, 20)
        self.gs.update(bottom = 0.05, top=0.95, left=0.05, right=0.95, wspace=5.0, hspace=5.0)
        plt.suptitle(str(os.path.basename(self.filepath)), size=12)
 
  
        ax1 = plt.subplot(self.gs[1:7,0:8])
        plt.tick_params(labelsize=12)
        #ax1.plot([1,2,3],[1,2,3])
        #plt.imshow(self.BG_frameSummed, origin='bottom', interpolation = 'None', cmap=plt.cm.jet)
        plt.imshow(self.BG_frameSummed, origin='upper', interpolation = 'None', cmap=plt.cm.jet)  # trackpy 0.3
        
        
        plt.colorbar()
        #ax1.plot(self.FIONA_centerIntensityAve)
        ax1.set_title(str(self.NBGpixel) + ' x ' + str(self.NBGpixel) + ', Frames: ' + str(self.aveFrameN1_begin) + ' - ' + str(self.aveFrameN1_end) + '\nIn Averaging Frames: (' + str(BG_x) + ', ' + str(BG_y) + ')' + '\nIn ImageJ: (' + str(self.BG_x_ImageJ) + ', ' + str(self.BG_y_ImageJ) + ')', size=12)
        #plt.xlabel('Frame #')
 
 
        plt.ion()
        plt.show()

         

        
        self.btn_ShowBackgroundIntensity.Show()
                
                        
        self.btn_ShowAnisotropyM.Hide()
        self.btn_ShowRawImages.Hide()
        self.btn_Save_data_figures.Hide()
          
        return
        
  
  
  



    def ShowBackgroundIntensity(self, event):
        print '\n ShowBackgroundIntensity starts '
        
        globalBackground_x_temp = np.array([int(self.BG_x.GetValue())])
        globalBackground_y_temp = np.array([int(self.BG_y.GetValue())])
        
        
        #globalBackground_x_temp_ImageJ = globalBackground_x_temp
        #globalBackground_y_temp_ImageJ = self.N_xPixel - globalBackground_y_temp - 1
        
        
        globalBackground_x_temp_ImageJ = globalBackground_x_temp # trackpy 0.3
        globalBackground_y_temp_ImageJ = globalBackground_y_temp # trackpy 0.3
        
        
        
        print 'globalBackground_x_temp_ImageJ ', globalBackground_x_temp_ImageJ
        print 'globalBackground_y_temp_ImageJ ', globalBackground_y_temp_ImageJ
        






        

        IntensityMethod = 1 # 1 is median
        centerPixelN = 1
        edgePixelN_temp = int(self.BG_pixel_size.GetValue())

        
        if edgePixelN_temp != self.NBGpixel:
            MessageBox(self, 'The background pixel number changed. Re-image the background, then try again', 'Error')
            return










        
        globalBackgroundCenterIntensityTrajectory, globalBackgroundEdgeIntensityTrajectory = IntensityTrajectoryPlotsFuncForM(self.framesLoaded, self.filepath, globalBackground_x_temp_ImageJ, globalBackground_y_temp_ImageJ, IntensityMethod, centerPixelN, edgePixelN_temp)


        self.GlobalBackgroundIntensityTrajectory = globalBackgroundEdgeIntensityTrajectory[0]



        print 'globalBackgroundCenterIntensityTrajectory ', globalBackgroundCenterIntensityTrajectory
        print 'globalBackgroundEdgeIntensityTrajectory ', globalBackgroundEdgeIntensityTrajectory



        '''
        plt.figure()
        plt.plot(globalBackgroundCenterIntensityTrajectory[0])
        plt.plot(globalBackgroundEdgeIntensityTrajectory[0])
        plt.ion()
        plt.show()
        '''




    
        IntensityTrajectoryGBG_fit_x = np.arange(len(self.GlobalBackgroundIntensityTrajectory))
        
        
        T_in_frame = float(self.fitInitialGuess_T.GetValue())



        A_guess = ( np.max(self.GlobalBackgroundIntensityTrajectory) + np.min(self.GlobalBackgroundIntensityTrajectory) )/2
        
        initial_guess = [A_guess, 0.2, T_in_frame, 0]
        
        try:
            popt, pcov = curve_fit(ModulationDepth_M_func, IntensityTrajectoryGBG_fit_x, self.GlobalBackgroundIntensityTrajectory, p0=initial_guess )
            self.GlobalBackgroundIntensityTrajectory_fit = ModulationDepth_M_func(IntensityTrajectoryGBG_fit_x, popt[0],popt[1],popt[2],popt[3])
            
            self.GlobalBackgroundIntensityTrajectory_fit_A = popt[0]
            self.GlobalBackgroundIntensityTrajectory_fit_M = popt[1]
            self.GlobalBackgroundIntensityTrajectory_fit_T = popt[2]
            self.GlobalBackgroundIntensityTrajectory_fit_phi = (popt[3]/(np.pi))*180
            
            self.GlobalBackgroundIntensityTrajectory_fit_A_err = np.sqrt(pcov[0,0])
            self.GlobalBackgroundIntensityTrajectory_fit_M_err = np.sqrt(pcov[1,1])
            self.GlobalBackgroundIntensityTrajectory_fit_T_err = np.sqrt(pcov[1,1])
            self.GlobalBackgroundIntensityTrajectory_fit_phi_err = np.sqrt(pcov[1,1])
            
            


            Fit_ave = int(self.GlobalBackgroundIntensityTrajectory_fit_A)
            Fit_err_amp = int(self.GlobalBackgroundIntensityTrajectory_fit_A * self.GlobalBackgroundIntensityTrajectory_fit_M)
            Fit_err_amp_per = (np.around(100.0*float(Fit_err_amp)/float(Fit_ave),1))
            Fit_T = int(self.GlobalBackgroundIntensityTrajectory_fit_T)
            
            DarkCount = float(self.DartCount.GetValue())
            
            Fit_ave_DarkCountCorrected = Fit_ave - DarkCount
            Fit_err_amp_per_DarkCountCorrected = (np.around(100.0*float(Fit_err_amp)/float(Fit_ave_DarkCountCorrected),1))
            

            
        except:
            
            print '\n\n fit failed.'
            
            self.GlobalBackgroundIntensityTrajectory_fit = ModulationDepth_M_func(IntensityTrajectoryGBG_fit_x, 0, 0, 10, 0)
            
            self.GlobalBackgroundIntensityTrajectory_fit_A = 0
            self.GlobalBackgroundIntensityTrajectory_fit_M = 0
            self.GlobalBackgroundIntensityTrajectory_fit_T = 0
            self.GlobalBackgroundIntensityTrajectory_fit_phi = 0
            
            self.GlobalBackgroundIntensityTrajectory_fit_A_err = 0
            self.GlobalBackgroundIntensityTrajectory_fit_M_err = 0
            self.GlobalBackgroundIntensityTrajectory_fit_T_err = 0
            self.GlobalBackgroundIntensityTrajectory_fit_phi_err = 0
            


            Fit_ave = int(self.GlobalBackgroundIntensityTrajectory_fit_A)
            Fit_err_amp = int(self.GlobalBackgroundIntensityTrajectory_fit_A * self.GlobalBackgroundIntensityTrajectory_fit_M)
            Fit_err_amp_per = np.nan
            Fit_T = int(self.GlobalBackgroundIntensityTrajectory_fit_T)
            
            DarkCount = float(self.DartCount.GetValue())
            
            Fit_ave_DarkCountCorrected = Fit_ave - DarkCount
            Fit_err_amp_per_DarkCountCorrected = np.nan
            
            
            pass
   
   








        ax2 = plt.subplot(self.gs[2:6,10:19])
        plt.tick_params(labelsize=12)
        plt.plot(self.GlobalBackgroundIntensityTrajectory)
        plt.plot(self.GlobalBackgroundIntensityTrajectory_fit, color='c')
        #ax2.set_title('Global Edge Background Median Intensity\nFit: ' + str(Fit_ave) + ' +- ' + str(Fit_err_amp) + ' (' + str(Fit_err_amp_per) + ' %)' + '   , T = ' + str(Fit_T) , size=12)
        ax2.set_title('Global Edge Background Median Intensity\nFit: ' + str(Fit_ave) + ' +- ' + str(Fit_err_amp) + ' (' + str(Fit_err_amp_per) + ' %)' + '   , T = ' + str(Fit_T) + '\nFit dark count corrected: ' + str(Fit_ave_DarkCountCorrected) + ' +- ' + str(Fit_err_amp) + ' (' + str(Fit_err_amp_per_DarkCountCorrected) + ' %)', size=12)
        plt.xlabel('Frame #')
 

        plt.ion()
        plt.show()

                                
        
        self.btn_ShowRawImages.Show()
        self.btn_ShowAnisotropyM.Hide()
        self.btn_Save_data_figures.Hide()







    def ShowRawImages(self, event):
        print '\nstarting ShowRawImages'
        
        
        

        try:
            for n in range(len(self.fig_AllDataCenterImagesOnly)):
                plt.close(self.fig_AllDataCenterImagesOnly[n])
        except: pass
            
        
        
        
        
        self.xposImageJFound = []
        self.yposImageJFound = []
                
        
        
            


        moreXtemp = []
        moreYtemp = []
        for n in range(len(self.MorePairs_x)):
            if self.MorePairs_x[n].GetValue() != '-' and self.MorePairs_y[n].GetValue() != '-':
                
                #moreXtemp.append(int(self.MorePairs_x[n].GetValue()) )  # 
                #moreYtemp.append(self.N_xPixel - int(self.MorePairs_y[n].GetValue()) -1 ) # -1 is for imageJ correction

                moreXtemp.append(int(self.MorePairs_x[n].GetValue()) )  # trackpy 0.3
                moreYtemp.append(int(self.MorePairs_y[n].GetValue()) )  # trackpy 0.3
            
        
        
        if len(moreXtemp) == 0:
            MessageBox(self, 'No features, check the feature (x,y)', 'Error')
            return



        
        
        
        
        xposTemp = np.append(np.around(self.xposImageJFound, 0) , moreXtemp)
        yposTemp = np.append(np.around(self.yposImageJFound, 0) , moreYtemp)
        



        xpos = xposTemp.astype(int)
        ypos = yposTemp.astype(int)

        self.xposAllTemp = xpos
        self.yposAllTemp = ypos


        
        print 'len(xpos) ', len(xpos)
        
        
        print 'xpos ', xpos
        print 'ypos ', ypos
        
        
        
        ################################################
        x_ImageJ, y_ImageJ = xpos, ypos
        
        #xc = y_ImageJ
        #yc = self.N_xPixel - x_ImageJ - 1
        
        xc = x_ImageJ           #trackpy 0.3
        yc = y_ImageJ        #trackpy 0.3
        
        
        ################################################
        
        
        
        frames = self.framesLoaded
        
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
        
        if frameNTemp > TotalNframes-1:
            MessageBox(self, 'Error: the frame # does not exist', 'Error')
            return
        
        
        
        for k in range(len(xpos)):            
            
            frameLargeTemp = frames[frameNTemp][yc[k]-NCenterSubPixel:yc[k]+NCenterSubPixel+1, xc[k]-NCenterSubPixel:xc[k]+NCenterSubPixel+1]
            #frameCenterTemp = frames[frameTemp][yc[k]-2:yc[k]+3, xc[k]-2:xc[k]+3]
            
            centerImage[k].append(frameLargeTemp) # it's not matching to the actual frame #, only first a few are saved



        self.centerImage = centerImage



        
        NFeatures = len(xpos)
  
       
        NCfig = int(  math.ceil((NFeatures/84.0)))
        self.fig_AllDataCenterImagesOnly = [[] for _ in xrange(NCfig)]
        
        
        
            



        x3 = np.linspace(0, 15, 16)
        y3 = np.linspace(0, 15, 16)
        x3, y3 = np.meshgrid(x3, y3)
        

        k = 0
        fn = 0
        for n in range(NCfig):
            #print 'n = ', n
            self.fig_AllDataCenterImagesOnly[n] = plt.figure('Raw_Data_Images'+ str(n), figsize = (18, 9))
            mngr = plt.get_current_fig_manager()
            mngr.window.setGeometry(820,230,1600, 800) 
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
                    #plt.imshow(np.rot90(centerImage[m+fn][0].reshape(15, 15),3), interpolation = 'None', cmap=plt.cm.jet, origin='upper', extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                    plt.imshow(centerImage[m+fn][0].reshape(15, 15), interpolation = 'None', cmap=plt.cm.jet, origin='upper',
                               extent=(x3.min(), x3.max(), y3.min(), y3.max()))   # trackpy 0.3
                except:
                    print '\nNo center image for ', 'M# = ', str(m + fn)


                
                k += 1
                if k%84 ==0:
                    fn += 84
                if k == NFeatures:
                    break




        self.btn_ShowAnisotropyM.Show()
        self.btn_Save_data_figures.Hide()




    

    def ShowAnisotropyM(self, event):
                
        print '\n Starting ShowAnisotropyM'
        


        
        firstFrameN = int(self.rawImageFN1.GetValue())
        secondFrameN = int(self.rawImageFN2.GetValue())

        
        if firstFrameN > self.TotalNframes-1 or secondFrameN > self.TotalNframes-1:
            MessageBox(self, 'Check the image frame numbers', 'Error')
            return
        
            
            


        if int(self.centerPixelN.GetValue())%2 == 0:
            MessageBox(self, 'Center pixel number must be an odd number', 'Error')
            return


        if int(self.edgePixelN.GetValue())%2 == 0:
            MessageBox(self, 'Edge pixel number must be an odd number', 'Error')
            return
        
        
        

        
        try: plt.close(self.fig_ExcitationM_Histogram)
        except: pass
    
        try:        
            for n in range(len(self.fig_AllDataCenterImages)):    
                plt.close(self.fig_AllDataCenterImages[n])
        except: pass




        '''
        self.xposImageJFound = []
        self.yposImageJFound = []
        




        moreXtemp = []
        moreYtemp = []
        for n in range(len(self.MorePairs_x)):
            if self.MorePairs_x[n].GetValue() != '-' and self.MorePairs_y[n].GetValue() != '-':
                moreXtemp.append(int(self.MorePairs_x[n].GetValue()) ) #
                moreYtemp.append(self.N_xPixel - int(self.MorePairs_y[n].GetValue()) -1) # -1 is for ImageJ correction.
             

  
        xposTemp = np.append(np.around(self.xposImageJFound, 0) , moreXtemp)
        yposTemp = np.append(np.around(self.yposImageJFound, 0) , moreYtemp)
        
        
        Allxpos = xposTemp.astype(int)
        Allypos = yposTemp.astype(int)

        '''




        Allxpos = self.xposAllTemp
        Allypos = self.yposAllTemp

        
        print 'Allxpos ', Allxpos
        print 'Allypos ', Allypos
        
        

        
        
    
        

        IntensityMethod = self.radio_IntensityMethod.GetSelection()
        
        centerPixelN = int(self.centerPixelN.GetValue())
        edgePixelN = int(self.edgePixelN.GetValue())
        
        featurePixelConditions = 'Center: ' + str(centerPixelN) + ' x ' + str(centerPixelN) + ',   Edge: ' + str(edgePixelN) + ' x ' + str(edgePixelN) + ',   Global BG: ' + str(self.NBGpixel) + ' x ' + str(self.NBGpixel) + ' at (' + str(self.BG_x_ImageJ) + ', ' + str(self.BG_y_ImageJ) + ')'
        
        centerIntensityTrajectory, edgeIntensityTrajectory = IntensityTrajectoryPlotsFuncForM(self.framesLoaded, self.filepath, Allxpos, Allypos, IntensityMethod, centerPixelN, edgePixelN)



        #print 'centerIntensityMaxTrajectory ', centerIntensityMaxTrajectory
        

 
        
        
        IntensityTrajectory = np.array(centerIntensityTrajectory)
        self.backgroundIntensityTrajectory = np.array(edgeIntensityTrajectory)
        
        
        







        
        

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
                IntensityTrajectoryBS = IntensityTrajectory - self.GlobalBackgroundIntensityTrajectory
                BG_ref_correction_text += 'Global Background Frame by Frame, '
                self.BGS_type = 'Global_BG_FrameByFrame'
                
            elif self.radio_BSMethod.GetSelection() == 1:
                IntensityTrajectoryBS = IntensityTrajectory - np.array(self.backgroundIntensityTrajectory)
                BG_ref_correction_text += 'Individual Background Frame by Frame, '
                self.BGS_type = 'Global_BG_FrameByFrame'
                
        else:
            IntensityTrajectoryBS = IntensityTrajectory
            #BG_ref_correction_text = 'Background Substraction 9x9 edge: Total Ave'
            BG_ref_correction_text = 'Background Subtraction: No,'
            self.BGS_type = 'noBGS'
            
        
        
        
        if IntensityMethod == 0:
            BG_ref_correction_text += '\nIntensity Method: Ave '
        elif IntensityMethod == 1:
            BG_ref_correction_text += '\nIntensity Method: Med '
        elif IntensityMethod == 2:
            BG_ref_correction_text += '\nIntensity Method: Max 1 '
            



            
        
        
        Pol_Cor_factor = float(self.Excitation_Pol_correction_factor.GetValue())
        print '\nPol_Cor_factor = ', Pol_Cor_factor
        
        
        if self.cb_Excitation_Pol_correction_yn.GetValue():
            Excitation_Pol_correction_text = 'Excitation Polarization Correction: YES by ' + str(Pol_Cor_factor)
        else:
            Excitation_Pol_correction_text = 'Excitation Polarization Correction: NO'
            



        
        
        self.AnalysisConditions = ' '
        
        
        
        
        xposChosen = []
        yposChosen = []
        
        self.backgroundIntensityTrajectoryChosen = []
        self.backgroundIntensityTrajectory_FitChosen = []
        
        self.IntensityTrajectoryChosen = []
        self.IntensityTrajectoryChosenBS = []
        self.IntensityTrajectoryChosenBS_fit = []
        
        self.IntensityTrajectoryChosenBS_fit_A = []
        self.IntensityTrajectoryChosenBS_fit_M = []
        self.IntensityTrajectoryChosenBS_fit_T = []
        self.IntensityTrajectoryChosenBS_fit_phi = []

        self.IntensityTrajectoryChosenBS_fit_A_err = []
        self.IntensityTrajectoryChosenBS_fit_M_err = []
        self.IntensityTrajectoryChosenBS_fit_T_err = []
        self.IntensityTrajectoryChosenBS_fit_phi_err = []



        
        
        self.IntensityTrajectoryChosenBS_fit_RsquareFitGoodness = []
        
        
        
        IntensityTrajectoryChosenBS_fit_x = np.arange(len(IntensityTrajectory[0]))
        
        
        T_in_frame = float(self.fitInitialGuess_T.GetValue())



        print 'IntensityTrajectoryBS ', IntensityTrajectoryBS
        
        for n in range(len(IntensityTrajectoryBS)):
           
            A_guess = ( np.max(IntensityTrajectoryBS[n]) + np.min(IntensityTrajectoryBS[n]) )/2
            
            initial_guess = [A_guess, 0.5, T_in_frame, 0]
            
            try:
                popt, pcov = curve_fit(ModulationDepth_M_func, IntensityTrajectoryChosenBS_fit_x, IntensityTrajectoryBS[n], p0=initial_guess )
                fit_temp = ModulationDepth_M_func(IntensityTrajectoryChosenBS_fit_x, popt[0],popt[1],popt[2],popt[3])
                
                A_temp = popt[0]
                #if A_temp < 0:
                #    A_temp *= -1
                    
                M_temp = popt[1]
                
                T_temp = popt[2]
                
                phi_temp = (popt[3]/(np.pi))*180
                
                A_temp_err = np.sqrt(pcov[0,0])
                M_temp_err = np.sqrt(pcov[1,1])
                T_temp_err = np.sqrt(pcov[1,1])
                phi_temp_err = np.sqrt(pcov[1,1])
                
                
                '''
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
                    
                '''
                
                print 'Fit: A, M, T, phi = ', A_temp, M_temp, T_temp, phi_temp
                print 'Fit Errors: A, M, T, Phi = ', A_temp_err, M_temp_err, T_temp_err, phi_temp_err
                                
                RsquareFitGoodness = R_SqureFitGoodness(IntensityTrajectoryBS[n], fit_temp)
                
                
            except:
                print 'fit failed for n: ', n
            
                fit_temp = IntensityTrajectoryChosenBS_fit_x*0 + A_guess
                
                A_temp = -999
                M_temp = -999
                T_temp = -999
                phi_temp = -999
                
                
                A_temp_err = -999
                M_temp_err = -999
                T_temp_err = -999
                phi_temp_err = -999
                
                
                
                RsquareFitGoodness = -999
               
                
                
            if self.cb_Excitation_Pol_correction_yn.GetValue():
                Excitation_Pol_correction_text = 'Excitation Polarization Correction: YES by ' + str(Pol_Cor_factor)
            else:
                Excitation_Pol_correction_text = 'Excitation Polarization Correction: NO'
            
     
     
     
     
            xposChosen.append(Allxpos[n])
            yposChosen.append(Allypos[n])
            self.IntensityTrajectoryChosen.append(IntensityTrajectory[n])
            self.IntensityTrajectoryChosenBS.append(IntensityTrajectoryBS[n])
            self.backgroundIntensityTrajectoryChosen.append(self.backgroundIntensityTrajectory[n])
            self.backgroundIntensityTrajectory_FitChosen.append(self.backgroundIntensityTrajectory_Fit[n])

            
            self.IntensityTrajectoryChosenBS_fit.append(fit_temp)
            
            self.IntensityTrajectoryChosenBS_fit_A.append(A_temp)
            self.IntensityTrajectoryChosenBS_fit_M.append(M_temp)
            self.IntensityTrajectoryChosenBS_fit_T.append(T_temp)
            self.IntensityTrajectoryChosenBS_fit_phi.append(phi_temp)
            
                      
            self.IntensityTrajectoryChosenBS_fit_A_err.append(A_temp_err)
            self.IntensityTrajectoryChosenBS_fit_M_err.append(M_temp_err)
            self.IntensityTrajectoryChosenBS_fit_T_err.append(T_temp_err)
            self.IntensityTrajectoryChosenBS_fit_phi_err.append(phi_temp_err)
            
                      
                      
            self.IntensityTrajectoryChosenBS_fit_RsquareFitGoodness.append(RsquareFitGoodness)

          

        Corrections = BG_ref_correction_text + ',       ' + Excitation_Pol_correction_text    
        Corrections2 = BG_ref_correction_text + ',  ' + Excitation_Pol_correction_text + '\nF#s: ' + str(self.rawImageFN1.GetValue())+ ', ' + str(self.rawImageFN2.GetValue()) + ',   ' + featurePixelConditions
        self.Corrections = Corrections        
        self.Corrections2 = Corrections2







        
        NFeatures = len(self.IntensityTrajectoryChosen)

        TotalFeatureN = 'Total # ' + str(NFeatures)






        print '\n For RawImages'
        
        xpos = np.array(xposChosen)
        ypos = np.array(yposChosen)       
        
        self.xposChosen = xpos
        self.yposChosen = ypos
      
        ################################################
        
        x_ImageJ, y_ImageJ = xpos, ypos
            
        
        #xc = y_ImageJ
        #yc = self.N_xPixel - x_ImageJ -1
        
                
        xc = x_ImageJ           #trackpy 0.3
        yc = y_ImageJ           #trackpy 0.3
        
        
        ################################################
        
      
        
        frames = self.framesLoaded
        
        TotalNframes = len(frames)
        print 'TotalNframes: ', TotalNframes

        
        NCenterSubPixel = 7  # number of pixels from the center
        
        
        
        #frameTemp0 = frame0[yc-NCenterSubPixel:yc+NCenterSubPixel + 1, xc-NCenterSubPixel:xc+NCenterSubPixel+1]
        x3 = np.linspace(0, 2 * NCenterSubPixel, 2* NCenterSubPixel + 1)
        y3 = np.linspace(0, 2 * NCenterSubPixel, 2* NCenterSubPixel + 1)
        x3, y3 = np.meshgrid(x3, y3)
        


        self.centerImage = [[] for _ in xrange(len(xc))] # centerImage[Molecule#][Frame#]
        

        for frameTemp in [firstFrameN,secondFrameN]: # from the 0th frame to the end frame.
            print 'Not rotated center image data for frame # ', frameTemp
            for k in range(len(xpos)):            
                frameLargeTemp = frames[frameTemp][yc[k]-NCenterSubPixel:yc[k]+NCenterSubPixel+1, xc[k]-NCenterSubPixel:xc[k]+NCenterSubPixel+1]
                #frameCenterTemp = frames[frameTemp][yc[k]-2:yc[k]+3, xc[k]-2:xc[k]+3]
                self.centerImage[k].append(frameLargeTemp) # it's not matching to the actual frame #, only first a few are saved
        







        
        NFeatures = len(xpos)
  
        NFeaturePerFig = 30.0
       
        NCfig = int(  math.ceil((NFeatures/float(NFeaturePerFig))))
        self.fig_AllDataCenterImages = [[] for _ in xrange(NCfig)]

        x3 = np.linspace(0, 15, 16)
        y3 = np.linspace(0, 15, 16)
        x3, y3 = np.meshgrid(x3, y3)
        

        k = 0
        fn = 0
        for n in range(NCfig):
            #print 'n = ', n
            self.fig_AllDataCenterImages[n] = plt.figure('M_Data_Images'+ str(n), figsize = (18, 9))
            self.fig_AllDataCenterImages[n].subplots_adjust(top = 0.9, bottom = 0.05, left = 0.02, right = 0.99, wspace = 0.36, hspace = 0.4)
            self.fig_AllDataCenterImages[n].text(0.02, 0.925, TotalFeatureN +'\n' +Corrections2, ha="left", va="bottom", size="small",color="black", weight='normal')
            
            plt.figtext(0.5, 0.94, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly+ '\n' + self.AnalysisConditions ), ha='center', color='black', weight='normal', size='small')
                


            for m in range(int(NFeaturePerFig)):

  
                titlecolorTemp = 'black'
                weight='normal'
                

                
                imageTemp = np.append(np.ravel(self.centerImage[m+fn][0]), np.ravel(self.centerImage[m+fn][1]))
                
                vminTemp = int(np.mean( np.sort(imageTemp, axis=None)[:5]) )
                vmaxTemp = int(np.mean( np.sort(imageTemp, axis=None)[-5:]) )
                
                
                plt.subplot(6,15,3*m+1)
                plt.title('# ' + str(m+fn) + ', (' + str(int(x_ImageJ[m+fn])) + ',' + str(int(y_ImageJ[m+fn]))+')'  , fontsize=9, color = titlecolorTemp, weight=weight)
                plt.tick_params(labelsize=7)
                #plt.imshow(np.rot90(self.centerImage[m+fn][0].reshape(15, 15),3), interpolation = 'None',vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='upper',
                #           extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                plt.imshow(self.centerImage[m+fn][0].reshape(15, 15), interpolation = 'None',vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='upper',
                           extent=(x3.min(), x3.max(), y3.min(), y3.max()))  # trackpy 0.3
                
                
                plt.subplot(6,15,3*m+2)
                plt.title('R2=' + str(np.around(self.IntensityTrajectoryChosenBS_fit_RsquareFitGoodness[m + fn],2)) + ',T=' + str(int(np.around(self.IntensityTrajectoryChosenBS_fit_T[m + fn],0)) )  , fontsize=9, color = titlecolorTemp, weight=weight)
                plt.tick_params(labelsize=7)
                #plt.imshow(np.rot90(self.centerImage[m+fn][1].reshape(15, 15),3), interpolation = 'None',vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='upper',
                #           extent=(x3.min(), x3.max(), y3.min(), y3.max()))
                plt.imshow(self.centerImage[m+fn][1].reshape(15, 15), interpolation = 'None',vmin = vminTemp, vmax=vmaxTemp, cmap=plt.cm.jet, origin='upper',
                           extent=(x3.min(), x3.max(), y3.min(), y3.max()))  # trackpy 0.3
                      
                           
                plt.subplot(6,15, 3*m+3)
                plt.tick_params(labelsize=7)
                maxFrameN = np.argmax(self.IntensityTrajectoryChosenBS_fit[m + fn])
                #print '\n:maxFrameN = ', maxFrameN 
                if maxFrameN >= 90:
                    maxFrameN = maxFrameN - 90
                #if maxFrameN >= 45:
                #    maxFrameN = maxFrameN - 45                    
                #print 'maxFrameN = ', maxFrameN 
                
                #plt.title( str(maxFrameN*2) + '$^\circ$'  + ', M=' + str(np.around(IntensityTrajectoryChosenBS_fit_M[m + fn],2) )  , fontsize=9, color = titlecolorTemp, weight=weight)
                plt.title( str(int(float(np.around(self.IntensityTrajectoryChosenBS_fit_phi[m + fn],0)))) + '$^\circ$'  + ', M=' + str(np.around(self.IntensityTrajectoryChosenBS_fit_M[m + fn],2) )  , fontsize=9, color = titlecolorTemp, weight=weight)

                
                plt.plot(self.IntensityTrajectoryChosenBS[m + fn], color = 'OrangeRed')
                plt.plot(self.backgroundIntensityTrajectoryChosen[m + fn], color = 'Gray')
                plt.plot(self.GlobalBackgroundIntensityTrajectory, color = 'Cyan')
                plt.plot(self.IntensityTrajectoryChosen[m + fn], color = 'Black')
                plt.xlim(0, len(self.GlobalBackgroundIntensityTrajectory)-1)

                
                
                                
                
                
                
                if self.cb_includeFit_yn.GetValue() == True:
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
        self.fig_ExcitationM_Histogram.text(0.02, 0.97, TotalFeatureN +',       ' + Corrections, ha="left", va="bottom", size=9,color="black", weight='normal')
        plt.figtext(0.5, 0.92, str(self.filedirectoryOnly) + '\n' + str(self.filenameOnly+ '\n' + self.AnalysisConditions ), ha='center', color='black', weight='normal', size=9)

        #plt.title('Total # ' + str(len(IntensityTrajectoryChosenBS_fit_M)))
        plt.subplot(2,2,1)
        plt.title('M', size=14)
        if len(self.IntensityTrajectoryChosenBS_fit_M) == 0:
            self.IntensityTrajectoryChosenBS_fit_M = [999,9999]            
        plt.hist(self.IntensityTrajectoryChosenBS_fit_M, bins=np.arange(0.0, 1.3, 0.1))
        plt.xlabel('Modulation Depth (M)')
        
        plt.xlim(0.0, 1.31)
        
        



        plt.subplot(2,2,3)
        plt.title('Phi')
        if len(self.IntensityTrajectoryChosenBS_fit_phi) == 0 or len(self.IntensityTrajectoryChosenBS_fit_phi) == 1:
            self.IntensityTrajectoryChosenBS_fit_phi = [11, 111]
        plt.hist(self.IntensityTrajectoryChosenBS_fit_phi)
        plt.xlabel('Fit phase angle phi')



       
        plt.show()

        self.btn_Save_data_figures.Show()
        
        self.text_NFeatures.SetLabel('# of features: ' + str(len(self.IntensityTrajectoryChosenBS)))
        
  
        







        
                
      
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
                    
                self.fig_AveragingFrames.savefig(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_2_AveragingFrames.png')
                plt.close(self.fig_AveragingFrames)
                print 'fig_AveragingFrames figure is saved'
            
            except:
                pass
            
    
            try:
                    
                self.fig_GlobalBackground.savefig(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_3_GlobalBackground.png')
                plt.close(self.fig_GlobalBackground)
                print 'fig_GlobalBackground figure is saved'
            
            except:
                save_data_error += 1
                pass
            
    
              
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
    
            
    


        
        if self.cb_SaveImageJ_xyData_yn.GetValue():
            print ''
            

            ff = open(self.filedirectoryOnly + '\\' + figFileName +'_' + todayDate +  '_6_xy_A_M_T(frames)_phi(deg).txt','w')
            
            for i in range(len(self.xpos_ImageJ)):
                
                ff.write(str(self.xposChosen[i]) + ' ' + str(self.yposChosen[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_A[i])
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_M[i]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_T[i]) 
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_phi[i]) + '\n'   )    
                
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
            
            ff.write('@agg_v1 , Aggregate_M_analysis , Background_Subtraction_Type: ' + self.BGS_type + '\n')
            
            
            ff.write('\nGlobal Background Data: ImageJ_x ImageJ_y pixel_by_pixel_size / BG_Intensity / BG_Image\n')
            ff.write(str(self.BG_x_ImageJ) + ' ' + str(self.BG_y_ImageJ) + ' , ' + str(self.NBGpixel) + '\n')
            for i in range(len(self.GlobalBackgroundIntensityTrajectory)):
                ff.write(str(int(self.GlobalBackgroundIntensityTrajectory[i])) + ' ')
            ff.write('\n')
            for i in range(len(self.BG_frameSummed)):
                for k in range(len(self.BG_frameSummed[i])):
                    ff.write(str(int(self.BG_frameSummed[i][k])) + ' ')
                ff.write(' , ')
            ff.write('\n\n')
            
            
            ff.write('x y A M T phi  / rawIntensity / intensityBS / fit / BG / BGfit / image0 / image1\n')



            for n in range(len(self.IntensityTrajectoryChosenBS)):
                

                ff.write('#' + str(n) + '\n')
                
                    

                    
                ff.write(str(self.xposChosen[n]) + ' ' + str(self.yposChosen[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_A[n])
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_M[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_T[n]) 
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_phi[n]) 
                + ' Err: ' + str(self.IntensityTrajectoryChosenBS_fit_A_err[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_M_err[n]) 
                + ' ' + str(self.IntensityTrajectoryChosenBS_fit_T_err[n]) + ' ' + str(self.IntensityTrajectoryChosenBS_fit_phi_err[n])
                + '\n')    
                



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
     
     
     
     












#----------------------------------------------------------------------
# Run the program
   
if __name__ == "__main__":
    app = wx.App(False)

    frame = Aggregate_M_analysis()
    frame.Show()
        
    app.MainLoop()
    
    
    
        
        

#MessageBox(self, 'Opening Data Analysis Window', 'Message')



       
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
