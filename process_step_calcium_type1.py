# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 19:44:20 2017

@author: jcoleman
"""
import numpy as np
from collections import OrderedDict
from matplotlib import pyplot as plt

import cPickle as pickle
import readroi as roizip

# Coleman lab modules - combine into a single class?
from calcium import OpenStepCalcium as OSC
from calcium import calcium_imaging_data_fast as cidfast
from calcium import open_calcium_data_fast as ocdfast
from calcium import StepCodeFile as SCF

#==============================================================================
#Get step code file (*.bin), TIF image, ROI.zip files
#==============================================================================
#fileDirectoryTS = 'F:/coleman lab/jasonc/thygcamp6s_test2/'
#fileDirectory = 'F:/coleman lab/jasonc/thygcamp6s_test2/'
fileDirectory = '/Users/jcoleman/Documents/--DATA--/in vivo gcamp analysis/'+ \
    'thygcamp6s_LT4(9-10-17)/'
#fileDirectory = '/Users/jcoleman/Documents/--DATA--/' + \
#    'in vivo gcamp analysis/thygcamp6s_D4 5Hz (9-30-17)/'

datafile = 'D2_001_Z1t0_hz05'

if datafile == 'D2_001_Z1t0_hz05':

    fileBin = 'mThy6s2_alldrift_D2_001_12_data.bin'
    imgname = 'STD_mThy6s2_alldrift_D2_001.tif'
    roizipname = 'mThy6s2_alldrift_D2_001_ROI.zip'
    

if datafile == 'D3_001_Z1t1_hz05':
    
    fileBin = 'mThy6s_001_D3_6_data.bin'
    imgname = 'STD_mThy6s2_alldrift_D3_001Z1.tif'
    roizipname = 'mThy6s2_alldrift_D3_001Z1_ROI.zip'
    

if datafile == 'D3_001_Z1t2_hz05':
    
    fileBin = 'mThy6s_001_D3Z1t2_13_data.bin'
    imgname = 'STD_mThy6s2_alldrift_D3_001Z1t2.tif'
    roizipname = 'mThy6s2_alldrift_D3_001Z1t2_ROI.zip'
    

if datafile == 'D4_001_Z1_hz1':

    fileBin = 'mThy6s_001_D4Z1_2_data.bin'
    imgname = 'STD_mThy6s2_alldrift_D4_001Z1.tif'
    roizipname = 'mThy6s2_alldrift_D4_001Z1_ROI.zip'
    

if datafile == 'D4_001_Z1_hz5':

    fileBin = 'mThy6s_001_D4Z1hz5_5_data.bin'
    imgname = 'STD_mThy6s2_alldrift_D4_001Z1hz5.tif'
    roizipname = 'mThy6s2_alldrift_D4_001Z1hz5_ROI.zip'
    
    
if datafile == 'D4_002_Z1_hz1':

    fileBin = 'mThy6s_002_D4Z1_3_data.bin'
    imgname = 'STD_mThy6s2_alldrift_D4_002Z1.tif'
    roizipname = 'mThy6s2_alldrift_D4_002Z_ROI.zip'
    
    
if datafile == 'D4_002_Z1_hz5':

    fileBin = 'mThy6s_002_D4Z1hz5_4_data.bin'
    imgname = 'STD_mThy6s2_alldrift_D4_002Z1hz5.tif'
    roizipname = 'mThy6s2_alldrift_D4_002Zhz5_ROI.zip'
    

# parameters for mulitplot function
a=roizip.read_roi_zip(fileDirectory + roizipname)

rows = len(a) # number of ROIs
cols = 8 # number of distinct stims

#==============================================================================

#==============================================================================
#PARAMETERS
#==============================================================================
daqRecordingFreq = 1000.0 #sampling rate - Hz

stimList = [0,45,90,135,180,225,270,315]

#stimList = [0,45,90,135,180] #designate flip/flop?

csvfiletype = 3 # 1= ; 2= ; 3=FIJI/ImageJ csv intensity file

"""
Set parameters for run_deltaf_ewma function below.
Recommended parameters for 30Hz scan/frame rate: 
    t_ 0= 0.2 #0.2 for 30fps; 0.8 for 7.6fps
    t_1 = 0.75
    t_2 = 3
"""
t_0 = 0.2 
t_1 = 0.75
t_2 = 3
mpmSamplingFreq = 30 # ~Hz/fps for 2p scan rate

# gray time (in between stims)
gray_offset= 7.0 #seconds

# set the window for averaging pre-stim gray (207 frames = 7s)
pre_timewindow = 207-30 #207-75 # e.g., 207-30 = 1s from 30fps sample rate

# delay time prior to first stim onset (in frames)
# used for spontaneous activity analysis
delayframes = 1700
#==============================================================================


# Begin decoding and data extraction
code_file = SCF.StepCodeFile(fileDirectory+fileBin,8,8,1000)

stims = code_file.get_stim_angles(stimList)

Bframes = [ts[0] for ts in code_file.timestamps[2]]


# Get intensity file
cellsraw, cellsmean, areas, xycoords, Data = ocdfast.openData(csvfiletype)

"""
Calculate deltaF/F for spontaneous activity
Run Konnerth lab deltaF/F and EWMA filter
From Konnerth lab Nature Protocols paper, for 30Hz scanning:
    t_0 = 0.2
    t_1 = 0.75
    t_2 = 3
    samplingfreq = 30
See paper for more details
"""   
dff, dff_ewma, dff_offset = (
    cidfast.run_deltaf_ewma(
        cellsmean,
        t_0,
        t_1,
        t_2,
        mpmSamplingFreq
        )
)

"""
Calculate deltaF/F for event-related activity
1) Filter raw data
2) Use filtered data for deltaF/F below (e.g., use 1s pre-stim gray as f0)

lo-pass filter function for raw data:
    butter_lowpass_filter(data,cutoff_freq,sampling_freq,order)
"""
cells = list()
for cell in range(len(cellsmean)):
    if mpmSamplingFreq == 30:
        lopass_temp = ( 
            OSC.butter_lowpass_filter
                (
                 cellsmean[cell],
                 1.0,
                 mpmSamplingFreq,
                 1
                 )
        )
        cells.append(lopass_temp)
        
    elif mpmSamplingFreq != 30:
        print('Enter correct sampling frequency (mpmSamplingFreq) for '+ \
                                                'butter_lowpass_filter')


handler = cidfast.FileHandler()

#makes dictionary of [cellnumber][orientation][block]
gray_offset *= daqRecordingFreq

response_rawdata = OrderedDict()

grays = OrderedDict()

response_indices = OrderedDict()

grays_indices = OrderedDict()

for cell in range(len(cells)):

    response_rawdata[cell] = OrderedDict()

    grays[cell] = OrderedDict()

    response_indices[cell] = OrderedDict()

    grays_indices[cell] = OrderedDict()

    for stim in stims:

        response_rawdata[cell][stim] = list()

        grays[cell][stim] = list()

        response_indices[cell][stim] = list()

        grays_indices[cell][stim] = list()

        for ts in stims[stim]:

            begin = float(ts[0])

            end = float(ts[1])

            begin_frame_time = handler.get_nearest(begin, Bframes)

            begin_gray_time = handler.get_nearest(begin - gray_offset, Bframes)

            

            end_frame_time = handler.get_nearest(end, Bframes)

            end_gray_time = handler.get_nearest(begin, Bframes)

            

            begin_frame = list(Bframes).index(begin_frame_time)

            begin_gray = list(Bframes).index(begin_gray_time)

            

            end_frame = list(Bframes).index(end_frame_time)

            end_gray = list(Bframes).index(end_gray_time)

            

            chunk = cells[cell][int(begin_frame):int(end_frame)]

            gray_chunk = cells[cell][int(begin_gray):int(end_gray)]

            

            response_rawdata[cell][stim].append(chunk)

            grays[cell][stim].append(gray_chunk)

            

            (response_indices[cell][stim].append(
                                                [int(begin_frame),
                                                 int(end_frame)])
                                                 )

            (grays_indices[cell][stim].append(
                                                [int(begin_gray),
                                                 int(end_gray)])
                                                 )



#example: plots all 5 block of degree 45 orientation for cell 0

#cell_0_45 = response_rawdata[0][45]
#
#for block in cell_0_45:
#
#    plt.plot(block)


"""
Calculate deltaF/F by block using mean of 1s pre-gray for each stim for f0 
(where deltaF = ft-f0, F = f0; ft = fluorescence at time t)
"""
# First, get 1s pre-stim gray data and store in dictionary: 
#   pre1sgray[cell][orientation][block]
pre1sgray_data = OrderedDict()

for cell in grays:
    
    pre1sgray_data[cell] = OrderedDict()
    
    for orientation in grays[cell]:
        
        pre1sgray_data[cell][orientation] = OrderedDict()
        
        for block in range(len(grays[cell][orientation])):
             
             pre1sgray_data[cell][orientation][block] = (
                 grays[cell][orientation][block][pre_timewindow:end])
             
# Second, store deltaF/F response data:
#   response_data[cell][orientation][block]
response_data = OrderedDict()

for cell in pre1sgray_data:
    
    response_data[cell] = OrderedDict()
    
    for orientation in pre1sgray_data[cell]:
            
        response_data[cell][orientation] = OrderedDict()
            
        for block in range(len(pre1sgray_data[cell][orientation])):
            
            tempmean = np.mean(pre1sgray_data[cell][orientation][block])
            
            response_data[cell][orientation][block] = (
                (response_rawdata[cell][orientation][block] - tempmean) / 
                    tempmean)
            
# Third, store deltaF/F response data for gray 1s:
#   pre1sgray_data[cell][orientation][block]
pre1sgray_dFoF_data = OrderedDict()

for cell in pre1sgray_data:
    
    pre1sgray_dFoF_data[cell] = OrderedDict()
    
    for orientation in pre1sgray_data[cell]:
            
        pre1sgray_dFoF_data[cell][orientation] = OrderedDict()
            
        for block in range(len(pre1sgray_data[cell][orientation])):
            
            tempmean = np.mean(pre1sgray_data[cell][orientation][block])
            
            pre1sgray_dFoF_data[cell][orientation][block] = (
                (pre1sgray_data[cell][orientation][block] - tempmean) / 
                    tempmean)

#"""
#How to access responses to a given orientation for a given cell:
#   response_data[cell][ori]
#"""
#for ori in response_data[26]:
#
#    cell_0_x = response_data[26][ori]
#    plt.subplots()
#    
#    for block in cell_0_x:
#    
#        plt.plot(cell_0_x[block])
#    tempvals = list()
#    for i in range(len(cell_0_x)):
#        tempvals.append(cell_0_x[i][0:147])
#    tempmean = np.mean(tempvals,axis=0)
#    plt.plot(tempmean, lw=3, c='k')

#Fourth, for plotting: concatenate pre1sgray_data:response_data:post3sgray_data                   
#response_prepost_data[cell][orientation][block] = (
#    np.concatenate((pre1sgray_data[cell][orientation][block], 
#                    response_data[cell][orientation][block]), axis=0))
pre_response_post_data = OrderedDict()

for cell in pre1sgray_data:
    
    pre_response_post_data[cell] = OrderedDict()
    
    for orientation in pre1sgray_data[cell]:
            
        pre_response_post_data[cell][orientation] = OrderedDict()
            
        for block in range(len(pre1sgray_data[cell][orientation])):
            
            pre_temp = pre1sgray_data[cell][orientation][block]
            response_temp = response_rawdata[cell][orientation][block]
            """
            Need to figure out how to get gray after each ORI
            (use stim indices to get 0:30 gray "chunk" after stim)
            post_temp = grays[cell][orientation][block][0:30:end]
            """
            
            concat_prestim_temp = (np.concatenate((pre_temp, response_temp),
                                                  axis=0))
                                                  
            #concat_prestimpost_temp = (np.concatenate(
            #                          (concat_prestim_temp, ???),axis=0)))
            
            tempmean = np.mean(pre1sgray_data[cell][orientation][block])
            
            pre_response_post_data[cell][orientation][block] = (
                (concat_prestim_temp - tempmean) / tempmean
            )

#for ori in pre_response_post_data[26]:
#
#    cell_0_x = pre_response_post_data[26][ori]
#    plt.subplots()
#    
#    for block in cell_0_x:
#    
#        plt.plot(cell_0_x[block])
#
#    tempvals = list()
#
#    for i in range(len(cell_0_x)):
#
#        tempvals.append(cell_0_x[i][0:147])


#Gets response averages and pre-stim response averages
response_avgs = OrderedDict()
pre_response_post_avgs = OrderedDict()
pregray1s_response_avgs = OrderedDict()

for cell in response_data:
    response_avgs[cell] = OrderedDict()
    pre_response_post_avgs[cell] = OrderedDict()
    pregray1s_response_avgs[cell] = OrderedDict()

    for orientation in response_data[cell]:

        cell_ori = response_data[cell][orientation].values()
        cell_ori_long = pre_response_post_data[cell][orientation].values()
        cell_ori_gray1s = pre1sgray_dFoF_data[cell][orientation].values()
        
        #trim signals down to shortest length for averaging
        min_chunk_length = min([len(x) for x in cell_ori])
        cell_ori = (
            [x[0:min_chunk_length] for x in cell_ori])
        
        min_chunk_length_long = min([len(x) for x in cell_ori_long])
        cell_ori_long = (
            [x[0:min_chunk_length_long] for x in cell_ori_long])
        
        min_chunk_length_gray1s = min([len(x) for x in cell_ori_gray1s])
        cell_ori_gray1s = (
            [x[0:min_chunk_length_gray1s] for x in cell_ori_gray1s])

        A = np.array(cell_ori)
        Along = np.array(cell_ori_long)
        Agray = np.array(cell_ori_gray1s)

        B = np.mean(A, axis = 0)
        Blong = np.mean(Along, axis = 0)
        Bgray = np.mean(Agray, axis = 0)

        response_avgs[cell][orientation] = B
        pre_response_post_avgs[cell][orientation] = Blong
        pregray1s_response_avgs[cell][orientation] = Bgray


#example: plots all 45 degree responses

#for cell in response_avgs:

#    plt.plot(response_avgs[cell][45])

#example: plots all orientations for cell 10

#for ori in response_avgs[10]:

#    plt.plot(response_avgs[0][ori])
    
# stim_onsetframes, grays_onsetframes = (
#    getOnsetIndices(response_indices, grays_indices))
    

# use for plotting t1 vs t2 responses over time, etc
responses_means = OSC.get_ori_responses(response_avgs)

"""   
Detect responsive cells: Use Mann-Whitney U to compare each cell's orientation
response (full mean trace) to corresponding pre-gray (full mean trace)

    stimwindow -> duration of stimulus in frames (e.g., 145 frames ~ 5s stim)
    
    pthresh -> minimum mean response to be considered for 'responsive' class
    
    togglePlot -> 0 = off or 1 = on; recommended ONLY for troubleshooting
                                    (creates #cells x 8 figures)

"""
#get all orienation response means for a given cell
#cell = 6
stimwindow = 145 
pthresh = 0.05
togglePlot = 0
all_response_indices = (
                        OSC.get_responseClass(
                        response_avgs,
                        pregray1s_response_avgs, 
                        pre_response_post_avgs, 
                        stimwindow, 
                        pthresh, 
                        togglePlot)
                        )

"""
Chunk out the first ~60sec (1700 frames) of signal during gray screen at the
beginning of session - raw signal and deltaF/F via EWMA-filter signal

delayframes -> the number of frames prior to first stim onset
    (ie during gray screen)
"""
# might be able to embed in graydff loop if len(dff_ewma) == len(cellsmean)
grayraw_frames = list()

for cell in range(len(cellsmean)):
    temp_gray = (cellsmean[cell][0:delayframes])
    grayraw_frames.append(temp_gray)
    
    
graydff_frames = list()

for cell in range(len(dff_ewma)):
    temp_gray = (dff_ewma[cell][0:delayframes])
    graydff_frames.append(temp_gray)

#==============================================================================  
# Run some functions
"""
1) Save all variables from the session; to load -->
    dill.load_session(fileDirectory+filenamepkl)
"""
filenamedill = roizipname.replace('ROI.zip','SESSION.pkl')
filenamepkl = roizipname.replace('ROI.zip','VARS.pickle')
#dill.dump_session(fileDirectory+'RAW_'+filenamepkl)


# Pickle specific variables for response t1 vs. t2, etc plotting
# obj0, obj1, obj2 are created here...
# Saving the objects:
with open(fileDirectory+filenamepkl, 'w') as f:  # Python 3: open(..., 'wb')
    pickle.dump({
                'spontaneous_raw': grayraw_frames,
                'spontaneous_dff': graydff_frames,
                'response_avgs':response_avgs,
                'pregray1s_response_avgs':pregray1s_response_avgs,
                'pre_response_post_avgs':pre_response_post_avgs, 
                'responses_means':responses_means,
                'stimwindow':stimwindow,
                'pthresh':pthresh,
                'all_response_indices':all_response_indices
                
                },
                 
                f)
#==============================================================================

# need tracePlot_raw function to increase space between traces
OSC.tracePlot(grayraw_frames,
              response_indices,
              grays_indices)
plt.title('Raw signal during gray screen delay prior to stim')


OSC.tracePlot(graydff_frames,
              response_indices,
              grays_indices)
plt.title('deltaF/F (EWMA) signal during gray screen delay prior to stim')


OSC.plotStack(
              fileDirectory,
              imgname,
              roizipname,
              all_response_indices
              )

OSC.multiPlot(
              pre_response_post_avgs,
              rows,
              cols,
              range(0,rows),
              [-0.5,1.0],
              207-pre_timewindow,
              all_response_indices
              )
                  
#==============================================================================

print ('Processed data for '+datafile)
