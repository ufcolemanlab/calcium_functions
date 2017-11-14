# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 19:44:20 2017

@author: jcoleman
"""
import numpy as np
from collections import OrderedDict
from matplotlib import pyplot as plt
import os
from scipy.io import matlab

import cPickle as pickle
import readroi as roizip

# Coleman lab modules - combine into a single class?
from calcium import OpenStepCalcium as OSC
from calcium import calcium_imaging_data_fast as cidfast
from calcium import open_calcium_data_fast as ocdfast
from calcium import StepCodeFile as SCF

#==============================================================================
#PART 1 - user input
#Get step code file (*.bin), TIF image, ROI.zip files
#==============================================================================  
datafile = 'DATA_mThy6f-2_alldrift_d1_002z1.csv' #unique key-string to ID data set
datafile = 'DATA_mThy6f-2_slowpr_d1_002z1.csv' #unique key-string to ID data set

datatype = 1 # 1 (alldrift) or 2 (slowpr)

#PARAMETERS
daqChannels = 8 # number of channels recorded on DAQ
daqRecordingFreq = 1000 #sampling rate - Hz

if datatype == 1: # all drift
    stimList = [0,45,90,135,180,225,270,315]

if datatype == 2: # (slow) phase-reversal
    stimList = [0,45,90,135,180] #designate flip/flop?

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

# Set the window for averaging pre-stim gray (207 frames = 7s)
#   This is index from which to START to chunk out the LAST 1s of pre-stim gray
#   (i.e. the 1s immediately preceeding the stim onset)
pre_timewindow = 207-30 #207-75 # e.g., 207-30 = 1s from 30fps sample rate

# delay time prior to first stim onset (in frames)
# used for spontaneous activity analysis
delayframes = 1700 #(60s~1700frames)

# Parameters for classifying responsive cells
# Used in OSC.get_responseClass function below
"""   
Detect responsive cells: Use Mann-Whitney U to compare each cell's orientation
response (full mean trace) to corresponding pre-gray (full mean trace)

    stimwindow -> duration of stimulus in frames (e.g., 145 frames ~ 5s stim)
    
    pthresh -> pval trheshold (ie p=0.05) for 'responsive' class
    
    dff_thresh -> minimum mean response to be considered for 'responsive' class
        e.g., dff_thresh = 0.1 for Thy1-GCaMP6s mice
    
    togglePlot -> 0 = off or 1 = on; recommended ONLY for troubleshooting
                                    (creates #cells x 8 figures)

"""
stimwindow = 145 # number of frames used to calculate response means
pthresh = 0.05
dff_thresh = 0.1 # 0.1 thygcamp6s
togglePlot = 0

user_parameters = {'datafile':datafile,
                   'datatype':datatype,
                   'stimList':stimList,
                   'csvfiletype':csvfiletype,
                   't_0':t_0,
                   't_1':t_1,
                   't_2':t_2,
                   'mpmSamplingFreq':mpmSamplingFreq,
                   'gray_offset':gray_offset,
                   'pre_timewindow':pre_timewindow,
                   'delayframes':delayframes,
                   'stimwindow':stimwindow,
                   'pthresh':pthresh,
                   'dff_thresh':dff_thresh,
                   'togglePlot':togglePlot}
#==============================================================================

#==============================================================================
#PART 2 - user input (check files if needed?)
#==============================================================================
# Select folder containing BIN, CSV, TIF, and ZIP files for ONLY one T-series
filenameInfo = OSC.load_datafiles('raw')

#data_index = 0
#
#column1 = [os.path.split(a[0][0])[1],
#           os.path.split(a[2][0])[1],
#           os.path.split(a[3][0])[1]]
#column2 = [os.path.split(a[0][1])[1],
#           os.path.split(a[2][1])[1],
#           os.path.split(a[3][1])[1]]
#padding = 9
#for c1, c2 in zip(column1, column2):
#    print "%s %s" % (c1.ljust(padding), c2)

# is this OS-dependent?    
fileBin = filenameInfo[0][0]
csvname = filenameInfo[1][0]
imgname = filenameInfo[2][0]
roizipname = filenameInfo[3][0]

#==============================================================================
# PART 3 - do not change code up to PART 4 (this should move to a module (e.g., OSC))
#==============================================================================
# Begin decoding and data extraction
nTimeStampCodes = len(stimList)
code_file = SCF.StepCodeFile(fileBin,
                             nTimeStampCodes,
                             daqChannels,
                             daqRecordingFreq)

stims = code_file.get_stim_angles(stimList)

Bframes = [ts[0] for ts in code_file.timestamps[2]]

# parameters for mulitplot function
a = roizip.read_roi_zip(roizipname)
rows = len(a) # number of ROIs
cols = nTimeStampCodes # number of distinct stims


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
dff, dff_ewma, dff_offset = (cidfast.run_deltaf_ewma(cellsmean,
                                                     t_0,
                                                     t_1,
                                                     t_2,
                                                     mpmSamplingFreq)
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
    
    pthresh -> pval trheshold (ie p=0.05) for 'responsive' class
    
    dff_thresh -> minimum mean response to be considered for 'responsive' class
        e.g., dff_thresh = 0.1 for Thy1-GCaMP6s mice
    
    togglePlot -> 0 = off or 1 = on; recommended ONLY for troubleshooting
                                    (creates #cells x 8 figures)

"""
#get all orienation response means for a given cell
#cell = 6
#stimwindow = 145 
#pthresh = 0.05
#dff_thresh = 0.01 
#togglePlot = 0
all_response_indices = (OSC.get_responseClass(
                        response_avgs,
                        pregray1s_response_avgs, 
                        pre_response_post_avgs, 
                        stimwindow, 
                        pthresh,
                        dff_thresh,
                        togglePlot))

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
## get xy centroid points or close
#centroids_x, centroids_y = (OSC.get_cellCentroids(
#    roizipname,
#    0))

# keeps geting error all of a suddnen - need to work for spatial info 
#Traceback (most recent call last):
#  File "<stdin>", line 3, in <module>
#  File "/Users/jcoleman/anaconda/lib/python2.7/site-packages/calcium/OpenStepCalcium.py", line 524, in get_cellCentroids
#    xlist = [a[j][1][i][1] for i in range(len(a[j][1]))]
#IndexError: invalid index to scalar variable.
   
centroids_x =[]
centroids_y =[]

#==============================================================================  
# Save processed data to files (PICKLE, MAT)
"""
1) Save all variables from the session; to load -->
    dill.load_session(filenamepkl)
2) Save selected variables to a pickle file
3) Save selected variables to a Matlab (mat) file
"""
# Create a sensible name
filenamedill = roizipname.replace('ROI','SESSION').replace('.zip','.pkl')
filenamepkl = roizipname.replace('ROI','VARS').replace('.zip','.pickle')
filenamemat = roizipname.replace('ROI','VARS').replace('.zip','.mat')

# 1) Save 'workspace' to a PKL file (re-loading is an issue, not working on Mac)
#dill.dump_session('RAW_'+filenamepkl)

# 2) Pickle specific variables for response t1 vs. t2, etc plotting
with open(filenamepkl, 'w') as f:  # Python 3: open(..., 'wb')
    pickle.dump({'user_parameters': user_parameters,
                 'grayraw_frames': grayraw_frames,
                 'graydff_frames': graydff_frames,
                 'centroids_x': centroids_x,
                 'centroids_y': centroids_y,
                 'response_avgs':response_avgs,
                 'pregray1s_response_avgs':pregray1s_response_avgs,
                 'pre_response_post_avgs':pre_response_post_avgs, 
                 'responses_means':responses_means,
                 'stimwindow':stimwindow,
                 'pthresh':pthresh,
                 'all_response_indices':all_response_indices
                 }, f)

# 3) Save initial "gray frames" to a MAT file
savedVariables = {'grayraw_frames': grayraw_frames,
                 'graydff_frames': graydff_frames,
                 'centroids_x': centroids_x,
                 'centroids_y': centroids_y,
                 'README': 'm = cell, n = time (i.e. frame#)'
                 }              
matlab.savemat(filenamemat, savedVariables)

# Need to figure out saving nested dict into MAT file - below works for cell#0
#matlab.savemat('keys.mat', {'response_avgs_keys':response_avgs.keys(),
#                                          'response_avgs_vals':response_avgs[0].values()
#                                          })

#==============================================================================
# PART 4 - Run some functions using the procesed data from PART 3
#==============================================================================
# need tracePlot_raw function to increase space between traces
OSC.tracePlot(grayraw_frames,
              response_indices,
              grays_indices,
              'False',
              'raw')
plt.xlabel('Frame (~'+ str(mpmSamplingFreq) +' fps)')
plt.ylabel('Mean gray value')
plt.title('Raw signal during gray screen delay prior to stim')


OSC.tracePlot(graydff_frames,
              response_indices,
              grays_indices,
              'True',
              'dFoF')
plt.xlabel('Frame (~'+ str(mpmSamplingFreq) +' fps)')
plt.ylabel('deltaF/F')
plt.title('Filtered/normalized signal during gray screen delay prior to stim')


OSC.plotStack(imgname,
              roizipname,
              all_response_indices)


plotcells = 10 # = len(pre_response_post_avgs) to plot all cells 
startcell = 0 # = 0 to plot all cells
rows = plotcells
OSC.multiPlot(pre_response_post_avgs,
              rows,
              cols,
              range(startcell,rows+startcell),
              [-0.5,1.0],
              207-pre_timewindow,
              all_response_indices)
              
                  
#==============================================================================

print ('Processed data for '+ datafile)
