# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 23:41:08 2017

@author: jcoleman
"""
import Tkinter as tk
import tkFileDialog
from glob import glob
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
import dill # need to download
import cPickle as pickle

# Coleman lab modules - combine into a single class?
from calcium import OpenStepCalcium as OSC
from calcium import calcium_imaging_data_fast as cidfast
from calcium import open_calcium_data_fast as ocdfast
from calcium import StepCodeFile as SCF

#%% ONLY RUN CODE BELOW AFTER EACH CSV FILE HAS BEEN PICKLED
# Getting back the objects:
    
#
responses_001D2Z1t0hz05 = OSC.load_responses()
responses_001D3Z1t1hz05 = OSC.load_responses()
#responses_001D3Z1t2hz05 = OSC.load_responses()
#responses_001D4hz1 = OSC.load_responses()
#responses_001D4hz5 = OSC.load_responses()
#responses_002D4hz1 = OSC.load_responses()
#responses_002D4hz5 = OSC.load_responses()
#    
#responses_6f2_002z1_d1 = OSC.load_responses() #means[0], indices[1]

#%% NORMALIZE RESPONSE DATA
from copy import deepcopy

def load_wholepickle():
    
    """
    ['user_parameters',
    'stimwindow',
    'pre_response_post_avgs',
    'centroids_y',
    'centroids_x',
    'all_response_indices',
    'grayraw_frames',
    'responses_means',
    'response_avgs',
    'pregray1s_response_avgs',
    'graydff_frames',
    'pthresh']
    
    """
    
    root = tk.Tk()
    root.withdraw()
    root.update()
    
    filepath = tkFileDialog.askopenfilename(parent=root,title='Choose a pickle file ...')

    with open(filepath) as f:  # Python 3: open(..., 'rb')
        alldata = pickle.load(f)
    
    return alldata

    
def normalize_data(data):
    
    data_norm = deepcopy(data)
    
    for cell in data_norm:
        
        for ori in data_norm[cell]:
            
            temp = data_norm[cell][ori]
            
            temp -= min(temp)
            temp /= max(temp)
            
            data_norm[cell][ori] = temp
        
    #    session_gray_avgs = np.array(session_gray_avgs).astype(np.float)
    return data_norm

# load data of interest from the pickle file        
all = load_wholepickle()

# normalize the data
normalized_data = normalize_data(all['response_avgs'])

# check some data
plt.plot(all['response_avgs'][6][45])
plt.plot(normalized_data[6][45])

#%% 

def plot_heatmap(data, stimList, cellsort):
    #Heat plots
    """
    data -> dictionary of data[cell][orientation]
    stimList -> list of integers (orientations)
    cellsort -> 'sort max' to rank by timing of max response
                'sort no' to leave unsorted (list by cell ID rank)
                ? 'sort 0' to sort based on sort-max of 0deg
    
    """
    
    # Sort cells by latency to max intensity value       
    # group all cells into one variable for one orientation
    for stim in range(len(stimList)):
        
        datatemp = []
    
        for cell in data:
            datatemp.append(data[cell][stimList[stim]])
        
        
        if cellsort == 'sort max':
            sorted_response_avgs = sorted(datatemp, key=lambda x: x.argmax())
            ## check sorted data:
            #plt.subplots()
            #for cell in range(len(sorted_response_avgs)):
            #    plt.plot(sorted_response_avgs[cell]+0.2*cell)
        
            sorted_response_avgs_keys = np.argsort(np.argmax(datatemp, axis=1))
            
        elif cellsort == 'sort no':
            sorted_response_avgs = datatemp
            sorted_response_avgs_keys = data.keys()                
        
        
#        if all_orientations==1:
#            np.concat()
#            
#        #apply hanning window smoothing?
#        
#        else:
        cidfast.plotHeatAvgs(sorted_response_avgs, sorted_response_avgs_keys, 30, 0, 1.0, 'deltaF/F (%)')
        plt.title('Sorted: '+str(stimList[stim]))


#user_parameters['stimList']
stimList = [0,45,90,135,180,225,270,315]

plot_heatmap(normalized_data, stimList, 'sort max')

#%%
## all cells
#plt.subplots()
#for cell in responses_means_D2t0:
#    for ori in responses_means_D2t0[cell]:
#        plt.plot(responses_means_D3t1[cell][ori], responses_means_D3t2[cell][ori],'o')
#        
#plt.title('D3t1 v D3t2')
#
#
#plt.ylim([-0.5, 1.0])
#plt.xlim([-0.5, 1.0])

# filter by responders

t1=list()
t2=list()

#t1,t2 = (OSC.plot_time_responses(responses_means_001D4hz1, responses_indices_001D4hz1, 
#                             responses_means_002D4hz5, responses_indices_002D4hz5, 0))
                             
t1,t2 = (OSC.plot_time_responses(responses_001D2Z1t0hz05[0], responses_001D2Z1t0hz05[1], 
                             responses_001D3Z1t1hz05[0], responses_001D3Z1t1hz05[1], 0))

#%% plot best fit
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from random import sample

# sample data
x=list()
y=list()

for i in range(len(t1)):
    
    if (t1[i] == 'NaN') is False:
        x.append(t1[i])#np.arange(10)
        y.append(t2[i]) 

# fit with np.polyfit
m, b = np.polyfit(x, y, 1)

#plt.plot(x, y, '.')
plt.plot(x, np.array(m)*x + b, '-')

plt.plot([0,0.5,1.0], [0,0.5,1.0], 'k-')

#axes
print 'min: ' + str(np.floor(min(y))) 
print 'max: ' + str(np.ceil(max(y))) 



f,p = stats.pearsonr(x,y)
# what about with min-max normalized? normalized = (x-min(x))/(max(x)-min(x))

plt.title('Pearson\'s p = '+str(p))

# shuffle x to see if random correlation (need xcorr)
x_shuf = sample(x,len(x))
f,p = stats.pearsonr(x_shuf,y)
np.correlate(x,y,'full')

#%%
import numpy as np
import matplotlib.pyplot as plt
b = []
for i in range(50):
    a = np.random.normal(5,i+1,10)
    b.append(a)
    
c = np.array(b)
cm = np.corrcoef(c)
plt.imshow(cm,interpolation='nearest')
#plt.imshow(cm,interpolation='hanning')
plt.colorbar()
plt.show()



#%%
# plot responses at orientation for a cell all 3 days

def plt_mean_oriR(data,colormark):
        
    templist_plot = list()
    
    for i in range(len(data)):
        
        tempkey = data.keys()
        tempkey_plot = [1,0,7,6,3,5,2,4] # tempkey indices for plotting order
        templist_plot.append(data[tempkey[tempkey_plot[i]]])
    
    plt.plot(templist_plot, colormark)       
    plt.title('cell# ' + str(cell))
    #plt.xlabel([225,0,315,90,45,180,135,270])
    plt.xlabel([135,90,45,0,315,270,225,180]) # GCaMP6s Nature paper

for cell in responses_001D2Z1t0hz05[0]: 
    plt.subplots()
    plt_mean_oriR(responses_001D2Z1t0hz05[0][cell], 'bo-')
    plt_mean_oriR(responses_001D3Z1t1hz05[0][cell], 'ro-')       
    
    plt.ylim([-0.5, 1.0])
    plt.xlim([-1.0, 8.0])
    
    plt.legend(['D1 1Hz', 'D3 1Hz'])
    

 #%%
#polar plots for normalized data - one celabs

for i in range(len(normalized_data)):
    r_avgs_cell = normalized_data[i]
    polar_data = []
    for ori in r_avgs_cell:
        polar_data.append(np.mean(r_avgs_cell[ori]))
    
    #ra = all['responses_means'] # not sure if should start with normalized vector data
    #ra_cell6 = ra[6]
    #aa=(ra_cell6.values())
    #aa-=min(aa)
    #aa/=max(aa)
    
    stimListRads = np.array(r_avgs_cell.keys())*(np.pi/180)
    
    #plt.figure()
    plt.title('cell '+ str(i))
    plt.polar(stimListRads,polar_data, 'o')



        


#Plot selected cells
#OSC.plotROIavgs(response_avgs,filenamepkl.replace('_SESSION.pkl',''),[6,7,9,10,11,13,14,15,20,25,26],[-0.5,2.0])
#OSC.plotROIavgs(pre_response_post_avgs,filenamepkl.replace('_SESSION.pkl',''),[6,7,9,10,11,13,14,15,20,25,26],[-0.5,3.0])


## find repsonsive cell (binary)
#c=list()
#for cell in all_response_indices:
#    b = list()
#    for ori in all_response_indices[cell]:
#        if all_response_indices[cell][ori][0] == 1:
#            b.append(1)
#        elif all_response_indices[cell][ori][0] != 1:
#            b.append(0)
#    if sum(b) > 0:
#        print cell
#        c.append(1)
#    elif sum(b) == 0 or sum(b) < 0:
#        c.append(0)
#j=0 #(cell)
#if c[j] > 0: # if cell is responsive (see loop above)
#    print 'do something'