# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 23:41:08 2017

@author: jcoleman
"""
# Coleman lab modules - combine into a single class?
from calcium import OpenStepCalcium as OSC
from calcium import calcium_imaging_data_fast as cidfast
from calcium import open_calcium_data_fast as ocdfast
from calcium import StepCodeFile as SCF

# Python modules
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt
import dill # need to download
import cPickle as pickle

#%% ONLY RUN CODE BELOW AFTER EACH CSV FILE HAS BEEN PICKLED
# Getting back the objects:
    
#
datafile = 'D4_001_Z1_hz5'

if datafile == 'D2_001_Z1t0_hz05':

    filedir = '/Users/jcoleman/Documents/--DATA--/in vivo gcamp analysis/thygcamp6s_LT4(9-10-17)/'
    picklefile = 'mThy6s2_alldrift_D2_001_VARS.pickle'
    responses_means_001D2Z1t0hz05, responses_indices_001D0Z1t0hz05 = OSC.load_responses(filedir, picklefile)
    

if datafile == 'D3_001_Z1t1_hz05':
    
    filedir = '/Users/jcoleman/Documents/--DATA--/in vivo gcamp analysis/thygcamp6s_LT4(9-10-17)/'
    picklefile = 'mThy6s2_alldrift_D3_001Z1_VARS.pickle'
    responses_means_001D3Z1t1hz05, responses_indices_001D3Z1t1hz05 = OSC.load_responses(filedir, picklefile)
    

if datafile == 'D3_001_Z1t2_hz05':
    
    filedir = '/Users/jcoleman/Documents/--DATA--/in vivo gcamp analysis/thygcamp6s_LT4(9-10-17)/'
    picklefile = 'mThy6s2_alldrift_D3_001Z1t2_VARS.pickle'
    responses_means_001D3Z1t2hz05, responses_indices_001D3Z1t2hz05 = OSC.load_responses(filedir, picklefile)
    

if datafile == 'D4_001_Z1_hz1':

    filedir = '/Users/jcoleman/Documents/--DATA--/in vivo gcamp analysis/thygcamp6s_D4 5Hz (9-30-17)/'
    picklefile = 'mThy6s2_alldrift_D4_001Z1_VARS.pickle'
    responses_means_001D4hz1, responses_indices_001D4hz1 = OSC.load_responses(filedir, picklefile)
    

if datafile == 'D4_001_Z1_hz5':

    filedir = '/Users/jcoleman/Documents/--DATA--/in vivo gcamp analysis/thygcamp6s_D4 5Hz (9-30-17)/'
    picklefile = 'mThy6s2_alldrift_D4_001Z1hz5_VARS.pickle'
    responses_means_001D4hz5, responses_indices_001D4hz5 = OSC.load_responses(filedir, picklefile)
    
    
if datafile == 'D4_002_Z1_hz1':

    filedir = '/Users/jcoleman/Documents/--DATA--/in vivo gcamp analysis/thygcamp6s_D4 5Hz (9-30-17)/'
    picklefile = 'mThy6s2_alldrift_D4_002Z1hz1_VARS.pickle'
    responses_means_002D4hz1, responses_indices_002D4hz1 = OSC.load_responses(filedir, picklefile)
    
    
if datafile == 'D4_002_Z1_hz5':

    filedir = '/Users/jcoleman/Documents/--DATA--/in vivo gcamp analysis/thygcamp6s_D4 5Hz (9-30-17)/'
    picklefile = 'mThy6s2_alldrift_D4_001Z1hz5_VARS.pickle'
    responses_means_002D4hz5, responses_indices_002D4hz5 = OSC.load_responses(filedir, picklefile)


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
                             
t1,t2 = (plot_time_responses(responses_means_001D4hz1, responses_indices_001D4hz1, 
                             responses_means_001D4hz5, responses_indices_001D4hz5, 0))

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

for cell in responses_means_002D4hz1: 
    plt.subplots()
    plt_mean_oriR(responses_means_002D4hz1[cell], 'bo-')
    plt_mean_oriR(responses_means_002D4hz5[cell], 'ro-')       
    
    plt.ylim([-0.5, 1.0])
    plt.xlim([-1.0, 8.0])
    
    plt.legend(['1Hz', '5Hz'])
    

        

        


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