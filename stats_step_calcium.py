# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 11:28:27 2017

@author: jcoleman
"""

import numpy as np
import matplotlib.pyplot as plt
import tkFileDialog
import csv
import Tkinter as tk
import scipy.stats as stats
from scipy import signal
from scipy.signal import butter, lfilter, freqz
from collections import OrderedDict


#test function for finding p-values
def get_cell_response(cell, window):
    flips, flops = Data.get_flip_flops(cell, Data.strobe_down, Data.stim_up, Data.Bframe_up)
    flip_vals = np.array(flips.values())
    flip_vals = flip_vals.transpose()
    for i in range(len(flip_vals)-window):
        groups = [flip_vals[j] for j in range(i, i + window)]
        if stats.f_oneway(*groups)[1] < 0.05:
            print "response detected"

# response_avgs[cell][orientation]
# pregray1s_response_avgs[cell][orientation]

def anova_response(data,pthresh):
    for cell in data:
        #f,p=stats.f_oneway(*pregray1s_response_avgs[cell].values())
        fstat,pval=stats.kruskal(*data[cell].values())
        if pval < pthresh:
            print "Response in cell "+str(cell)+" p="+str(pval)
        elif pval > pthresh:
            print "NO response in cell "+str(cell)
            
    return fstat, pval
    
    
def is_responsive(stim_stats, cells):
    responsive = list()
    nonresponsive = list()
    for i in stim_stats:
        if stim_stats[i][1] < 0.01:
            #print str(i) + " responsive"
            responsive.append(i)
        elif stim_stats[i][1] >= 0.01:
            nonresponsive.append(i)
    
    return responsive, nonresponsive


def get_response_graph(delta_f, flips, flops, figure):
    plt.figure(figure)
    for i in range(len(delta_f)):
        plt.plot(delta_f[i] + float(i))
    plt.axis('tight')
    for ts in flips[0].keys():
        plt.axvline(x = ts[0], color = 'red')
        plt.axvline(x = ts[1], color = 'green', ls = 'dashed')
    for ts in flops[0].keys():
        plt.axvline(x = ts[1], color = 'black')


def get_prestim_values(cells, stim, i):
    values = dict()
    for ts in stim[i]:
        values[(ts[0] - (ts[1] - ts[0]), ts[0])] = cells[i][ts[0] - (ts[1] - ts[0]): ts[0]]
    
    return values
    

def is_thresh_responsive(stim, cells):
    pass


def normalize_signals(delta_f):
    norms = list()
    for i in range(len(delta_f)):
        norms.append(delta_f[i] / np.max(delta_f[i]))
    
    return norms
    
    
def butter_lowpass(cutoff, fs, order):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    
    return y
    
#    #Example:
#    cells_lopass = list()
#    cutoff = 0.8 # in Hz Rose et al. used 0.8Hz
#    fs = 30 # sampling freq in Hz
#    order = 2 # not sure what this should be exactly?
#    for i in range(len(cells)):
#        data = cells[i][:]
#        y = butter_lowpass_filter(data, cutoff, fs, order)
#        cells_lopass.append(y)
#        
#    delta_f = calc_delta_f(cells_lopass)
#    Example: 'cells' (raw instensity data)
#    plt.plot(butter_lowpass_filter(a,1.0,30,1),color='cyan')
    

def cumplot(x):
    #% plots cdf function as a solid line connecting points.

    #% plot(sort(abs(x)),((1:length(x))./length(x)),varargin{:})
    #plot(sort(x),((1:length(x))./length(x)),varargin{:})

    #plt.subplots()
    xvals = np.sort(x)
    nvals = np.float(len(x)) # need this as integer to avoid 'truncation' during division below
    ytemp = np.array(range(len(x)))
    yvals = list()
    for i in range(len(x)):
        yvals.append(ytemp[i]/nvals)
    plt.plot(xvals,yvals)
    plt.ylabel('Cumulative probability')

    
    

## ANOVA comparing grays and ori responses for one cell
#pthresh=0.05
#for cell in pregray1s_response_avgs:
#    tempmerge=list()
#    for ori in pregray1s_response_avgs[cell]:
#        tempmerge.append(pregray1s_response_avgs[cell][ori])
#        tempmerge.append(response_avgs[cell][ori])
#    
#    fstat,pval=stats.kruskal(*tempmerge)
#    if pval < pthresh:
#            print "Response in cell "+str(cell)+" p="+str(pval)
#    elif pval > pthresh:
#        print "NO response in cell "+str(cell)
        
        
## ANOVA comparing mean grays and mean ori responses for one cell
#pthresh=0.001
#postwindow = 30
#for cell in pregray1s_response_avgs:
#    tempmerge=list()
#    for ori in pregray1s_response_avgs[cell]:
#        tempmerge.append(pregray1s_response_avgs[cell][ori])
#        tempmerge.append(response_avgs[cell][ori][0:postwindow])
#    
#    fstat,pval=stats.kruskal(*tempmerge)
#    if pval < pthresh:
#            print "Response in cell "+str(cell)+" p="+str(pval)
#    elif pval > pthresh:
#        print "NO response in cell "+str(cell)
        
        
## ANOVA comparing all cells for one orientation response
#pthresh=0.05
#tempmerge=list()
#ori = 0
#for cell in response_avgs:
#
#    tempmerge.append(response_avgs[cell][ori])
#    
#fstat,pval=stats.kruskal(*tempmerge)
#if pval < pthresh:
#    print "Response in cell "+str(cell)+" p="+str(pval)
#elif pval > pthresh:
#    print "NO response in cell "+str(cell)
#    
#anova_response(response_avgs,.001)
#
#
#response_avgs_means = get_ori_responsemeans(response_avgs)
#response_gray1savgs_means = get_ori_responsemeans(pregray1s_response_avgs)

#-----------------------------------------------------------------------------------------
## use K-S to compare each cell to cell0 (means of all ori)
## consider using only 1-2s window of post-stim data (currently entire 5s window used for mean, 1s preceding for gray mean)
#baseline_cell = 0
#pthresh = 0.01
#postwindow = 30 #window for analysis, e.g., 1sec after onset in frames, fps = 30
##cell = 9
##ori = 0
#for cell in response_data:
#    for ori in response_data[cell]:
#        aa=response_data[cell][ori] # dFoF stim = 5 trials of [cell][orientation]
#        bb=pre1sgray_dFoF_data[cell][ori] # dFoF gray = 5 trials of [cell][orientation]
#        
#        # attempt 2-sample KS test for pre1sgray_dFoF_data[cell][ori] vs response_data[cell][ori]
#        aanew=list() #response data
#        for trial in aa:
#            aanew.append(aa[trial][0:postwindow])
#        aanew = np.concatenate(aanew,axis=0)
#        
#        bbnew=list() #gray data
#        for trial in bb:
#            bbnew.append(bb[trial])
#        bbnew = np.concatenate(bbnew,axis=0)
#        
#        #KSstat,pval = stats.ks_2samp(aanew-min(aanew), bbnew-min(bbnew))
#        KSstat,pval = stats.ks_2samp(aanew, bbnew)
#        
#        if pval < pthresh:
#            print "Response in cell "+str(cell)+" for "+str(ori)+"deg"+" p="+str(pval)
#        elif pval > pthresh:
#            print "NO response in cell "+str(cell)+" for "+str(ori)+"deg"
#        
#plt.subplots()
#cumplot(aanew) #*args
#cumplot(bbnew)
#plt.title('Cell# '+str(cell)+'ori='+str(ori)+'deg'+'; p='+str(pval))  
#-----------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
# use anova on response average/trial (response_data[cell][ori][trial])
##baseline_cell = 0
##pthresh = 0.01
##postwindow = 30 #window for analysis, e.g., 1sec after onset in frames, fps = 30
##cell = 9
##ori = 0
#for cell in response_data:
#    for ori in response_data[cell]:
#        stimresponse_trials = response_data[cell][ori] # dFoF stim = 5 trials of [cell][orientation]
#        grayresponse_trials = pre1sgray_dFoF_data[cell][ori] # dFoF gray = 5 trials of [cell][orientation]
#        
#        # attempt 2-sample KS test for pre1sgray_dFoF_data[cell][ori] vs response_data[cell][ori]
#        stimresponse_trialmean=list() #response data
#        for trial in stimresponse_trials:
#            stimresponse_trialmean.append(np.mean(stimresponse_trials[trial][0:postwindow]))
#        #stimresponse_trialmean = np.concatenate(stimresponse_trialmean,axis=0)
#        
#        grayresponse_trialmean=list() #gray data
#        for trial in grayresponse_trials:
#            grayresponse_trialmean.append(np.mean(grayresponse_trials[trial]))
#        #grayresponse_trialmean = np.concatenate(grayresponse_trialmean,axis=0)
#        
#        fstat,pval=stats.kruskal(stimresponse_trialmean,grayresponse_trialmean)
#        
#        if pval < pthresh:
#            print "Response in cell "+str(cell)+" for "+str(ori)+"deg"+" p="+str(pval)
#        elif pval > pthresh:
#            print "NO response in cell "+str(cell)+" for "+str(ori)+"deg"  

#-----------------------------------------------------------------------------------------
            