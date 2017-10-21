# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 15:10:12 2016

@author: jcoleman
@author: jtrinity
"""

#%%
import matplotlib.pyplot as plt

import matplotlib
matplotlib.use("wx")
from pylab import *

def topfig(userX,userY):
    figmgr = get_current_fig_manager()
    figmgr.canvas.manager.window.raise_()
    geom = figmgr.window.geometry()
    x,y,dx,dy = geom.getRect()
    figmgr.window.setGeometry(userX, userY, dx, dy)

def figAllResponsiveANOVA(keyarray):
    plt.subplots()
    r_avgs_sorted_keys = keyarray #sorted(r_both_avgs, key = lambda key: np.argmax(r_both_avgs[key]))
    r_avgs_sorted = list()
    for key in r_avgs_sorted_keys:
        r_avgs_sorted.append(r_both_avgs[key])
    get_colormaps(r_avgs_sorted,r_avgs_sorted_keys, 0, 0.6)
    plt.title('Responsive Cell Averages (ANOVA)').set_family('helvetica')
    plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'limegreen', ls = 'dashed')
    plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]) + len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'green', ls = 'dashed')
    plt.axis('tight')
    #plt.savefig('Responsive_Cell_Averages_(ANOVA).eps', format = 'eps', dpi = 1200)

def figAllThresh(keyarray): 
    plt.subplots()
    all_avgs_sorted_keys = keyarray #sorted(all_avgs, key = lambda key: np.argmax(all_avgs[key]))
    all_avgs_sorted = list()
    for key in all_avgs_sorted_keys:
        all_avgs_sorted.append(all_avgs[key])
    get_colormaps(all_avgs_sorted, all_avgs_sorted_keys, 0, 0.6)
    plt.title('All Cell Averages - Thresholded')
    plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'green', ls = 'dashed')
    plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]) + len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'green', ls = 'dashed')
    plt.axis('tight')

#==============================================================================
# for i in range(len(r_flips)):
#     plt.plot(r_flips[i]+(i*0.25))
# 
# for i in range(len(nr_flips)):
#     tempa = nr_flips[i]
#     a = tempa[tempa >0.5]
# 
# # plot trace on new figure window    
# plt.figure(11)
# plt.plot(cells[100])    
#==============================================================================

#%% PLOTS TRACES AND HEAT - slow stim (0.2Hz)

# plot all cell averages
data = all_avgs_sorted
threshkey_all=list()
plt.subplots()
#axhline(y=0.25,linewidth=1, color='r')
for i in range(len(all_avgs_sorted_keys)):   
    tempkey = all_avgs_sorted_keys[i]
    if max(data[tempkey]) > 0.25:
        #plt.pause(0.8)
        threshkey_all.append(tempkey)
        currPlot=plt.plot(data[tempkey])
#topfig(10,10)
plt.title('Responsive Cell Averages (th)').set_family('helvetica')
plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'limegreen', ls = 'dashed')
plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]) + len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'green', ls = 'dashed')
plt.axis('tight')

all_avgs_sorted_keys = threshkey_all #sorted(all_avgs, key = lambda key: np.argmax(all_avgs[key]))
all_avgs_sorted_thresholded = list()
for key in all_avgs_sorted_keys:
    all_avgs_sorted_thresholded.append(all_avgs[key])
get_colormaps(all_avgs_sorted_thresholded, all_avgs_sorted_keys, 0, 0.6)
plt.title('Responsive Cell Averages (th)')
plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'green', ls = 'dashed')
plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]) + len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'green', ls = 'dashed')
plt.axis('tight')
#topfig(800)

#%% PLOTS TRACES - slow stim (005 = 0.2Hz)
#data = all_flip_avgs
# plot responsive avgs
data = all_avgs_sorted
datakeys = all_avgs_sorted_keys

for avg in data:
    avg -= min(avg)
    avg /= max(avg)
    data = np.array(data).astype(np.float)
get_colormaps(data, datakeys, 0, 1.0)
plt.title('All Cell Averages')
plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'green', ls = 'dashed')
plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]) + len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'green', ls = 'dashed')
plt.axis('tight')
        
d=list()
plt.subplots()
#axhline(y=0.2,linewidth=1, color='r')
for i in range(len(datakeys)):   
    tempkey = datakeys[i]
    if max(data[tempkey]) > 0.25:
        #plt.pause(0.8)
        d.append(tempkey)
        currPlot=plt.plot(data[tempkey])

#topfig(10,10)
plt.title('Responsive Cell Averages (ANOVA + th)').set_family('helvetica')
plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'limegreen', ls = 'dashed')
plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]) + len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'green', ls = 'dashed')
plt.axis('tight')

# plot non-responsive avgs
data = nr_both_avgs
d=list()
plt.subplots()
axhline(y=0.2,linewidth=1, color='r')
for i in range(len(nr_avgs_sorted_keys)):   
    tempkey = nr_avgs_sorted_keys[i]
    if max(data[tempkey]) > 0.25:
        #plt.pause(0.8)
        d.append(tempkey)
        currPlot=plt.plot(data[tempkey])

#topfig(10,10)
plt.title('NonResponsive Cell Averages (ANOVA + th)').set_family('helvetica')
plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'limegreen', ls = 'dashed')
plt.axvline(x = len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]) + len(r_flip_both_avgs[r_flip_both_avgs.keys()[0]]), color = 'green', ls = 'dashed')
plt.axis('tight')


# plot responsive traces
data = cells
d=list()
plt.subplots()
axhline(y=0.2,linewidth=1, color='r')
for i in range(len(r_avgs_sorted_keys)):   
    tempkey = r_avgs_sorted_keys[i]
    #plt.pause(0.8)
    d.append(tempkey)
    currPlot=plt.plot(data[tempkey]+(i*.4))
plt.title('Responsive Cell Averages (ANOVA)').set_family('helvetica')
    
# plot non-responsive traces
data = cells
d=list()
plt.subplots()
axhline(y=0.2,linewidth=1, color='r')
for i in range(len(nr_avgs_sorted_keys)):   
    tempkey = nr_avgs_sorted_keys[i]
    #plt.pause(0.8)
    d.append(tempkey)
    currPlot=plt.plot(data[tempkey]+(i*.4))
plt.title('NonResponsive Cell Averages (ANOVA)').set_family('helvetica')    

#==============================================================================
# plt.figure(12)
# for i in range(len(a)):
#     plt.plot(a[i])
#==============================================================================
#==============================================================================
#%% PLOT ROIs over image
     
# Create figure from TIF    
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm

import numpy as np 
imgdir = '/Users/jcoleman/Documents/PYTHON/calcium imaging/sample data/'

img = mpimg.imread(imgdir+'test.tif')
plt.subplots()
imgplot = plt.imshow(img)
imgplot = plt.imshow(img, cmap="gray")#, origin='lower')
#==============================================================================

#%% 
#read ROI zip file from Fiji/ImageJ (rectangle ROI objects)
import readroi as roizip
imgdir = '/Users/jcoleman/Documents/PYTHON/calcium imaging/sample data/'
a=roizip.read_roi_zip(imgdir+'RoiSetsoma.zip')

# Plot ROI shapes from Fiji ROI zip file
#plt.subplots() 
for j in range(len(a)):
    xy=list()
    for i in range(len(a[j][1])):
        xy.append([a[j][1][i][1],a[j][1][i][0]])
        plt.plot(xy[i][0], xy[i][1], 'm.')

#%% 
## Define centroid from square ROIs
#roiCentroids = list()
#roiXY = list()
#    
#for i in range(len(a)):
#    y1 = a[i][1][0][0]
#    y3 = a[i][1][2][0]  
#    x1 = a[i][1][0][1]
#    x2 = a[i][1][1][1]   
#     
#    roiTempCentroidxy = [x1+(0.5*(x2-x1)), y1+(0.5*(y3-y1))];
#    roiCentroids.append(roiTempCentroidxy)
#    
#    roiTempxy = [x1, y1];
#    roiXY.append(roiTempxy)


#%%
#import matplotlib
#
#nx, ny = 10, 10
#poly_verts = [(1,1), (5,1), (5,9),(3,2),(1,1)]
#
## Create vertex coordinates for each grid cell...
## (<0,0> is at the top left of the grid in this system)
#x, y = np.meshgrid(np.arange(nx), np.arange(ny))
#x, y = x.flatten(), y.flatten()
#
#points = np.vstack((x,y)).T
#
#grid = matplotlib.path.Path.contains(points, poly_verts)
#grid = grid.reshape((ny,nx))
#
#print grid
#
def outline_to_mask(line, x, y):
    """
    Create mask from outline contour

    Parameters
    ----------
    line: array-like (N, 2)
    x, y: 1-D grid coordinates (input for meshgrid)

    Returns
    -------
    mask : 2-D boolean array (True inside)
    """
    import matplotlib.path as mplp
    mpath = mplp.Path(line)
    X, Y = np.meshgrid(x, y)
    points = np.array((X.flatten(), Y.flatten())).T
    mask = mpath.contains_points(points).reshape(X.shape)
    return mask

   
#==============================================================================
# plot ROIs based on keys
#plt.subplots()
myarray = np.asarray(roiCentroids)

allROIs = myarray[threshkey_all]
catROIs = myarray[r_avgs_sorted_keys]
#flipROIs = myarray[flipKey]
#flopROIs = myarray[flopKey]
#grayROIs = myarray[grayKey]

plt.plot(myarray[:,0],myarray[:,1], 'o', markerfacecolor='none', markeredgecolor='b')
#plt.plot(allROIs[:,0],allROIs[:,1], 'o', markerfacecolor='none', markeredgecolor='r')
#plt.plot(catROIs[:,0],catROIs[:,1], 'y,')
plt.axis('tight')
plt.axis('equal')


# plot ROIs as rainbow colors
colors = cm.rainbow(np.linspace(0, 1, len(myarray)))
colors2 = ( np.random.permutation ( colors ))

for j in range(len(myarray)):
    colortemp = colors2[j]
    plt.plot(myarray[j,0],myarray[j,1], 'o', markerfacecolor='none', markeredgecolor=colortemp)


# Data from fast stim

#one ROI, entire trace (6500)
def figFastStimTrace(data,step_factor):
    fig, ax = plt.subplots()
    #topfig(10,10)
    for k in range(len(data)):
        #roi_fast1 = cells[2][:]
        plt.plot(data[k][:]+k*step_factor) 

#%%        
import numpy as np
import cv2


imgdir = '/Users/jcoleman/Documents/PYTHON/calcium imaging/coutours/' 
im = cv2.imread(imgdir+'test22.jpg')
imgray = cv2.cvtColor(im,cv2.COLOR_BGR2GRAY)
ret,thresh = cv2.threshold(imgray,127,255,0)
im2, contours, hierarchy = cv2.findContours(thresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)

img=cv2.imshow('im',im)
cv2.drawContours(img, contours[0], -1, (0,255,0), 3)

#cnt = contours[0]
#cv2.drawContours(img, [cnt], 0, (0,255,0), 3)

#%%
cv2.imshow('im',im)
cv2.destroyAllWindows()
        
#==============================================================================

#%%
#====== plottting session timestamps
plt.subplots()
a=np.asarray(flipTimeStamps)
b=np.asarray(flopTimeStamps)
plt.plot(a[1:148,1],np.zeros(147),'gx')
plt.plot(b[2:149,1],np.zeros(147),'r+')

#get means of truncated lsist
session_avgs = list()
plt.subplots()
for cell in sessions:
    temp_session_avg = np.nanmean(np.array(sessions[cell].values()), axis = 0)
    session_avgs.append(temp_session_avg)
    plt.plot(temp_session_avg+(cell*0.1))
    
session_avgs_sorted = sorted(session_avgs, key=lambda x: x.argmax())
plt.subplots()
for i in range(len(session_avgs_sorted)):        
    plt.plot(session_avgs_sorted[i]+(i*0.1))


# ****************************************************************************************
#%% SESSION DATA - fast stim (003 - 2Hz)
#==============================================================================       
# get GRAY SESSION means and adjust lengths of sessions so equal
zeromin = min([len(grays[0][session]['data']) for session in grays[0]])              
gray_avgs = list()
for cell in grays:
    temp_gray_avg = np.nanmean(np.array([grays[cell][s]["data"][-zeromin:] for s in grays[cell]]), axis = 0)
    gray_avgs.append(temp_gray_avg)
#==============================================================================
  
#==============================================================================
# get STIM SESSION means of truncated lsist
session_avgs = list()
for cell in sessions:
    temp_session_avg = np.nanmean(np.array([sessions[cell][s]["data"] for s in sessions[cell]]), axis = 0)
    session_avgs.append(temp_session_avg)
#==============================================================================
#==============================================================================            
# Combine GRAY and STIM SESSION DATA 
session_gray_avgs = [np.concatenate((session_avgs[i], gray_avgs[i]), axis = 0) for i in range(len(session_avgs))]
#==============================================================================        
        
#%% PLOT - fast stim (003 - 2Hz) - Sessions
#==============================================================================        
# DATA VISUALIZATION

# Plot STIM SESSIONs and mean by cell of truncated lsist
session_cell = list()
for cell in sessions:
    temp_session = np.array([sessions[cell][s]["data"] for s in sessions[cell]])
    session_cell.append(temp_session)
  
# Plot each trace and mean (JESSE: need to work into subplot?)
nrows = 17
ncols = 4
fig, axes = plt.subplots(nrows, ncols, sharex='all', sharey='all')      
  
for i in range(nrows):
    for j in range(ncols):
        for k in range(len(session_cell[i*ncols+j])):
            axes[i,j].plot(session_cell[i*ncols+j][k])
            #temptitle = 'Cell '+str(i*ncols+j)
            #axes[i,j].title(temptitle)
            
        axes[i,j].plot(session_avgs[i*ncols+j], linewidth=2, color='k')

# Normalize session_avgs
# Can we create a new var/copy session_gray_avgs called session_gray_avgsNORM?
for avg in session_gray_avgs:
    avg -= min(avg)
    avg /= max(avg)
session_gray_avgs = np.array(session_gray_avgs).astype(np.float)        


# Create a list of origincal-ordered ROI#s 1 to x
cell_axis_labels = [x+1 for x in range(len(cells))]        

# Sort cells by latency to max intensity value        
sorted_both = sorted(session_gray_avgs, key=lambda x: x.argmax())
sorted_both_keys = np.argsort(np.argmax(session_gray_avgs, axis=1))
sorted_both_axis_labels = [x+1 for x in sorted_both_keys]

## Normalize session_avgs
#gray_avgs_norm = gray_avgs
#for avg in gray_avgs_norm:
#    avg -= min(avg)
#    avg /= max(avg)
#gray_avgs_norm = np.array(gray_avgs_norm).astype(np.float)
#
## Sort cells by latency to max intensity value
#sorted_avgs = sorted(session_avgs_norm, key=lambda x: x.argmax())
#sorted_avgs_keys = np.argsort(np.argmax(session_avgs_norm, axis=1))
#sorted_avgs_axis_labels = [x+1 for x in sorted_avgs_keys]
#
#sorted_gray = sorted(gray_avgs_norm, key=lambda x: x.argmax())
#sorted_gray_keys = np.argsort(np.argmax(gray_avgs_norm, axis=1))
#sorted_gray_axis_labels = [x+1 for x in sorted_gray_keys]
#
#
#both_norm = session_gray_avgs
#for avg in both_norm:
#    avg -= min(avg)
#    avg /= max(avg)
#both_norm = np.array(both_norm).astype(np.float)


# Generate trace plots
figFastStimTrace(session_gray_avgs,0.25)
for i in range(len(flip_i)):
    templine = flip_i[i][0]
    plt.axvline(x=templine, ymin=0, ymax = len(session_gray_avgs), linewidth=0.5, color='red')
    templine = flop_i[i][0]
    plt.axvline(x=templine, ymin=0, ymax = len(session_gray_avgs), linewidth=0.5, color='limegreen')
plt.title('All Cell Averages').set_family('helvetica')

figFastStimTrace(sorted_both,0.25)
for i in range(len(flip_i)):
    templine = flip_i[i][0]
    plt.axvline(x=templine, ymin=0, ymax = len(sorted_both), linewidth=0.5, color='red')
    templine = flop_i[i][0]
    plt.axvline(x=templine, ymin=0, ymax = len(sorted_both), linewidth=0.5, color='limegreen')
plt.title('All Cell Averages (sorted)').set_family('helvetica')        

# Generate heat plots        
#        get_colormaps(gray_avgs_norm,cell_axis_labels, 0, 1.0)
#        plt.title('All Cell GRAY Averages').set_family('helvetica')
#        plt.axis('tight')
#
#        get_colormaps(sorted_gray,sorted_gray_axis_labels, 0, 1.0)
#        plt.title('All Cell GRAY Averages').set_family('helvetica')
#        plt.axis('tight')

get_colormaps(session_gray_avgs,cell_axis_labels, 0, 1.0)
for i in range(len(flip_i)):
    templine = flip_i[i][0]
    plt.axvline(x=templine, ymin=0, ymax = len(session_avgs), linewidth=0.5, color='red')
    templine = flop_i[i][0]
    plt.axvline(x=templine, ymin=0, ymax = len(session_avgs), linewidth=0.5, color='limegreen')
plt.title('All Cell Averages').set_family('helvetica')
plt.axis('tight')        
      
get_colormaps(sorted_both,sorted_both_axis_labels, 0, 1.0)
for i in range(len(flip_i)):
    templine = flip_i[i][0]
    plt.axvline(x=templine, ymin=0, ymax = len(session_avgs), linewidth=0.5, color='red')
    templine = flop_i[i][0]
    plt.axvline(x=templine, ymin=0, ymax = len(session_avgs), linewidth=0.5, color='limegreen')
plt.title('All Cell Averages (sorted)').set_family('helvetica')
plt.axis('tight')
plt.savefig('All_Cell_Averages_(sorted_fast135).eps', format = 'eps', dpi = 1200)

        
        
#        Normalize to pre-stim gray
#        Divide cell0 session2 by the mean of the last 30-60 frames (~1-2s) of gray1
#        ...f0_cell0_session2 = grays[cell0][s1]["data"][-60:]
#        ...f0_cell0_session3 = grays[cell0][s2]["data"][-60:]
#==============================================================================        
        
#==============================================================================  
#%%
import process_function_jc as pf

t_0 = 0.2
t_1 = 0.75
t_2 = 3
samplingfreq = 30

#R_0 = pf.process_function(cells45[0], t_0, t_1, t_2, samplingfreq)
R_0, R0_EWMA, diffSignal_R0 = pf.process_function(cells45[0], t_0, t_1, t_2, samplingfreq)        
        
        
        
        
        
        