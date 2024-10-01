#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 13:40:50 2024

@author: jcoleman
"""

from matplotlib import pyplot as plt
import pandas as pd


#Calcium imaging data from Jason Coleman
# 

datadirectory = '/Users/jcoleman/Documents/--LARGE DATA--/#Febo/thy1gcamp6s data set/'
datafile = 'mThy6s2_alldrift_D4_001Z1_ROIdata.csv'

df = pd.read_csv(datadirectory + datafile)
# df = pd.DataFrame(data=data, index=row_labels)

roimap_png = 'MAP mThy6s2_alldrift_D4_001Z1_ROIdata.png'
roimap_eps = 'MAP mThy6s2_alldrift_D4_001Z1_ROIdata.eps'

csv_timestamps_frames = {'startDelay': [0,1800], 'stim1': [1801,1950], 'gray1': [1951,2160]}
csv_timestamps_sec = {'startDelay': [0,60], 'stim1': [61,65], 'gray1': [66,72]} #... 'stim2', 'gray2', etc...
# need to def() into tuples for all stims/grays
#example code plot - raw data
# dont know exact stim, just that each epoch was different, although likely the same order as below
# Orientation of stims (relative to horizon = 0deg)
# stim1=0, stim2=45, stim3=90,stim4=135, stim5=180, stim6=225,stim7=270, stim8=315
# gray (gray1 throuhg gray8; gray = 50% contrast)

"""
STORED OPTIONS, FILE = mThy6s_001_D4Z1
Time of Recording (recording stopped)
Tue Sep 26 15:12:21 2017
Session Name
mThy6s_001_D4Z1
Recording Frequency
1000
Total Channels: 8
Channel 1 Name: 
Channel 2 Name: 
Channel 3 Name: 
Channel 4 Name: 
Channel 5 Name: 
Channel 6 Name: 
Channel 7 Name: 
Channel 8 Name: 
Experimental Notes
8days Lt4/2p

0,45,90,135,180,225,270,315

duration = 2.5 (true = 5s)
inter-session = 3.5 (true = 7s)
start delay = 30 (true = 60sec)
16200 frames (9m7s on bscope streaming)
drift rate = 3Hz? (true = 1.5?); try 6Hz, true = 3Hz?

mon width  = 52 cm
distance = 14 cm
hoiz = 1920; vert = 1080 pix

5s, 7s gray
0,45,90,135,180,225,270,315
5 sessions
60s delay
= (16200 frames = 9m7s)
"""