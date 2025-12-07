WORKFLOW:

1) run process_step_calcium_type1.py for "alldrift" (select CSV)

Use *_step_calcium_type*.py files as follows:

type1

Drifting gratings @ ~30fps
rows = len(a) # number of ROIs
cols = 8 # number of distinct stims
stimList = [0,45,90,135,180,225,270,315]
daqRecordingFreq = 1000.0 #sampling rate - Hz

#For classifying responsive cells with Mann-W U
stimwindow = 145 
pthresh = 0.05


type2

Phase-reversing gratings @ ~30fps