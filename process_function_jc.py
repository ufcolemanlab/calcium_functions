import math
import EWMA
import numpy as np

def runningMeanFast(x, N):
    """
    The running mean is equivalent to convolving x with a vector that is N long, with all members equal to 1/N.
    The numpy implementation of convolve includes the starting transient, so you have to remove the first N-1 points.
    
    Note that convolve does include a 'same' mode which seems like it should address the starting transient issue,
    but it splits it between the beginning and end.
    """
    return np.convolve(x, np.ones((N,))/N)[(N-1):]

def smooth(signal, span):
    
    final = []

    l = len(signal)
    neighbors = span // 2

    for i in range(l):
        sum = 0
        if (i < neighbors):
            sp = i

        elif ((neighbors+i) >  l-1):
            sp = (l-1)-i

        else:
            sp = neighbors

        for j in range(int(i-sp), int(i+sp+1)):
            sum+=signal[j]

        final.append(sum/((sp*2)+1))

    return final

def process_function(signalin, t_0, t_1, t_2, Fs):
    """
    Ported from original MATLAB function to Python by Z. Royston (Coleman lab):
    Implementation of Nature protocol
    Hongbo Jia, Nathalie L Rochefort1, Xiaowei Chen & Arthur Konnerth1 "In
    vivo two-photon imaging of sensory-evoked dendritic calcium signals in cortical neurons"

    Implementation copyright Petros Xanthopoulos 2013-2014
    usage: signalout=process_function(signalin,t_0,t_1,t_2)
    where input: signalin is the raw signal 
    t_0,t_1,t_2 are the parameters described in Nature protocol paper
    comments: for a 30Hz (Fs) imaging systems the following parameter setup is
    recommended (empirical note on Nature paper): 
    Fs = 30
    t_0= 0.2;
    t_1=0.75;
    t_2=3;
    
    9/24/16 - validated by JEC
    """

    F_0 = []

    t_0_s = math.floor(t_0 * Fs)
    t_1_s = math.floor(t_1 * Fs)
    t_2_s = int(math.floor(t_2 * Fs))

    F_sm = smooth(signalin, t_1_s)

    for i in range((t_2_s), len(signalin)):
        F_0.append(min(F_sm[i-t_2_s:i]))

    R_0 = np.divide((signalin[t_2_s:] - F_0), F_0)
    # R_0 is reduced by 90 elements, zeroed and magnitude reduced by ~10-fold
    R_0_sm = np.divide((signalin[t_2_s:] - F_0), F_0)

    R_0_sm = EWMA.smoothed_z(R_0_sm, t_0_s)
    
    diffSignal_R0 = len(signalin) - len(R_0)

    return R_0, R_0_sm, diffSignal_R0

# smooth() test
#F = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
#G = smooth(F, 3)
#print(G)
