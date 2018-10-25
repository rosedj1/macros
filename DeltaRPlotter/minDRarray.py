import numpy as np
from lhetoarray import lheToArray
from dRCalc import *

def minDRarray(inputFile,inputTree,numEvents,PRINT):
    eta_arr = lheToArray(inputFile,inputTree,"eta",numEvents,PRINT)
    phi_arr = lheToArray(inputFile,inputTree,"phi",numEvents,PRINT)

    if PRINT == 1:
        print "eta array1\n:",eta_arr,"\n",eta_arr.shape
        print "phi array2\n:",phi_arr,"\n",phi_arr.shape

    ## Initialize empty dR arrays
    if numEvents == -1:
        numEvents = eta_arr.shape[0]

    dRarr = np.zeros((numEvents,6))
    min_dRarr = np.zeros(numEvents)
    ## For each event...
    for event in range(dRarr.shape[0]):
        dRarr[event][0] = dRCalc(eta_arr[event][0], phi_arr[event][0], eta_arr[event][1], phi_arr[event][1])
        dRarr[event][1] = dRCalc(eta_arr[event][0], phi_arr[event][0], eta_arr[event][2], phi_arr[event][2])
        dRarr[event][2] = dRCalc(eta_arr[event][0], phi_arr[event][0], eta_arr[event][3], phi_arr[event][3])
        dRarr[event][3] = dRCalc(eta_arr[event][1], phi_arr[event][1], eta_arr[event][2], phi_arr[event][2])
        dRarr[event][4] = dRCalc(eta_arr[event][1], phi_arr[event][1], eta_arr[event][3], phi_arr[event][3])
        dRarr[event][5] = dRCalc(eta_arr[event][2], phi_arr[event][2], eta_arr[event][3], phi_arr[event][3])

    ## Fill min_dRarr
    for event in range(dRarr.shape[0]):
        min_dRarr[event] = min(dRarr[event])

    if PRINT == 1: 
        print "dRarr is:\n",dRarr,"\n",dRarr.shape
        print "min_dRarr is:\n",min_dRarr,"\n",min_dRarr.shape

    return min_dRarr
