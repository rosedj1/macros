#######################################
##### READ NTUPLE AND MAKE HISTOS #####
#######################################

import ROOT as root
import numpy as np
import root_numpy
import matplotlib.pyplot as plt
from tdrstyle import *

## enter batch mode in root (so python can access displays)
root.gROOT.SetBatch(True)
    
def Plotter(x_array, y_array, x1_array=None, y1_array=None):
    plt.scatter(x_array, y_array, c='b')
    plt.scatter(x1_array, y1_array, c='r')
    plt.title('H-->ZdZ-->4l')
    plt.xlabel('eta')
    plt.ylabel('phi')
    #plt.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    plt.show()

## Return just one value: minDeltaRvalue of a single lepton
def deltaRMin(lepindex, eta_array, phi_array, showInfo=True):
    ## Convert to lists for easy concatenation
    eta_list = list(eta_array)
    phi_list = list(phi_array)
    
    if showInfo==True:
        print "Initial eta_list:", eta_list
        print "Initial phi_list:", phi_list
    
    ## Compare this lepton to its neighbors
    etaDiff_list = []
    phiDiff_list = []

    ## Get etaDiff_list
    for etaNeighbor in eta_list[:lepindex]+eta_list[(lepindex+1):]:
        etadiff = etaNeighbor - eta_list[lepindex]
        etaDiff_list.append(etadiff)
    if showInfo==True: print "etaDiff_list:\n", etaDiff_list

    ## Get phiDiff_list
    for phiNeighbor in phi_list[:lepindex]+phi_list[(lepindex+1):]:            
        phidiff = deltaPhi(phiNeighbor, phi_list[lepindex])
        phiDiff_list.append(phidiff)
    if showInfo==True: print "phiDiff_list:\n", phiDiff_list

    ## Calculate all deltaRvals for this lepton
    ## Choose min
    deltaR_array = np.sqrt( np.array(etaDiff_list)**2 + np.array(phiDiff_list)**2 )
    if showInfo==True: 
        print "deltaR_array for lep %i" % (lepindex), "\n", deltaR_array
        print "deltaR min value:", min(deltaR_array)
    return min(deltaR_array)

## User sets these parameters:
## LHE file
inputFile = "/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/gridpack/cmsgrid_final_eps1e-2_mzd20_lhaid306000.lhesyntax2.root"
inputTree = "lheEvents_tchan"
numEvents = 20
showInfo = False
bins = 50
## Signal ROOT fie
#inputFile = "/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180706/ZD_UpTo0j_MZD15_Eps1e-2_klo.root"
#inputTree = "Ana/passedEvents"

## read data from file into a numpy array
lep_eta_branch  = root_numpy.root2array(inputFile, inputTree, "eta3")
#lep_phi_branch  = root_numpy.root2array(inputFile, inputTree, "lep_phi")
lep_eta_array   = analyzeBranch("lep_eta", lep_eta_branch, -1,  4, showInfo=False) 
#lep_phi_array   = analyzeBranch("lep_phi", lep_phi_branch, -1,  4, showInfo=False) 

print "initial lep_eta_array/deltaRmin_array:\n", lep_eta_array

##################################
##### Create deltaRmin_array #####
##################################

## Initialize deltaRmin_array to have same dim as eta_array and phi_array
deltaRmin_array = lep_eta_array

## Fill up deltaRmin_array
numLeptons = np.shape(lep_eta_array)[1]
#numEvents = np.shape(lep_eta_array)[0]
numEvents = 1
eventCount = 0

for event in range(numEvents):
    print "Event", eventCount
    eventCount += 1
    for lepindex in range(numLeptons):
         deltaRmin_array[event][lepindex] = deltaRMin(lepindex, lep_eta_array[event], lep_phi_array[event], showInfo=True)
print deltaRmin_array
