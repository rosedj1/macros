import numpy as np
import root_numpy
import ROOT

from lookatbranch import lookAtBranch
from lhetoarray import lheToArray
from dRCalc import dRCalc, dPhiCalc
## User sets these parameters:                                                                               
## Possible kinemVar names in LHE file:
    #weight
    #id[3-6]
    #eta[3-6]
    #pt[3-6]
    #pto[1-4]   This is the important pT! o = ordered
    #mass4l
    #pT4l
    #eta4l
    #rapidity4l
    #massZ1
    #massZ2
    #costheta1
    #costheta2
    #costhetastar
    #phi
    #phi1
inputFile = "/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/gridpack/cmsgrid_final_eps1e-2_mzd20_lhaid306000.lhesyntax2.root" 
inputTree = "lheEvents_tchan"                                                                                
kinemVar1 = "eta"
kinemVar2 = "phi"

## Use -1 to print ALL events
numEvents = -1                                                                                               
showInfo = False                                                                                             
lep_eta_arr = lheToArray(inputFile,inputTree,"eta",numEvents,showInfo)
#lep_phi_arr = lheToArray(inputFile,inputTree,"phi",numEvents,showInfo)

print lep_eta_arr
#print lep_phi_arr
