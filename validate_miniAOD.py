##################################################################################
##### Purpose:
## Analyze a miniAOD file by following Lucien's ORIGINALvalidate_GEN-SIM.py script.
## Figure out how events and handles work.
## Eventually find the ZZ-->4l events and analyze them.
## Specifically, plot the histogram of the minDeltaR's of the leps in each event.
##################################################################################

import ROOT
from DataFormats.FWLite import Events, Handle
import numpy as np
import sys
#import matplotlib.pyplot as plt

PRINT = 1

events = Events('HIG-RunIIFall17MiniAOD-00079_99.root')
handle = Handle('std::vector<pat::Electron>')
label = ("slimmedElectrons","","PAT") # a tuple

numProcess = 0
maxEvents = 20
if PRINT == 1:
    print "events looks like:\n", events
    print "handle looks like:\n", handle

### Each item in this for loop is an entire gg-->H-->ZZd-->4lep process ###
while numProcess < maxEvents:
    for count,e in enumerate(events, 1):
    # Only look at first "count" events
        print "e.getByLabel:", e.getByLabel(label,handle)
        print "handle.product():", handle.product()
        #patParticles = handle.product()
        numProcess+=1
sys.exit("End of program.")
#print "The first %i events in events are:\n", 
### Each item in this for loop is a single particle in the process above ###
for p in patParticles:
    #if not p.isHardProcess(): continue
    print p.pdgId(),p.pt(),p.eta()
# Potentially interesting quantities:
# charge(), p4(), p(), energy(), et?, et2?,
# mass(), massSqr(), pt(), phi(), theta(), eta()
# rapidity(), vertex(), status(), longLived()
# isElectron(), isMuon(), is StandAloneMuon()
# isGlobalMuon(), isTrackerMuon(), isCaloMuon()
# is Photon(), isConvertedPhoton(), isJet()
if p.pdgId() == 1023: 
        #print "-"*20 
        print "Process %d: %f" % (count, p.mass() )
print "\n", "Number of events analyzed:", numEvents  
#print "Total number of isHardProcess events in root file:", count, "\n"
