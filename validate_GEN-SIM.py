import ROOT
from DataFormats.FWLite import Events, Handle
import numpy as np
#import matplotlib.pyplot as plt

events = Events('zd0j_MZd30_eps1e-3_LHE-GEN-SIM_57.root')
handle = Handle('std::vector<reco::GenParticle>')
label = ("genParticles","","SIM") 
numProcess = 0

### Each item in this for loop is an entire gg-->H-->ZZd-->4lep process ###
for count,e in enumerate(events, 1):
# Only look at first "count" events
    if count > 100: continue
    e.getByLabel(label,handle)
    genParticles = handle.product()
    numProcess+=1
### Each item in this for loop is a single particle in the process above ###
    for p in genParticles:
        if not p.isHardProcess(): continue
        #print p.pdgId(),p.pt(),p.eta()
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
print "\n", "Number of events analyzed:", numProcess  
print "Total number of isHardProcess events in root file:", count, "\n"
