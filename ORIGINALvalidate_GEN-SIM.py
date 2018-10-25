import ROOT
from DataFormats.FWLite import Events, Handle

# Change this line
events = Events('step1_GEN-SIM.root')			# Events takes .root file.
# and these
handle = Handle('std::vector<reco::GenParticle>')	# Handle takes objects.
label = ("genParticles","","GEN") 
for count,e in enumerate(events):
    if count > 10: continue			# Turn this to while loop?
    e.getByLabel(label,handle)
    genParticles = handle.product()		# No longer genParticles
    print "-"*20 
    for p in genParticles:
        if not p.isHardProcess(): continue
        print p.pdgId(),p.pt(),p.eta()
