from ROOT import TFile
import sys
import numpy as np
import matplotlib.pyplot as plt

pT_min = 5.0
pT_max = 7.0
eta_min = 0.0
eta_max = 2.0

infile = sys.argv[1]

f = TFile.Open(infile)
t = f.Get("passedEvents")

total_evts = t.GetEntries()
qd0_ls = [] 
qd0_bins = np.linspace(-0.01, 0.01, 200)

for evt in range(total_evts):
    t.GetEntry(evt)

    # Make q*d0 dist.
    pT1 = getattr(t, "pT1")
    pT2 = getattr(t, "pT2")
    eta1 = getattr(t, "eta1")
    eta2 = getattr(t, "eta2")

    def calc_qd0(Id_str, d0BS_str):
        Id = getattr(t, Id_str)
        d0BS = getattr(t, d0BS_str)
        qd0 = d0BS * float(Id) / -13.
        return qd0

    if (pT1 > pT_min) and (pT1 < pT_max):
        if (abs(eta1) > eta_min) and (abs(eta1) < eta_max):
            # Passed cuts.
            qd0 = calc_qd0("Id1", "d0BS1")
            qd0_ls.append(qd0)

    if (pT2 > pT_min) and (pT2 < pT_max):
        if (abs(eta2) > eta_min) and (abs(eta2) < eta_max):
            # Passed cuts.
            qd0 = calc_qd0("Id2", "d0BS2")
            qd0_ls.append(qd0)

print "Number of leptons in low pT region, in central barrel:\n{}".format(len(qd0_ls))
#fig, ax = plt.subplots()
#ax.hist(qd0_ls, bins=qd0_bins)
#plt.show()
#plt.savefig("/home/rosedj1/HiggsMeasurement/CMSSW_8_0_32/src/jake_root_skims/test1.pdf")
