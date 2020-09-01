from ROOT import *

f = TFile("/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180706/ZD_UpTo0j_MZD20_Eps1e-2_klo.root","READ")
t = f.Get("Ana/passedEvents")

n_4e = t.GetEntries("Sum$(abs(GENlep_id[])==11)==4")

n_acc_4e = t.GetEntries("Sum$(abs(GENlep_id[])==11)==4 && passedFiducialSelection==1")
n_acceff_4e = t.GetEntries("Sum$(abs(GENlep_id[])==11)==4 && passedFullSelection==1")

print "acc. 4e:",float(n_acc_4e)/float(n_4e),"acc*eff 4e:",float(n_acceff_4e)/float(n_4e)

h = {}

h["massZ2"] = TH1D("h_massZ2","h_massZ2",48,4,52)
h["massZ2"].Sumw2()

t.Draw("massZ2>>h_massZ2","passedFullSelection==1","goff")

from tdrStyle import *
setTDRStyle()

c1 = TCanvas("c1","c1",800,800)
c1.cd()
h["massZ2"].SetLineColor(2)
h["massZ2"].Draw("hist")

c1.SaveAs("test.pdf")
