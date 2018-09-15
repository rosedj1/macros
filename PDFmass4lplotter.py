## PURPOSE: Plot overlaid histograms of, e.g., mass4l, pT4l, etc., from 
##          lhe files (NTuples) which have been processed by ZZD_lhe.C
## SYNTAX:  python <this_exe>
## AUTHOR:  Jake Rosenzweig
## DATE:    2018-09-02
## NOTES:   Make sure to point to the correct NTuples!

from ROOT import *
import numpy as np
from tdrStyle import *
#import root_numpy
#import matplotlib.pyplot as plt
#import pandas as pd
#import seaborn as sns

## enter batch mode in root (so python can access displays)
#root.gROOT.SetBatch(True)
ROOT.gROOT.IsBatch()

## Parameters
bins = 50      #mass4l=100 // pT4l=50
xmin = 0.0      #mass4l=124.9 // pT4l=0.0
xmax = 2E-9    #mass4l=125.1 // pT4l=2E-9
pdf1 = "262000"
pdf2 = "263000"
## Branch in NTuple you want to plot
## Try something like: "mass4l" or "pT4l"
kinem = "pT4l"

f1 = TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/gridpack/cmsgrid_final_eps1e-2_mzd20_lhaid"+pdf1+".lhesyntax2.root","READ")
f2 = TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/gridpack/cmsgrid_final_eps1e-2_mzd20_lhaid"+pdf2+".lhesyntax2.root","READ")
t1 = f1.Get("lheEvents_tchan")
t2 = f2.Get("lheEvents_tchan")

binwidth = (xmax-xmin)/double(bins)

## Create histograms
h1 = TH1D("h1_%s" % kinem,"h1_%s" % kinem, bins, xmin, xmax)
h1.Sumw2()                              
h2 = TH1D("h2_%s" % kinem,"h2_%s" % kinem, bins, xmin, xmax)
h2.Sumw2()                              

#t.Draw("mass4l>>h_mass4l","passedFullSelection==1","goff")
t1.Draw("%s>>h1_%s" % (kinem,kinem),"","goff")
t2.Draw("%s>>h2_%s" % (kinem,kinem),"","goff")

from tdrStyle import *                                    
setTDRStyle()                                             
                                                          
c1 = TCanvas("c1","c1",1000,1000)                           
c1.cd()                                                   

#h2.SetTitle("H-->ZdZ-->4l") #Title not working
h1.Scale(1./h1.Integral())
h2.Scale(1./h2.Integral())

## For setting the range on axes
## SetAxisRange creates weird spike behavior
h1.SetAxisRange(0.0, 0.1, "Y")
h2.SetAxisRange(0.0, 0.1, "Y")

## Try GetXaxis instead if you need to set range
#h1.GetXaxis().SetRangeUser(xmin, xmax)
#h2.GetXaxis().SetRangeUser(xmin, xmax)
h1.SetLineColor(1)                               
h2.SetLineColor(2)                               
h1.SetXTitle("%s (GeV)" % kinem)

h1.SetLineWidth(2)
h2.SetLineWidth(2)

## For LO pdfs, pT4l is nearly 0 GeV
## So use eV scale
if binwidth <= 1E-5:#GeV
    h1.SetYTitle("Fraction of Events / %4.4f (eV)" % (binwidth * 1E9))
else:
    h1.SetYTitle("Fraction of Events / %4.4f (GeV)" % binwidth)

h1.SetLabelSize(0.025, "X")
h1.SetLabelSize(0.03, "Y")
h2.SetLabelSize(0.025, "X")
h2.SetLabelSize(0.03, "Y")

h1.Draw("hist")                                  
h2.Draw("hist same")                                  

legend = TLegend(0.60,0.7,0.8,0.9)
legend.AddEntry(h1,"lhaid = %s" % pdf1,"lpf")
legend.AddEntry(h2,"lhaid = %s" % pdf2,"lpf")
legend.SetLineWidth(3)
legend.SetBorderSize(0)
legend.SetTextSize(0.03)
legend.Draw("same")
## Wait for user input to close the canvas when ready
wait = raw_input()
