## PURPOSE: Plot overlaid histograms of MadGraph5 vs. JHUGen kinematics
##          using lhe files (NTuples) which have been processed by ZZD_lhe.C
## SYNTAX:  python <script.py>
## AUTHOR:  Jake Rosenzweig
## DATE:    2018-10-03
## NOTES:   Make sure to point to the correct NTuples!

import numpy as np
import os
from shutil import copy2
from ROOT import *
from tdrStyle import *

## enter batch mode in root (so python can access displays)
ROOT.gROOT.SetBatch(kTRUE) #kTRUE = will NOT draw plots!
setTDRStyle()                                             

##### PARAMETERS
cuts = "WITH"
outDirectory = "/home/rosedj1/public_html/DarkZ/MadGraphvsJHUGenPlotsTEST_"+cuts+"cuts/"
## Branch in NTuple you want to plot
## Try something like: "mass4l" or "pT4l"
#kinemlist = [   pto1,pto2,pto3,pto4,
#                eta3,eta4,eta5,eta6,
#                phi3,phi4,phi5,phi6,
#                massZ1,massZ2,mass4l,pT4l
#                costheta1,costheta2,costhetastar,
#                phi,phi1 ]
#kinemlist = ["mass4l","eta3","eta4","eta5","eta6","phi3","phi4","phi5","phi6","pto1","pto2","pto3","pto4","massZ1","massZ2","pT4l","costheta1","costheta2","costhetastar","phi","phi1"]
kinemlist = ["mass4l","eta3","phi3","pto1","massZ1","massZ2","costheta1","costheta2","costhetastar","phi","phi1","pT4l"]
#kinemlist = ["mass4l"]
#masslist = [1,2,3,4,7,10,15,20,25,30,35]
#masslist = [1,7,20,35]
masslist = [4,25,30,35]

if not os.path.exists(outDirectory): 
    os.makedirs(outDirectory)
copy2('/home/rosedj1/index.php',outDirectory)

for mass in masslist:
    ##### User Parameters
    ## Convert masses to strings for naming purposes
    mass     = str(mass)
    mass_int = int(mass)

    ## The file below has great kinematics!
    #f1 = TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/ggHZZd4l_fromJHUGen_MZd30_15000events_test20181018_2025.lhe_15000events_WITHcuts.root","READ")
    #f1 = TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/new_ggHZZd4l_MZd"+mass+".lhe_15000events_"+cuts+"cuts.root","READ")
    f1 = TFile.Open("root://cmsio5.rc.ufl.edu//store/user/drosenzw/JHUGenStudies/ggHZZd4l_fromJHUGen_MZd"+mass+"_15000events.lhe_"+cuts+"cuts.root","READ")
    f2 = TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/gridpack/gridpack_HAHM_variablesw_v3_eps1e-2_MZd"+mass+"_lhaid306000/cmsgrid_final.lheMZd"+mass+"_15000events_lhaid306000"+cuts+"cuts.root","READ")
    t1 = f1.Get("lheEvents_tchan")
    t2 = f2.Get("lheEvents_tchan")

    for kinem in kinemlist:
        ## Setting plotting parameters
        if kinem in ["pto1","pto2","pto3","pto4"]:
            xmin = 0.0
            xmax = 120.0
            bins = 120
            binwidth = (xmax-xmin)/float(bins)
            xlabel = "p_{T} (GeV)"
            ylabel = str("%.5f" % binwidth)+" GeV"
            xminrange = xmin - 1*binwidth
            xmaxrange = xmax + 1*binwidth
        elif kinem in ["eta3","eta4","eta5","eta6"]:
            xmin = -2.8
            xmax = 2.8
            bins = 112
            binwidth = (xmax-xmin)/float(bins)
            xlabel = "#eta"
            ylabel = str("%.5f" % binwidth)
            xminrange = xmin
            xmaxrange = xmax
        elif kinem in ["phi3","phi4","phi5","phi6","phi","phi1"]:
            xmin = -np.pi
            xmax = np.pi
            bins = 40
            binwidth = (xmax-xmin)/float(bins)
            xlabel = kinem + " [radians]"
            ylabel = str("%.5f" % binwidth)+" radians"
            xminrange = -3.5
            xmaxrange = 3.5
        elif kinem in ["costheta1","costheta2","costhetastar"]:
            xmin = -1.0
            xmax = 1.0
            bins = 40
            binwidth = (xmax-xmin)/float(bins)
            xlabel = kinem
            ylabel = str("%.5f" % binwidth)
            xminrange = -1.2
            xmaxrange = 1.2
        elif kinem in ["massZ1"]:
            xmin = 75.0
            xmax = 100.0
            bins = 200
            binwidth = (xmax-xmin)/float(bins)
            xlabel = "m_{Z1} (GeV)"
            ylabel = str("%.5f" % binwidth)+" GeV"
            #xminrange = 65.0
            #xmaxrange = 120.0
            xminrange = xmin - 50*binwidth
            xmaxrange = xmax + 50*binwidth
        elif kinem in ["massZ2"]:
            #xmin = float(mass_int)-0.001
            #xmax = float(mass_int)+0.001
            xmin = float(mass_int)-1.0
            xmax = float(mass_int)+1.0
            bins = 100
            binwidth = (xmax-xmin)/float(bins)
            xlabel =  "m_{Z2} (GeV)"
            ylabel = str("%.5f" % binwidth)+" GeV"
            xminrange = xmin
            xmaxrange = xmax
        elif kinem in ["mass4l"]:
            xmin = 124.8
            xmax = 125.2
            bins = 500
            binwidth = (xmax-xmin)/float(bins)
            xlabel = "m_{4l} (GeV)"
            ylabel = str("%.5f" % binwidth)+" GeV"
            xminrange = xmin
            xmaxrange = xmax
        elif kinem in ["pT4l"]:
            xmin = 0.0
            xmax = 150.0
            bins = 150
            binwidth = (xmax-xmin)/float(bins)
            xlabel = "p_{T,4l} (GeV)"
            ylabel = str("%.5f" % binwidth)+" GeV"
            xminrange = xmin
            xmaxrange = xmax

        ## Create histograms
        h1 = TH1D("h1_mZd%s_%s" % (mass,kinem),"h1_mZd%s_%s" % (mass,kinem), int(bins), xmin, xmax)
        h1.Sumw2()                              
        h2 = TH1D("h2_mZd%s_%s" % (mass,kinem),"h2_mZd%s_%s" % (mass,kinem), int(bins), xmin, xmax)
        h2.Sumw2()                              

        c1 = TCanvas("c_mZd%s_%s" % (mass,kinem),"c_mZd%s_%s" % (mass,kinem),800,800)
        c1.cd()                                                   

        ## Fill histograms
        #t.Draw("mass4l>>h_mass4l","passedFullSelection==1","goff")
        t1.Draw("%s>>h1_mZd%s_%s" % (kinem,mass,kinem), "", "goff")
        t2.Draw("%s>>h2_mZd%s_%s" % (kinem,mass,kinem), "", "goff")

        ## Normalize the histograms
        if h1.Integral() != 0: 
            #print("Scaling "+kinem)
            h1.Scale(1./h1.Integral())
        else: print("WARNING!: "+kinem+" for mass "+mass+" has Integral() = 0!")
        if h2.Integral() != 0: 
            #print("Scaling "+kinem)
            h2.Scale(1./h2.Integral())
        else: print("WARNING!: "+kinem+" for mass "+mass+" has Integral() = 0!")

        ## Try GetXaxis instead if you need to set range
        #h1.GetXaxis().SetRangeUser(xmin, xmax)
        #h2.GetXaxis().SetRangeUser(xmin, xmax)

        h1.SetLineWidth(3)
        h2.SetLineWidth(2)
        h1.SetLineColor(1)                               
        h2.SetLineColor(2)                               

        #h1.SetTitle("MadGraph5 vs. JHUGen, %s" % kinem)
        h1.SetTitle("H#rightarrowZZd#rightarrow4l, mZd = %d GeV" % mass_int)
        h1.SetXTitle("%s" % xlabel)
        #h1.SetTitleSize(0.03)
        #h1.SetTitleSize(0.03, "XY")

        h1.SetYTitle("Fraction of Events / %s" % ylabel )
        h1.SetLabelSize(0.025, "XY")

        h1.SetAxisRange(xminrange, xmaxrange, "X")
        h2.SetAxisRange(xminrange, xmaxrange, "X")
        #h1.SetAxisRange(yminrange, ymaxrange, "Y")
        #h2.SetAxisRange(yminrange, ymaxrange, "Y")

        h1.Draw("hist goff")                                  
        h2.Draw("hist same goff")                                  

        ## Make sure legend doesn't block histogram lines
        if kinem in ["phi","phi1","phi3","phi4","phi5","phi6","costheta1","costheta2","costhetastar"]:
            legend = TLegend(0.75,0.20,0.95,0.40)
        else:
            legend = TLegend(0.75,0.7,0.95,0.9)
        legend.AddEntry(h1, "JHUGen", "lpf")
        legend.AddEntry(h2, "MadGraph5", "lpf")
        legend.SetLineWidth(2)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.025)
        legend.Draw("same goff")

        ## For plot-ordering purposes on tier2 website, change 01 --> 01, 2-->02, etc.
        if len(mass) == 1: mass = "0"+mass
        outputFile = outDirectory + kinem+"_mZd"+mass+".png"
        c1.SaveAs( outputFile )
    ## Wait for user input to close the canvas when ready
    #wait = raw_input()
    #//#h1.GetXaxis().SetRangeUser(xmin, xmax)
    #//#h2.GetXaxis().SetRangeUser(xmin, xmax)
    #//h6.SetLineColor(kBlue)                               
    #//h7.SetLineColor(kCyan+1)                               
    #//
    #//h1.SetLineWidth(3)
    #//
    #//h1.SetYTitle("Events")
    #//
    #//h1.SetLabelSize(0.025, "X")
    #//h1.SetLabelSize(0.025, "Y")
    #//#h1.SetLabelOffset(0.003, "X")
    #//#h1.SetLabelOffset(0.003, "Y")
    #//
    #//h1.Draw("hist")                                  
    #//#c1.SetLogx()
    #//#c1.SetLogy()
