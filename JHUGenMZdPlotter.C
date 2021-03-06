#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "tdrstyle.C"

using namespace std;

void JHUGenMZdPlotter() {

setTDRStyle();

//string outputFile = "";
//string branch = "massZ2";

TCanvas * c1 = new TCanvas("c1","c1",900,900);

TFile* f1 = new TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/ggHZZd4l_MZd1.lhe_2000events_WITHcuts.root");
TTree* t1 = (TTree*)f1->Get("lheEvents_tchan");
//t1->SetFillColor(kBlack);
t1->SetFillColor(kBlack);
t1->Draw("massZ2>>h1(120,0.0,40.0)","");

TFile* f2 = new TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/ggHZZd4l_MZd2.lhe_2000events_WITHcuts.root");
TTree* t2 = (TTree*)f2->Get("lheEvents_tchan");
t2->SetFillColor(kBlue);
t2->Draw("massZ2>>h2(120,0.0,40.0)","","same");

TFile* f3 = new TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/ggHZZd4l_MZd3.lhe_2000events_WITHcuts.root");
TTree* t3 = (TTree*)f3->Get("lheEvents_tchan");
t3->SetFillColor(kGreen);
t3->Draw("massZ2>>h3(120,0.0,40.0)","","same");

TFile* f4 = new TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/ggHZZd4l_MZd4.lhe_2000events_WITHcuts.root");
TTree* t4 = (TTree*)f4->Get("lheEvents_tchan");
t4->SetFillColor(kViolet);
t4->Draw("massZ2>>h4(120,0.0,40.0)","","same");

TFile* f7 = new TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/ggHZZd4l_MZd7.lhe_2000events_WITHcuts.root");
TTree* t7 = (TTree*)f7->Get("lheEvents_tchan");
t7->SetFillColor(kRed);
t7->Draw("massZ2>>h7(120,0.0,40.0)","","same");

TFile* f10 = new TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/ggHZZd4l_MZd10.lhe_2000events_WITHcuts.root");
TTree* t10 = (TTree*)f10->Get("lheEvents_tchan");
t10->SetFillColor(kRed);
t10->Draw("massZ2>>h10(120,0.0,40.0)","","same");

TFile* f15 = new TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/ggHZZd4l_MZd15.lhe_2000events_WITHcuts.root");
TTree* t15 = (TTree*)f15->Get("lheEvents_tchan");
t15->SetFillColor(kRed);
t15->Draw("massZ2>>h15(120,0.0,40.0)","","same");

TFile* f20 = new TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/ggHZZd4l_MZd20.lhe_2000events_WITHcuts.root");
TTree* t20 = (TTree*)f20->Get("lheEvents_tchan");
t20->SetFillColor(kRed);
t20->Draw("massZ2>>h20(120,0.0,40.0)","","same");

TFile* f25 = new TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/ggHZZd4l_MZd25.lhe_2000events_WITHcuts.root");
TTree* t25 = (TTree*)f25->Get("lheEvents_tchan");
t25->SetFillColor(kRed);
t25->Draw("massZ2>>h25(120,0.0,40.0)","","same");

TFile* f30 = new TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/ggHZZd4l_MZd30.lhe_2000events_WITHcuts.root");
TTree* t30 = (TTree*)f30->Get("lheEvents_tchan");
t30->SetFillColor(kRed);
t30->Draw("massZ2>>h30(120,0.0,40.0)","","same");

TFile* f35 = new TFile("/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/genproductions/bin/JHUGen/ggHZZd4l_MZd35.lhe_2000events_WITHcuts.root");
TTree* t35 = (TTree*)f35->Get("lheEvents_tchan");
t35->SetFillColor(kRed);
t35->Draw("massZ2>>h35(120,0.0,40.0)","","same");

TLegend * leg = new TLegend(0.7, 0.7, 0.9, 0.9, "Zd mass");
//leg->SetFillColor(0);
leg->AddEntry(t1, "1 GeV");
leg->AddEntry(t2, "2 GeV");
leg->AddEntry(t3, "3 GeV");
leg->AddEntry(t4, "4 GeV");
leg->Draw();

// JAKE MAKE SURE YOU CHANGE THE PLOT NAME!
//c1->SaveAs("/home/rosedj1/public_html/DarkZ/JHUGenkinematics/MadGraphvsJHUGen.png");

//TFile* f2 = new TFile("cmsgrid_final_eps1e-3.lhesyntax2.root");
//TTree* t2 = (TTree*)f2->Get("lheEvents_tchan");
//t2->SetLineColor(kRed);
//t2->Draw("massZ2>>h2(90,19.8,20.2)","","same");
//
//TFile* f3 = new TFile("cmsgrid_final_eps1e-4.lhesyntax2.root");
//TTree* t3 = (TTree*)f3->Get("lheEvents_tchan");
//t3->SetLineColor(kGreen);
//t3->Draw("massZ2>>h3(90,19.8,20.2)","","same");
//
//TFile* f4 = new TFile("cmsgrid_final_eps1p5e-2.lhesyntax2.root");
//TTree* t4 = (TTree*)f4->Get("lheEvents_tchan");
//t4->SetLineColor(kBlue);
//t4->Draw("massZ2>>h4(90,19.8,20.2)","","same");
//
//TFile* f5 = new TFile("cmsgrid_final_eps5e-3.lhesyntax2.root");
//TTree* t5 = (TTree*)f5->Get("lheEvents_tchan");
//t5->SetLineColor(kViolet);
//t5->Draw("massZ2>>h5(90,19.8,20.2)","","same");


// Build a legend
//TLegend leg(0.1, 0.7, 0.3, 0.9, "Epsilon Values");
//leg.SetFillColor(0);
//leg.AddEntry(t1, "1e-2");
//leg.Draw();
}
//import numpy as np
//import root_numpy
//from ROOT import *
//
//from lookatbranch import lookAtBranch
//from lhetoarray import lheToArray
//from dRCalc import dRCalc, dPhiCalc
//from minDRarray import minDRarray
//
//import matplotlib.pyplot as plt
//from rootpy.plotting import histogram
//from tdrStyle import *
//
//## Enter batch mode in ROOT (so python can access displays)
//ROOT.gROOT.IsBatch()
//setTDRStyle()                                             
//## User sets these parameters:                                                                               
//## Possible kinemVar names in LHE file:
//    #weight
//    #id[3-6]
//    #eta[3-6]
//    #phi[3-6]
//    #pt[3-6]
//    #pto[1-4]   This is the important pT! o = ordered
//    #mass4l
//    #pT4l
//    #eta4l
//    #rapidity4l
//    #massZ1
//    #massZ2
//    #costheta1
//    #costheta2
//    #costhetastar
//    #phi
//    #phi1
//
//inputDir = "/home/rosedj1/DarkZ-EvtGeneration/CMSSW_9_4_2/src/DarkZ-EvtGeneration/gridpack/"
//outputFile = "/home/rosedj1/public_html/DarkZ/LHEplotsdiffMZdmasses/LHEplot200bins.png"
//
//inputFile1 = inputDir+"gridpack_HAHM_variablesw_v3_eps1e-2_MZd1_lhaid306000/cmsgrid_final.lheMZd1_15000events_lhaid306000WITHcuts.root" 
//inputFile2 = inputDir+"gridpack_HAHM_variablesw_v3_eps1e-2_MZd2_lhaid306000/cmsgrid_final.lheMZd2_15000events_lhaid306000WITHcuts.root" 
//inputFile3 = inputDir+"gridpack_HAHM_variablesw_v3_eps1e-2_MZd3_lhaid306000/cmsgrid_final.lheMZd3_15000events_lhaid306000WITHcuts.root" 
//inputFile4 = inputDir+"gridpack_HAHM_variablesw_v3_eps1e-2_MZd4_lhaid306000/cmsgrid_final.lheMZd4_15000events_lhaid306000WITHcuts.root" 
//inputFile5 = inputDir+"gridpack_HAHM_variablesw_v3_eps1e-2_MZd8_lhaid306000/cmsgrid_final.lheMZd8_15000events_lhaid306000WITHcuts.root" 
//inputFile6 = inputDir+"gridpack_HAHM_variablesw_v3_eps1e-2_MZd16_lhaid306000/cmsgrid_final.lheMZd16_15000events_lhaid306000WITHcuts.root" 
//inputFile7 = inputDir+"gridpack_HAHM_variablesw_v3_eps1e-2_MZd32_lhaid306000/cmsgrid_final.lheMZd32_15000events_lhaid306000WITHcuts.root" 
//
//inputTree = "lheEvents_tchan"
//
//## Use -1 to print ALL events
//numEvents = -1                                                                                               
//bins = 200
//PRINT = 0                                                                                             
//kinemVar1 = "eta"
//kinemVar2 = "phi"
//
//
//min_dRarr1 = minDRarray(inputFile1,inputTree,numEvents,PRINT)
//min_dRarr2 = minDRarray(inputFile2,inputTree,numEvents,PRINT)
//min_dRarr3 = minDRarray(inputFile3,inputTree,numEvents,PRINT)
//min_dRarr4 = minDRarray(inputFile4,inputTree,numEvents,PRINT)
//min_dRarr5 = minDRarray(inputFile5,inputTree,numEvents,PRINT)
//min_dRarr6 = minDRarray(inputFile6,inputTree,numEvents,PRINT)
//min_dRarr7 = minDRarray(inputFile7,inputTree,numEvents,PRINT)
//
//setTDRStyle()
//
//## Parameters
//#xmin = 0.025  
//xmin = 0.0
//xmax = 2.5 
//
//binwidth = (xmax-xmin)/float(bins)
//
//## Create histograms
//#bins, h1 = histogram(min_dRarr1, bins, xmin, xmax, drawstyle='hist')
//#bins, h2 = histogram(min_dRarr2, drawstyle='hist')
//h1 = TH1F("h1","h1", bins, xmin, xmax)
//for event in range(len(min_dRarr1)):
//    h1.Fill(min_dRarr1[event])
//h2 = TH1F("h2","h2", bins, xmin, xmax)
//for event in range(len(min_dRarr2)):
//    h2.Fill(min_dRarr2[event])
//h3 = TH1F("h3","h3", bins, xmin, xmax)
//for event in range(len(min_dRarr3)):
//    h3.Fill(min_dRarr3[event])
//h4 = TH1F("h4","h4", bins, xmin, xmax)
//for event in range(len(min_dRarr4)):
//    h4.Fill(min_dRarr4[event])
//h5 = TH1F("h5","h5", bins, xmin, xmax)
//for event in range(len(min_dRarr5)):
//    h5.Fill(min_dRarr5[event])
//h6 = TH1F("h6","h6", bins, xmin, xmax)
//for event in range(len(min_dRarr6)):
//    h6.Fill(min_dRarr6[event])
//h7 = TH1F("h7","h7", bins, xmin, xmax)
//for event in range(len(min_dRarr3)):
//    h7.Fill(min_dRarr7[event])
//
//c1 = TCanvas("c1","c1",1000,1000)                           
//c1.cd()                                                   
//
//h1.SetTitle("Rando Title") #Title not working
//h1.Scale(1./h1.Integral())
//h2.Scale(1./h2.Integral())
//h3.Scale(1./h3.Integral())
//h4.Scale(1./h4.Integral())
//h5.Scale(1./h5.Integral())
//h6.Scale(1./h6.Integral())
//h7.Scale(1./h7.Integral())
//
//## For setting the range on axes
//## SetAxisRange creates weird spike behavior
//#h1.SetAxisRange(0.0, 0.1, "Y")
//#h2.SetAxisRange(0.0, 0.1, "Y")
//
//## Try GetXaxis instead if you need to set range
//#h1.GetXaxis().SetRangeUser(xmin, xmax)
//#h2.GetXaxis().SetRangeUser(xmin, xmax)
//h1.SetLineColor(kRed)                               
//h2.SetLineColor(kRed+2)
//h3.SetLineColor(kMagenta)                               
//h4.SetLineColor(kGreen)                               
//h5.SetLineColor(kBlack)                               
//h6.SetLineColor(kBlue)                               
//h7.SetLineColor(kCyan+1)                               
//
//h1.SetXTitle("deltaR min (4 leptons/event)")
//
//h1.SetLineWidth(3)
//h2.SetLineWidth(3)
//h3.SetLineWidth(3)
//h4.SetLineWidth(3)
//h5.SetLineWidth(2)
//h6.SetLineWidth(2)
//h7.SetLineWidth(2)
//
//h1.SetYTitle("Events")
//
//h1.SetLabelSize(0.025, "X")
//h1.SetLabelSize(0.025, "Y")
//#h1.SetLabelOffset(0.003, "X")
//#h1.SetLabelOffset(0.003, "Y")
//
//h1.Draw("hist")                                  
//h2.Draw("hist same")                                  
//h3.Draw("hist same")                                  
//h4.Draw("hist same")                                  
//h5.Draw("hist same")                                  
//h6.Draw("hist same")                                  
//h7.Draw("hist same")                                  
//
//legend = TLegend(0.75,0.70,0.95,0.90)
//legend.AddEntry(h1,"mZd = 1 GeV" ,"lpf")
//legend.AddEntry(h2,"mZd = 2 GeV" ,"lpf")
//legend.AddEntry(h3,"mZd = 3 GeV" ,"lpf")
//legend.AddEntry(h4,"mZd = 4 GeV" ,"lpf")
//legend.AddEntry(h5,"mZd = 8 GeV" ,"lpf")
//legend.AddEntry(h6,"mZd = 16 GeV" ,"lpf")
//legend.AddEntry(h7,"mZd = 32 GeV" ,"lpf")
//#legend.AddEntry(h2,"lhaid = %s" % pdf2,"lpf")
//legend.SetLineWidth(3)
//legend.SetBorderSize(0)
//#legend.SetTextSize(0.03)
//legend.Draw("same")
//#c1.SetLogx()
//#c1.SetLogy()
//c1.SaveAs(outputFile)
//## Wait for user input to close the canvas when ready
//wait = raw_input()
