from ROOT import *
import pandas as pd
#import seaborn as sns

f = TFile("/raid/raid7/lucien/Higgs/DarkZ-NTuple/20180706/ZD_UpTo0j_MZD30_Eps1e-2_klo.root","READ")
t = f.Get("Ana/passedEvents")

## Get total number of 4e events and 4mu events
n_4e 	= t.GetEntries("Sum$(abs(GENlep_id[])==11)==4")
n_4mu 	= t.GetEntries("Sum$(abs(GENlep_id[])==13)==4")
n_2e2mu = t.GetEntries("Sum$(abs(GENlep_id[])==11)==2 && Sum$(abs(GENlep_id[])==13)==2")

n_acc_4e 	= t.GetEntries("Sum$(abs(GENlep_id[])==11)==4 && passedFiducialSelection==1")
n_acc_4mu 	= t.GetEntries("Sum$(abs(GENlep_id[])==13)==4 && passedFiducialSelection==1")
n_acc_2e2mu 	= t.GetEntries("Sum$(abs(GENlep_id[])==11)==2 && Sum$(abs(GENlep_id[])==13)==2 && passedFiducialSelection==1")
n_acceff_4e 	= t.GetEntries("Sum$(abs(GENlep_id[])==11)==4 && passedFullSelection==1")
n_acceff_4mu 	= t.GetEntries("Sum$(abs(GENlep_id[])==13)==4 && passedFullSelection==1")
n_acceff_2e2mu 	= t.GetEntries("Sum$(abs(GENlep_id[])==11)==2 && Sum$(abs(GENlep_id[])==13)==2 && passedFullSelection==1")

A_4e 		= float(n_acc_4e)/float(n_4e)
A_4mu 		= float(n_acc_4mu)/float(n_4mu)
A_2e2mu 	= float(n_acc_2e2mu)/float(n_2e2mu)
A_eff_4e 	= float(n_acceff_4e)/float(n_4e)
A_eff_4mu 	= float(n_acceff_4mu)/float(n_4mu)
A_eff_2e2mu = float(n_acceff_2e2mu)/float(n_2e2mu)
eff_4e		= A_eff_4e / A_4e
eff_4mu		= A_eff_4mu / A_4mu
eff_2e2mu	= A_eff_2e2mu / A_2e2mu

## Python 2.7.X must use float division since by default it is integer division
#print "mZd30"
#print "acc. 4e = ", 		A_4e
#print "acc*eff 4e = ",		A_eff_4e
#print "eff 4e = ",		eff_4e, "\n"
#print "acc. 4mu = ", 		A_4mu
#print "acc*eff 4mu = ",		A_eff_4mu
#print "eff 4mu = ",		eff_4mu, "\n"	 
#print "acc. 2e2mu = ", 		A_2e2mu
#print "acc*eff 2e2mu = ",	A_eff_2e2mu
#print "eff 2e2mu = ",		eff_2e2mu, "\n"

#h = {}

#h["massZ2"] = TH1D("h_massZ2","h_massZ2",48,4,52)
#h["massZ2"].Sumw2()

#t.Draw("massZ2>>h_massZ2","passedFullSelection==1","goff")

#from tdrStyle import *
#setTDRStyle()

#c1 = TCanvas("c1","c1",800,800)
#c1.cd()
#h["massZ2"].SetLineColor(2)
#h["massZ2"].Draw("hist")

#c1.SaveAs("test.pdf")
