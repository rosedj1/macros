import ROOT, sys, os, string, re
from ROOT import *
from array import array
import math
from math import *

from tdrStyle import *
setTDRStyle()

from LoadData import *

def plot_m4l(channel, var, mh, bin, low, high, xlabel, xunits, prelim, setLog):

    print ''
    print '=========================================='
    print channel,var,low,'-',high,xunits

    save = var+'_'+str(low).replace('.','pt')+'to'+str(high).replace('.','pt')+'_'+channel+'_mh'+mh

    lumi2015 = 10000
    lumiplot2015 = '10.0'

    List = ['ZZTo4L']

    if (channel == '4l'):      
        cut = 'passedFullSelection==1 && mass4l>'+str(low)+' && mass4l<'+str(high)
    if (channel == '2e2mu'):      
        cut = 'passedFullSelection==1 && mass2e2mu>'+str(low)+' && mass2e2mu<'+str(high)
    if (channel == '4mu'):      
        cut = 'passedFullSelection==1 && mass4mu>'+str(low)+' && mass4mu<'+str(high)
    if (channel == '4e'):      
        cut = 'passedFullSelection==1 && mass4e>'+str(low)+' && mass4e<'+str(high)

    Variable = {}
    
    stack = THStack('a', 'a')
    added = TH1D('a', 'a',bin,low,high)

    for Sample in List:

        if (not Sample in TreesPassedEvents): continue

        histName = Sample+channel
        Variable[histName] = TH1D(histName, histName, bin, low, high)        

        TreesPassedEvents[Sample].Draw(var + " >> " + histName, "(" + cut + ")", 'goff')
            
        if nEvents[Sample] != 0:
            Variable[histName].Scale(lumi2015) 
            
    c1 = TCanvas("c1","c1", 600, 633)


    ################################
    #### Irreducible Backgrounds ###
    ################################
    
    Variable['ZZTo4L'+channel].SetLineColor(4)
    Variable['ZZTo4L'+channel].SetLineWidth(2)
    Variable['ZZTo4L'+channel].SetFillColor(ROOT.kAzure+6)

    stack.Add(Variable['ZZTo4L'+channel])

    stack.SetMinimum(0.0)
    stack.Draw()
    stack.GetXaxis().SetTitle(xlabel+' ['+xunits+']')
    binsize = str(int((high-low)/bin))
    ylabel = 'Events / '+binsize+' '+xunits
    stack.GetYaxis().SetTitle(ylabel)
        
    legend = TLegend(.62,.60,.92,.92)
    legend.AddEntry(Variable['ZZTo4L'+channel], 'Z#gamma*, ZZ', "f")
    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw("same")  

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.03)   
    latex2.SetTextFont(42);
    latex2.SetTextAlign(31) # align right
    if (not prelim):
        latex2.DrawLatex(0.2, 0.95, "CMS"),
    else:
        save = save+'_prelim'
        latex2.DrawLatex(0.2, 0.95, "CMS Simulation"),

    latex2.DrawLatex(0.87,0.95, "#sqrt{s} = 13 TeV, L="+lumiplot2015+" fb^{-1}");

    print ''
    #print 'ZZ Bkg:  ',irr_bkg.Integral(1,bin)
    print ''

    c1.SaveAs('Histo_' + save + '.pdf')
      

    print '=========================================='
    print ''

    del c1
    del stack
    del added
    del Variable

#plot_m4l('4e', 'mass4l', '125', 95, 50.0, 1000.0, 'm_{4e}', 'GeV', false, false)
plot_m4l('4l', 'mass4l', '125', 30, 50.0, 140.0, 'm_{4l}', 'GeV', true, false)
