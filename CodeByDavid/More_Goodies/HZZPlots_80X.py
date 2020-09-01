import ROOT, sys, os, string, re
from ROOT import *
from array import array
import math
from math import *

from tdrStyle import *
setTDRStyle()

from LoadData import *
LoadData()

def plot_m4l(channel, var, mh, bin, low, high, m4llow, m4lhigh, xlabel, xunits, prelim, setLogX, setLogY,EventCat):

    print ''
    print '=========================================='
    print channel,var,low,'-',high,xunits

    if (var=='abs((Z_phi[0]-Z_phi[1])>TMath::Pi()?(Z_phi[0]-Z_phi[1]-2*TMath::Pi()):(Z_phi[0]-Z_phi[1]+2*TMath::Pi()))'): savevar = 'dPhiZZ'
    elif (var=='@jet_pt.size()'): savevar = 'njets_pt30_eta4p7'
    elif (var=='-1.0*TMath::Log(me_qqZZ_MCFM)'): savevar = 'nl_me_qqZZ_MCFM'
    elif (var=='-1.0*TMath::Log(me_0plus_JHU)'): savevar = 'nl_me_0plus_JHU'
    elif (var=='MyMath::Rank(me_qqZZ_MCFM)'): savevar = 'rank_me_qqZZ_MCFM'
    elif (var=="mass4lErr/mass4l"): savevar = "relmasserr"
    else: savevar = var
    
    save = savevar+'_'+str(low).replace('.','pt')+'to'+str(high).replace('.','pt')+'_'+channel+'_mh'+mh

    doratio = False

    #lumi2016 = 7649.0 
    #lumiplot2016 = '7.65 fb^{-1}'
    
    #lumi2016 = 12900.0 
    #lumiplot2016 = '12.9 fb^{-1}'
    #save = save+'_12pt9fb'

    lumi2016 = 33590.0 
    lumiplot2016 = '33.6 fb^{-1}'
    save = save+'_33pt6fb'

    List = [
        'ZZTo4L_13TeV_powheg_pythia8_RunIISpring16MiniAODv2',
        'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2',
        'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2',
        'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2',
        'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2',
        'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2',
        'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2',
        'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2',
        'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2',
        'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_RunIISpring16MiniAODv2',
        'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_RunIISpring16MiniAODv2',
        'ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8_RunIISpring16MiniAODv2',
        'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2',
        #'Data_2016B_4lskim',
        #'Data_2016C_4lskim',
        #'Data_2016_13fb_4lskim',
        'Data_Run2016_PromptReco' 
        ]

    if (EventCat=="-1"):
        passedSelection  = "passedZ4lSelection==1"
    else:
        passedSelection  = "passedZ4lSelection==1 && EventCat=="+EventCat
    print "EventCat==",EventCat 
 
    if (channel == '4l'):      
        cut = passedSelection+' && mass4l>0 && '+var+'>='+str(low)+' && '+var+'<'+str(high)+' && mass4l>='+str(m4llow)+' && mass4l<'+str(m4lhigh)
    if (channel == '2e2mu'):      
        cut = passedSelection+' && mass2e2mu>0 && '+var+'>='+str(low)+' && '+var+'<'+str(high)+' && mass4l>='+str(m4llow)+' && mass4l<'+str(m4lhigh)
    if (channel == '4mu'):      
        cut = passedSelection+' && mass4mu>0 && '+var+'>='+str(low)+' && '+var+'<'+str(high)+' && mass4l>='+str(m4llow)+' && mass4l<'+str(m4lhigh)
    if (channel == '4e'):      
        cut = passedSelection+' && mass4e>0 && '+var+'>='+str(low)+' && '+var+'<'+str(high)+' && mass4l>='+str(m4llow)+' && mass4l<'+str(m4lhigh)

    Variable = {}
    
    stack = THStack('stack', 'stack')
    added = TH1D('a', 'a',bin,low,high)
    red_bkg = TH1D('red_bkg', 'red_bkg',bin,low,high)
    irr_bkg = TH1D('irr_bkg', 'irr_bkg',bin,low,high)
    rare_bkg = TH1D('rare_bkg', 'rare_bkg',bin,low,high)
    qqzz = TH1D('qqzz', 'qqzz',bin,low,high)
    ggzz = TH1D('ggzz', 'ggzz',bin,low,high)
    sig = TH1D('sig', 'sig',bin,low,high)
    data = TH1D('data', 'data',bin,low,high)

    for Sample in List:
        print Sample
        if (not Sample in Tree): continue

        if (Sample.startswith('Data')): 
            weight  = '(1.0)'
            #cut = cut + ' && !(mass4l>110 && mass4l<150) && !(mass4l>500.0)'
            cut = cut + ' && !(mass4l>110 && mass4l<150) '
        else: 
            #weight = '(genWeight*crossSection*'+str(lumi2016)+'/'+str(sumw[Sample])+')'
            #weight = '(genWeight*crossSection*dataMCWeight*'+str(lumi2016)+'/'+str(sumw[Sample])+')'
            weight = '(genWeight*crossSection*dataMCWeight*pileupWeight*'+str(lumi2016)+'/'+str(sumw[Sample])+')'
            #datamcweight = "passedZ4lSelection*lep_dataMC[passedZ4lSelection*lep_Hindex[0]]*lep_dataMC[passedZ4lSelection*lep_Hindex[1]]*lep_dataMC[passedZ4lSelection*lep_Hindex[2]]*lep_dataMC[passedZ4lSelection*lep_Hindex[3]]"
            #weight = '(genWeight*crossSection*'+datamcweight+'*'+str(lumi2016)+'/'+str(sumw[Sample])+')'
            
        if (Sample.startswith('ZZTo4L')):
            weight = weight + '*(k_qqZZ_qcd_M*k_qqZZ_ewk)'
        if ('MCFM' in Sample): 
            weight = weight + '*(k_ggZZ)'

        histName = Sample+channel
        Variable[histName] = TH1D(histName, histName, bin, low, high)        
        if ('Data' in Sample):
            cut = cut.replace("passedZ4lSelection","passedFullSelection")
            Tree[Sample].Draw(var + " >> " + histName, weight+"*("+ cut + ")", 'goff')
        else:
            if ('125' in Sample): cut = cut.replace("passedZ4lSelection","passedFullSelection")
            Tree[Sample].Draw(var + " >> " + histName, weight+"*("+ cut + ")", 'goff')

        print Sample,Variable[histName].Integral()
        
        if (Sample.startswith('Data')):
            data.Add(Variable[histName])
        elif (Sample.startswith('ZZTo4L')):
            irr_bkg.Add(Variable[histName])
            qqzz.Add(Variable[histName])
        elif ('MCFM' in Sample):
            irr_bkg.Add(Variable[histName])
            ggzz.Add(Variable[histName])
        elif (('M125' in Sample) or ('Seesaw' in Sample)):
            sig.Add(Variable[histName])
        elif ('ttZ' in Sample):
            rare_bkg.Add(Variable[histName])
#        else:
#            red_bkg.Add(Variable[histName])

        if (not 'Data' in Sample):
            added.Add(Variable[histName])

    #if (var=='mass'+channel): c1 = TCanvas("c1","c1", 1200, 800)
    c1 = TCanvas("c1","c1", 800, 800)
    if (setLogX): c1.SetLogx()
    if (setLogY): c1.SetLogy()
    if (doratio): c1.SetBottomMargin(0.3)
    c1.SetRightMargin(0.03);

    ################################
    #### Reducible Background ###
    ################################

    #red_bkg.Smooth(1)
    red_bkg.SetLineColor(TColor.GetColor("#003300"))
    red_bkg.SetLineWidth(2)
    red_bkg.SetFillColor(TColor.GetColor("#669966"))

    ### red_bkg ###
    if (channel=="4e"):
        fsum = TF1("fsum","landau(0)")
        fsum.SetParameter(0,1.0)
        fsum.SetParameter(1,141.9)
        fsum.SetParameter(2,21.3)
    if (channel=="4mu"):
        fsum = TF1("fsum","landau(0)")
        fsum.SetParameter(0,1.0)
        fsum.SetParameter(1,130.4)
        fsum.SetParameter(2,15.6)
    if (channel=="2e2mu"):
        fsum = TF1("fsum","landau(0)+landau(3)")
        fsum.SetParameter(0,0.55)
        fsum.SetParameter(1,131.1)
        fsum.SetParameter(2,18.1)
        fsum.SetParameter(3,0.45)
        fsum.SetParameter(4,133.8)
        fsum.SetParameter(5,18.9)
    if (channel=="4l"):
        fsum = TF1("fsum","landau(0)+landau(3)+landau(6)+landau(9)")
        fsum.SetParameter(0,9.8)
        fsum.SetParameter(1,141.9)
        fsum.SetParameter(2,21.3)
        fsum.SetParameter(3,10.2)
        fsum.SetParameter(4,130.4)
        fsum.SetParameter(5,15.6)
        fsum.SetParameter(6,0.55*20.4)
        fsum.SetParameter(7,131.1)
        fsum.SetParameter(8,18.1)
        fsum.SetParameter(9,0.45*20.4)
        fsum.SetParameter(10,133.8)
        fsum.SetParameter(11,18.9)
    
    nfill=100000
    htemp = TH1F("htemp","htemp",1300,70,2000);
    htemp.FillRandom("fsum",nfill);
    if (channel=="4e"): htemp.Scale(9.8/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
    if (channel=="4mu"): htemp.Scale(10.2/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
    if (channel=="2e2mu"): htemp.Scale(20.4/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
    if (channel=="4l"): htemp.Scale((9.8+10.2+20.4)/htemp.Integral(htemp.FindBin(70),htemp.FindBin(2000)))
    red_bkg.FillRandom("fsum",nfill);
    red_bkg.Scale(htemp.Integral(htemp.FindBin(low),htemp.FindBin(high))/red_bkg.Integral(red_bkg.FindBin(low),red_bkg.FindBin(high)))

    #red_bkg.Scale(lumi2016/(1000.0*7.65)) # lumi

    ################################
    #### Rare Backgrounds ###
    ################################

    rare_bkg.SetFillColor(kOrange+2)
    rare_bkg.SetLineColor(kOrange+3)
    rare_bkg.SetLineWidth(2)


    ################################
    #### Irreducible Backgrounds ###
    ################################

    irr_bkg.SetLineColor(TColor.GetColor("#000099"))
    irr_bkg.SetLineWidth(2)
    irr_bkg.SetFillColor(TColor.GetColor("#3366ff"))
        
    qqzz.SetLineColor(TColor.GetColor("#000099"))
    qqzz.SetLineWidth(2)
    qqzz.SetFillColor(TColor.GetColor("#99ccff"))

    ggzz.SetFillColor(TColor.GetColor("#3366ff"))
    ggzz.SetLineWidth(2)
    ggzz.SetLineColor(TColor.GetColor("#000099"))

    sig.SetFillColor(TColor.GetColor("#ffafaf"))
    sig.SetLineColor(TColor.GetColor("#cc0000"))
    sig.SetLineWidth(2)



    #stack.Add(rare_bkg)
    if (var=='mass'+channel and EventCat=="-1"): stack.Add(red_bkg)
    stack.Add(ggzz)
    stack.Add(qqzz)
    stack.Add(sig)    

    stack.Draw("hist")
    stack.GetXaxis().SetMoreLogLabels(kTRUE)
    stack.GetXaxis().SetNoExponent(kTRUE)
    stack.SetMinimum(0.01)
    if (setLogY): stack.SetMaximum(14*max(stack.GetMaximum(),data.GetMaximum()+1))
    else: stack.SetMaximum(2.0*max(stack.GetMaximum(),data.GetMaximum()+1))
    if (xunits==''): stack.GetXaxis().SetTitle(xlabel)
    else: stack.GetXaxis().SetTitle(xlabel+' ['+xunits+']')
    if (doratio): stack.GetXaxis().SetTitleSize(0)
    if (doratio): stack.GetXaxis().SetLabelSize(0)
    binsize = str(round(float((high-low)/bin),2))
    if binsize.endswith('.0'): binsize=binsize.rstrip('.0')
    ylabel = 'Events / '+binsize+' '+xunits
    stack.GetYaxis().SetTitle(ylabel)           
    
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetBinErrorOption(TH1.kPoisson)
    data.Draw("ex0psame")

    legend = TLegend(.62,.60,.90,.92)
    #legend = TLegend(.2,.30,.52,.60)
    legend.AddEntry(data, 'Data', "ep")
    legend.AddEntry(sig, 'm_{H} = 125 GeV', "f")
    #legend.AddEntry(irr_bkg, 'Z#gamma*, ZZ', "f")
    legend.AddEntry(qqzz, 'qq#rightarrowZZ', "f")
    legend.AddEntry(ggzz, 'gg#rightarrowZZ', "f")
    legend.AddEntry(red_bkg, 'Z+X', "f")
    legend.SetShadowColor(0);
    legend.SetFillColor(0);
    legend.SetLineColor(0);
    legend.Draw("same")  

    latex2 = TLatex()
    latex2.SetNDC()
    latex2.SetTextSize(0.6*c1.GetTopMargin())
    latex2.SetTextFont(42)
    latex2.SetTextAlign(31) # align right
    latex2.DrawLatex(0.92, 0.94,lumiplot2016+" (13 TeV)")
    latex2.SetTextSize(1.0*c1.GetTopMargin())
    latex2.SetTextFont(62)
    latex2.SetTextAlign(11) # align right
    latex2.DrawLatex(0.2, 0.84, "CMS")
    latex2.SetTextSize(0.8*c1.GetTopMargin())
    latex2.SetTextFont(52)
    latex2.SetTextAlign(11)
    #latex2.DrawLatex(0.28, 0.945, "Unpublished")    
    latex2.DrawLatex(0.2, 0.78, "Preliminary")

    lastbin=bin
    if (var=='mass4l' and high>700.0): lastbin=bin+1
    print ''
    print 'ZZ Bkg:  ',irr_bkg.Integral(1,lastbin)
    print 'ggZZ Bkg:  ',ggzz.Integral(1,lastbin)
    print 'qqZZ Bkg:  ',qqzz.Integral(1,lastbin)
    print 'Red Bkg:  ',red_bkg.Integral(1,lastbin)
    print 'Rare Bkg:  ',rare_bkg.Integral(1,lastbin)
    print 'Signal:  ',sig.Integral(1,lastbin)
    print 'Data  :  ',data.Integral(1,lastbin)
    print ''

    print "KS:",data.KolmogorovTest(added)
    print "Data mean:",data.GetMean()
    print "bkg Mean:",added.GetMean()
    print "sig mean:",sig.GetMean()

    gPad.RedrawAxis()

    if (doratio):

        ratio = data.Clone('ratio')
        ratio.Divide(added)
        ratio.GetXaxis().SetMoreLogLabels(kTRUE)
        ratio.GetXaxis().SetNoExponent(kTRUE)

        pad = TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
        if (setLogX): pad.SetLogx()
        pad.SetTopMargin(0.7);
        pad.SetRightMargin(0.03);
        pad.SetFillColor(0);
        pad.SetGridy(1);
        pad.SetFillStyle(0);
        pad.Draw();        
        pad.cd(0);
        
        #ratio.SetLineColor(1);
        #ratio.SetMarkerColor(1);
        if (xunits==''):
            ratio.GetXaxis().SetTitle(xlabel)
        else:    
            ratio.GetXaxis().SetTitle(xlabel+' ['+xunits+']')
        ratio.GetXaxis().SetTitle
        ratio.GetYaxis().SetTitleSize(0.04);
        ratio.GetYaxis().SetTitleOffset(1.8);
        ratio.GetYaxis().SetTitle("Data/Bkg.");
        ratio.GetYaxis().CenterTitle();
        ratio.GetYaxis().SetLabelSize(0.03);
        #ratio.SetMarkerStyle(20);
        #ratio.SetMarkerSize(1.2);
        ratio.SetLineWidth(2)
        ratio.SetLineColor(1)
        ratio.SetMarkerColor(1)
        ratio.SetLineStyle(1)
        ratio.SetMinimum(0.61);
        ratio.SetMaximum(1.39);
        ratio.Draw("hist");

    c1.SaveAs('Histo_' + save + '.pdf')
    c1.SaveAs('Histo_' + save + '.png')
      

    print '=========================================='
    print ''

    del c1
    del stack
    del added
    del Variable


#plot_m4l('4l', 'mass4l', '125', 204, 70.0, 886.0, 70.0, 886.0, 'm_{4l}', 'GeV', True, True, False,"-1")
#plot_m4l('4l', 'mass4l', '125', 204, 70.0, 2000.0, 70.0, 2000.0, 'm_{4l}', 'GeV', True, True, False,"-1")
#plot_m4l('4l', 'mass4l', '125', 12, 118.0, 130.0, 118, 130.0, 'm_{4l}', 'GeV', True, False, False,"-1")
plot_m4l('4l', 'mass4l', '125', 156,  105.0, 1055.0, 105.0, 1055.0, 'm_{4l}', 'GeV', True, True, False,"-1")
#plot_m4l('4l', 'mass4l', '125', 15, 105.0, 140.0, 105.0, 140.0, 'm_{4l}', 'GeV', True, False, False,"-1")
#plot_m4l('4l', 'mass4l', '125', 35, 105.0, 175.0, 105.0, 175.0, 'm_{4l}', 'GeV', True, False, False,"-1")
#plot_m4l('4l', 'mass4l', '125', 100, 200.0, 1000.0, 200.0, 1000.0, 'm_{4l}', 'GeV', True, False, True,"-1")
#plot_m4l('4l', 'mass4lREFIT', '125', 15, 105.0, 140.0, 105.0, 140.0, "m_{4l}^{'}", 'GeV', True, False, False,"-1")
#plot_m4l('4l', 'mass4lErr/mass4l', '125', 30, 0.0, 0.03, 118.0, 130.0, "#sigma(m_{4l})/m_{4l}", '', True, False, False,"-1")

#plot_m4l('4mu', 'mass4mu', '125', 204, 70.0, 886.0, 70.0, 886.0, 'm_{4#mu}', 'GeV', True, True, False,"-1")
#plot_m4l('4mu', 'mass4mu', '125', 204, 70.0, 2000.0, 70.0, 2000.0, 'm_{4#mu}', 'GeV', True, True, False,"-1")
#plot_m4l('4mu', 'mass4mu', '125', 12, 118.0, 130.0, 118, 130.0, 'm_{4#mu}', 'GeV', True, False, False,"-1")
plot_m4l('4mu', 'mass4mu', '125', 156,  105.0, 1055.0, 105.0, 1055.0, 'm_{4#mu}', 'GeV', True, True, False,"-1")
#plot_m4l('4mu', 'mass4mu', '125', 15, 105.0, 140.0, 105.0, 140.0, 'm_{4#mu}', 'GeV', True, False, False,"-1")
#plot_m4l('4mu', 'mass4mu', '125', 30, 105.0, 175.0, 105.0, 175.0, 'm_{4#mu}', 'GeV', True, False, False,"-1")
#plot_m4l('4mu', 'mass4lREFIT', '125', 15, 105.0, 140.0, 105.0, 140.0, "m_{4#mu}^{'}", 'GeV', True, False, False,"-1")
#plot_m4l('4mu', 'mass4lErr/mass4l', '125', 30, 0.0, 0.03, 118.0, 130.0, "#sigma(m_{4l})/m_{4l}", '', True, False, False,"-1")

#plot_m4l('4e', 'mass4e', '125', 204, 70.0, 886.0, 70.0, 886.0, 'm_{4e}', 'GeV', True, True, False,"-1")
#plot_m4l('4e', 'mass4e', '125', 204, 70.0, 2000.0, 70.0, 2000.0, 'm_{4e}', 'GeV', True, True, False,"-1")
#plot_m4l('4e', 'mass4e', '125', 12, 118.0, 130.0, 118, 130.0, 'm_{4e}', 'GeV', True, False, False,"-1")
plot_m4l('4e', 'mass4e', '125', 156,  105.0, 1055.0, 105.0, 1055.0, 'm_{4e}', 'GeV', True, True, False,"-1")
#plot_m4l('4e', 'mass4e', '125', 15, 105.0, 140.0, 105.0, 140.0, 'm_{4e}', 'GeV', True, False, False,"-1")
#plot_m4l('4e', 'mass4e', '125', 30, 105.0, 175.0, 105.0, 175.0, 'm_{4e}', 'GeV', True, False, False,"-1")
#plot_m4l('4e', 'mass4e', '125', 20, 70.0, 110.0, 70.0, 110.0, 'm_{4e}', 'GeV', True, True, False,"-1")
#plot_m4l('4e', 'mass4e', '125', 54, 130.0, 400.0, 130.0, 400.0, 'm_{4e}', 'GeV', True, True, False,"-1")
#plot_m4l('4e', 'mass4lREFIT', '125', 15, 105.0, 140.0, 105.0, 140.0, "m_{4e}^{'}", 'GeV', True, False, False,"-1")
#plot_m4l('4e', 'mass4lErr/mass4l', '125', 30, 0.0, 0.03, 118.0, 130.0, "#sigma(m_{4l})/m_{4l}", '', True, False, False,"-1")

#plot_m4l('2e2mu', 'mass2e2mu', '125', 204, 70.0, 886.0, 70.0, 886.0, 'm_{2e2#mu}', 'GeV', True, True, False,"-1")
#plot_m4l('2e2mu', 'mass2e2mu', '125', 204, 70.0, 2000.0, 70.0, 2000.0, 'm_{2e2#mu}', 'GeV', True, True, False,"-1")
#plot_m4l('2e2mu', 'mass2e2mu', '125', 12, 118.0, 130.0, 118, 130.0, 'm_{2e2#mu}', 'GeV', True, False, False,"-1")
plot_m4l('2e2mu', 'mass2e2mu', '125', 156,  105.0, 1055.0, 105.0, 1055.0, 'm_{2e2#mu}', 'GeV', True, True, False,"-1")
#plot_m4l('2e2mu', 'mass2e2mu', '125', 15, 105.0, 140.0, 105.0, 140.0, 'm_{2e2#mu}', 'GeV', True, False, False,"-1")
#plot_m4l('2e2mu', 'mass2e2mu', '125', 30, 105.0, 175.0, 105.0, 175.0, 'm_{2e2#mu}', 'GeV', True, False, False,"-1")
#plot_m4l('2e2mu', 'mass4lREFIT', '125', 15, 105.0, 140.0, 105.0, 140.0, "m_{2e2#mu}^{'}", 'GeV', True, False, False,"-1")
#plot_m4l('2e2mu', 'mass4lErr/mass4l', '125', 30, 0.0, 0.03, 118.0, 130.0, "#sigma(m_{4l})/m_{4l}", '', True, False, False,"-1")


#for i in range(0,6):
#    plot_m4l('4l', 'mass4l', '125', 12, 118.0, 130.0, 118, 130.0, 'm_{4l}', 'GeV', True, False, False,str(i))
#    plot_m4l('4mu', 'mass4mu', '125', 12, 118.0, 130.0, 118, 130.0, 'm_{4#mu}', 'GeV', True, False, False,str(i))
#    plot_m4l('4e', 'mass4e', '125', 12, 118.0, 130.0, 118, 130.0, 'm_{4e}', 'GeV', True, False, False,str(i))
#    plot_m4l('2e2mu', 'mass2e2mu', '125', 12, 118.0, 130.0, 118, 130.0, 'm_{2e2#mu}', 'GeV', True, False, False,str(i))

#plot_m4l('4l', 'mass4lErr', '125', 20, 0.0, 10.0, 150.0, 2000.0, "#sigma(m_{4l}) (reco.)", 'GeV', True, False, False,"-1")
#plot_m4l('4l', 'mass4lErrREFIT', '125', 20, 0.0, 10.0, 150.0, 2000.0, "#sigma(m_{4l}) (refit.)", 'GeV', True, False, False,"-1")

#plot_m4l('4l', 'mass4lErr', '125', 20, 0.0, 10.0, 70.0, 105.0, "#sigma(m_{4l}) (reco.)", 'GeV', True, False, False,"-1")
#plot_m4l('4l', 'mass4lErrREFIT', '125', 20, 0.0, 10.0, 70.0, 105.0, "#sigma(m_{4l}) (refit.)", 'GeV', True, False, False,"-1")


