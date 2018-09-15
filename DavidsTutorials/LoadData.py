from ROOT import *
from array import array
import os

# 80X samples
dirMC_80 = '/cms/data/store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_MC80X_20160725/'
#dirData_80 = '/cms/data/store/user/t2/users/dsperka/rootfiles_Data80X_20160725/'
dirData_80 = '/cms/data/store/user/t2/users/dsperka/Run2/HZZ4l/SubmitArea_13TeV/rootfiles_Data80X_20161031_4lskim/'

SamplesMC_80 = [
#'DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2.root',
#'DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2.root',
#'DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2.root',
#'DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2.root',
#'DYBJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2.root',
#'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring16MiniAODv2.root',
#'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2.root',
#'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring16MiniAODv2.root',
#'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring16MiniAODv2.root',
#'GluGluHToZZTo4L_M120_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
#'GluGluHToZZTo4L_M124_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
'GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
#'GluGluHToZZTo4L_M126_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv1.root',
#'GluGluHToZZTo4L_M130_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
'GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2.root',
'GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2.root',
'GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2.root',
'GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2.root',
'GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2.root',
'GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2.root',
#'GluGluToHiggs0PMContinToZZTo2e2mu_M125_GaSM_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2.root',
#'GluGluToHiggs0PMContinToZZTo2e2tau_M125_GaSM_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2.root',
#'GluGluToHiggs0PMContinToZZTo2mu2tau_M125_GaSM_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2.root',
#'GluGluToHiggs0PMContinToZZTo4tau_M125_GaSM_13TeV_MCFM701_pythia8_RunIISpring16MiniAODv2.root',
#'VBF_HToZZTo4L_M120_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
#'VBF_HToZZTo4L_M124_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
'VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
#'VBF_HToZZTo4L_M126_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
#'VBF_HToZZTo4L_M130_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
#'WWTo2L2Nu_13TeV-powheg_RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1.root',
#'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_RunIISpring16MiniAODv2.root',
#'WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring16MiniAODv2.root',
'WminusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
'WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
'ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
#'ZZTo2L2Nu_13TeV_powheg_pythia8_RunIISpring16MiniAODv2.root',
#'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_RunIISpring16MiniAODv2.root',
'ZZTo4L_13TeV_powheg_pythia8_RunIISpring16MiniAODv2.root',
'ttH_HToZZ_4LFilter_M125_13TeV_powheg2_JHUgenV6_pythia8_RunIISpring16MiniAODv2.root',
#'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIISpring16MiniAODv2_4lSkim.root',
#'TTTo2L2Nu_13TeV-powheg_RunIISpring16MiniAODv2_4lSkim.root',
]

SamplesData_80 = [
#'Data_2016B_4lskim.root',
#'Data_2016C_4lskim.root',
#'Data_2016_13fb_4lskim.root'
'Data_Run2016_PromptReco.root'
]


###################################################### 
RootFile = {} 
Tree = {} 
nEvents = {} 
sumw = {}

def LoadData():
    # 80X MC
    for i in range(0,len(SamplesMC_80)):
        
        sample = SamplesMC_80[i].rstrip('.root')
        
        RootFile[sample] = TFile(dirMC_80+'/'+sample+'.root',"READ")
        if ('4lskim' in sample):
            Tree[sample]  = RootFile[sample].Get("passedEvents")
        else:
            Tree[sample]  = RootFile[sample].Get("Ana/passedEvents")
        
        h_nevents = RootFile[sample].Get("Ana/nEvents")
        h_sumw = RootFile[sample].Get("Ana/sumWeights")
        
        if (h_nevents): nEvents[sample] = h_nevents.Integral()
        else: nEvents[sample] = 0.

        if (h_sumw): sumw[sample] = h_sumw.Integral()
        else: sumw[sample] = 0.
        
        if (not Tree[sample]): print sample+' has no passedEvents tree'
        else:
            print sample,"nevents",nEvents[sample],"sumw",sumw[sample]
            
            
    for i in range(0,len(SamplesData_80)):
                
        sample = SamplesData_80[i].replace('.root','')
        
        RootFile[sample] = TFile(dirData_80+'/'+sample+'.root',"READ")
        Tree[sample]  = RootFile[sample].Get("passedEvents")
        if (not Tree[sample]): Tree[sample]  = RootFile[sample].Get("Ana/passedEvents")

        h_nevents = RootFile[sample].Get("nEvents")
        h_sumw = RootFile[sample].Get("sumWeights")
        
        if (h_nevents):
            nEvents[sample] = h_nevents.Integral()
            sumw[sample] = h_sumw.Integral()
        else:
            nEvents[sample] = 0.
            sumw[sample] = 0.
            
        if (not Tree[sample]): print sample+' has no passedEvents tree'
        else:
            print sample,"nevents",nEvents[sample],"sumw",sumw[sample]
            
#LoadData()

