# import ROOT in batch mode
import sys
oldargv = sys.argv[:]
sys.argv = [ '-b-' ]
from ROOT import *
from array import *


#setenv LHAPDF_BASE `scram tool info lhapdf | grep LHAPDF_BASE | sed -e s%LHAPDF_BASE=%%`
#setenv PATH $LHAPDF_BASE/bin/:$PATH
#setenv LHAPDF_DATA_PATH `$LHAPDF_BASE/bin/lhapdf-config --datadir`
#setenv PYTHONPATH $LHAPDF_BASE/lib/python2.7/site-packages:$PYTHONPATH
#import lhapdf

gROOT.SetBatch(True)
sys.argv = oldargv
gStyle.SetOptStat(0)
gStyle.SetOptFit(1)

# load FWLite C++ libraries
gSystem.Load("libFWCoreFWLite.so");
AutoLibraryLoader.enable()
gSystem.Load("libDataFormatsFWLite.so");
gSystem.Load("libSimDataFormatsHZZFiducial.so")

# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

gen, genLabel = Handle("GenEventInfoProduct"), ("generator")
lhe, lheLabel = Handle("LHEEventProduct"), ("externalLHEProducer")

#p1 = lhapdf.mkPDF("MSTW2008nnlo68cl", 0)
#p2 = lhapdf.mkPDF(260000)

#handlePruned, prunedLabel  = Handle ("std::vector<reco::GenParticle>"), ("mergedGenParticles")
handlePruned, prunedLabel  = Handle ("std::vector<reco::GenParticle>"), ("prunedGenParticles")
#fid, fidLabel  = Handle ("HZZFid::FiducialSummary"), ("rivetProducerHZZ","FiducialSummary","runRivetAnalysis")

bins = 100; lo = -5; hi = 5.

h1 = TH1F("h1","",bins,lo,hi); h1.Sumw2()
h2 = TH1F("h2","",bins,lo,hi); h2.Sumw2()
h3 = TH1F("h3","",2000,-100,100); h3.Sumw2()
h4 = TH1F("h4","",2000,-100,100); h4.Sumw2()

#events = Events(["testHZZFidRivet_ggH4l_NNLOPS.root"])
#events = Events(["testHZZFidRivet_ggH4l.root"])
events = Events([
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/08914500-C0F0-E611-9948-02163E01417A.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1048526D-BEF0-E611-BDC0-02163E012AEB.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/12404DEA-6ADC-E611-A3F5-0242AC130003.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/148787A1-36DB-E611-88CC-0CC47A4D7630.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/148913CC-C7F0-E611-99EB-02163E01344E.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/16AB26AE-72DC-E611-BE0D-FA163E93EDBE.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1E1422EF-39DB-E611-8DE1-001E675A6AB8.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/1E410583-2ADD-E611-8D4A-0019B9CB01E8.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/264493C4-3DDB-E611-A820-0CC47A4D7600.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/34952E8D-6BDB-E611-BAA0-001E674DA347.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/34FB1712-1DDD-E611-8956-0CC47A6C1056.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/3607D3D1-95DB-E611-919A-ECF4BBE1CEB0.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/364A5E51-39DB-E611-8E6E-008CFA165F5C.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/402ACC27-24DB-E611-B355-001E67A3EC2D.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/407E9F5E-CEF0-E611-A2BA-02163E011837.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/48696757-26DB-E611-B98F-008CFA1113F4.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/48F02F66-0EDB-E611-8430-008CFA111334.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/4C5A0215-C3F0-E611-A988-02163E01A5FE.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/50F795A0-36DB-E611-89A6-0025905B858A.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/52E5B7FE-BFF0-E611-9AC8-02163E019E49.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/54B8D9E8-DEF0-E611-9B45-02163E01A491.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5AC2379E-C5F0-E611-AA22-02163E01419C.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5C6071B3-D6DC-E611-9259-24BE05C48831.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5CB13FC3-C0F0-E611-B32B-02163E019C9B.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5E79909D-36DB-E611-9727-0CC47A4D7634.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/5EBE3DC2-3DDB-E611-B18E-0CC47A4C8E0E.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/646754C9-6ADC-E611-AC66-0025905A608E.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/66A9CB17-52DC-E611-9CA5-008CFA111170.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/66B9DCFA-CFF0-E611-B5F4-02163E019DA2.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/6E1CC376-CBF0-E611-81B9-02163E012B28.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/767C7D76-9CDA-E611-B1BF-24BE05CEECD1.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/7AAB0919-26DB-E611-A426-A4BF01025C16.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/84E62A82-4CDB-E611-9017-0CC47A7C340E.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/86BDFDEF-4BDB-E611-A702-0CC47A4D7670.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/86F0CFB4-BDDB-E611-B63A-90B11C12EA74.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/88FB99A3-6BDB-E611-AFE8-001E67457E9F.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/8C3B2DCF-8BDB-E611-A185-0CC47A4D766C.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/8CBF15D0-32DB-E611-9A61-001E67D80528.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/948D64B1-DCDA-E611-898A-E0071B7AC760.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/961BB3D4-2FDB-E611-B586-0025905A60A8.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/982CAD04-5CDB-E611-98A0-0025905A60CA.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/984CA849-CDF0-E611-ADC9-02163E0134A9.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/AC06AA43-B0DA-E611-A47A-A0369F310374.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/B29E5AD7-8BDB-E611-AD77-0025905B859A.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/B42B3DFB-EFDC-E611-990C-0025907D24F0.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/B8490553-0FDB-E611-9D88-001E675A681F.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/BA3F5450-2ADB-E611-984A-001E67A3FBAA.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/BA4F36C3-3DDB-E611-8E3E-0CC47A4C8F08.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/C0100039-A0DB-E611-96D5-0CC47A13D09C.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/C4C1E9A1-15DC-E611-8D9E-002590FD5A48.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/C85254B4-BDDB-E611-8116-24BE05C4D801.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/CEB24329-CFF0-E611-8887-02163E019E5C.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/D009CA42-A0DB-E611-B62C-001E67444EAC.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/D283E727-6DDB-E611-9DA7-001E67397CB5.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/D4394090-44DB-E611-98E8-0CC47A7C34B0.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/D497112F-61DC-E611-A473-0025905B85EC.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E66413C3-3DDB-E611-8EA6-0CC47A745250.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/E6A7C707-5CDB-E611-911D-0025905B8610.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/EAB15212-2ADD-E611-9C58-A4BF0102621F.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/F0998864-3CDC-E611-909D-FA163EF4B995.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/F601FF1C-C2F0-E611-9A3C-02163E019C74.root",
"eos/cms/store/mc/RunIISummer16MiniAODv2/HJ_MiNLO_NNLOPS_HToWWTo2L2Nu_M125_NNPDF30_13TeV_JHUGen702_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/FE693535-B0DA-E611-B6D5-B8CA3A709028.root"])


fileW=TFile("medianWeights.root")

median0=fileW.Get("hMedians0")
median9=fileW.Get("hMedians9")

ntot=0.0; 
nfid=0.0;
n2p5=0.0;

outfile = TFile("outtree.root", 'RECREATE')

tree = TTree('tree', 'tree')

rap = array( 'd', [ 0.0 ] )
wnlo = array( 'd', [ 0.0 ] )
wnnlo = array( 'd', [ 0.0 ] )
worig = array( 'd', [ 0.0 ] )
tree.Branch('rap', rap, 'rap/D')
tree.Branch('w_nlo', wnlo, 'w_nlo/D')
tree.Branch('w_nnlo', wnnlo, 'w_nnlo/D')
tree.Branch('w_orig', worig, 'w_orig/D')

for i,event in enumerate(events):

  #if i>10000: break
  print i

  event.getByLabel(prunedLabel, handlePruned)
  event.getByLabel(genLabel, gen)
  event.getByLabel(lheLabel, lhe)
  #event.getByLabel(fidLabel, fid)

  #id1=lhe.product().pdf().id.first
  #id2=lhe.product().pdf().id.second
  #x1=lhe.product().pdf().x.first
  #x2=lhe.product().pdf().x.second
  #Q=lhe.product().pdf().scalePDF
  #origPDFWeight=p1.xfxQ(id1, x1, Q)*p1.xfxQ(id2, x2, Q)
  #newPDFWeight =p2.xfxQ(id1, x1, Q)*p2.xfxQ(id2, x2, Q) 

  #strange cases
  #if origPDFWeight*newPDFWeight < 0:
  #  print "id1", id1, "id2", id2, "x1", x1, "x2", x2, "Q", Q, "originalPDF", origPDFWeight, "newPDF", newPDFWeight 

  w_nlo = lhe.product().weights()[0].wgt
  w_nnlo = lhe.product().weights()[9].wgt

  pruned = handlePruned.product()

  higgs=TLorentzVector(0,0,0,0)
  for p in pruned :
    if (p.pdgId()==25 and p.status()==62):
      higgs=TLorentzVector(p.px(),p.py(),p.pz(),p.energy())

  #badw=False
  #medianWeight0=median0.GetBinContent(median0.FindBin(higgs.Rapidity()));
  #medianWeight9=median9.GetBinContent(median9.FindBin(higgs.Rapidity()));
  #if abs((w_nlo-medianWeight0)/medianWeight0)>50:
  #  w_nlo=medianWeight0
  #if abs((w_nnlo-medianWeight9)/medianWeight9)>50:
  #  w_nnlo=medianWeight9

  h3.Fill(w_nnlo)
  h4.Fill(w_nlo)

  rap[0] = higgs.Rapidity()
  wnlo[0] = w_nlo
  wnnlo[0] = w_nnlo
  worig[0] = lhe.product().originalXWGTUP()

  w_nlo /= lhe.product().originalXWGTUP()
  w_nnlo /= lhe.product().originalXWGTUP()

  #protection for some strange cases  
  #if newPDFWeight/origPDFWeight>0 and abs(newPDFWeight/origPDFWeight) < 100 and abs(newPDFWeight/origPDFWeight)>0.01:
  #  w=w*newPDFWeight/origPDFWeight

  ntot+=w_nnlo
  #if (fid.product().passedFiducial): nfid+=w_nnlo
  if (abs(higgs.Rapidity())<2.5): n2p5+=w_nnlo

  h1.Fill(higgs.Rapidity(), w_nlo)
  if (abs(w_nnlo/w_nlo)<300.0):
    h2.Fill(higgs.Rapidity(), w_nnlo)

  tree.Fill()

outfile.cd();

h1.Write();
h2.Write();
h3.Write();
h4.Write();

outfile.Write()
outfile.Close()

    
#h1.SaveAs("rapiditynlo.root")
#h2.SaveAs("rapiditynnlo.root")
#h3.SaveAs("nnloweight.root")
#h4.SaveAs("nloweight.root")
print ntot,nfid,nfid/ntot,n2p5,n2p5/ntot
