from ROOT import *

f = TFile('DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1.root',"READ")
t = f.Get("passedEvents")

nummc = TH1D("nummc","nummc",50,0.0,50.0)
denmc = TH1D("denmc","denmc",50,0.0,50.0)

numtnp = TH1D("numtnp","numtnp",50,0.0,50.0)
dentnp = TH1D("dentnp","dentnp",50,0.0,50.0)

t.Draw("pTl2>>denmc","(VlepId==11 && pTl1>50.0)","goff")
t.Draw("pTl2>>nummc","(VlepId==11 && pTl1>50.0 && passTrig==1)","goff")
effmc = TEfficiency(nummc,denmc)
effmc.SetStatisticOption(TEfficiency.kBUniform)

t.Draw("pTl2>>dentnp","(VlepId==11 && pTl1>50.0)","goff")
t.Draw("pTl2>>numtnp","mceff*(VlepId==11 && pTl1>50.0)","goff")
efftnp = TEfficiency(numtnp,dentnp)
efftnp.SetStatisticOption(TEfficiency.kBUniform)

effmc.Draw("AP")
efftnp.SetMarkerColor(4)
efftnp.Draw("PSame")

c1.SaveAs("test.pdf")
