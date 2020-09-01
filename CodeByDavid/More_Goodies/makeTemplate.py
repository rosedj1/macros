from ROOT import *

f = TFile('ggH_2015MC_mH125.root','READ')
t = f.Get("Ana/passedEvents")

h1D = TH1D("h1D","h1D",160,40,120)
h1D.Sumw2()
h2D = TH2D("h2D","h2D",120,12,72,160,40,120)
h2D.Sumw2()

t.Draw("GENmassZ1>>h1D","passedFiducialSelection==1","goff")
t.Draw("GENmassZ1:GENmassZ2>>h2D","passedFiducialSelection==1","goff")

h1D.Scale(1.0/h1D.Integral())
h2D.Scale(1.0/h2D.Integral())

for x in range(h2D.GetNbinsX()):
  for y in range(h2D.GetNbinsY()):
    if h2D.GetBinContent(x,y)==0.0:
      print "here"
      h2D.SetBinContent(x,y,0.0001)

h2D.Scale(1.0/h2D.Integral())

var_mZ1 = RooRealVar("var_mZ1","var_mZ1",40,120);
var_mZ2 = RooRealVar("var_mZ2","var_mZ2",12,72);

plot = var_mZ1.frame(40,120)

GENmZ1DataHist = RooDataHist("GENmZ1DataHist","GENmZ1DataHist",RooArgList(var_mZ1),h1D);
GENmZ1DataHist.plotOn(plot,RooFit.LineColor(1))

GENmZ1vsmZ2DataHist = RooDataHist("GENmZ1vsmZ2DataHist","GENmZ1vsmZ2DataHist",RooArgList(var_mZ1,var_mZ2),h2D);
GENmZ1vsmZ2 = RooHistPdf("GENmZ1vsmZ2","GENmZ1vsmZ2",RooArgSet(var_mZ1,var_mZ2),GENmZ1vsmZ2DataHist,2);

GENmZ1vsmZ2.plotOn(plot,RooFit.LineColor(2))
plot.Draw()
c1.SaveAs("shapecompare.pdf")

print h2D.Integral()

h2D.SaveAs("mz1_vs_mz2_125.root")
