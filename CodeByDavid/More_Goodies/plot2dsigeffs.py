from ROOT import *

from tdrStyle import *
setTDRStyle()

obsName = 'pT4l'
label = 'p_{T}(H)'
#obsName = 'massZ2'
#label = 'm(Z_{2})'

_temp = __import__('inputs_sig_'+obsName, globals(), locals(), ['eff'], -1)
eff = _temp.eff

modelNames = ['ggH_powheg15_JHUgen_125','VBF_powheg_125','WH_pythia_125','ZH_pythia_125','ttH_pythia_125']
fStates = ['4e','4mu','2e2mu']

for model in modelNames:
    for fState in fStates:
        eff2d = TH2D("eff2d", label, 4, 0, 4, 4, 0, 4)
        for x in range(0,4):
            for y in range(0,4):
                eff2d.GetXaxis().SetBinLabel(x+1,str(x))
                eff2d.GetYaxis().SetBinLabel(y+1,str(y))
                eff2d.Fill(x,y,eff[model+'_'+fState+'_'+obsName+'_genbin'+str(x)+'_recobin'+str(y)])
        c=TCanvas("c","c",1000,800)
        c.cd()
        c.SetTopMargin(0.10)
        c.SetRightMargin(0.20)
        eff2d.GetXaxis().SetTitle('gen. bin '+label+' '+fState)
        eff2d.GetYaxis().SetTitle('reco. bin '+label+' '+fState)
        eff2d.GetZaxis().SetTitle('#epsilon^{ij}')
        eff2d.GetZaxis().SetRangeUser(0.0,1.0)
        eff2d.Draw("colzTEXT0")
        latex2 = TLatex()
        latex2.SetNDC()
        latex2.SetTextSize(0.5*c.GetTopMargin())
        latex2.SetTextFont(42)
        latex2.SetTextAlign(31) # align right   
        latex2.SetTextSize(0.5*c.GetTopMargin())
        latex2.SetTextFont(62)
        latex2.SetTextAlign(11) # align right  
        latex2.DrawLatex(0.2, 0.92, "CMS")
        latex2.SetTextSize(0.4*c.GetTopMargin())
        latex2.SetTextFont(52)
        latex2.SetTextAlign(11)
        latex2.DrawLatex(0.29, 0.92, "Simulation")
        c.SaveAs("eff2d_"+model+"_"+obsName+"_"+fState+".png")
        
