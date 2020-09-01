from ROOT import *
from tdrStyle import *
setTDRStyle()

def plotAsimov_simultaneous(asimovDataModel, asimovPhysicalModel, obsName, fstate, recobin):

    if (fstate=="4mu"): channel = "1"
    if (fstate=="4e"): channel = "2"
    if (fstate=="2e2mu"): channel = "3"

    # Load some libraries                                
    ROOT.gSystem.AddIncludePath("-I$CMSSW_BASE/src/ ")
    ROOT.gSystem.Load("$CMSSW_BASE/lib/slc5_amd64_gcc472/libHiggsAnalysisCombinedLimit.so")
    ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include")
    ROOT.gSystem.AddIncludePath("-Iinclude/")
    
    RooMsgService.instance().setGlobalKillBelow(RooFit.WARNING)
    
    f = TFile(asimovDataModel+'_all_'+obsName+'_8TeV_Asimov_'+asimovPhysicalModel+'.root','READ')
    
    data = f.Get("toys/toy_asimov");
    #data->Print("");

    w = f.Get("w")

    print "loading snap shot clean"

    w.loadSnapshot("clean");

    sim_orig = w.pdf("model_s")
    pdfi_orig = sim_orig.getPdf("ch"+channel+"_ch"+recobin);

    trueHBin0_orig = w.function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_trueH"+fstate+"Bin0");
    trueHBin1_orig = w.function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_trueH"+fstate+"Bin1");
    trueHBin2_orig = w.function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_trueH"+fstate+"Bin2");
    trueHBin3_orig = w.function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_trueH"+fstate+"Bin3");
    out_trueH_orig = w.function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_out_trueH");
    fakeH_orig = w.function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_fakeH");
    ggzz_orig = w.function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_bkg_ggzz")
    qqzz_orig = w.function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_bkg_qqzz")
    zjets_orig = w.function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_bkg_zjets")

    n_trueH_orig = trueHBin0_orig->getVal()+trueHBin1_orig->getVal()+trueHBin2_orig->getVal()+trueHBin3_orig->getVal();
    n_out_trueH_orig = out_trueH_orig->getVal();
    n_fakeH_orig = fakeH_orig->getVal();
    n_zz_orig = ggzz_orig->getVal()+qqzz_orig->getVal();
    n_zjets_orig = zjets_orig->getVal();

    print "load snap shot MultiDimFit"

    w->loadSnapshot("MultiDimFit");

    cout<<"load snap shot done"<<endl;

    w->Print("");

    const RooSimultaneous *sim  = dynamic_cast<const RooSimultaneous *>(w->pdf("model_s"));
    pdfi  = sim->getPdf("ch"+channel+"_ch"+str(recobin));

    trueHBin0 = w->function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_trueH"+tchannel+"Bin0");
    trueHBin1 = w->function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_trueH"+tchannel+"Bin1");
    trueHBin2 = w->function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_trueH"+tchannel+"Bin2");
    trueHBin3 = w->function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_trueH"+tchannel+"Bin3");

    out_trueH = w->function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_out_trueH");

    fakeH = w->function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_fakeH");

    ggzz = w->function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_bkg_ggzz");
    qqzz = w->function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_bkg_qqzz");
    zjets = w->function("n_exp_final_binch"+channel+"_ch"+str(recobin)+"_proc_bkg_zjets");

    n_trueH = trueHBin0.getVal()+trueHBin1.getVal()+trueHBin2.getVal()+trueHBin3.getVal();
    n_out_trueH = out_trueH.getVal();
    n_fakeH = fakeH.getVal();
    n_zz = ggzz.getVal()+qqzz.getVal();
    n_zjets = zjets.getVal();

    CMS_zz4l_mass = w.var("CMS_zz4l_mass");
    #CMS_zz4l_mass.Print("");

    mass = w.var("CMS_zz4l_mass").frame(RooFit.Bins(17));
    data.plotOn(mass);

    sbin = "ch"+channel+"_ch"+str(recobin)

    pdfi.plotOn(mass, LineColor(kRed))
    pdfi.plotOn(mass, LineColor(kMagenta), Components("shapeBkg_bkg_zjets_"+sbin+",shapeBkg_bkg_ggzz_"+sbin+",shapeBkg_bkg_qqzz_"+sbin+",shapeBkg_fakeH_"+sbin+",shapeBkg_out_trueH_"+bin))
    pdfi.plotOn(mass, LineColor(kOrange-3), Components("shapeBkg_bkg_zjets_"+sbin+",shapeBkg_bkg_ggzz_"+sbin+",shapeBkg_bkg_qqzz_"+sbin+",shapeBkg_fakeH_"+bin))
    pdfi.plotOn(mass, LineColor(kAzure-3), Components("shapeBkg_bkg_zjets_"+sbin+",shapeBkg_bkg_ggzz_"+sbin+",shapeBkg_bkg_qqzz_"+bin))
    pdfi.plotOn(mass, LineColor(kGreen-3), Components("shapeBkg_bkg_zjets_"+bin))

    gStyle.SetOptStat(0);

    c = TCanvas("c","c",1000,800);
    c.cd();

    dummy = TH1D("","",1,105.6,140.6);
    dummy.SetBinContent(1,2);
    dummy.SetFillColor(0)
    dummy.SetLineColor(0)
    dummy.SetLineWidth(0)
    dummy.SetMarkerSize(0)
    dummy.SetMarkerColor(0); 
    dummy.GetYaxis().SetTitle("Events / (2 GeV)");
    dummy.GetXaxis()->SetTitle("m_{"+fstate.replace("mu","#mu")+"} [GeV]");
    dummy.Draw();

    dummy_trueH = TH1D("dummy_trueH","",1,105.6,140.6); 
    dummy_trueH.SetLineColor(kRed);
    dummy_outH = TH1D("dummy_outH","",1,105.6,140.6); 
    dummy_outH.SetLineColor(kMagenta);
    dummy_fakeH = TH1D("dummy_fakeH","",1,105.6,140.6); 
    dummy_fakeH.SetLineColor(kOrange-3);

    dummy_zz = TH1D("dummy_zz","",1,105.6,140.6)
    dummy_zz.SetLineColor(kAzure-3);
    dummy_zjets = TH1D("dummy_zjets","",1,105.6,140.6)
    dummy_zjets.SetLineColor(kGreen-3);

    leg = TLegend(0.7,0.3,0.9,0.7);
    leg.SetBorderSize(0);
    leg.SetFillColor(0);
    leg.AddEntry(dummy_trueH, "trueH","L");
    leg.AddEntry(dummy_outH, "outH","L");
    leg.AddEntry(dummy_fakeH, "fakeH","L");
    leg.AddEntry(dummy_zz, "ZZ","L");
    leg.AddEntry(dummy_zjets, "Z+Jets","L");

    leg.Draw("same");

    mass.Draw("same");

    TLatex *latex = new TLatex();
    latex.SetNDC();
    latex.SetTextSize(0.6*c.getTopMargin());
    latex.SetTextFont(42);
    latex.SetTextAlign(11);
    latex.DrawLatex(0.14, 0.92, "CMS Preliminary");
    latex.SetTextSize(0.6*c.getTopMargin());
    latex.SetTextFont(42);
    latex.DrawLatex(0.22, 0.85, "Asimov Data");
    latex.SetTextSize(0.4*c.getTopMargin());
    latex.SetTextFont(42);

    latex.DrawLatex(0.22, 0.70, "N_{right sig.}^{fit} = %.2f (%.2f)"%(n_trueH,n_trueH_orig))
    latex.DrawLatex(0.22, 0.65, "N_{out sig.}^{fit} = %.2f (%.2f)"%(n_out_trueH, n_out_trueH_orig))
    latex.DrawLatex(0.22, 0.60, "N_{ZZ}^{fit} = %.2f (%.2f)"%(n_zz, n_zz_orig))
    latex.DrawLatex(0.22, 0.50, "N_{Z+X}^{fit} = %.2f (%.2f)"%(n_zjets, n_zjets_orig))
    latex.DrawLatex(0.22, 0.40, "N_{wrong sig.}^{fit} = %.2f (%.2f)" %(n_fakeH, n_fakeH_orig))

    latex.SetTextSize(0.035);
    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.DrawLatex(0.87,0.92, "#sqrt{s} = 8 TeV, L=19.7 fb^{-1}");

    c.SaveAs("plots/asimovdata_"+asimovDataPhyiscalModel+"_"+asimovDataModelName+"_"+fstate+"_recobin"+str(recobin)+".pdf");

plotAsimov_simultaneous('ggH_powheg15_JHUgen_125', 'v2', 'pT4l', '4mu', 0)
