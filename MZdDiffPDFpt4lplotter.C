#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"

using namespace std;

void MZdDiffPDFpt4lplotter() {

//string branch = "mass4l_tchan";

TFile* f1 = new TFile("cmsgrid_final_eps1e-2_mzd20_lhaid262000.lhesyntax2.root");
TTree* t1 = (TTree*)f1->Get("lheEvents_tchan");
t1->SetLineColor(kBlack);
t1->Draw("pT4l>>pT4l(200, 124.8, 125.2)","","");

TFile* f2 = new TFile("cmsgrid_final_eps1e-2_mzd20_lhaid306000.lhesyntax2.root");
TTree* t2 = (TTree*)f2->Get("lheEvents_tchan");
t2->SetLineColor(kRed);
t2->Draw("pT4l>>h2(200, 124.8, 125.2)","", "same");

// Build a legend
TLegend leg(0.1, 0.7, 0.3, 0.9, "Epsilon Values");
leg.SetFillColor(0);
leg.AddEntry(t1, "1e-2");
leg.Draw();
}
