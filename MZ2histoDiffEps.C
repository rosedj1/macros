#include "TTree.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"

using namespace std;

void MZ2histoDiffEps() {

//string branch = "massZ2";

TFile* f1 = new TFile("cmsgrid_final_eps1e-2.lhesyntax2.root");
TTree* t1 = (TTree*)f1->Get("lheEvents_tchan");
t1->SetLineColor(kBlack);
t1->Draw("massZ2>>h1(90,19.8,20.2)","");

TFile* f2 = new TFile("cmsgrid_final_eps1e-3.lhesyntax2.root");
TTree* t2 = (TTree*)f2->Get("lheEvents_tchan");
t2->SetLineColor(kRed);
t2->Draw("massZ2>>h2(90,19.8,20.2)","","same");

TFile* f3 = new TFile("cmsgrid_final_eps1e-4.lhesyntax2.root");
TTree* t3 = (TTree*)f3->Get("lheEvents_tchan");
t3->SetLineColor(kGreen);
t3->Draw("massZ2>>h3(90,19.8,20.2)","","same");

TFile* f4 = new TFile("cmsgrid_final_eps1p5e-2.lhesyntax2.root");
TTree* t4 = (TTree*)f4->Get("lheEvents_tchan");
t4->SetLineColor(kBlue);
t4->Draw("massZ2>>h4(90,19.8,20.2)","","same");

TFile* f5 = new TFile("cmsgrid_final_eps5e-3.lhesyntax2.root");
TTree* t5 = (TTree*)f5->Get("lheEvents_tchan");
t5->SetLineColor(kViolet);
t5->Draw("massZ2>>h5(90,19.8,20.2)","","same");

// Build a legend
TLegend leg(0.1, 0.7, 0.3, 0.9, "Epsilon Values");
leg.SetFillColor(0);
leg.AddEntry(t1, "1e-2");
leg.Draw();
}
