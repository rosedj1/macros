#include "TFile.h"
#include <iostream>
#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include "TPad.h"

using namespace RooFit;

void MergeHisto(TString fs, bool DATA, TString reconstruction, TString year){

	  if(fs!="2e" && fs!="2mu"){
  		cout<<"fs has to be 2e, or 2mu"<<endl;
	  	return;
	  }	
	
	TString fs_dir;

	std::vector<Double_t> pT_bins;
	if(fs=="2mu"){
		pT_bins.push_back(5);
		pT_bins.push_back(20);
		pT_bins.push_back(30);
		pT_bins.push_back(40);
		pT_bins.push_back(50);
		pT_bins.push_back(60);
		pT_bins.push_back(100);
	}
	else{
		pT_bins.push_back(7);
		pT_bins.push_back(20);
		pT_bins.push_back(30);
		pT_bins.push_back(40);
		pT_bins.push_back(50);
		pT_bins.push_back(60);
		pT_bins.push_back(100);
	}

	std::vector<Double_t> eta_bins;
	eta_bins.push_back(0);
	eta_bins.push_back(0.5);
	if(fs=="2mu"){
		eta_bins.push_back(0.9);	
		eta_bins.push_back(1.4);
		eta_bins.push_back(1.8);
		eta_bins.push_back(2.1);
		eta_bins.push_back(2.4);
	}
	else{
		eta_bins.push_back(0.8);	
		eta_bins.push_back(1.2);
		eta_bins.push_back(1.5);
		eta_bins.push_back(2.2);
		eta_bins.push_back(2.5);
	}


	std::vector<TString> eta_bins_name;
	eta_bins_name.push_back("0");
	eta_bins_name.push_back("0p5");
	if(fs=="2mu"){
		eta_bins_name.push_back("0p9");
		eta_bins_name.push_back("1p4");
		eta_bins_name.push_back("1p8");
		eta_bins_name.push_back("2p1");
		eta_bins_name.push_back("2p4");
	}
	else{
		eta_bins_name.push_back("0p8");
		eta_bins_name.push_back("1p2");
		eta_bins_name.push_back("1p5");
		eta_bins_name.push_back("2p2");
		eta_bins_name.push_back("2p5");
	}
	
	TString histo_name;


	TH1F *full = new TH1F("full", "full", 1000, 60, 120);
	TH1F *full_pt[pT_bins.size()-1];
	TH1F *full_eta[eta_bins.size()-1];	


    TH1F *h_mass[pT_bins.size()-1][3];
    for(int pt = 0; pt < pT_bins.size()-1; pt ++){
	    	histo_name = Form("Mass_pT_%d_%d_eta_0_0p9", (int)pT_bins[pt], (int)pT_bins[pt+1]);
	    	h_mass[pt][0] = new TH1F(histo_name, histo_name, 1000, 60, 120);
	    	histo_name = Form("Mass_pT_%d_%d_eta_0p9_1p8", (int)pT_bins[pt], (int)pT_bins[pt+1]);
	    	h_mass[pt][1] = new TH1F(histo_name, histo_name, 1000, 60, 120);
	    	histo_name = Form("Mass_pT_%d_%d_eta_1p8_2p4", (int)pT_bins[pt], (int)pT_bins[pt+1]);
	    	h_mass[pt][2] = new TH1F(histo_name, histo_name, 1000, 60, 120);
    }	

	for(int pt = 0; pt < pT_bins.size()-1; pt++){
		histo_name = Form("Mass_pT_%d_%d", (int)pT_bins[pt], (int)pT_bins[pt+1]);
		full_pt[pt] = new TH1F(histo_name, histo_name, 1000, 60, 120);			
 	}

	for(int eta = 0; eta < eta_bins.size()-1; eta++){
		histo_name = Form("Mass_eta_%s_%s", eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
		full_eta[eta] = new TH1F(histo_name, histo_name, 1000, 60, 120);
 	}
	
	
	TFile* infile;
	TFile* outfile;
	
	if(fs=="2mu")  fs_dir = "Muon";
	else fs_dir = "Electron";
	infile = new TFile("./Diffential_Scale_Root/" + year + "/" + fs_dir + "/OutFile_m" + fs + "_" + year + "_" + reconstruction + "_Differential_v3.root");

	if(!infile) std::cout<<"Problem in opening the file!"<<std::endl;
		
	for(int pt = 0; pt < pT_bins.size()-1; pt++){
		for(int eta = 0; eta < eta_bins.size()-1; eta++){

				histo_name = Form("Mass_pT_%d_%d_eta_%s_%s", (int)pT_bins[pt], (int)pT_bins[pt+1], eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
				TH1F* histo = (TH1F*)gDirectory->Get(histo_name);
				full->Add(histo,1);
				full_pt[pt]->Add(histo,1);
				full_eta[eta]->Add(histo,1);
				if(eta == 0 || eta == 1)
					h_mass[pt][0]->Add(histo,1);
				if(eta == 2 || eta == 3)
					h_mass[pt][1]->Add(histo,1);
				if(eta == 4 || eta == 5)
					h_mass[pt][2]->Add(histo,1);
		}
	}


	outfile = new TFile("./Final_Scale_Root/" + year + "/" + fs_dir + "/OutFile_m" + fs + "_" + year + "_" + reconstruction + "_Final.root", "RECREATE");

	full->Write();	

	for(int pt = 0; pt < pT_bins.size()-1; pt++){
		full_pt[pt]->Write();
		for(int eta = 0; eta < 3; eta++)	h_mass[pt][eta]->Write();
	}
			
	for(int eta = 0; eta < eta_bins.size()-1; eta++)
		full_eta[eta]->Write();	
				
	infile->Close();
	

}