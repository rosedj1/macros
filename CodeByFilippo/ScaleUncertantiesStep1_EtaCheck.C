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

bool muon;
bool DATA;
TString fs;
TString directory;
bool verbose;


void ScaleUncertantiesStep1_EtaCheck(){

	gStyle->SetOptStat(111);
	TRandom3 rand;
	double u1;  
	TString nome_histo; 
	TFile* outfile;
	TFile* infile;
	float chosen_eta;
	float chosen_pT;
	unsigned int ketaBin, kpTBin;
	           
// 	verbose = true;
	verbose = false;
	                                                                                         
	
	muon = true;
// 	muon = false;

// 	DATA = true;
	DATA = false;
	
	Double_t pT1, eta1, phi1;
	Double_t genLep_pt1;
	Double_t pterr1;
	Double_t genLep_eta1;
	Double_t genLep_phi1;
	Double_t pT2, eta2, phi2;
	Double_t genLep_pt2;
	Double_t pterr2;
	Double_t genLep_eta2;
	Double_t genLep_phi2;
	Double_t weight;
	Int_t lep1_ecalDriven;
	Int_t lep2_ecalDriven;
	Double_t massZ;
	Double_t genzm;
	Double_t GENmass2l;
	
	std::vector<Double_t> pT_bins;
	pT_bins.push_back(5);
// 	pT_bins.push_back(20);
	pT_bins.push_back(30);
// 	pT_bins.push_back(40);
// 	pT_bins.push_back(50);
	pT_bins.push_back(60);
	pT_bins.push_back(100);

	std::vector<Double_t> eta_bins;
	eta_bins.push_back(0);
	eta_bins.push_back(0.5);
	eta_bins.push_back(0.9);
	eta_bins.push_back(1.2);
	eta_bins.push_back(1.8);
	eta_bins.push_back(2.1);
	eta_bins.push_back(2.4);

	std::vector<TString> eta_bins_name;
   	eta_bins_name.push_back("0");
	eta_bins_name.push_back("0p5");
	eta_bins_name.push_back("0p9");
	eta_bins_name.push_back("1p2");
	eta_bins_name.push_back("1p8");
	eta_bins_name.push_back("2p1");
	eta_bins_name.push_back("2p4");

    TH1F *h_mass[pT_bins.size()-1][eta_bins.size()-1];
    TH1F *h_pT[pT_bins.size()-1];
    TH1F *h_pT_tot[pT_bins.size()-1];
    TH1F *h_eta[eta_bins.size()-1];
    TH1F *h_eta_tot[eta_bins.size()-1];
    
	TLorentzVector lep_1, lep_2;
	TLorentzVector ZPrime; 
	
	TH2F* h_pt_vs_eta[4];



    for(int pt = 0; pt < pT_bins.size()-1; pt ++){
    	nome_histo = Form("pT_%d_%d",(int)pT_bins[pt], (int)pT_bins[pt+1]);
    	h_pT[pt] = new TH1F(nome_histo, nome_histo, pT_bins.size()-1, &pT_bins[0]);
    	nome_histo = Form("pT_Tot_%d_%d", (int)pT_bins[pt], (int)pT_bins[pt+1]);
    	h_pT_tot[pt] = new TH1F(nome_histo, nome_histo, pT_bins.size()-1, &pT_bins[0]);

	    for(int eta = 0; eta <  eta_bins.size()-1; eta ++){
	    	if(pt == 0){
    			nome_histo = Form("eta_%s_%s", eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
		    	h_eta[eta] = new TH1F(nome_histo, nome_histo, eta_bins.size()-1, &eta_bins[0]);
    			nome_histo = Form("eta_Tot_%s_%s", eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
		    	h_eta_tot[eta] = new TH1F(nome_histo, nome_histo,  eta_bins.size()-1, &eta_bins[0]);
	    	}
	    	nome_histo = Form("Mass_pT_%d_%d_eta_%s_%s", (int)pT_bins[pt], (int)pT_bins[pt+1], eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
	    	if(verbose) std::cout<<nome_histo<<std::endl;
	    	h_mass[pt][eta] = new TH1F(nome_histo, nome_histo, 50, 60, 110);
	    }
    }
    
//    	TH1F* HISTO_MASS;
//     TH1F* HISTO_MASS_LEP;
    
    for(int mode = 0; mode < 4; mode++){

//     	HISTO_MASS = new TH1F("HISTO_MASS", "HISTO_MASS", 50, 60, 110);
// 	    HISTO_MASS_LEP = new TH1F("HISTO_MASS_LEP", "HISTO_MASS_LEP", 50, 60, 110);
    
	    if(mode < 2) muon = true;
    	else muon = false;
	    if(mode == 0 || mode == 2) DATA = false;
	    else DATA = true;
	    
	    std::cout<<muon<<"\t"<<DATA<<std::endl;
	    	
		if(muon)
			fs = "mu";
		else
			fs = "e";
		

	    if(DATA){
			h_pt_vs_eta[mode] = new TH2F("pt_vs_eta_" + fs + "DATA", "pt_vs_eta_" + fs + "DATA", pT_bins.size()-1, &pT_bins[0], eta_bins.size()-1, &eta_bins[0]);
		    infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/makeSlimTree/Data_m2" + fs + ".root");
		}
		else{
			h_pt_vs_eta[mode] = new TH2F("pt_vs_eta_" + fs + "MC", "pt_vs_eta_" + fs + "MC", pT_bins.size()-1, &pT_bins[0], eta_bins.size()-1, &eta_bins[0]);
		    infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/makeSlimTree/DYJetsToLL_M-50_Rochester_v4_m2" + fs + ".root");
		}

		TTree* tree;                                               
		if(infile) tree = (TTree*) infile->Get("passedEvents");		
		else 
			std::cout<<"ERROR could not find the file"<<std::endl;

		tree->SetBranchAddress("massZ", &massZ);       	
		tree->SetBranchAddress("genzm", &genzm);       	
		tree->SetBranchAddress("GENmass2l", &GENmass2l);   
		tree->SetBranchAddress("pT1", &pT1); 
		tree->SetBranchAddress("eta1", &eta1);   
		tree->SetBranchAddress("phi1", &phi1);   
		tree->SetBranchAddress("pterr1", &pterr1);   
		tree->SetBranchAddress("genLep_pt1", &genLep_pt1);   
		tree->SetBranchAddress("genLep_eta1", &genLep_eta1);       		
		tree->SetBranchAddress("genLep_phi1", &genLep_phi1);       		
		tree->SetBranchAddress("pT2", &pT2);   
		tree->SetBranchAddress("eta2", &eta2);   
		tree->SetBranchAddress("phi2", &phi2);   
		tree->SetBranchAddress("pterr2", &pterr2);   
		tree->SetBranchAddress("genLep_pt2", &genLep_pt2);   
		tree->SetBranchAddress("genLep_eta2", &genLep_eta2);  
		tree->SetBranchAddress("genLep_phi2", &genLep_phi2);  
		tree->SetBranchAddress("weight", &weight);
		
		Long64_t nentries = tree->GetEntries();
		std::cout<<nentries<<std::endl;		

		for(int i = 0; i < nentries; i++){
// 		for(int i = 0; i < (int)nentries/10; i++){
// 		for(int i = 0; i < 1000; i++){

			tree->GetEntry(i);                                   
			if(i % 1000000 == 0){  
				if(DATA) 	std::cout<<i<<" --- Dentro il TREE 2018 --- DATA -------------- "<<std::endl;        
				if(muon)	std::cout<<i<<" --- Dentro il TREE 2018 --- MUON"<<std::endl; 
				else 
					std::cout<<i<<" --- Dentro il TREE 2018 --- ELECTRON"<<std::endl; 
			}
		
			if(massZ < 60 || massZ > 110) continue;
				
			rand.SetSeed(abs(static_cast<int>(sin((phi1+phi2)/2)*100000)));                                                               
			u1 = rand.Uniform(1.);
		
			if(u1 < 0.5){
				chosen_eta = eta1;
				chosen_pT = pT1;
			}
			else{
				chosen_eta = eta2;
				chosen_pT = pT2;
			}
		
			if(verbose) std::cout<<u1<<"\t"<<"\teta = "<<chosen_eta<<"\teta1= = "<<eta1<<"\teta2="<<eta2<<
									"\tpt = "<<chosen_pT<<"\tpt1 = "<<pT1<<"\tpt2 = "<<pT2<<std::endl;
		
			for (unsigned int kbin=1; kbin < eta_bins.size(); ++kbin) {
				if (fabs(chosen_eta) < eta_bins.at(kbin)) {
					ketaBin = kbin-1;
					break;
				}
			}		

			for (unsigned int kbin=1; kbin < pT_bins.size(); ++kbin) {
				if (chosen_pT < pT_bins.at(kbin)) {
					kpTBin = kbin-1;
					break;
				}
			}	
			
			h_pt_vs_eta[mode]->Fill(chosen_pT, chosen_eta);
				
		   	if(verbose) std::cout<<"Bins: pT= "<<kpTBin<<"\teta = "<<ketaBin<<"\t\tpt = "<<pT_bins[kpTBin]<<"\teta = "<<eta_bins[ketaBin]<<std::endl;
	   		
		   	h_mass[kpTBin][ketaBin]->Fill(massZ);
		   	h_pT[kpTBin]->Fill(chosen_pT);
	   		h_pT_tot[kpTBin]->Fill(pT1);
		   	h_pT_tot[kpTBin]->Fill(pT2);
		   	h_eta[ketaBin]->Fill(chosen_eta);
		   	h_eta_tot[ketaBin]->Fill(eta1);
		   	h_eta_tot[ketaBin]->Fill(eta2);

		   	if(verbose) std::cout<<" ----------------- "<<std::endl;	
	   	
// 		lep_1.SetPtEtaPhiM(pT1, eta1, phi1, 0.105);
// 		lep_2.SetPtEtaPhiM(pT2, eta2, phi2, 0.105);
// 		ZPrime = lep_1 + lep_2;
// 		
// 		HISTO_MASS->Fill(massZ);
// 		HISTO_MASS_LEP->Fill(ZPrime.M());
		
		   	if(verbose) std::cout<<massZ<<"\t"<<ZPrime.M()<<std::endl;
	

		}

	    if(DATA)
		    outfile = new TFile("OutFile_m2" + fs + "_data_EtaCheck.root", "RECREATE");
		else
		    outfile = new TFile("OutFile_m2" + fs + "_MC_EtaCheck.root", "RECREATE");

		if ( outfile->IsOpen() ) printf("File opened successfully\n");
		else {
			std::cout<<"Problem in opening output file."<<std::endl;
			return;
		}

	    for(int pt = 0; pt < pT_bins.size()-1; pt ++){
// 	    	nome_histo = Form("pT_%d_%d", (int)pT_bins[pt], (int)pT_bins[pt+1]);
//  	   	h_pT[pt]->Write(nome_histo);
    		nome_histo = Form("pT_Tot_%d_%d", (int)pT_bins[pt], (int)pT_bins[pt+1]);
	    	h_pT_tot[pt]->Write(nome_histo);
        
		    for(int eta = 0;eta <  eta_bins.size()-1; eta ++){
	    	if(pt == 0){
//     			nome_histo = Form("eta_%s_%s", eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
// 		    	h_eta[eta]->Write(nome_histo);
    			nome_histo = Form("eta_ToT_a%s_%s", eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
		    	h_eta_tot[eta]->Write(nome_histo);
		    }
		    	nome_histo = Form("Mass_pT_%d_%d_eta_%s_%s", (int)pT_bins[pt], (int)pT_bins[pt+1], eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
		    	h_mass[pt][eta]->Write(nome_histo);
		    	h_mass[pt][eta]->Reset("ICESM");
		    }
	    }
    
	    outfile->Close();
    
    } // for on fs and DATA
    

//     for(int h = 0; h < 4; h++){
//     	TCanvas* c1 = new TCanvas("pt_vs_eta", "pt_vs_eta", 700, 500);
//     	h_pt_vs_eta[h]->Draw("COLZ");
//     	if(h == 0) c1->Print("pt_vs_eta_EtaCheck.pdf[");
//     	c1->Print("pt_vs_eta_EtaCheck.pdf");
//     	if(h == 3) c1->Print("pt_vs_eta_EtaCheck.pdf]");
//     }


}