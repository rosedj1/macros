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

void Norm(TH2F* h1);
void Save(TH2F* h1, TString nome);

TFile* infile_MU;
TFile* infile_ELE;
TTree* tree_MU;                                               
TTree* tree_ELE;                                               

// std::vector<TString> year_name;
TString year_name;

TString directory_year;
TString reconstruction = "madgraph";
//TString reconstruction = "amcatnlo";
TString directory = "2D_plot_" + reconstruction;

void TwoDMap() {

	TString execute = "mkdir " + directory;
    std::cout<<execute<<std::endl;
	gSystem->Exec(execute);
	directory += "/";

	gStyle->SetOptStat(0);
	
	TH1F *h_mass = new TH1F("massZ", "massZ", 60, 60, 120);

	TH1F *h_0_09 = new TH1F("0_09", "0_09", 48, -2.4, 2.4);
	TH1F *h_09_18 = new TH1F("09_18", "09_18", 48, -2.4, 2.4);
	TH1F *h_18_24 = new TH1F("18_24", "18_24", 48, -2.4, 2.4);

	TH1F *h_ELE_0_08 = new TH1F("ELE_0_08", "ELE_0_08", 48, -2.4, 2.4);
	TH1F *h_ELE_08_1 = new TH1F("ELE_08_1", "ELE_08_1", 48, -2.4, 2.4);
	TH1F *h_ELE_1_1p44 = new TH1F("ELE_1_1p44", "ELE_1_1p44", 48, -2.4, 2.4);
	TH1F *h_ELE_1p44_1p57 = new TH1F("ELE_1p44_1p57", "ELE_1p44_1p57", 48, -2.4, 2.4);
	TH1F *h_ELE_1p57_2 = new TH1F("ELE_1p57_2", "ELE_1p57_2", 48, -2.4, 2.4);
	TH1F *h_ELE_2_25 = new TH1F("ELE_2_25", "ELE_2_25", 48, -2.4, 2.4);

	TH2F *ErrPt_vs_eta_MU = new TH2F("2D_Map_vs_eta", "2D_Map_vs_eta", 100, -2.4, 2.4, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_0_09_MU = new TH2F("2D_Map_vs_pt_0_09", "2D_Map_vs_pt_0_09", 100, 5, 100, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_09_18_MU = new TH2F("2D_Map_vs_pt_09_18", "2D_Map_vs_pt_09_18", 100, 5, 100, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_18_24_MU = new TH2F("2D_Map_vs_pt_18_24", "2D_Map_vs_pt_18_24", 100, 5, 100, 100, 0.0, 0.1);

	TH2F *ErrPt_vs_pt_0_08_ELE = new TH2F("2D_Map_vs_pt_0_08", "2D_Map_vs_pt_0_08", 100, 0.0, 100, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_08_1_ELE = new TH2F("2D_Map_vs_pt_08_1", "2D_Map_vs_pt_08_1", 100, 0.0, 100, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_1_1p44_ELE = new TH2F("2D_Map_vs_pt_1p44", "2D_Map_vs_pt_1p44", 100, 0.0, 100, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_1p44_1p57_ELE = new TH2F("2D_Map_vs_pt_1p44_1p57", "2D_Map_vs_pt_1p44_1p57", 100, 0.0, 100, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_1p57_2_ELE = new TH2F("2D_Map_vs_pt_1p57_2", "2D_Map_vs_pt_1p57_2", 100, 0.0, 100, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_2_25_ELE = new TH2F("2D_Map_vs_pt_2_25", "2D_Map_vs_pt_2_25", 100, 0.0, 100, 100, 0.0, 0.1);

	TH2F *ErrPt_vs_eta_ECAL = new TH2F("2D_Map_vs_eta", "2D_Map_vs_eta", 100, -2.4, 2.4, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_Less1_ECAL = new TH2F("ErrPt_vs_pt_Less1_ECAL", "ErrPt_vs_pt_Less1_ECAL", 100, 0.0, 100, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_Greater1_ECAL = new TH2F("ErrPt_vs_pt_Greater1_ECAL", "ErrPt_vs_pt_Greater1_ECAL", 100, 0.0, 100, 100, 0.0, 0.1);

	TH2F *ErrPt_vs_eta_TRACKER = new TH2F("2D_Map_vs_eta", "2D_Map_vs_eta", 100, -2.4, 2.4, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_0_14_TRACKER = new TH2F("ErrPt_vs_pt_0_14_TRACKER", "ErrPt_vs_pt_0_14_TRACKER", 100, 0.0, 100, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_14_16_TRACKER = new TH2F("ErrPt_vs_pt_14_16_TRACKER", "ErrPt_vs_pt_14_16_TRACKER", 100, 0.0, 100, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_16_2_TRACKER = new TH2F("ErrPt_vs_pt_16_2_TRACKER", "ErrPt_vs_pt_16_2_TRACKER", 100, 0.0, 100, 100, 0.0, 0.1);
	TH2F *ErrPt_vs_pt_2_25_TRACKER = new TH2F("ErrPt_vs_pt_2_25_TRACKER", "ErrPt_vs_pt_2_25_TRACKER", 100, 0.0, 100, 100, 0.0, 0.1);



// 	ErrPt_vs_eta_MU->GetYaxis()->SetTitle("p_{T}Err/pT_{GEN}");
// 	ErrPt_vs_eta_MU->GetXaxis()->SetTitle("eta_{GEN}");
// 	ErrPt_vs_pt_MU->GetYaxis()->SetTitle("p_{T}Err/pT_{GEN}");
// 	ErrPt_vs_pt_MU->GetXaxis()->SetTitle("pt_{GEN}");


	Double_t pT1, eta1;
	Double_t genLep_pt1;
	Double_t pterr1;
	Double_t genLep_eta1;
	Double_t pT2, eta2;
	Double_t genLep_pt2;
	Double_t pterr2;
	Double_t genLep_eta2;
	Double_t weight;
	Int_t lep1_ecalDriven;
	Int_t lep2_ecalDriven;
	Double_t massZ;
	
// 	year_name.push_back("2016");
// 	year_name.push_back("2017");
// 	year_name.push_back("2018");
	
	for(int year = 0; year < 3; year++){
	
		if(year == 0){
			infile_MU = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/makeSlimTree/DYJetsToLL_M-50_Rochester_v4_m2mu.root");
			infile_ELE = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/makeSlimTree/DYJetsToLL_M-50_Rochester_v4_m2e.root");
			year_name = "2016";
		}
		else if(year == 1){
			infile_MU = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/" + reconstruction + "/DYJetsToLL_M-50_Full_RunII_" + reconstruction + "_m2mu_2017.root");
			infile_ELE = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/" + reconstruction + "/DYJetsToLL_M-50_Full_RunII_" + reconstruction + "_m2e_2017.root");
			year_name = "2017";
		}
		else{
            infile_MU = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/" + reconstruction + "/DYJetsToLL_M-50_Full_RunII_" + reconstruction + "_m2mu_2018.root");
                        infile_ELE = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/" + reconstruction + "/DYJetsToLL_M-50_Full_RunII_" + reconstruction + "_m2e_2018.root");
            year_name = "2018";
		}

		if(infile_MU) tree_MU = (TTree*) infile_MU->Get("passedEvents");		
		else std::cout<<"ERROR could not find the MUON file"<<std::endl;
		if(infile_ELE) tree_ELE = (TTree*) infile_ELE->Get("passedEvents");		
		else std::cout<<"ERROR could not find the ELECTRON file"<<std::endl;

		tree_MU->SetBranchAddress("pT1", &pT1);   
		tree_MU->SetBranchAddress("eta1", &eta1);   
		tree_MU->SetBranchAddress("genLep_pt1", &genLep_pt1);   
		tree_MU->SetBranchAddress("pterr1", &pterr1);   
		tree_MU->SetBranchAddress("genLep_eta1", &genLep_eta1);       	
		tree_MU->SetBranchAddress("pT2", &pT2);   
		tree_MU->SetBranchAddress("eta2", &eta2);   
		tree_MU->SetBranchAddress("genLep_pt2", &genLep_pt2);   
		tree_MU->SetBranchAddress("pterr2", &pterr2);   
		tree_MU->SetBranchAddress("genLep_eta2", &genLep_eta2);  
		tree_MU->SetBranchAddress("massZ", &massZ);       	
		tree_MU->SetBranchAddress("weight", &weight);       	

		tree_ELE->SetBranchAddress("pT1", &pT1);   
		tree_ELE->SetBranchAddress("eta1", &eta1);   
		tree_ELE->SetBranchAddress("genLep_pt1", &genLep_pt1);   
		tree_ELE->SetBranchAddress("pterr1", &pterr1);   
		tree_ELE->SetBranchAddress("genLep_eta1", &genLep_eta1);       	
		tree_ELE->SetBranchAddress("pT2", &pT2);   
		tree_ELE->SetBranchAddress("eta2", &eta2);   
		tree_ELE->SetBranchAddress("genLep_pt2", &genLep_pt2);   
		tree_ELE->SetBranchAddress("pterr2", &pterr2);   
		tree_ELE->SetBranchAddress("genLep_eta2", &genLep_eta2);  
		tree_ELE->SetBranchAddress("massZ", &massZ);       	
		tree_ELE->SetBranchAddress("lep1_ecalDriven", &lep1_ecalDriven);       	
		tree_ELE->SetBranchAddress("lep2_ecalDriven", &lep2_ecalDriven);       	
		tree_ELE->SetBranchAddress("weight", &weight);       	
	
		Long64_t nentries_MU = tree_MU->GetEntries();
		std::cout<<nentries_MU<<std::endl;	

        if(year == 0) continue;

 		for(int i = 0; i < nentries_MU; i++){
//		for(int i = 0; i < (int)nentries_MU/100; i++){
		
			tree_MU->GetEntry(i);                                   
			if(i % 1000000 == 0)         
				std::cout<<i<<" --- Dentro il TREE --- MUON "<<year<<std::endl; 
			
// 		if(pT1 < 5 || pT2 < 5) std::cout<<"OK APPLICA il taglio ai reco"<<std::endl;
// 		 || genLep_pt2 < 5) std::cout<<"OK APPLICA il taglio ai GEN"<<std::endl;
		
			if(massZ < 80 || massZ > 100) continue;
		
// 		h_mass->Fill(massZ);

			ErrPt_vs_eta_MU->Fill(eta1, pterr1/pT1, weight);
			ErrPt_vs_eta_MU->Fill(eta2, pterr2/pT2, weight);
// 		if(pterr1/pT1 < 0.004) std::cout<<"Mu_1 = "<<eta1<<std::endl;
// 		if(pterr2/pT2 < 0.004) std::cout<<"Mu_2 = "<<eta2<<std::endl;

			if(abs(eta1) < 0.9 && abs(eta2) < 0.9){		
// 		if(fabs(eta1) < 0.9)
					ErrPt_vs_pt_0_09_MU->Fill(pT1, pterr1/pT1, weight);
// 		if(fabs(eta2) < 0.9)
					ErrPt_vs_pt_0_09_MU->Fill(pT2, pterr2/pT2, weight);
// 				if(pT1 < 10) std::cout<<"First muon = "<<pterr1/genLep_pt1<<"\t"<<fabs(eta1)<<std::endl;
// 				if(pT2 < 10) std::cout<<"Second muon = "<<pterr2/genLep_pt2<<"\t"<<fabs(eta2)<<std::endl;
// 			
// 			h_0_09->Fill(genLep_eta1);
// 			h_0_09->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 0.9 && fabs(genLep_eta2) > 0.9 && 
					fabs(genLep_eta1) < 1.8 && fabs(genLep_eta2) < 1.8){

				ErrPt_vs_pt_09_18_MU->Fill(pT1, pterr1/pT1, weight);
				ErrPt_vs_pt_09_18_MU->Fill(pT2, pterr2/pT2, weight);

				h_09_18->Fill(genLep_eta1);
				h_09_18->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 1.8 && fabs(genLep_eta2) > 1.8 && 
				fabs(genLep_eta1) < 2.4 && fabs(genLep_eta2) < 2.4){

				ErrPt_vs_pt_18_24_MU->Fill(pT1, pterr1/pT1, weight);
				ErrPt_vs_pt_18_24_MU->Fill(pT2, pterr2/pT2, weight);

				h_18_24->Fill(genLep_eta1);
				h_18_24->Fill(genLep_eta2);
			}
		
		}
	
// 	TCanvas *c1 = new TCanvas("c1", "c1", 500 ,500);
// 	h_mass->Draw();
// 	
// 	return;


		Long64_t nentries_ELE = tree_ELE->GetEntries();
		std::cout<<nentries_ELE<<std::endl;	
 		for(int i = 0; i < nentries_ELE; i++){
//   		for(int i = 0; i < (int)nentries_ELE/100; i++){

			tree_ELE->GetEntry(i);                                   
			if(i % 1000000 == 0)         
				std::cout<<i<<" --- Dentro il TREE --- ELECTRON "<<year<<std::endl; 

			if(massZ < 60 || massZ > 120) continue;
		
			if(lep1_ecalDriven == 1)
				ErrPt_vs_eta_ECAL->Fill(eta1, pterr1/pT1, weight);
			else
				ErrPt_vs_eta_TRACKER->Fill(eta1, pterr1/pT1, weight);
		
			if(lep2_ecalDriven == 1)
				ErrPt_vs_eta_ECAL->Fill(eta2, pterr2/pT2, weight);
			else
				ErrPt_vs_eta_TRACKER->Fill(eta2, pterr2/pT2, weight);

			if(fabs(genLep_eta1) > 0.0 && fabs(genLep_eta2) > 0.0 && 
					fabs(genLep_eta1) < 0.8 && fabs(genLep_eta2) < 0.8){
									ErrPt_vs_pt_0_08_ELE->Fill(pT1, pterr1/pT1, weight);
									ErrPt_vs_pt_0_08_ELE->Fill(pT2, pterr1/pT1, weight);
									h_ELE_0_08->Fill(genLep_eta1);
									h_ELE_0_08->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 0.8 && fabs(genLep_eta2) > 0.8 && 
					fabs(genLep_eta1) < 1. && fabs(genLep_eta2) < 1.){
									ErrPt_vs_pt_08_1_ELE->Fill(pT1, pterr1/pT1, weight);
									ErrPt_vs_pt_08_1_ELE->Fill(pT2, pterr1/pT1, weight);
									h_ELE_08_1->Fill(genLep_eta1);
									h_ELE_08_1->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 1. && fabs(genLep_eta2) > 1. && 
					fabs(genLep_eta1) < 1.44 && fabs(genLep_eta2) < 1.44){
									ErrPt_vs_pt_1_1p44_ELE->Fill(pT1, pterr1/pT1, weight);
									ErrPt_vs_pt_1_1p44_ELE->Fill(pT2, pterr1/pT1, weight);
									h_ELE_1_1p44->Fill(genLep_eta1);
									h_ELE_1_1p44->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 1.44 && fabs(genLep_eta2) > 1.44 && 
					fabs(genLep_eta1) < 1.57 && fabs(genLep_eta2) < 1.57){
									ErrPt_vs_pt_1p44_1p57_ELE->Fill(pT1, pterr1/pT1, weight);
									ErrPt_vs_pt_1p44_1p57_ELE->Fill(pT2, pterr1/pT1, weight);
									h_ELE_1p44_1p57->Fill(genLep_eta1);
									h_ELE_1p44_1p57->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 1.57 && fabs(genLep_eta2) > 1.57 && 
					fabs(genLep_eta1) < 2. && fabs(genLep_eta2) < 2){
									ErrPt_vs_pt_1p57_2_ELE->Fill(pT1, pterr1/pT1, weight);
									ErrPt_vs_pt_1p57_2_ELE->Fill(pT2, pterr1/pT1, weight);
									h_ELE_1p57_2->Fill(genLep_eta1);
									h_ELE_1p57_2->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 2. && fabs(genLep_eta2) > 2. && 
					fabs(genLep_eta1) < 2.5 && fabs(genLep_eta2) < 2.5){
									ErrPt_vs_pt_2_25_ELE->Fill(pT1, pterr1/pT1, weight);
									ErrPt_vs_pt_2_25_ELE->Fill(pT2, pterr1/pT1, weight);
									h_ELE_2_25->Fill(genLep_eta1);
									h_ELE_2_25->Fill(genLep_eta2);
			}
		
			if(massZ > 60 && massZ < 120){
				if(!lep1_ecalDriven && !lep2_ecalDriven){
					if(fabs(eta1) < 1.44 && fabs(eta2) < 1.44){
						if(pT1 < 100 && pT2 < 100){
							if(pT1 < 100 && pT2 < 100)
								h_mass->Fill(massZ);
						}
					}
				}
			}




		}

// 		TCanvas *c_pt_0_09_MU = new TCanvas("c_pt_0_09_MU", "c_pt_0_09_MU", 500, 500);
// 		h_mass->Draw();
// 	
// 		return;
// 		TCanvas *c_pt_09_18_MU = new TCanvas("c_pt_09_18_MU", "c_pt_09_18_MU", 500, 500);
// 		h_09_18->Draw();
// 		TCanvas *c_pt_18_24_MU = new TCanvas("c_pt_18_24_MU", "c_pt_18_24_MU", 500, 500);
// 		h_18_24->Draw();
// 		TCanvas *c_eta_0_08_MU = new TCanvas("c_eta_0_08_MU", "c_eta_0_08_MU", 500, 500);
// 		h_ELE_0_08->Draw();
// 		TCanvas *c_eta_08_1_MU = new TCanvas("c_eta_08_1_MU", "c_eta_08_1_MU", 500, 500);
// 		h_ELE_08_1->Draw();
// 		TCanvas *c_eta_1_1p44_MU = new TCanvas("c_eta_1_1p44_MU", "c_eta_1_1p44_MU", 500, 500);
// 		h_ELE_1_1p44->Draw();
// 		TCanvas *c_eta_1p44_1p57_MU = new TCanvas("c_eta_1p44_1p57_MU", "c_eta_1p44_1p57_MU", 500, 500);
// 		h_ELE_1p44_1p57->Draw();
// 		TCanvas *c_eta_1p57_2_MU = new TCanvas("c_eta_1p57_2_MU", "c_eta_1p57_2_MU", 500, 500);
// 		h_ELE_1p57_2->Draw();
// 		TCanvas *c_eta_2_25_MU = new TCanvas("c_eta_2_25_MU", "c_eta_2_25_MU", 500, 500);
// 		h_ELE_2_25->Draw();
	
// 		return;

		gSystem->Exec("cd " + directory);
		directory_year = year_name;
		std::cout<<directory_year<<std::endl;
		execute = "mkdir " + directory + directory_year;
		std::cout<<execute<<std::endl;
		gSystem->Exec(execute);
		directory_year += "/";
		gSystem->Exec("cd ..");
		

// 	
		Norm(ErrPt_vs_eta_MU);
		Norm(ErrPt_vs_pt_0_09_MU);
		Norm(ErrPt_vs_pt_09_18_MU);
		Norm(ErrPt_vs_pt_18_24_MU);
		Norm(ErrPt_vs_eta_ECAL);
		Norm(ErrPt_vs_eta_TRACKER);
		Norm(ErrPt_vs_pt_0_08_ELE);
		Norm(ErrPt_vs_pt_08_1_ELE);
		Norm(ErrPt_vs_pt_1_1p44_ELE);
		Norm(ErrPt_vs_pt_1p44_1p57_ELE);
		Norm(ErrPt_vs_pt_1p57_2_ELE);
		Norm(ErrPt_vs_pt_2_25_ELE);

        //gROOT->Reset();
        gROOT->SetBatch();

		Save(ErrPt_vs_eta_MU, "ErrPt_vs_eta_MU");
		Save(ErrPt_vs_pt_0_09_MU, "ErrPt_vs_pt_0_09_MU");
		Save(ErrPt_vs_pt_09_18_MU, "ErrPt_vs_pt_09_18_MU");
		Save(ErrPt_vs_pt_18_24_MU, "ErrPt_vs_pt_18_24_MU");
		Save(ErrPt_vs_eta_ECAL, "ErrPt_vs_eta_ECAL");
		Save(ErrPt_vs_eta_TRACKER, "ErrPt_vs_eta_TRACKER");
		Save(ErrPt_vs_pt_0_08_ELE, "ErrPt_vs_pt_0_08_ELE");
		Save(ErrPt_vs_pt_08_1_ELE, "ErrPt_vs_pt_08_1_ELE");
		Save(ErrPt_vs_pt_1_1p44_ELE, "ErrPt_vs_pt_1_1p44_ELE");
		Save(ErrPt_vs_pt_1p44_1p57_ELE, "ErrPt_vs_pt_1p44_1p57_ELE");
		Save(ErrPt_vs_pt_1p57_2_ELE, "ErrPt_vs_pt_1p57_2_ELE");
		Save(ErrPt_vs_pt_2_25_ELE, "ErrPt_vs_pt_2_25_ELE");
	}	
                                                                                            
}                                                   




void Norm(TH2F* h1){
	Double_t norm = h1->GetEntries();
	h1->Scale(1/norm);
}

void Save(TH2F* h1, TString name){
	TCanvas *c1 = new TCanvas("c1", "c1", 750, 750);
	c1->SetLeftMargin(0.15);
	c1->SetRightMargin(0.15);
	h1->Draw("COLZ1");
// 	h1->Draw();
// 	h1->SetMinimum(0.1);
// 	h1->Draw();
	TString save_name_pdf = directory + directory_year + name + ".pdf";
	TString save_name_png = directory + directory_year + name + ".png";
	c1->Print(save_name_pdf);
	c1->Print(save_name_png);	
// 	TString html_save_name_pdf = "/home/ferrico/public_html/lepCorr_2017/" + name + ".pdf";
// 	TString html_save_name_png = "/home/ferrico/public_html/lepCorr_2017/" + name + ".png";
// 	c1->Print(html_save_name_pdf);
// 	c1->Print(html_save_name_png);	

}
