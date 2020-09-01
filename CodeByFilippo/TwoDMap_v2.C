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
float EvalMax(TH2F* h1, TH2F* h2, TH2F* h3);
void  Save(TH2F* h1, TString name, float maximum);
void DrawSame(TH2F* h1, TH2F* h2,TString name);

TFile* infile_MU;
TFile* infile_ELE;
TTree* tree_MU;                                               
TTree* tree_ELE;                                               

// std::vector<TString> year_name;
TString year_name;

TString directory_year;
// TString reconstruction = "madgraph";
//TString reconstruction = "amcatnlo";
TString directory;// = "2D_plot_" + reconstruction;

void TwoDMap_v2(TString reconstruction) {

	 directory = "2D_plot_" + reconstruction;
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

	TH2F *ErrPt_vs_eta_MU[3];
	TH2F *ErrPt_vs_pt_0_09_MU[3];
	TH2F *ErrPt_vs_pt_09_18_MU[3];
	TH2F *ErrPt_vs_pt_18_24_MU[3];

	TH2F *ErrPt_vs_pt_0_08_ELE[3];
	TH2F *ErrPt_vs_pt_08_1_ELE[3];
	TH2F *ErrPt_vs_pt_1_1p44_ELE[3];
	TH2F *ErrPt_vs_pt_1p44_1p57_ELE[3];
	TH2F *ErrPt_vs_pt_1p57_2_ELE[3];
	TH2F *ErrPt_vs_pt_2_25_ELE[3];

	TH2F *ErrPt_vs_eta_ECAL[3];
	TH2F *ErrPt_vs_pt_Less1_ECAL[3];
	TH2F *ErrPt_vs_pt_Greater1_ECAL[3];

	TH2F *ErrPt_vs_eta_TRACKER[3];
	TH2F *ErrPt_vs_pt_0_14_TRACKER[3];
	TH2F *ErrPt_vs_pt_14_16_TRACKER[3];
	TH2F *ErrPt_vs_pt_16_2_TRACKER[3];
	TH2F *ErrPt_vs_pt_2_25_TRACKER[3];



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
	
	float count_entry[3] = {0};
	
// 	year_name.push_back("2016");
// 	year_name.push_back("2017");
// 	year_name.push_back("2018");
	
	for(int year = 0; year < 3; year++){
	
		if(year == 0){
			infile_MU = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/" + reconstruction + "/DYJetsToLL_M-50_Full_RunII_" + reconstruction + "_m2mu_2016.root");
			infile_ELE = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/" + reconstruction + "/DYJetsToLL_M-50_Full_RunII_" + reconstruction + "_m2e_2016.root");
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


	ErrPt_vs_eta_MU[year] = new TH2F("2D_Map_vs_eta_muon_" + year_name, "2D_Map_vs_eta_muon_" + year_name, 100, -2.4, 2.4, 100, 0.0, 0.1);
	ErrPt_vs_pt_0_09_MU[year] = new TH2F("2D_Map_vs_pt_0_09_" + year_name, "2D_Map_vs_pt_0_09_" + year_name, 100, 5, 100, 100, 0.0, 0.1);
	ErrPt_vs_pt_09_18_MU[year] = new TH2F("2D_Map_vs_pt_09_18_" + year_name, "2D_Map_vs_pt_09_18_" + year_name, 100, 5, 100, 100, 0.0, 0.1);
	ErrPt_vs_pt_18_24_MU[year] = new TH2F("2D_Map_vs_pt_18_24_" + year_name, "2D_Map_vs_pt_18_24_" + year_name, 100, 5, 100, 100, 0.0, 0.1);

	ErrPt_vs_pt_0_08_ELE[year] = new TH2F("2D_Map_vs_pt_0_08_" + year_name, "2D_Map_vs_pt_0_08_" + year_name, 100, 0.0, 100, 100, 0.0, 0.1);
	ErrPt_vs_pt_08_1_ELE[year] = new TH2F("2D_Map_vs_pt_08_1_" + year_name, "2D_Map_vs_pt_08_1_" + year_name, 100, 0.0, 100, 100, 0.0, 0.1);
	ErrPt_vs_pt_1_1p44_ELE[year] = new TH2F("2D_Map_vs_pt_1p44_" + year_name, "2D_Map_vs_pt_1p44_" + year_name, 100, 0.0, 100, 100, 0.0, 0.1);
	ErrPt_vs_pt_1p44_1p57_ELE[year] = new TH2F("2D_Map_vs_pt_1p44_1p57_" + year_name, "2D_Map_vs_pt_1p44_1p57_" + year_name, 100, 0.0, 100, 100, 0.0, 0.1);
	ErrPt_vs_pt_1p57_2_ELE[year] = new TH2F("2D_Map_vs_pt_1p57_2_" + year_name, "2D_Map_vs_pt_1p57_2_" + year_name, 100, 0.0, 100, 100, 0.0, 0.1);
	ErrPt_vs_pt_2_25_ELE[year] = new TH2F("2D_Map_vs_pt_2_25_", "2D_Map_vs_pt_2_25_", 100, 0.0, 100, 100, 0.0, 0.1);

	ErrPt_vs_eta_ECAL[year] =new TH2F("2D_Map_vs_eta_ECAL_" + year_name, "2D_Map_vs_eta_ECAL_" + year_name, 100, -2.4, 2.4, 100, 0.0, 0.1);
	ErrPt_vs_pt_Less1_ECAL[year] =new TH2F("ErrPt_vs_pt_Less1_ECAL_", "ErrPt_vs_pt_Less1_ECAL_" + year_name, 100, 0.0, 100, 100, 0.0, 0.1);
	ErrPt_vs_pt_Greater1_ECAL[year] =new TH2F("ErrPt_vs_pt_Greater1_ECAL_" + year_name, "ErrPt_vs_pt_Greater1_ECAL_" + year_name, 100, 0.0, 100, 100, 0.0, 0.1);

	ErrPt_vs_eta_TRACKER[year] =new TH2F("2D_Map_vs_eta_TRACKER_" + year_name, "2D_Map_vs_eta_TRACKER_" + year_name, 100, -2.4, 2.4, 100, 0.0, 0.1);
	ErrPt_vs_pt_0_14_TRACKER[year] =new TH2F("ErrPt_vs_pt_0_14_TRACKER_" + year_name, "ErrPt_vs_pt_0_14_TRACKER_" + year_name, 100, 0.0, 100, 100, 0.0, 0.1);
	ErrPt_vs_pt_14_16_TRACKER[year] =new TH2F("ErrPt_vs_pt_14_16_TRACKER_" + year_name, "ErrPt_vs_pt_14_16_TRACKER_" + year_name, 100, 0.0, 100, 100, 0.0, 0.1);
	ErrPt_vs_pt_16_2_TRACKER[year] =new TH2F("ErrPt_vs_pt_16_2_TRACKER_" + year_name, "ErrPt_vs_pt_16_2_TRACKER_" + year_name, 100, 0.0, 100, 100, 0.0, 0.1);
	ErrPt_vs_pt_2_25_TRACKER[year] =new TH2F("ErrPt_vs_pt_2_25_TRACKER_" + year_name, "ErrPt_vs_pt_2_25_TRACKER_" + year_name, 100, 0.0, 100, 100, 0.0, 0.1);







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
      		
		Long64_t nentries_MU = tree_MU->GetEntries();
		std::cout<<nentries_MU<<std::endl;	

 		for(int i = 0; i < nentries_MU; i++){
//		for(int i = 0; i < (int)nentries_MU/100; i++){
//   		for(int i = 0; i < 100; i++){
		
			tree_MU->GetEntry(i);                                   
			if(i % 1000000 == 0)         
				std::cout<<i<<" --- Dentro il TREE --- MUON "<<year<<std::endl; 
			
// 		if(pT1 < 5 || pT2 < 5) std::cout<<"OK APPLICA il taglio ai reco"<<std::endl;
// 		 || genLep_pt2 < 5) std::cout<<"OK APPLICA il taglio ai GEN"<<std::endl;
		
			if(massZ < 60 || massZ > 120) continue;
		
// 		h_mass->Fill(massZ);

			ErrPt_vs_eta_MU[year]->Fill(eta1, pterr1/pT1);//, weight);
			ErrPt_vs_eta_MU[year]->Fill(eta2, pterr2/pT2);//, weight);
// 		if(pterr1/pT1 < 0.004) std::cout<<"Mu_1 = "<<eta1<<std::endl;
// 		if(pterr2/pT2 < 0.004) std::cout<<"Mu_2 = "<<eta2<<std::endl;

			if(abs(eta1) < 0.9 && abs(eta2) < 0.9){		
// 		if(fabs(eta1) < 0.9)
					ErrPt_vs_pt_0_09_MU[year]->Fill(pT1, pterr1/pT1);//, weight);
// 		if(fabs(eta2) < 0.9)
					ErrPt_vs_pt_0_09_MU[year]->Fill(pT2, pterr2/pT2);//, weight);
// 				if(pT1 < 10) std::cout<<"First muon = "<<pterr1/genLep_pt1<<"\t"<<fabs(eta1)<<std::endl;
// 				if(pT2 < 10) std::cout<<"Second muon = "<<pterr2/genLep_pt2<<"\t"<<fabs(eta2)<<std::endl;
// 			
// 			h_0_09->Fill(genLep_eta1);
// 			h_0_09->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 0.9 && fabs(genLep_eta2) > 0.9 && 
					fabs(genLep_eta1) < 1.8 && fabs(genLep_eta2) < 1.8){

				ErrPt_vs_pt_09_18_MU[year]->Fill(pT1, pterr1/pT1);//, weight);
				ErrPt_vs_pt_09_18_MU[year]->Fill(pT2, pterr2/pT2);//, weight);

				h_09_18->Fill(genLep_eta1);
				h_09_18->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 1.8 && fabs(genLep_eta2) > 1.8 && 
				fabs(genLep_eta1) < 2.4 && fabs(genLep_eta2) < 2.4){

				ErrPt_vs_pt_18_24_MU[year]->Fill(pT1, pterr1/pT1);//, weight);
				ErrPt_vs_pt_18_24_MU[year]->Fill(pT2, pterr2/pT2);//, weight);

				h_18_24->Fill(genLep_eta1);
				h_18_24->Fill(genLep_eta2);
			}
		
		}
	
// 	TCanvas *c1 = new TCanvas("c1", "c1", 500 ,500);
// 	h_mass->Draw();
// 	
// 	return;

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
	
		Long64_t nentries_ELE = tree_ELE->GetEntries();
		std::cout<<nentries_ELE<<std::endl;	
 		for(int i = 0; i < nentries_ELE; i++){
//   		for(int i = 0; i < (int)nentries_ELE/100; i++){
// 	  		for(int i = 0; i < 100; i++){

			tree_ELE->GetEntry(i);                                   
			if(i % 1000000 == 0)         
				std::cout<<i<<" --- Dentro il TREE --- ELECTRON "<<year<<std::endl; 

			if(massZ < 60 || massZ > 120) continue;
		
			if(lep1_ecalDriven == 1){
// 				if(pterr1/pT1 > 0.09 && pterr1/pT1 < 0.1 && eta1 > -0.48 && eta1 < 0.48){
// 					std::cout<<"Entry 1 = "<<pterr1/pT1<<"\t"<<year_name<<std::endl;
// 					count_entry[year]+=1/nentries_ELE;
// 				}
				ErrPt_vs_eta_ECAL[year]->Fill(eta1, pterr1/pT1);//, weight);
			}
			else
				ErrPt_vs_eta_TRACKER[year]->Fill(eta1, pterr1/pT1);//, weight);
		
			if(lep2_ecalDriven == 1){
// 				if(pterr2/pT2 > 0.09 && pterr2/pT2 < 0.1 && eta2 > -0.48 && eta2 < 0.48){
// 					std::cout<<"Entry 2 = "<<pterr2/pT2<<"\t"<<year_name<<std::endl;
// 					count_entry[year]+=1/nentries_ELE;
// 				}
				ErrPt_vs_eta_ECAL[year]->Fill(eta2, pterr2/pT2);//, weight);
			}
			else
				ErrPt_vs_eta_TRACKER[year]->Fill(eta2, pterr2/pT2);//, weight);

			if(fabs(genLep_eta1) > 0.0 && fabs(genLep_eta2) > 0.0 && 
					fabs(genLep_eta1) < 0.8 && fabs(genLep_eta2) < 0.8){
									ErrPt_vs_pt_0_08_ELE[year]->Fill(pT1, pterr1/pT1);//, weight);
									ErrPt_vs_pt_0_08_ELE[year]->Fill(pT2, pterr1/pT1);//, weight);
									h_ELE_0_08->Fill(genLep_eta1);
									h_ELE_0_08->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 0.8 && fabs(genLep_eta2) > 0.8 && 
					fabs(genLep_eta1) < 1. && fabs(genLep_eta2) < 1.){
									ErrPt_vs_pt_08_1_ELE[year]->Fill(pT1, pterr1/pT1);//, weight);
									ErrPt_vs_pt_08_1_ELE[year]->Fill(pT2, pterr1/pT1);//, weight);
									h_ELE_08_1->Fill(genLep_eta1);
									h_ELE_08_1->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 1. && fabs(genLep_eta2) > 1. && 
					fabs(genLep_eta1) < 1.44 && fabs(genLep_eta2) < 1.44){
									ErrPt_vs_pt_1_1p44_ELE[year]->Fill(pT1, pterr1/pT1);//, weight);
									ErrPt_vs_pt_1_1p44_ELE[year]->Fill(pT2, pterr1/pT1);//, weight);
									h_ELE_1_1p44->Fill(genLep_eta1);
									h_ELE_1_1p44->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 1.44 && fabs(genLep_eta2) > 1.44 && 
					fabs(genLep_eta1) < 1.57 && fabs(genLep_eta2) < 1.57){
									ErrPt_vs_pt_1p44_1p57_ELE[year]->Fill(pT1, pterr1/pT1);//, weight);
									ErrPt_vs_pt_1p44_1p57_ELE[year]->Fill(pT2, pterr1/pT1);//, weight);
									h_ELE_1p44_1p57->Fill(genLep_eta1);
									h_ELE_1p44_1p57->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 1.57 && fabs(genLep_eta2) > 1.57 && 
					fabs(genLep_eta1) < 2. && fabs(genLep_eta2) < 2){
									ErrPt_vs_pt_1p57_2_ELE[year]->Fill(pT1, pterr1/pT1);//, weight);
									ErrPt_vs_pt_1p57_2_ELE[year]->Fill(pT2, pterr1/pT1);//, weight);
									h_ELE_1p57_2->Fill(genLep_eta1);
									h_ELE_1p57_2->Fill(genLep_eta2);
			}
			else if(fabs(genLep_eta1) > 2. && fabs(genLep_eta2) > 2. && 
					fabs(genLep_eta1) < 2.5 && fabs(genLep_eta2) < 2.5){
									ErrPt_vs_pt_2_25_ELE[year]->Fill(pT1, pterr1/pT1);//, weight);
									ErrPt_vs_pt_2_25_ELE[year]->Fill(pT2, pterr1/pT1);//, weight);
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

		
// 	
// 		Norm(ErrPt_vs_eta_MU[year]);
// 		Norm(ErrPt_vs_pt_0_09_MU[year]);
// 		Norm(ErrPt_vs_pt_09_18_MU[year]);
// 		Norm(ErrPt_vs_pt_18_24_MU[year]);
// 		Norm(ErrPt_vs_eta_ECAL[year]);
// 		Norm(ErrPt_vs_eta_TRACKER[year]);
// 		Norm(ErrPt_vs_pt_0_08_ELE[year]);
// 		Norm(ErrPt_vs_pt_08_1_ELE[year]);
// 		Norm(ErrPt_vs_pt_1_1p44_ELE[year]);
// 		Norm(ErrPt_vs_pt_1p44_1p57_ELE[year]);
// 		Norm(ErrPt_vs_pt_1p57_2_ELE[year]);
// 		Norm(ErrPt_vs_pt_2_25_ELE[year]);

        //gROOT->Reset();
        gROOT->SetBatch();
    }


		
	for(int year = 0; year < 3; year++){

		if(year == 0)
			year_name = "2016";
		else if(year == 1)
			year_name = "2017";
		else
            year_name = "2018";

		gSystem->Exec("cd " + directory);
		directory_year = year_name;
		std::cout<<directory_year<<std::endl;
		execute = "mkdir " + directory + directory_year;
		std::cout<<execute<<std::endl;
		gSystem->Exec(execute);
		directory_year += "/";
		gSystem->Exec("cd ..");
		
		Save(ErrPt_vs_eta_MU[year], "ErrPt_vs_eta_MU_" + year_name, EvalMax(ErrPt_vs_eta_MU[0], ErrPt_vs_eta_MU[1], ErrPt_vs_eta_MU[2]));
		Save(ErrPt_vs_pt_0_09_MU[year], "ErrPt_vs_pt_0_09_MU_" + year_name, EvalMax(ErrPt_vs_pt_0_09_MU[0], ErrPt_vs_pt_0_09_MU[1], ErrPt_vs_pt_0_09_MU[2]));
		Save(ErrPt_vs_pt_09_18_MU[year], "ErrPt_vs_pt_09_18_MU_" + year_name, EvalMax(ErrPt_vs_pt_09_18_MU[0], ErrPt_vs_pt_09_18_MU[1], ErrPt_vs_pt_09_18_MU[2]));
		Save(ErrPt_vs_pt_18_24_MU[year], "ErrPt_vs_pt_18_24_MU_" + year_name, EvalMax(ErrPt_vs_pt_18_24_MU[0], ErrPt_vs_pt_18_24_MU[1], ErrPt_vs_pt_18_24_MU[2]));
		Save(ErrPt_vs_eta_ECAL[year], "ErrPt_vs_eta_ECAL_" + year_name, EvalMax(ErrPt_vs_eta_ECAL[0], ErrPt_vs_eta_ECAL[1], ErrPt_vs_eta_ECAL[2]));
		Save(ErrPt_vs_eta_TRACKER[year], "ErrPt_vs_eta_TRACKER_" + year_name, EvalMax(ErrPt_vs_eta_TRACKER[0], ErrPt_vs_eta_TRACKER[1], ErrPt_vs_eta_TRACKER[2]));
		Save(ErrPt_vs_pt_0_08_ELE[year], "ErrPt_vs_pt_0_08_ELE_" + year_name, EvalMax(ErrPt_vs_pt_0_08_ELE[0], ErrPt_vs_pt_0_08_ELE[1], ErrPt_vs_pt_0_08_ELE[2]));
		Save(ErrPt_vs_pt_08_1_ELE[year], "ErrPt_vs_pt_08_1_ELE_" + year_name, EvalMax(ErrPt_vs_pt_08_1_ELE[0], ErrPt_vs_pt_08_1_ELE[1], ErrPt_vs_pt_08_1_ELE[2]));
		Save(ErrPt_vs_pt_1_1p44_ELE[year], "ErrPt_vs_pt_1_1p44_ELE_" + year_name, EvalMax(ErrPt_vs_pt_1_1p44_ELE[0], ErrPt_vs_pt_1_1p44_ELE[1], ErrPt_vs_pt_1_1p44_ELE[2]));
		Save(ErrPt_vs_pt_1p44_1p57_ELE[year], "ErrPt_vs_pt_1p44_1p57_ELE_" + year_name, EvalMax(ErrPt_vs_pt_1p44_1p57_ELE[0], ErrPt_vs_pt_1p44_1p57_ELE[1], ErrPt_vs_pt_1p44_1p57_ELE[2]));
		Save(ErrPt_vs_pt_1p57_2_ELE[year], "ErrPt_vs_pt_1p57_2_ELE_" + year_name, EvalMax(ErrPt_vs_eta_MU[0], ErrPt_vs_pt_1p57_2_ELE[1], ErrPt_vs_pt_1p57_2_ELE[2]));
		Save(ErrPt_vs_pt_2_25_ELE[year], "ErrPt_vs_pt_2_25_ELE_" + year_name, EvalMax(ErrPt_vs_pt_2_25_ELE[0], ErrPt_vs_pt_2_25_ELE[1], ErrPt_vs_pt_2_25_ELE[2]));
	}


	DrawSame(ErrPt_vs_eta_MU[0], ErrPt_vs_eta_MU[1], "ErrPt_vs_eta_MU_2016_vs_2017");
	DrawSame(ErrPt_vs_pt_0_09_MU[0], ErrPt_vs_pt_0_09_MU[1], "ErrPt_vs_pt_0_09_MU_2016_vs_2017");
	DrawSame(ErrPt_vs_pt_09_18_MU[0], ErrPt_vs_pt_09_18_MU[1], "ErrPt_vs_pt_09_18_MU_2016_vs_2017");
	DrawSame(ErrPt_vs_pt_18_24_MU[0], ErrPt_vs_pt_18_24_MU[1], "ErrPt_vs_pt_18_24_MU_2016_vs_2017");
	DrawSame(ErrPt_vs_eta_ECAL[0], ErrPt_vs_eta_ECAL[1], "ErrPt_vs_eta_ECAL_2016_vs_2017");
	DrawSame(ErrPt_vs_eta_TRACKER[0], ErrPt_vs_eta_TRACKER[1], "ErrPt_vs_eta_TRACKER_2016_vs_2017");
	DrawSame(ErrPt_vs_pt_0_08_ELE[0], ErrPt_vs_pt_0_08_ELE[1], "ErrPt_vs_pt_0_08_ELE_2016_vs_2017");
	DrawSame(ErrPt_vs_pt_08_1_ELE[0], ErrPt_vs_pt_08_1_ELE[1], "ErrPt_vs_pt_08_1_ELE_2016_vs_2017");
	DrawSame(ErrPt_vs_pt_1_1p44_ELE[0], ErrPt_vs_pt_1_1p44_ELE[1], "ErrPt_vs_pt_1_1p44_ELE_2016_vs_2017");
	DrawSame(ErrPt_vs_pt_1p44_1p57_ELE[0], ErrPt_vs_pt_1p44_1p57_ELE[1], "ErrPt_vs_pt_1p44_1p57_ELE_2016_vs_2017");
	DrawSame(ErrPt_vs_pt_1p57_2_ELE[0], ErrPt_vs_pt_1p57_2_ELE[1], "ErrPt_vs_pt_1p57_2_ELE_2016_vs_2017");
	DrawSame(ErrPt_vs_pt_2_25_ELE[0], ErrPt_vs_pt_2_25_ELE[1], "ErrPt_vs_pt_2_25_ELE_2016_vs_2017");

	DrawSame(ErrPt_vs_eta_MU[1], ErrPt_vs_eta_MU[2], "ErrPt_vs_eta_MU_2017_vs_2018");
	DrawSame(ErrPt_vs_pt_0_09_MU[1], ErrPt_vs_pt_0_09_MU[2], "ErrPt_vs_pt_0_09_MU_2017_vs_2018");
	DrawSame(ErrPt_vs_pt_09_18_MU[1], ErrPt_vs_pt_09_18_MU[2], "ErrPt_vs_pt_09_18_MU_2017_vs_2018");
	DrawSame(ErrPt_vs_pt_18_24_MU[1], ErrPt_vs_pt_18_24_MU[2], "ErrPt_vs_pt_18_24_MU_2017_vs_2018");
	DrawSame(ErrPt_vs_eta_ECAL[1], ErrPt_vs_eta_ECAL[2], "ErrPt_vs_eta_ECAL_2017_vs_2018");
	DrawSame(ErrPt_vs_eta_TRACKER[1], ErrPt_vs_eta_TRACKER[2], "ErrPt_vs_eta_TRACKER_2017_vs_2018");
	DrawSame(ErrPt_vs_pt_0_08_ELE[1], ErrPt_vs_pt_0_08_ELE[2], "ErrPt_vs_pt_0_08_ELE_2017_vs_2018");
	DrawSame(ErrPt_vs_pt_08_1_ELE[1], ErrPt_vs_pt_08_1_ELE[2], "ErrPt_vs_pt_08_1_ELE_2017_vs_2018");
	DrawSame(ErrPt_vs_pt_1_1p44_ELE[1], ErrPt_vs_pt_1_1p44_ELE[2], "ErrPt_vs_pt_1_1p44_ELE_2017_vs_2018");
	DrawSame(ErrPt_vs_pt_1p44_1p57_ELE[1], ErrPt_vs_pt_1p44_1p57_ELE[2], "ErrPt_vs_pt_1p44_1p57_ELE_2017_vs_2018");
	DrawSame(ErrPt_vs_pt_1p57_2_ELE[1], ErrPt_vs_pt_1p57_2_ELE[2], "ErrPt_vs_pt_1p57_2_ELE_2017_vs_2018");
	DrawSame(ErrPt_vs_pt_2_25_ELE[1], ErrPt_vs_pt_2_25_ELE[2], "ErrPt_vs_pt_2_25_ELE_2017_vs_2018");



// 	std::cout<<count_entry[0]<<"\t"<<count_entry[1]<<"\t"<<count_entry[2]<<std::endl;
                                                                                            
}                                                   




void Norm(TH2F* h1){
	Double_t norm = h1->Integral();
	h1->Scale(1/norm);
}

void Save(TH2F* h1, TString name, float maximum){

	Norm(h1);
// 	std::cout<<"Final max value = "<<maximum<<std::endl;
	h1->GetZaxis()->SetRangeUser(0, maximum);
	TCanvas *c1 = new TCanvas("c1", "c1", 750, 750);
	c1->SetLeftMargin(0.15);
	c1->SetRightMargin(0.15);
	h1->Draw("COLZ1");
    h1->GetXaxis()->SetTile("#eta");
    h1->GetYaxis()->SetTile("p_{T}err/p_{T}");
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

void DrawSame(TH2F* h1, TH2F* h2,TString name){

	Norm(h1);
	Norm(h2);
	
	TH2F *h22 = (TH2F*)h2->Clone();
	
	h1->RebinX(2);
	h1->RebinY(2);
    h1->SetTitle(name);
	h22->RebinX(2);
	h22->RebinY(2);
    h22->SetTitle(name);
		
	TCanvas *c1 = new TCanvas("c1", "c1", 750, 750);
	c1->SetLeftMargin(0.15);
	c1->SetRightMargin(0.15);
	h1->Divide(h22);
	h1->Draw("COLZ1");
// 	h1->SetMinimum(0.1);
// 	h1->Draw();
	TString save_name_pdf = directory + name + ".pdf";
	TString save_name_png = directory + name + ".png";
	c1->Print(save_name_pdf);
	c1->Print(save_name_png);	
	

}

float EvalMax(TH2F* h1, TH2F* h2, TH2F* h3){

	float maximum = -999;
	
	Norm(h1);
	Norm(h2);
	Norm(h3);
		
	float h1_max = h1->GetBinContent(h1->GetMaximumBin());
	float h2_max = h2->GetBinContent(h2->GetMaximumBin());
	float h3_max = h3->GetBinContent(h3->GetMaximumBin());
	
	std::cout<<h1_max<<"\t"<<h2_max<<"\t"<<h3_max<<std::endl;
	
	if(h1_max > h2_max){
		if(h1_max > h3_max) maximum =  h1_max;
		else maximum =  h3_max;
	}
	else{
		if(h2_max > h3_max) maximum =  h2_max;
		else maximum =  h3_max;
	}
	
	return maximum;
	
}
