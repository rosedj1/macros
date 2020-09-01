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

void Draw(TH1F *h1, TH1F *h2, TString nome_canvas, TString save, TLegend *legend);
void Draw_TH2F(TH2F* h1, TH2F* h2, TString nome_canvas, TString save, TLegend *legend);

TLegend *legend_reco_gen = new TLegend(0.75,0.75,0.9,0.9);
TLegend *legend_year = new TLegend(0.75,0.75,0.9,0.9);
TH1F *h_blank = new TH1F("h_blank", "h_blank", 10, -0.5, 9.5);

TFile* infile;
TTree* tree_MU;                                               

bool muon;
TString fs;
TString directory;
TString execute;
TString save_nome;

void YearComparison_check_ele_noFSR(TString reconstruction){

	gStyle->SetOptStat(0);
		
	Double_t pT1, eta1, phi1, Iso1;
	Double_t genLep_pt1;
	Double_t pterr1;
	Double_t genLep_eta1;
	Double_t genLep_phi1;
	Double_t pT2, eta2, phi2, Iso2;
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
	
	TString reco_gen[2] = {"reco_", "gen_"};
	std::vector<TString> variables;
	variables.push_back("mass2l");
	variables.push_back("massZ");
	variables.push_back("LeptonPt1");
	variables.push_back("LeptonPt2");
	variables.push_back("LeptonPtErr1");
	variables.push_back("LeptonPtErr2");
	variables.push_back("LeptonEta");
	variables.push_back("LeptonPhi");
	variables.push_back("LeptonPtErr1OverPt1");
	variables.push_back("LeptonPtErr2OverPt2");
	variables.push_back("RelIso1");
	variables.push_back("RelIso2");
	TString years[3] = {"2016", "2017", "2018"};
	
	TH1F* histograms[2][(int)variables.size()][3];
	TH2F* DeltaPt_vs_eta[2][3];
	TH2F* DeltaPt_vs_pT_0_0p9[2][3];
	TH2F* DeltaPt_vs_pT_0p9_1p8[2][3];
	TH2F* DeltaPt_vs_pT_1p8_2p4[2][3];
	
	gROOT->Reset();
	gROOT->SetBatch();

//     TString reconstruction = "amcatnlo";
//     TString reconstruction = "madgraph";
	
	for(int mu = 0; mu < 2; mu++){

		if(mu == 0){	
			muon = true;
			fs = "mu";
		}
		else{
			muon = false;
			fs = "e";
		}
	
		for(int i = 0; i < 2; i++){
			for(int var = 0; var < variables.size(); var++){
				for(int year = 0; year < 3; year++){
					TString nome_histo = reco_gen[i] + variables[var] + years[year];
					if(var < 2){
						histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 100, 60, 120);
        	            histograms[i][var][year]->GetXaxis()->SetTitle("m_{ll} [GeV]");
            	    }
					else if(var == 2 || var == 3){
						histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 100, 0, 200);
        	            histograms[i][var][year]->GetXaxis()->SetTitle("p_{T} [GeV]");
            	    }
					else if(var == 4 || var == 5){
						histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 200, 0, 10);
        	            histograms[i][var][year]->GetXaxis()->SetTitle("#sigma_{p_{T}}");
	                }
					else if(var == 6){
						histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 24, -2.4, 2.4);
            	        histograms[i][var][year]->GetXaxis()->SetTitle("#eta");
	                }
					else if(var == 7){
						histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 40, -TMath::Pi(), TMath::Pi());
            	        histograms[i][var][year]->GetXaxis()->SetTitle("#phi");
	                }
	                else if(var == 8 || var == 9){
						histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 100, 0, 0.1);
            	        histograms[i][var][year]->GetXaxis()->SetTitle("#sigma_{p_{T}}/p_{T}");
	                }
	                else{
						histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 100, 0, 0.1);
            	        histograms[i][var][year]->GetXaxis()->SetTitle("Rel. Iso.");
	                }
						
					if(var == 0){
						TString nome_histo = reco_gen[i] + "DeltaPt_vs_eta_" + years[year];
						DeltaPt_vs_eta[i][year] = new TH2F(nome_histo, nome_histo, 50, -2.4, 2.4, 50, 0.0, 0.1);
						nome_histo = reco_gen[i] + "DeltaPt_vs_pT_0_0p9_" + years[i];
						DeltaPt_vs_pT_0_0p9[i][year] = new TH2F(nome_histo, nome_histo, 50, 0, 100, 50, 0.0, 0.1);
						nome_histo = reco_gen[i] + "DeltaPt_vs_pT_0p9_1p8_" + years[i];
						DeltaPt_vs_pT_0p9_1p8[i][year] = new TH2F(nome_histo, nome_histo, 50, 0, 100, 50, 0.0, 0.1);
						nome_histo = reco_gen[i] + "DeltaPt_vs_pT_1p8_2p4_" + years[i];
						DeltaPt_vs_pT_1p8_2p4[i][year] = new TH2F(nome_histo, nome_histo, 50, 0, 100, 50, 0.0, 0.1);	
					}
				} // for on year
			} // for on variable
		} // for on Reco Gen
			
		for(int year = 0; year < 2; year++){
	
			if(year == 0){
			    if(muon) infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/" + reconstruction + "/DYJetsToLL_M-50_Full_RunII__m2mu_2017check_Rochester_amcatnlo_noFSR.root");                                                                               
        	    else infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/" + reconstruction + "/DYJetsToLL_M-50_Full_RunII__m2e_2017check_Rochester_amcatnlo_noFSR.root");                                                                                          
	        }		
	        else if(year == 1){                                                                                                                 
			    if(muon) infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/amcatnlo/DYJetsToLL_M-50_Full_RunII_amcatnlo_m2mu_2017check_noFSR.root");                                                                                  
        	    else infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/amcatnlo/DYJetsToLL_M-50_Full_RunII_amcatnlo_m2e_2017check_noFSR.root");                                                                                                                                         
			}
			else{
			    if(muon) infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/" + reconstruction + "/DYJetsToLL_M-50_Full_RunII_" + reconstruction + "_m2mu_2018.root");                                                                                                                                         
            	else infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/" + reconstruction + "/DYJetsToLL_M-50_Full_RunII_" + reconstruction + "_m2e_2018.root");                                                                                                                                         
	        }
	
			if(infile) tree_MU = (TTree*) infile->Get("passedEvents");		
			else std::cout<<"ERROR could not find the file"<<std::endl;

			tree_MU->SetBranchAddress("massZ", &massZ);       	
			tree_MU->SetBranchAddress("genzm", &genzm);       	
			tree_MU->SetBranchAddress("GENmass2l", &GENmass2l);   
			tree_MU->SetBranchAddress("pT1", &pT1); 
			tree_MU->SetBranchAddress("eta1", &eta1);   
			tree_MU->SetBranchAddress("phi1", &phi1);   
			tree_MU->SetBranchAddress("pterr1", &pterr1);   
			tree_MU->SetBranchAddress("Iso1", &Iso1);   
			tree_MU->SetBranchAddress("genLep_pt1", &genLep_pt1);   
			tree_MU->SetBranchAddress("genLep_eta1", &genLep_eta1);       	
			tree_MU->SetBranchAddress("genLep_phi1", &genLep_phi1);       	
			tree_MU->SetBranchAddress("pT2", &pT2);   
			tree_MU->SetBranchAddress("eta2", &eta2);   
			tree_MU->SetBranchAddress("phi2", &phi2);   
			tree_MU->SetBranchAddress("pterr2", &pterr2);   
			tree_MU->SetBranchAddress("Iso2", &Iso2);   
			tree_MU->SetBranchAddress("genLep_pt2", &genLep_pt2);   
			tree_MU->SetBranchAddress("genLep_eta2", &genLep_eta2);  
			tree_MU->SetBranchAddress("genLep_phi2", &genLep_phi2);  
			tree_MU->SetBranchAddress("weight", &weight);
            tree_MU->SetBranchAddress("lep1_ecalDriven", &lep1_ecalDriven);
            tree_MU->SetBranchAddress("lep2_ecalDriven", &lep2_ecalDriven);
	
			Long64_t nentries_MU = tree_MU->GetEntries();
			std::cout<<nentries_MU<<std::endl;	
	 		for(int entry = 0; entry < nentries_MU; entry++){
//			for(int entry = 0; entry < (int)nentries_MU/100; entry++){
// 	 		for(int entry = 0; entry < 1000; entry++){

				tree_MU->GetEntry(entry);                                   
				if(entry % 1000000 == 0){         
					if(muon)	std::cout<<entry<<" --- Dentro il TREE "<<year<<" --- MUON "<<year<<std::endl; 
					else 	std::cout<<entry<<" --- Dentro il TREE "<<year<<" --- ELECTRON "<<year<<std::endl; 
				}
				

				
// 				if(pterr1 < 0.5 || pterr1 > 1.5) continue;
// 				if(pterr2 < 0.5 || pterr2 > 1.5) continue;
// 				if(pterr1 > 0.5) continue;
// 				if(pterr2 > 0.5) continue;

			
//		 			mass_reco->Fill(massZ);
// 					reco_1.SetPtEtaPhiM(pT1, eta1, phi1, 0.105);
// 					reco_2.SetPtEtaPhiM(pT2, eta2, phi2, 0.105);
// 				float reco_mass = (reco_1+reco_2).M();
					histograms[0][1][year]->Fill(massZ);
					histograms[0][2][year]->Fill(pT1);
					histograms[0][3][year]->Fill(pT2);
                    if(lep1_ecalDriven)                             
                        histograms[0][4][year]->Fill(pterr1);       
                    else{                                   
                        if(fs!="mu")                        
                            histograms[0][4][year]->Fill(pterr1+8);                     
                        else                                                            
                            histograms[0][4][year]->Fill(pterr1);                                                       
                    }                                                                                                   
                    if(lep2_ecalDriven)                                                                         
                        histograms[0][5][year]->Fill(pterr2);                                                   
                    else{                                                                                                           
                        if(fs!="mu")                                                                                                            
                            histograms[0][5][year]->Fill(pterr2+8);                                                 
                        else
                            histograms[0][5][year]->Fill(pterr2);
                    }   
                    histograms[0][6][year]->Fill(eta1);
					histograms[0][6][year]->Fill(eta2);
					histograms[0][7][year]->Fill(phi1);
					histograms[0][7][year]->Fill(phi2);
					histograms[0][8][year]->Fill(pterr1/pT1);
					histograms[0][9][year]->Fill(pterr2/pT2);
					histograms[0][10][year]->Fill(Iso1);
					histograms[0][11][year]->Fill(Iso2);
		
// 					mass_gen->Fill(genzm);
					histograms[1][1][year]->Fill(GENmass2l);
					histograms[1][2][year]->Fill(genLep_pt1);
					histograms[1][3][year]->Fill(genLep_pt2);
					histograms[1][6][year]->Fill(genLep_eta1);
					histograms[1][6][year]->Fill(genLep_eta2);
					histograms[1][7][year]->Fill(genLep_phi1);
					histograms[1][7][year]->Fill(genLep_phi2);
			
// 					DeltaPt_vs_eta[0][year]->Fill(eta1, pterr1/pT1);//, weight);
// 					DeltaPt_vs_eta[0][year]->Fill(eta2, pterr2/pT2);//, weight);
// 	
// 					if(fabs(eta1) < 0.9)
// 						DeltaPt_vs_pT_0_0p9[0][year]->Fill(pT1, pterr1/pT1);
// 					else if(fabs(eta1) > 0.9 && fabs(eta1) <  1.8)
// 						DeltaPt_vs_pT_0p9_1p8[0][year]->Fill(pT1, pterr1/pT1);//, weight);
// 					else
// 						DeltaPt_vs_pT_1p8_2p4[0][year]->Fill(pT1, pterr1/pT1);//, weight);
// 	
// 					if(fabs(eta2) < 0.9)
// 						DeltaPt_vs_pT_0_0p9[0][year]->Fill(pT2, pterr2/pT2);
// 					else if(fabs(eta2) > 0.9 && fabs(eta2) <  1.8)
// 					DeltaPt_vs_pT_0p9_1p8[0][year]->Fill(pT2, pterr2/pT2);//, weight);
// 					else
// 						DeltaPt_vs_pT_1p8_2p4[0][year]->Fill(pT2, pterr2/pT2);//, weight);
			} // for on entry
		} // for on year

	
// 	TCanvas *c1 = new TCanvas("c1", "c1", 750, 500);
// 	histograms[0][5][0]->Draw();
// 	TCanvas *c2 = new TCanvas("c2", "c2", 750, 500);
// 	histograms[0][5][1]->Draw();
// 	return;
	
		if(mu == 0){
			legend_reco_gen->AddEntry(histograms[0][1][0], "Reco");
			legend_reco_gen->AddEntry(histograms[1][1][0], "Gen");
	
			legend_year->AddEntry(histograms[0][1][0], "2017_94X");
			legend_year->AddEntry(histograms[0][1][1], "2017_102");
		}

	
// 		directory = "./PlotYearComparison_pTerr_0p5_1p5_" + reconstruction;;
// 		directory = "./PlotYearComparison_pTerr_below0p5_" + reconstruction;;
		directory = "./PlotYearComparison_checkEle_noFSR" + reconstruction;;
		execute = "mkdir " + directory;
		gSystem->Exec(execute);
		directory += "/";
		if(muon){ 
            gSystem->Exec("mkdir " + directory + "MUON");
            directory += "MUON/";
        }
		else{
            gSystem->Exec("mkdir " + directory + "ELE");
            directory += "ELE/";
        }

		save_nome = "Reco_vs_gen_" + fs;
// 	Draw(h_blank, h_blank, "h_blank", save_nome + ".pdf[", legend_reco_gen);
// 
// 	for(int var = 0; var < variables.size(); var++){
// 			TString nome_canvas = variables[var] + years[0];
// 			Draw(histograms[0][var][0], histograms[1][var][0], nome_canvas, save_nome + ".pdf", legend_reco_gen);
// 	}
// 	
// 	for(int var = 0; var < variables.size(); var++){
// 			TString nome_canvas = variables[var] + years[1];
// 			Draw(histograms[0][var][1], histograms[1][var][1], nome_canvas, save_nome + ".pdf", legend_reco_gen);
// 	}
// 	Draw(h_blank, h_blank, "h_blank", save_nome + ".pdf]", legend_reco_gen);


		save_nome = "YearComparison_" + fs;	
		Draw(h_blank, h_blank, "h_blank", directory + save_nome + ".pdf[", legend_reco_gen);

		for(int i = 0; i < 2; i++){
			for(int var = 0; var < variables.size(); var++){
				TString nome_canvas = reco_gen[i] + variables[var];
				std:;cout<<nome_canvas<<std::endl;
// 				if(var == 0 || var == 4 || var == 5) continue; 
				if(i == 0 && var == 0) continue; 
				if(i == 1 && (var == 0 ||  var == 4 || var == 5 || var > 7)) continue;
				Draw(histograms[i][var][0], histograms[i][var][1], nome_canvas, directory + save_nome + ".pdf", legend_year);
				Draw(histograms[i][var][0], histograms[i][var][1], nome_canvas, directory + nome_canvas + ".png", legend_year);
			}
		}
	
// 	DeltaPt_vs_eta[0]->Scale(1/DeltaPt_vs_eta[0]->Integral());
// 	DeltaPt_vs_eta[1]->Scale(1/DeltaPt_vs_eta[1]->Integral());	
// 	DeltaPt_vs_eta[0]->Divide(DeltaPt_vs_eta[1]);
// 	TCanvas *c_DeltaPt_vs_eta = new TCanvas("c_DeltaPt_vs_eta", "c_DeltaPt_vs_eta", 750, 750);
// 	DeltaPt_vs_eta[0]->Draw("COLZ");
// 	DeltaPt_vs_eta[0]->GetZaxis()->SetRangeUser(0, 2);
// 	legend_year->Draw();
// 	c_DeltaPt_vs_eta->Print(directory + save_nome + ".pdf");
// 	c_DeltaPt_vs_eta->Print(directory + save_nome + ".png");
	
// 		Draw_TH2F(DeltaPt_vs_eta[0][0], DeltaPt_vs_eta[0][1], "DeltaPt_vs_eta_2016_2017", directory + save_nome + ".pdf", legend_year);
// 		Draw_TH2F(DeltaPt_vs_pT_0_0p9[0][0], DeltaPt_vs_pT_0_0p9[0][1], "DeltaPt_vs_pT_0_0p9_2016_2017", directory + save_nome + ".pdf", legend_year);
// 		Draw_TH2F(DeltaPt_vs_pT_0p9_1p8[0][0], DeltaPt_vs_pT_0p9_1p8[0][1], "DeltaPt_vs_pT_0p9_1p8_2016_2017", directory + save_nome + ".pdf", legend_year);
// 		Draw_TH2F(DeltaPt_vs_pT_1p8_2p4[0][0], DeltaPt_vs_pT_1p8_2p4[0][1], "DeltaPt_vs_pT_1p8_2p4_2016_2017", directory + save_nome + ".pdf", legend_year);
// 		Draw_TH2F(DeltaPt_vs_eta[0][0], DeltaPt_vs_eta[0][2], "DeltaPt_vs_eta_2016_2018", directory + save_nome + ".pdf", legend_year);
// 		Draw_TH2F(DeltaPt_vs_pT_0_0p9[0][0], DeltaPt_vs_pT_0_0p9[0][2], "DeltaPt_vs_pT_0_0p9_2016_2018", directory + save_nome + ".pdf", legend_year);
// 		Draw_TH2F(DeltaPt_vs_pT_0p9_1p8[0][0], DeltaPt_vs_pT_0p9_1p8[0][2], "DeltaPt_vs_pT_0p9_1p8_2016_2018", directory + save_nome + ".pdf", legend_year);
// 		Draw_TH2F(DeltaPt_vs_pT_1p8_2p4[0][0], DeltaPt_vs_pT_1p8_2p4[0][2], "DeltaPt_vs_pT_1p8_2p4_2016_2018", directory + save_nome + ".pdf", legend_year);

		Draw(h_blank, h_blank, "h_blank", directory + save_nome + ".pdf]", legend_reco_gen);
	
// 		Draw_TH2F(DeltaPt_vs_eta[0][0], DeltaPt_vs_eta[0][1], "DeltaPt_vs_eta_2016_2017", directory + "DeltaPt_vs_eta_2016_2017.png", legend_year);
// 		Draw_TH2F(DeltaPt_vs_pT_0_0p9[0][0], DeltaPt_vs_pT_0_0p9[0][1], "DeltaPt_vs_pT_0_0p9_2016_2017", directory + "DeltaPt_vs_pT_0_0p9_2016_2017.png", legend_year);
// 		Draw_TH2F(DeltaPt_vs_pT_0p9_1p8[0][0], DeltaPt_vs_pT_0p9_1p8[0][1], "DeltaPt_vs_pT_0p9_1p8_2016_2017", directory + "DeltaPt_vs_pT_0p9_1p8_2016_2017.png", legend_year);
// 		Draw_TH2F(DeltaPt_vs_pT_1p8_2p4[0][0], DeltaPt_vs_pT_1p8_2p4[0][1], "DeltaPt_vs_pT_1p8_2p4_2016_2017", directory + "DeltaPt_vs_pT_0p9_1p8_2016_2017.png", legend_year);
// 		Draw_TH2F(DeltaPt_vs_eta[0][0], DeltaPt_vs_eta[0][2], "DeltaPt_vs_eta_2016_2018", directory + "DeltaPt_vs_eta_2016_2018.png", legend_year);
// 		Draw_TH2F(DeltaPt_vs_pT_0_0p9[0][0], DeltaPt_vs_pT_0_0p9[0][2], "DeltaPt_vs_pT_0_0p9_2016_2018", directory + "DeltaPt_vs_pT_0_0p9_2016_2018.png", legend_year);
// 		Draw_TH2F(DeltaPt_vs_pT_0p9_1p8[0][0], DeltaPt_vs_pT_0p9_1p8[0][2], "DeltaPt_vs_pT_0p9_1p8_2016_2018", directory + "DeltaPt_vs_pT_0p9_1p8_2016_2018.png", legend_year);
// 		Draw_TH2F(DeltaPt_vs_pT_1p8_2p4[0][0], DeltaPt_vs_pT_1p8_2p4[0][2], "DeltaPt_vs_pT_1p8_2p4_2016_2018", directory + "DeltaPt_vs_pT_1p8_2p4_2016_2018.png", legend_year);
		
	} // for on mu = muon or electron



// 	Draw(mass_reco, mass_gen, "mass");
// 	Draw(pT1_reco, pT1_gen, "pT1");
// 	Draw(pT2_reco, pT2_gen, "pT2");
// 	Draw(eta_reco, eta_gen, "eta");
// 	Draw(phi_reco, phi_gen, "phi");

}


void Draw_TH2F(TH2F* h1, TH2F* h2, TString nome_canvas, TString save, TLegend *legend){

	h1->Scale(1/h1->Integral());
	h2->Scale(1/h2->Integral());	
	h1->Divide(h2);
	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
	h1->Draw("COLZ");
	h1->GetZaxis()->SetRangeUser(0, 2);
// 	legend_year->Draw();

	canvas->Print(save);

}

void Draw(TH1F *h1, TH1F *h2, TString nome_canvas, TString save, TLegend *legend){

	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();

	Double_t min, max;
	
	h1->Scale(1/h1->Integral());
	h2->Scale(1/h2->Integral());
	
	if(h1->GetMaximum() > h2->GetMaximum())
			max = h1->GetMaximum();
	else
			max = h2->GetMaximum();
	
	max = 1.1 * max;
	
	if(max == 0) max = 1;
	
	min = 0;
	
	h1->SetTitle(nome_canvas);
	h2->SetTitle(nome_canvas);
	h1->SetLineColor(kGreen+2);
	h2->SetLineColor(kRed+2);


	h1->GetYaxis()->SetRangeUser(0, max);//SetMaximum(max);
	h2->GetYaxis()->SetRangeUser(0, max);//SetMaximum(max);
	
	canvas->Update();

	
	h1->Draw();
	h2->Draw("Same");
	legend->Draw();	
	
	TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
	pt->SetBorderSize(0);
	pt->SetTextAlign(12);
	pt->SetFillStyle(0);
	pt->SetTextFont(42);
	pt->SetTextSize(0.05);
	TText *text;
	if(muon)
		text = pt->AddText(0.01,0.5,"Muon");
	else
		text = pt->AddText(0.01,0.5,"Electron");
	pt->Draw();   
// 	TLatex *fs_label = new TLatex();
//     fs_label->SetTextFont(42);
//     fs_label->SetTextSize(0.075);
// 	TString longstring;
// 	if(muon)
// 		longstring	= "Muon";
// 	else
// 		longstring	= "Electron";
// // 	fs_label->SetTextColor(kGreen+2);
//    	fs_label->DrawLatex(-0.85, yPos, longstring);
	TH1F *ratio_2017 = (TH1F*) h2->Clone();

	canvas->Update();
	canvas->cd();

	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.25);
// 	pad22->SetTopMargin(0);
// 	pad22->SetTopMargin(0.95);
// 	pad22->SetBottomMargin(0.25);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	
	ratio_2017->Divide(h1);
	ratio_2017->GetYaxis()->SetRangeUser(0.5, 1.5);
	ratio_2017->SetTitle("");
	ratio_2017->SetStats(0);
	ratio_2017->GetYaxis()->SetTitleSize(0.1);
// 	ratio->GetYaxis()->SetTitleFont(43);
	ratio_2017->GetYaxis()->SetLabelSize(0.14);
	ratio_2017->GetYaxis()->SetTitle("2017_102 / 2017_94X");//Reco / Gen");
	ratio_2017->GetYaxis()->SetTitleOffset(0.50);
	ratio_2017->GetYaxis()->SetNdivisions(506); 
// 	ratio->SetLineColor(kBlack);	
	ratio_2017->SetLineColor(kRed+2);
	ratio_2017->Draw();
	canvas->Update();
	pad22->Update();
	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(1);
	line->Draw();

	canvas->Print(save);
	
}






