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

void Draw(TH1F *h1, TH1F *h2, TH1F *h3, TString nome_canvas, TString save, TLegend *legend);
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


void YearComparison_EleNoScale(){
	
	float sigmaPt;
	float etaMin; 
	float etaMax;

	gStyle->SetOptStat(0);

   	std::vector<float> eta_bins;
   	eta_bins.clear();
	eta_bins.push_back(0);
	eta_bins.push_back(0.8);
	eta_bins.push_back(1);
	eta_bins.push_back(1.2);
	eta_bins.push_back(1.44);
	eta_bins.push_back(1.57);
	eta_bins.push_back(2);
	eta_bins.push_back(2.5);

   	std::vector<TString> eta_bins_name;
   	eta_bins_name.clear();
	eta_bins_name.push_back("0");
	eta_bins_name.push_back("0p8");
	eta_bins_name.push_back("1");
	eta_bins_name.push_back("1p2");
	eta_bins_name.push_back("1p44");
	eta_bins_name.push_back("1p57");
	eta_bins_name.push_back("2p0");
	eta_bins_name.push_back("2p5");

	Double_t pT1, eta1, phi1, Iso1;
	Double_t genLep_pt1;
	Double_t pterr1, pterr1old;
	Double_t genLep_eta1;
	Double_t genLep_phi1;
	Double_t pT2, eta2, phi2, Iso2;
	Double_t genLep_pt2;
	Double_t pterr2, pterr2old;
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
    variables.push_back("pterr1old");
    variables.push_back("pterr2old");
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
	for(int etaB = 0; etaB <eta_bins.size()-1; etaB++){
		for(int mu = 1; mu < 2; mu++){

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
						else if(var == 4 || var == 5 || var == 12 || var == 13){
							histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 100, 0, 5);
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
					} // for on year
				} // for on variable
			} // for on Reco Gen
			
			for(int year = 0; year < 3; year++){
	    	                
	
	            if(year == 0){                                                                                                                                   
        	        if(muon) infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/amcatnlo/DYJetsToLL_M-50_Full_RunII__m2mu_2017_check_Rochester_amcatnlo.root");                                                                                                    
//     	            else infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/2017_94X_check/DYJetsToLL_M-50_Full_RunII__m2e_2017_check_Rochester_amcatnlo.root");                                                                                                    
    	            else infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Electron_2017_amcatnlo_rochester.root");                                                                                                    
	            }                                                                                                                                    
        	    else if(year == 2 ){                                                                                                                 
    	            if(muon) infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/amcatnlo/DYJetsToLL_M-50_Full_RunII_amcatnlo_m2mu_2017.root");                                                                                                                      
	                else infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src//DYJetsToLL_M-50_Full_RunII_amcatnlo_2017_ElecNoScale.root");                                                                                                                                      
            	}
        	    else{
    	    	      if(muon) infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII//DYJetsToLL_M-50_Full_RunII_m2mu_2018.root");                                                                                          
	    	          else infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/amcatnlo/DYJetsToLL_M-50_Full_RunII_amcatnlo_m2e_2017.root");                                                                                                             
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
         	   tree_MU->SetBranchAddress("pterr1old", &pterr1old);   
				tree_MU->SetBranchAddress("Iso1", &Iso1);   
				tree_MU->SetBranchAddress("genLep_pt1", &genLep_pt1);   
				tree_MU->SetBranchAddress("genLep_eta1", &genLep_eta1);       	
				tree_MU->SetBranchAddress("genLep_phi1", &genLep_phi1);       	
				tree_MU->SetBranchAddress("pT2", &pT2);   
				tree_MU->SetBranchAddress("eta2", &eta2);   
				tree_MU->SetBranchAddress("phi2", &phi2);   
				tree_MU->SetBranchAddress("pterr2", &pterr2);
        	    tree_MU->SetBranchAddress("pterr2old", &pterr2old);   
				tree_MU->SetBranchAddress("Iso2", &Iso2);   
				tree_MU->SetBranchAddress("genLep_pt2", &genLep_pt2);   
				tree_MU->SetBranchAddress("genLep_eta2", &genLep_eta2);  
				tree_MU->SetBranchAddress("genLep_phi2", &genLep_phi2);  
				tree_MU->SetBranchAddress("weight", &weight);
	
				Long64_t nentries_MU = tree_MU->GetEntries();
				std::cout<<nentries_MU<<std::endl;	
	 			for(int entry = 0; entry < nentries_MU; entry++){
//				for(int entry = 0; entry < (int)nentries_MU/100; entry++){
//	 	 		for(int entry = 0; entry < 1000; entry++){

					tree_MU->GetEntry(entry);                                   
					if(entry % 1000000 == 0){         
					if(muon)	std::cout<<entry<<" --- Dentro il TREE "<<year<<" --- MUON "<<year<<std::endl; 
						else 	std::cout<<entry<<" --- Dentro il TREE "<<year<<" --- ELECTRON "<<year<<std::endl; 
					}
				
						if(etaB < 2) sigmaPt = 0.03;
						else sigmaPt = 0.07;
						
						if(fabs(eta1) > eta_bins.at(etaB) && fabs(eta1) < eta_bins.at(etaB+1) && 
							fabs(eta2) > eta_bins.at(etaB) && fabs(eta2) < eta_bins.at(etaB+1) && 
							pterr1/pT1 < sigmaPt && pterr2/pT2 < sigmaPt){
//		 			mass_reco->Fill(massZ);
// 					reco_1.SetPtEtaPhiM(pT1, eta1, phi1, 0.105);
// 					reco_2.SetPtEtaPhiM(pT2, eta2, phi2, 0.105);
// 				float reco_mass = (reco_1+reco_2).M();
							histograms[0][1][year]->Fill(massZ);
							histograms[0][2][year]->Fill(pT1);
							histograms[0][3][year]->Fill(pT2);
							histograms[0][4][year]->Fill(pterr1);
							histograms[0][5][year]->Fill(pterr2);
							histograms[0][6][year]->Fill(eta1);
							histograms[0][6][year]->Fill(eta2);
							histograms[0][7][year]->Fill(phi1);
							histograms[0][7][year]->Fill(phi2);
							histograms[0][8][year]->Fill(pterr1/pT1);
							histograms[0][9][year]->Fill(pterr2/pT2);
							histograms[0][10][year]->Fill(Iso1);
							histograms[0][11][year]->Fill(Iso2);
	            	        histograms[0][12][year]->Fill(pterr1old);
            	    	    histograms[0][13][year]->Fill(pterr2old);
			
// 						mass_gen->Fill(genzm);
							histograms[1][1][year]->Fill(GENmass2l);
							histograms[1][2][year]->Fill(genLep_pt1);
							histograms[1][3][year]->Fill(genLep_pt2);
							histograms[1][6][year]->Fill(genLep_eta1);
							histograms[1][6][year]->Fill(genLep_eta2);
							histograms[1][7][year]->Fill(genLep_phi1);
							histograms[1][7][year]->Fill(genLep_phi2);
						} //for on selection
			
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
	
// 				legend_year->AddEntry(histograms[0][1][0], "2017_942");
				legend_year->AddEntry(histograms[0][1][2], "2017_102_15_NoScale");
				legend_year->AddEntry(histograms[0][1][1], "2017_102_15");
			}

	
// 			directory = "./PlotYearComparison_pTerr_0p5_1p5_" + reconstruction;;
// 			directory = "./PlotYearComparison_pTerr_below0p5_" + reconstruction;;
			directory = "./PlotYearComparison_2017_94X_check";
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


			save_nome = "YearComparison_" + eta_bins_name.at(etaB) + "_" + eta_bins_name.at(etaB+1) + "_";
			if(etaB < 2)
				save_nome+="0p03";
			else	
				save_nome+="0p07";
			Draw(h_blank, h_blank, h_blank, "h_blank", directory + save_nome + ".pdf[", legend_reco_gen);

			for(int i = 0; i < 2; i++){
				for(int var = 0; var < variables.size(); var++){
					TString nome_canvas = reco_gen[i] + variables[var];
					std:;cout<<nome_canvas<<std::endl;
// 					if(var == 0 || var == 4 || var == 5) continue; 
					if(i == 0 && var == 0) continue; 
					if(i == 1 && (var == 0 ||  var == 4 || var == 5 || var > 7)) continue;
					Draw(histograms[i][var][0], histograms[i][var][1], histograms[i][var][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
					Draw(histograms[i][var][0], histograms[i][var][1], histograms[i][var][2], nome_canvas, directory + nome_canvas + ".png", legend_year);
				}
			}
	
			Draw(h_blank, h_blank, h_blank, "h_blank", directory + save_nome + ".pdf]", legend_reco_gen);
	
		} // for on mu = muon or electron
	}// for on sigmaPt
	



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

void Draw(TH1F *h1, TH1F *h2, TH1F *h3, TString nome_canvas, TString save, TLegend *legend){

	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
//    	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.05, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();

	Double_t min, max;
	
// 	h1->Scale(1/h1->Integral());
// 	h2->Scale(1/h2->Integral());
// 	h3->Scale(1/h3->Integral());
	
	if(h1->GetMaximum() > h2->GetMaximum())
		if(h1->GetMaximum() > h3->GetMaximum())
			max = h1->GetMaximum();
		else
			max = h3->GetMaximum();
	else
		if(h2->GetMaximum() > h3->GetMaximum())
			max = h2->GetMaximum();
		else
			max = h3->GetMaximum();
	
	max = 1.1 * max;
	
	if(max == 0) max = 1;
	
	min = 0;
	
	h1->SetTitle(nome_canvas);
	h2->SetTitle(nome_canvas);
	h3->SetTitle(nome_canvas);
	h1->SetLineColor(kGreen+2);
	h2->SetLineColor(kRed+2);
	h3->SetLineColor(kBlue);


	h1->GetYaxis()->SetRangeUser(0, max);//SetMaximum(max);
	h2->GetYaxis()->SetRangeUser(0, max);//SetMaximum(max);
	h3->GetYaxis()->SetRangeUser(0, max);//SetMaximum(max);
	
	canvas->Update();

	
	h1->Draw();
	h2->Draw("Same");
	h3->Draw("Same");
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
	TH1F *ratio_2018 = (TH1F*) h3->Clone();

	legend->Draw();
/*
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
	ratio_2017->GetYaxis()->SetTitle("X / 2016");//Reco / Gen");
	ratio_2017->GetYaxis()->SetTitleOffset(0.50);
	ratio_2017->GetYaxis()->SetNdivisions(506); 
	ratio_2018->Divide(h1);
	ratio_2018->GetYaxis()->SetRangeUser(0.5, 1.5);
	ratio_2018->SetTitle("");
	ratio_2018->SetStats(0);
	ratio_2018->GetYaxis()->SetTitleSize(0.2);
// 	ratio->GetYaxis()->SetTitleFont(43);
	ratio_2018->GetYaxis()->SetLabelSize(0.14);
	ratio_2018->GetYaxis()->SetTitle("X / 2016");//Reco / Gen");
	ratio_2018->GetYaxis()->SetTitleOffset(0.50);
	ratio_2018->GetYaxis()->SetNdivisions(506); 
// 	ratio->GetXaxis()->SetTitle(name_axis);
// 	ratio->GetXaxis()->SetLabelSize(0.15);
// 	ratio->GetXaxis()->SetTitleSize(0.15);
// 	ratio->GetXaxis()->SetTitleOffset(0.75);
// 	ratio->SetLineColor(kBlack);	
	ratio_2017->SetLineColor(kRed+2);
	ratio_2018->SetLineColor(kBlue);
	ratio_2017->Draw();
	ratio_2018->Draw("same");
	canvas->Update();
	pad22->Update();
	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(1);
	line->Draw();
// 	TLegend *legend_comparison = new TLegend(0.75,0.75,0.9,0.9);
// 	legend_comparison->AddEntry(ratio_2017, "2017");
// 	legend_comparison->AddEntry(ratio_2018, "2018");
// 	legend_comparison->Draw();
	*/
	canvas->Print(save);
	
}






