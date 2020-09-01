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
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

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

void YearComparison_ggH(){

	gStyle->SetOptStat(0);
	
	TString file_dir;

	Float_t pTL1, etaL1, phiL1, mL1, pTErrL1;
	Float_t pTL2, etaL2, phiL2, mL2, pTErrL2;
	Float_t pTL3, etaL3, phiL3, mL3, pTErrL3;
	Float_t pTL4, etaL4, phiL4, mL4, pTErrL4;
	Int_t idL1, idL2, idL3, idL4;
	Bool_t passedFullSelection;
	Int_t finalState;
	Float_t mass4l, mass4lErr, mass4lREFIT, mass4lErrREFIT;
	Float_t mass4mu, mass4e, mass2e2mu;
		
	TString reco_gen[2] = {"reco_", "gen_"};
	TString lep_num[4] = {"_1", "_2", "_3", "_4"};
	std::vector<TString> lep_variables;
	lep_variables.push_back("Lep_pT");
	lep_variables.push_back("Lep_eta");
	lep_variables.push_back("Lep_phi");
	lep_variables.push_back("Lep_m");
	lep_variables.push_back("Lep_id");
	lep_variables.push_back("Lep_pTErr");
		
	std::vector<TString> variables;
	variables.push_back("mass4l");
	variables.push_back("mass4lREFIT");
	variables.push_back("mass4mu");
	variables.push_back("mass4e");
	variables.push_back("mass2e2mu");
	variables.push_back("mass4lErr");
	variables.push_back("mass4lErrREFIT");
	for(int j = 0; j < lep_variables.size(); j++){
		for(int i = 0; i < 4; i++)
			variables.push_back(lep_variables.at(j) + lep_num[i]);
	}
	TString years[3] = {"2016", "2017", "2018"};
	
	TH1F* histograms[2][(int)variables.size()][3];
	
	gROOT->Reset();
	gROOT->SetBatch();
	
	for(int i = 0; i < 2; i++){
		for(int var = 0; var < variables.size(); var++){
			for(int year = 0; year < 3; year++){
				TString nome_histo = reco_gen[i] + variables[var] + years[year];
				if(var < 5){
					histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 200, 100, 150);
       	            histograms[i][var][year]->GetXaxis()->SetTitle("m_{ll} [GeV]");
           	    }
				else if(var == 5 || var == 6){
					histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 100, 0, 5);
       	            histograms[i][var][year]->GetXaxis()->SetTitle("#sigma_{m_{4l}} [GeV]");
           	    }
				else if(var >= 7 && var <= 10){
				histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 60, 0, 300);
       	            histograms[i][var][year]->GetXaxis()->SetTitle("p_{T}");
                }
				else if(var >= 11 && var <= 14){
					histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 24, -2.4, 2.4);
           	        histograms[i][var][year]->GetXaxis()->SetTitle("#eta");
                }
				else if(var >= 15 && var <= 18){
					histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 40, -TMath::Pi(), TMath::Pi());
           	        histograms[i][var][year]->GetXaxis()->SetTitle("#phi");
                }
				else if(var >= 19 && var <= 22){
					histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 200, -1, 1);
           	        histograms[i][var][year]->GetXaxis()->SetTitle("m_{ll} [GeV]");
                }
				else if(var >= 23 && var <= 26){
					histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 40, -20, 20);
           	        histograms[i][var][year]->GetXaxis()->SetTitle("ID");
                }
                else{
					histograms[i][var][year] = new TH1F(nome_histo, nome_histo, 100, 0, 5);
           	        histograms[i][var][year]->GetXaxis()->SetTitle("#sigma_{p_{T}}");
				}
			} // for on year
		} // for on variable
	} // for on Reco Gen
			
	for(int year = 0; year < 3; year++){
		
		if(year == 0)
		    infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/ggH/GluGluHToZZTo4L_M125_2016.root");                                                                                                                                         
        else if(year == 1)
		    infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/ggH/GluGluHToZZTo4L_M125_2017.root");                                                                                                                                         
		else
		    infile = new TFile("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/ggH/GluGluHToZZTo4L_M125_2018.root");                                                                                                                                         
	
		if(infile){ 
// 			infile->ls();
			infile->cd("Ana");
// 			infile->ls();
// 			tree_MU = (TTree*) infile->Get("passedEvents");		
			tree_MU = (TTree*)gDirectory->Get("passedEvents");		
		}
		else std::cout<<"ERROR could not find the file"<<std::endl;
	
		Long64_t nentries_MU = tree_MU->GetEntries();
		std::cout<<nentries_MU<<std::endl;	
 		for(int entry = 0; entry < nentries_MU; entry++){
// 		for(int entry = 0; entry < (int)nentries_MU/100; entry++){
// 		for(int entry = 0; entry < 1000; entry++){

			tree_MU->GetEntry(entry);                                   
			if(entry % 100000 == 0)       
				std::cout<<entry<<" --- Dentro il TREE --- "<<year<<std::endl; 
				
				tree_MU->SetBranchAddress("passedFullSelection", &passedFullSelection);       	
				tree_MU->SetBranchAddress("finalState", &finalState);       	
				tree_MU->SetBranchAddress("mass4l", &mass4l);       	
				tree_MU->SetBranchAddress("mass4lErr", &mass4lErr);       	
				tree_MU->SetBranchAddress("mass4lREFIT", &mass4lREFIT);       	
				tree_MU->SetBranchAddress("mass4lErrREFIT", &mass4lErrREFIT);       	
				tree_MU->SetBranchAddress("mass4mu", &mass4mu);       	
				tree_MU->SetBranchAddress("mass4e", &mass4e);       	
				tree_MU->SetBranchAddress("mass2e2mu", &mass2e2mu);       	
				tree_MU->SetBranchAddress("pTL1", &pTL1);       	
				tree_MU->SetBranchAddress("etaL1", &etaL1);       	
				tree_MU->SetBranchAddress("phiL1", &phiL1);       	
				tree_MU->SetBranchAddress("mL1", &mL1);       	
				tree_MU->SetBranchAddress("idL1", &idL1);       	
				tree_MU->SetBranchAddress("pTErrL1", &pTErrL1);       	
				tree_MU->SetBranchAddress("pTL2", &pTL2);       	
				tree_MU->SetBranchAddress("etaL2", &etaL2);       	
				tree_MU->SetBranchAddress("phiL2", &phiL2);       	
				tree_MU->SetBranchAddress("mL2", &mL2);       	
				tree_MU->SetBranchAddress("idL2", &idL2);       	
				tree_MU->SetBranchAddress("pTErrL2", &pTErrL2);       	
				tree_MU->SetBranchAddress("pTL3", &pTL3);       	
				tree_MU->SetBranchAddress("etaL3", &etaL3);       	
				tree_MU->SetBranchAddress("phiL3", &phiL3);       	
				tree_MU->SetBranchAddress("mL3", &mL3);       	
				tree_MU->SetBranchAddress("idL3", &idL3);       	
				tree_MU->SetBranchAddress("pTErrL3", &pTErrL3);       	
				tree_MU->SetBranchAddress("pTL4", &pTL4);       	
				tree_MU->SetBranchAddress("etaL4", &etaL4);       	
				tree_MU->SetBranchAddress("phiL4", &phiL4);       	
				tree_MU->SetBranchAddress("mL4", &mL4);       	
				tree_MU->SetBranchAddress("idL4", &idL4);       	
				tree_MU->SetBranchAddress("pTErrL4", &pTErrL4);       	

				if(!passedFullSelection) continue;

				histograms[0][0][year]->Fill(mass4l);	
				histograms[0][1][year]->Fill(mass4lREFIT);
				if(finalState == 1)
					histograms[0][2][year]->Fill(mass4mu);
				if(finalState == 2)
					histograms[0][3][year]->Fill(mass4e);
				if(finalState > 2)
					histograms[0][4][year]->Fill(mass2e2mu);
				histograms[0][5][year]->Fill(mass4lErr);
				histograms[0][6][year]->Fill(mass4lErrREFIT);
				histograms[0][7][year]->Fill(pTL1);
				histograms[0][8][year]->Fill(pTL2);
				histograms[0][9][year]->Fill(pTL3);
				histograms[0][10][year]->Fill(pTL4);
				histograms[0][11][year]->Fill(etaL1);
				histograms[0][12][year]->Fill(etaL2);
				histograms[0][13][year]->Fill(etaL3);
				histograms[0][14][year]->Fill(etaL4);
				histograms[0][15][year]->Fill(phiL1);
				histograms[0][16][year]->Fill(phiL2);
				histograms[0][17][year]->Fill(phiL3);
				histograms[0][18][year]->Fill(phiL4);
				histograms[0][19][year]->Fill(mL1);
				histograms[0][20][year]->Fill(mL2);
				histograms[0][21][year]->Fill(mL3);
				histograms[0][22][year]->Fill(mL4);
				histograms[0][23][year]->Fill(idL1);
				histograms[0][24][year]->Fill(idL2);
				histograms[0][25][year]->Fill(idL3);
				histograms[0][26][year]->Fill(idL4);
				histograms[0][27][year]->Fill(pTErrL1);
				histograms[0][28][year]->Fill(pTErrL2);
				histograms[0][29][year]->Fill(pTErrL3);
				histograms[0][30][year]->Fill(pTErrL4);					
		
// 					mass_gen->Fill(genzm);
// 					histograms[1][1][year]->Fill(GENmass2l);
// 					histograms[1][2][year]->Fill(genLep_pt1);
// 					histograms[1][3][year]->Fill(genLep_pt2);
// 					histograms[1][6][year]->Fill(genLep_eta1);
// 					histograms[1][6][year]->Fill(genLep_eta2);
// 					histograms[1][7][year]->Fill(genLep_phi1);
// 					histograms[1][7][year]->Fill(genLep_phi2);
			

		} // for on entry
	} // for on year

	
// 	TCanvas *c1 = new TCanvas("c1", "c1", 750, 500);
// 	histograms[0][5][0]->Draw();
// 	TCanvas *c2 = new TCanvas("c2", "c2", 750, 500);
// 	histograms[0][5][1]->Draw();
// 	return;
	
// 	legend_reco_gen->AddEntry(histograms[0][1][0], "Reco");
// 	legend_reco_gen->AddEntry(histograms[1][1][0], "Gen");
	
	legend_year->AddEntry(histograms[0][0][0], "2016");
	legend_year->AddEntry(histograms[0][0][1], "2017");
	legend_year->AddEntry(histograms[0][0][2], "2018");

	
// 		directory = "./PlotYearComparison_pTerr_0p5_1p5_" + reconstruction;
// 		directory = "./PlotYearComparison_pTerr_below0p5_" + reconstruction;
	directory = "./PlotYearComparison_ggH";
	execute = "mkdir " + directory;
	gSystem->Exec(execute);
	directory += "/";

// 	save_nome = "Reco_vs_gen_" + fs;
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


	save_nome = "YearComparison_ggH";
	Draw(h_blank, h_blank, h_blank, "h_blank", directory + save_nome + ".pdf[", legend_reco_gen);

	for(int i = 0; i < 1; i++){
		for(int var = 0; var < variables.size(); var++){
			TString nome_canvas = reco_gen[i] + variables[var];
			std:;cout<<nome_canvas<<std::endl;
// 			if(var == 0 || var == 4 || var == 5) continue; 
// 			if(i == 1) continue; 
// 			if(i == 1 && (var == 0 ||  var == 4 || var == 5 || var > 7)) continue;
			Draw(histograms[i][var][0], histograms[i][var][1], histograms[i][var][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
			Draw(histograms[i][var][0], histograms[i][var][1], histograms[i][var][2], nome_canvas, directory + nome_canvas + ".png", legend_year);
		}

		Draw(h_blank, h_blank, h_blank, "h_blank", directory + save_nome + ".pdf]", legend_reco_gen);
	
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

void Draw(TH1F *h1, TH1F *h2, TH1F *h3, TString nome_canvas, TString save, TLegend *legend){

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
	h3->Scale(1/h3->Integral());
	
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
	TLegend *legend_comparison = new TLegend(0.75,0.75,0.9,0.9);
	legend_comparison->AddEntry(ratio_2017, "2017");
	legend_comparison->AddEntry(ratio_2018, "2018");
	legend_comparison->Draw();

	canvas->Print(save);
	
}