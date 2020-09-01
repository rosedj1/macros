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
TLegend *legend_year = new TLegend(0.7,0.75,0.9,0.9);
TH1F *h_blank = new TH1F("h_blank", "h_blank", 10, -0.5, 9.5);

TFile* infile;
TTree* tree_MU;                                               

bool muon;
TString fs;
TString directory;
TString execute;
TString save_nome;

using namespace std;

void VtxContraint_Studies(){

	gStyle->SetOptStat(0);
	
	TString file_dir;

	std::vector<float> *lep_pt = 0; 
	std::vector<float> *lep_eta = 0;
	std::vector<float> *lep_phi = 0;
	std::vector<float> *lep_mass = 0;
	std::vector<float> *vtxLep_pt = 0;
	std::vector<float> *vtxLep_eta = 0;
	std::vector<float> *vtxLep_phi = 0;
	std::vector<float> *vtxLep_mass = 0;
// 	std::vector<float> *vtxLep_BS_pt = 0;
// 	std::vector<float> *vtxLep_BS_eta = 0;
// 	std::vector<float> *vtxLep_BS_phi = 0;
// 	std::vector<float> *vtxLep_BS_mass = 0;
// 	std::vector<float> *GENlep_pt = 0;
// 	std::vector<float> *GENlep_eta = 0;
// 	std::vector<float> *GENlep_phi = 0;
// 	std::vector<float> *GENlep_mass = 0;

	Bool_t passedFullSelection;
	Int_t finalState;
	Float_t mass4mu, mass4e, mass2e2mu;

	Float_t mass4l, mass4lErr, mass4lREFIT, mass4lErrREFIT;
	Float_t mass4l_vtx, mass4lErr_vtx, mass4lREFIT_vtx, mass4lErrREFIT_vtx;
	Float_t mass4l_vtx_BS, mass4lErr_vtx_BS, mass4lREFIT_vtx_BS, mass4lErrREFIT_vtx_BS;
	
	float GENMH;

	TH1F* One_over_pt[5][3];
	TH1F* Eta[5][3];		
	TH1F* Phi[5][3];
	TH1F* mass[3];
	TH1F* massREFIT[3];

	TH1F* PtDistribution[5][3];		
	TH1F* EtaDistribution[5][3];		
	TH1F* PhiDistribution[5][3];
	TH1F* massDistribution[3];
	TH1F* massREFITDistribution[3];
	

	TString years[3] = {"2018", "2017", "2018"};
	
	gROOT->Reset();
	gROOT->SetBatch();

	TString histo_name;	

	for(int year = 0; year < 3; year++){	

		for(int i = 0; i < 5; i++){
			for(int j = 0; j < 3; j++){
			
				if(j == 0) histo_name = "Base_";
				if(j == 1) histo_name = "Base_vtx_";
				if(j == 2) histo_name = "Base_vtx_BS_";
				
				if(i == 0){
					histo_name += "One_over_pt";					
					One_over_pt[i][j] = new TH1F(histo_name, histo_name, 200, -1.5, 0.5);
					histo_name += "Eta";
					Eta[i][j] = new TH1F(histo_name, histo_name, 100, -5, 5);
					histo_name += "Phi";
					Phi[i][j] = new TH1F(histo_name, histo_name, 100, -5, 5);
					histo_name += "mass";
					mass[j] = new TH1F(histo_name, histo_name, 100, -1, 1);
					histo_name += "mass4lREFIT";
					massREFIT[j] = new TH1F(histo_name, histo_name, 100, -1, 1);

					histo_name += "PtDistribution";					
					PtDistribution[i][j] = new TH1F(histo_name, histo_name, 200, 0, 500);
					histo_name += "EtaDistribution";
					EtaDistribution[i][j] = new TH1F(histo_name, histo_name, 48, -2.4, 2.4);
					histo_name += "PhiDistribution";
					PhiDistribution[i][j] = new TH1F(histo_name, histo_name, 62, -3.14, 3.14);
					histo_name += "massDistribution";
					massDistribution[j] = new TH1F(histo_name, histo_name, 50, 100, 150);
					histo_name += "massREFITDistribution";
					massREFITDistribution[j] = new TH1F(histo_name, histo_name, 50, 100, 150);
					
				}
				else{
					histo_name += Form("One_over_pt_Lep%d", i);
					One_over_pt[i][j] = new TH1F(histo_name, histo_name, 200, -1.5, 0.5);
					histo_name += Form("Eta_Lep%d", i);
					Eta[i][j] = new TH1F(histo_name, histo_name, 100, -5, 5);
					histo_name += Form("Phi_Lep%d", i);
					Phi[i][j] = new TH1F(histo_name, histo_name, 100, -5, 5);

					histo_name += Form("PtDistribution_Lep%d", i);
					PtDistribution[i][j] = new TH1F(histo_name, histo_name, 200, 0, 500);
					histo_name += Form("EtaDistribution_Lep%d", i);
					EtaDistribution[i][j] = new TH1F(histo_name, histo_name, 48, -2.4, 2.4);
					histo_name += Form("PhiDistribution_Lep%d", i);
					PhiDistribution[i][j] = new TH1F(histo_name, histo_name, 62, -3.14, 3.14);


				}		
			}
		}
		
	    TString filename = "Skimmed_" + years[year] + "_skimmed_v2.root";
	    infile = new TFile(filename);
	
		if(infile){
               infile->cd("Ana");
			tree_MU = (TTree*)gDirectory->Get("passedEvents");
        }        
		else std::cout<<"ERROR could not find the file"<<std::endl;
		
		tree_MU->SetBranchAddress("passedFullSelection", &passedFullSelection);       	
		tree_MU->SetBranchAddress("finalState", &finalState);       	
		tree_MU->SetBranchAddress("mass4l", &mass4l);       	
		tree_MU->SetBranchAddress("mass4lErr", &mass4lErr);       	
		tree_MU->SetBranchAddress("mass4lREFIT", &mass4lREFIT);       	
		tree_MU->SetBranchAddress("mass4lErrREFIT", &mass4lErrREFIT);       	
		tree_MU->SetBranchAddress("mass4mu", &mass4mu);       	
		tree_MU->SetBranchAddress("mass4e", &mass4e);       	
		tree_MU->SetBranchAddress("mass2e2mu", &mass2e2mu);       	

		tree_MU->SetBranchAddress("mass4l_vtx", &mass4l_vtx);       	
		tree_MU->SetBranchAddress("mass4lREFIT_vtx", &mass4lREFIT_vtx);       	

		tree_MU->SetBranchAddress("mass4l_vtx_BS", &mass4l_vtx_BS);       	
		tree_MU->SetBranchAddress("mass4lREFIT_vtx_BS", &mass4lREFIT_vtx_BS);       	
      	
   		tree_MU->SetBranchAddress("lep_pt", &lep_pt);       	
   		tree_MU->SetBranchAddress("lep_eta", &lep_eta);       	
   		tree_MU->SetBranchAddress("lep_phi", &lep_phi);       	
   		tree_MU->SetBranchAddress("lep_mass", &lep_mass);       	

   		tree_MU->SetBranchAddress("vtxLep_pt", &vtxLep_pt);
   		tree_MU->SetBranchAddress("vtxLep_eta", &vtxLep_eta);       	
   		tree_MU->SetBranchAddress("vtxLep_phi", &vtxLep_phi);       	
   		tree_MU->SetBranchAddress("vtxLep_mass", &vtxLep_mass);       	

// 		tree_MU->SetBranchAddress("vtxLep_BS_pt", &vtxLep_BS_pt);       	
// 		tree_MU->SetBranchAddress("vtxLep_BS_eta", &vtxLep_BS_eta);       	
// 		tree_MU->SetBranchAddress("vtxLep_BS_phi", &vtxLep_BS_phi);       	
// 		tree_MU->SetBranchAddress("vtxLep_BS_mass", &vtxLep_BS_mass);       	
// 
// 		tree_MU->SetBranchAddress("GENMH", &GENMH);       	
// 		tree_MU->SetBranchAddress("GENlep_pt", &GENlep_pt);       	
// 		tree_MU->SetBranchAddress("GENlep_eta", &GENlep_eta);       	
// 		tree_MU->SetBranchAddress("GENlep_phi", &GENlep_phi);       	
// 		tree_MU->SetBranchAddress("GENlep_mass", &GENlep_mass);  

		Long64_t nentries_MU = tree_MU->GetEntries();

		std::cout<<nentries_MU<<std::endl;			
		
 		for(int entry = 0; entry < nentries_MU; entry++){
// 		for(int entry = 0; entry < (int)nentries_MU/100; entry++){
// 		for(int entry = 0; entry < 1000; entry++){

			tree_MU->GetEntry(entry);  

			if(entry % 1000 == 0)       
				std::cout<<entry<<" --- Dentro il TREE --- "<<year<<std::endl;      	
   					
			if(!passedFullSelection || finalState != 1) continue;
			
				mass[0]->Fill((mass4l-GENMH)/GENMH);
				mass[1]->Fill((mass4l_vtx-GENMH)/GENMH);
				mass[2]->Fill((mass4l_vtx_BS-GENMH)/GENMH);
												
				massREFIT[0]->Fill((mass4lREFIT-GENMH)/GENMH);
				massREFIT[1]->Fill((mass4lREFIT_vtx-GENMH)/GENMH);
				massREFIT[2]->Fill((mass4lREFIT_vtx_BS-GENMH)/GENMH);

				massDistribution[0]->Fill(mass4l);
				massDistribution[1]->Fill(mass4l_vtx);
				massDistribution[2]->Fill(mass4l_vtx_BS);
												
				massREFITDistribution[0]->Fill(mass4lREFIT);
				massREFITDistribution[1]->Fill(mass4lREFIT_vtx);
				massREFITDistribution[2]->Fill(mass4lREFIT_vtx_BS);

				std::cout<<"lep\t"<<(*lep_pt)[0]<<"\t"<<lep_eta->at(0)<<"\t"<<lep_phi->at(0)<<std::endl;
				std::cout<<"\t\t1 vtxlep\t"<<float((*vtxLep_pt)[0])<<"\t"<<vtxLep_eta->at(0)<<"\t"<<vtxLep_phi->at(0)<<std::endl;
// 				std::cout<<"2 vtxlep\t"<<(*vtxLep_pt)[1]<<"\t"<<vtxLep_eta->at(1)<<"\t"<<vtxLep_phi->at(1)<<std::endl;
// 				std::cout<<"3 vtxlep\t"<<(*vtxLep_pt)[2]<<"\t"<<vtxLep_eta->at(2)<<"\t"<<vtxLep_phi->at(2)<<std::endl;
// 				std::cout<<"4 vtxlep\t"<<(*vtxLep_pt)[3]<<"\t"<<vtxLep_eta->at(3)<<"\t"<<vtxLep_phi->at(3)<<std::endl;

// 				for(int lep = 0; lep < 4; lep++){
// 					std::cout<<"lep\t"<<(*lep_pt)[lep]<<"\t"<<lep_eta->at(lep)<<"\t"<<lep_phi->at(lep)<<std::endl;
// 					std::cout<<"vtxlep\t"<<(*vtxLep_pt)[lep]<<"\t"<<vtxLep_eta->at(lep)<<"\t"<<vtxLep_phi->at(lep)<<std::endl;
// 					std::cout<<"vtxlepBS\t"<<(*vtxLep_BS_pt)[lep]<<"\t"<<vtxLep_BS_eta->at(lep)<<"\t"<<vtxLep_BS_phi->at(lep)<<std::endl;
				
/*					PtDistribution[0][0]->Fill(lep_pt->at(lep));
					PtDistribution[lep+1][0]->Fill(lep_pt->at(lep));
					PtDistribution[0][1]->Fill(vtxLep_pt->at(lep));
					PtDistribution[lep+1][1]->Fill(vtxLep_pt->at(lep));
					PtDistribution[0][2]->Fill(vtxLep_BS_pt->at(lep));
					PtDistribution[lep+1][2]->Fill(vtxLep_BS_pt->at(lep));

					EtaDistribution[0][0]->Fill(lep_eta->at(lep));
					EtaDistribution[lep+1][0]->Fill(lep_eta->at(lep));
					EtaDistribution[0][1]->Fill(vtxLep_eta->at(lep));
					EtaDistribution[lep+1][1]->Fill(vtxLep_eta->at(lep));
					EtaDistribution[0][2]->Fill(vtxLep_BS_eta->at(lep));
					EtaDistribution[lep+1][2]->Fill(vtxLep_BS_eta->at(lep));

					PhiDistribution[0][0]->Fill(lep_phi->at(lep));
					PhiDistribution[lep+1][0]->Fill(lep_phi->at(lep));
					PhiDistribution[0][1]->Fill(vtxLep_phi->at(lep));
					PhiDistribution[lep+1][1]->Fill(vtxLep_phi->at(lep));
					PhiDistribution[0][2]->Fill(vtxLep_BS_phi->at(lep));
					PhiDistribution[lep+1][2]->Fill(vtxLep_BS_phi->at(lep));


// 					One_over_pt[0][0]->Fill((1/lep_pt->at(lep) - 1/GENlep_pt->at(lep))/(1/GENlep_pt->at(lep)));			
// 					One_over_pt[lep+1][0]->Fill((1/lep_pt->at(lep)) - (1/GENlep_pt->at(lep))/(1/GENlep_pt->at(lep)));			
					One_over_pt[0][0]->Fill((lep_pt->at(lep) - GENlep_pt->at(lep))/GENlep_pt->at(lep));			
					One_over_pt[lep+1][0]->Fill((lep_pt->at(lep) - GENlep_pt->at(lep))/GENlep_pt->at(lep));			

					Eta[0][0]->Fill((lep_eta->at(lep) - GENlep_eta->at(lep))/GENlep_eta->at(lep));
					Eta[lep+1][0]->Fill((lep_eta->at(lep) - GENlep_eta->at(lep))/GENlep_eta->at(lep));

					Phi[0][0]->Fill((lep_phi->at(lep) - GENlep_phi->at(lep))/GENlep_phi->at(lep));
					Phi[lep+1][0]->Fill((lep_phi->at(lep) - GENlep_phi->at(lep))/GENlep_phi->at(lep));



// 					One_over_pt[0][1]->Fill((1/vtxLep_pt->at(lep)) - (1/GENlep_pt->at(lep))/(1/GENlep_pt->at(lep)));			
// 					One_over_pt[lep+1][1]->Fill((1/vtxLep_pt->at(lep)) - (1/GENlep_pt->at(lep))/(1/GENlep_pt->at(lep)));			
					One_over_pt[0][0]->Fill((vtxLep_pt->at(lep) - GENlep_pt->at(lep))/GENlep_pt->at(lep));			
					One_over_pt[lep+1][0]->Fill((vtxLep_pt->at(lep) - GENlep_pt->at(lep))/GENlep_pt->at(lep));			

					Eta[0][1]->Fill((vtxLep_eta->at(lep) - GENlep_eta->at(lep))/GENlep_eta->at(lep));
					Eta[lep+1][1]->Fill((vtxLep_eta->at(lep) - GENlep_eta->at(lep))/GENlep_eta->at(lep));

					Phi[0][1]->Fill((vtxLep_phi->at(lep) - GENlep_phi->at(lep))/GENlep_phi->at(lep));
					Phi[lep+1][1]->Fill((vtxLep_phi->at(lep) - GENlep_phi->at(lep))/GENlep_phi->at(lep));
					

				
// 					One_over_pt[0][2]->Fill((1/vtxLep_BS_pt->at(lep)) - (1/GENlep_pt->at(lep))/(1/GENlep_pt->at(lep)));			
// 					One_over_pt[lep+1][2]->Fill((1/vtxLep_BS_pt->at(lep)) - (1/GENlep_pt->at(lep))/(1/GENlep_pt->at(lep)));			
					One_over_pt[0][0]->Fill((vtxLep_BS_pt->at(lep) - GENlep_pt->at(lep))/GENlep_pt->at(lep));			
					One_over_pt[lep+1][0]->Fill((vtxLep_BS_pt->at(lep) - GENlep_pt->at(lep))/GENlep_pt->at(lep));			

					Eta[0][2]->Fill((vtxLep_BS_eta->at(lep) - GENlep_eta->at(lep))/GENlep_eta->at(lep));
					Eta[lep+1][2]->Fill((vtxLep_BS_eta->at(lep) - GENlep_eta->at(lep))/GENlep_eta->at(lep));

					Phi[0][2]->Fill((vtxLep_BS_phi->at(lep) - GENlep_phi->at(lep))/GENlep_phi->at(lep));
					Phi[lep+1][2]->Fill((vtxLep_BS_phi->at(lep) - GENlep_phi->at(lep))/GENlep_phi->at(lep));
				}			
*/
		} // for on entry


		if(year == 0){
			legend_year->AddEntry(One_over_pt[0][0], "Baseline");
			legend_year->AddEntry(One_over_pt[0][1], "Common Vertex");
			legend_year->AddEntry(One_over_pt[0][2], "Common Vertex: BS");
		}

		directory = "./VtxContraint_Studies";
		execute = "mkdir " + directory;
		gSystem->Exec(execute);
		directory += "/";
		execute = "mkdir " + directory + years[year];
		gSystem->Exec(execute);
		directory += years[year];
		directory += "/";	

		save_nome = "VtxContraint_Studies";
		TString nome_canvas;
		Draw(h_blank, h_blank, h_blank, "h_blank", directory + save_nome + ".pdf[", legend_reco_gen);
		for(int i = 0; i < 5; i++){
			std::cout<<i<<std::endl;
			if(i == 0){
				nome_canvas = "One_over_pt";					
				Draw(One_over_pt[0][0], One_over_pt[0][1], One_over_pt[0][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = "Eta";
				Draw(Eta[0][0], Eta[0][1], Eta[0][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = "Phi";
				Draw(Phi[0][0], Phi[0][1], Phi[0][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = "Mass";
				Draw(massDistribution[0], massDistribution[1], massDistribution[2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = "MassREFITTED";
				Draw(massREFIT[0], massREFIT[1], massREFIT[2], nome_canvas, directory + save_nome + ".pdf", legend_year);

				nome_canvas = "PtDistribution";					
				Draw(PtDistribution[0][0], PtDistribution[0][1], PtDistribution[0][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = "EtaDistribution";
				Draw(EtaDistribution[0][0], EtaDistribution[0][1], EtaDistribution[0][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = "PhiDistribution";
				Draw(PhiDistribution[0][0], PhiDistribution[0][1], PhiDistribution[0][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = "MassDistribution";
				Draw(mass[0], mass[1], mass[2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = "massREFITDistribution";
				Draw(massREFITDistribution[0], massREFITDistribution[1], massREFITDistribution[2], nome_canvas, directory + save_nome + ".pdf", legend_year);
			}
			else{
				nome_canvas = Form("One_over_pt_Lep%d", i);
				Draw(One_over_pt[i][0], One_over_pt[i][1], One_over_pt[i][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = Form("Eta_Lep%d", i);
				Draw(Eta[i][0], Eta[i][1], Eta[i][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = Form("Phi_Lep%d", i);
				Draw(Phi[i][0], Phi[i][1], Phi[i][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = Form("PtDistribution_Lep%d", i);
				Draw(PtDistribution[i][0], PtDistribution[i][1], PtDistribution[i][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = Form("EtaDistribution_Lep%d", i);
				Draw(EtaDistribution[i][0], EtaDistribution[i][1], EtaDistribution[i][2], nome_canvas, directory + save_nome + ".pdf", legend_year);
				nome_canvas = Form("PhiDistribution_Lep%d", i);
				Draw(PhiDistribution[i][0], PhiDistribution[i][1], PhiDistribution[i][2], nome_canvas, directory + save_nome + ".pdf", legend_year);

			}		
		}
		Draw(h_blank, h_blank, h_blank, "h_blank", directory + save_nome + ".pdf]", legend_reco_gen);


	} // for on year
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

	TH1F *ratio_2017 = (TH1F*) h2->Clone();
	TH1F *ratio_2018 = (TH1F*) h3->Clone();

	canvas->Update();
	canvas->cd();

	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.25);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	
	ratio_2017->Divide(h1);
	ratio_2017->GetYaxis()->SetRangeUser(0.5, 1.5);
	ratio_2017->SetTitle("");
	ratio_2017->SetStats(0);
	ratio_2017->GetYaxis()->SetTitleSize(0.1);
	ratio_2017->GetYaxis()->SetLabelSize(0.14);
	ratio_2017->GetYaxis()->SetTitle("X / 2016");//Reco / Gen");
	ratio_2017->GetYaxis()->SetTitleOffset(0.50);
	ratio_2017->GetYaxis()->SetNdivisions(506); 
	ratio_2018->Divide(h1);
	ratio_2018->GetYaxis()->SetRangeUser(0.5, 1.5);
	ratio_2018->SetTitle("");
	ratio_2018->SetStats(0);
	ratio_2018->GetYaxis()->SetTitleSize(0.2);
	ratio_2018->GetYaxis()->SetLabelSize(0.14);
	ratio_2018->GetYaxis()->SetTitle("X / 2016");//Reco / Gen");
	ratio_2018->GetYaxis()->SetTitleOffset(0.50);
	ratio_2018->GetYaxis()->SetNdivisions(506); 
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
