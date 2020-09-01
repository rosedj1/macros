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

TString directory;
bool verbose;
float pt_cut;


void ScaleUncertantiesStep1(TString fs, bool DATA, TString reconstruction, TString year){

	  if(fs!="2e" && fs!="2mu"){
  		cout<<"fs has to be 2e, or 2mu"<<endl;
	  	return;
	  }	

	TString fs_dir;
	gStyle->SetOptStat(111);
	TRandom3 rand;
	double u1;  
	TString nome_histo; 
	TFile* outfile;
	TFile* infile;
	float chosen_eta;
	float chosen_pT;
	unsigned int ketaBin, kpTBin;
	TString infile_name;
	           
// 	verbose = true;
	verbose = false;
	                                                                                         	
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
    TH1F *h_mass_tot = new TH1F("Mass tot", "Mass tot", 1000, 60 ,120);
    TH1F *h_mass[pT_bins.size()-1][eta_bins.size()-1];
    TH1F *h_pT[pT_bins.size()-1];
    TH1F *h_pT_tot[pT_bins.size()-1];
    TH1F *h_eta[eta_bins.size()-1];
    TH1F *h_eta_tot[eta_bins.size()-1];
    
	TLorentzVector lep_1, lep_2;
	TLorentzVector ZPrime; 



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
	    	h_mass[pt][eta] = new TH1F(nome_histo, nome_histo, 1000, 60, 120);
	    }
    }
    
//    	TH1F* HISTO_MASS;
//     TH1F* HISTO_MASS_LEP;

//     	HISTO_MASS = new TH1F("HISTO_MASS", "HISTO_MASS", 50, 60, 110);
// 	    HISTO_MASS_LEP = new TH1F("HISTO_MASS_LEP", "HISTO_MASS_LEP", 50, 60, 110);
        
    std::cout<<fs<<"\t"<<DATA<<std::endl;
    	
	if(fs=="2mu")
		pt_cut = 5;
	else
		pt_cut = 7;

	
    if(DATA)
		infile_name = "/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/" + reconstruction + "/" + year + "/Data_m" + fs + ".root";
	else
		infile_name = "/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/" + reconstruction + "/DYJetsToLL_M-50_Full_RunII_" + reconstruction + "_m" + fs + "_" + year + ".root";

    infile = new TFile(infile_name);

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
	std::cout<<nentries<<"\t"<<infile_name<<std::endl;	

	for(int i = 0; i < nentries; i++){
// 	for(int i = 0; i < (int)nentries/10; i++){
// 	for(int i = 0; i < 1000; i++){




		tree->GetEntry(i);                                   
		if(i % 1000000 == 0){  
			if(DATA) 	std::cout<<i<<" --- Dentro il TREE 2018 --- DATA -------------- "<<std::endl;        
			if(fs=="2mu")	std::cout<<i<<" --- Dentro il TREE 2018 --- MUON --- "<<pt_cut<<std::endl; 
			else 
					std::cout<<i<<" --- Dentro il TREE 2018 --- ELECTRON --- "<<pt_cut<<std::endl; 
		}
		
		if(massZ < 60 || massZ > 120) continue;

		if(pT1 < pt_cut || pT1 > 100 || pT2 < pt_cut || pT2 > 100) continue;
		
		h_mass_tot->Fill(massZ);
					
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
	
/*	
	if(!DATA){
		TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 750, 500);
   		TPad* pad11 = new TPad("pad1", "pad1", 0, 0.20, 1, 1);
	   	pad11->SetGrid();
   		pad11->SetBottomMargin(0.1);
	   	pad11->Draw();
   		pad11->cd();
	   	pad11->SetTicks();

		Double_t max;
		Double_t min = 0;
		if(HISTO_MASS->GetMaximum() > HISTO_MASS_LEP->GetMaximum())
			max = HISTO_MASS->GetMaximum();
		else
			max = HISTO_MASS_LEP->GetMaximum();	
	
		max = 1.1 * max;

		HISTO_MASS->SetLineColor(kRed);
		HISTO_MASS_LEP->SetLineColor(kBlue);

		HISTO_MASS->Draw();
		HISTO_MASS_LEP->Draw("same");
		
		c_MC->Update();
		c_MC->cd();

		TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.20);
		pad22->SetGrid();
		pad22->Draw();
		pad22->cd();
		pad22->SetTicks();

		TH1F *ratio = (TH1F*) HISTO_MASS->Clone();	
		ratio->Divide(HISTO_MASS_LEP);
		ratio->GetYaxis()->SetRangeUser(0, 2);
		ratio->SetTitle("");
		ratio->SetStats(0);
		ratio->GetYaxis()->SetTitleSize(0.1);
		ratio->GetYaxis()->SetLabelSize(0.14);
		ratio->GetYaxis()->SetTitleOffset(0.50);
		ratio->GetYaxis()->SetNdivisions(506); 
		ratio->SetLineColor(kBlack);
		ratio->Draw();
		c_MC->Update();
		pad22->Update();
	
		TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
		line->SetLineColor(kGreen);
		line->SetLineWidth(1);
		line->Draw();
	}
	else{
		TCanvas *c_DATA = new TCanvas("c_DATA", "c_DATA", 750, 500);
   		TPad* pad11 = new TPad("pad1", "pad1", 0, 0.20, 1, 1);
	   	pad11->SetGrid();
   		pad11->SetBottomMargin(0.1);
	   	pad11->Draw();
   		pad11->cd();
	   	pad11->SetTicks();

		Double_t max;
		Double_t min = 0;
		if(HISTO_MASS->GetMaximum() > HISTO_MASS_LEP->GetMaximum())
			max = HISTO_MASS->GetMaximum();
		else
			max = HISTO_MASS_LEP->GetMaximum();	
	
		max = 1.1 * max;

		HISTO_MASS->SetLineColor(kRed);
		HISTO_MASS_LEP->SetLineColor(kBlue);

		HISTO_MASS->Draw();
		HISTO_MASS_LEP->Draw("same");
		
		c_DATA->Update();
		c_DATA->cd();

		TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.20);
		pad22->SetGrid();
		pad22->Draw();
		pad22->cd();
		pad22->SetTicks();

		TH1F *ratio = (TH1F*) HISTO_MASS->Clone();	
		ratio->Divide(HISTO_MASS_LEP);
		ratio->GetYaxis()->SetRangeUser(0, 2);
		ratio->SetTitle("");
		ratio->SetStats(0);
		ratio->GetYaxis()->SetTitleSize(0.1);
		ratio->GetYaxis()->SetLabelSize(0.14);
		ratio->GetYaxis()->SetTitleOffset(0.50);
		ratio->GetYaxis()->SetNdivisions(506); 
		ratio->SetLineColor(kBlack);
		ratio->Draw();
		c_DATA->Update();
		pad22->Update();
	
		TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
		line->SetLineColor(kGreen);
		line->SetLineWidth(1);
		line->Draw();
	}	
*/	
	

//     if(DATA)
// 	    outfile = new TFile("OutFile_m2" + fs + "_data_inclusivePt.root", "RECREATE");
// 	else
// 	    outfile = new TFile("OutFile_m2" + fs + "_MC_inclusivePt.root", "RECREATE");
//     if(DATA)
// 	    outfile = new TFile("OutFile_m2" + fs + "_data_Differential.root", "RECREATE");
// 	else
// 	    outfile = new TFile("OutFile_m2" + fs + "_MC_Differential.root", "RECREATE");
	}	

	if(fs=="2mu")  fs_dir = "Muon";
	else fs_dir = "Electron";
    outfile = new TFile("Diffential_Scale_Root/" + year + "/" + fs_dir + "/OutFile_m" + fs + "_" + year + "_" + reconstruction + "_Differential_v3.root", "RECREATE");
	if ( outfile->IsOpen() ) printf("File opened successfully\n");
	else {
		std::cout<<"Problem in opening output file."<<std::endl;
		return;
	}

    for(int pt = 0; pt < pT_bins.size()-1; pt ++){
    	nome_histo = Form("pT_%d_%d", (int)pT_bins[pt], (int)pT_bins[pt+1]);
    	h_pT[pt]->Write(nome_histo);
    	nome_histo = Form("pT_Tot_%d_%d", (int)pT_bins[pt], (int)pT_bins[pt+1]);
    	h_pT_tot[pt]->Write(nome_histo);
        
	    for(int eta = 0;eta <  eta_bins.size()-1; eta ++){
	    	if(pt == 0){
    			nome_histo = Form("eta_%s_%s", eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
		    	h_eta[eta]->Write(nome_histo);
    			nome_histo = Form("eta_ToT_a%s_%s", eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
		    	h_eta_tot[eta]->Write(nome_histo);
		    }
	    	nome_histo = Form("Mass_pT_%d_%d_eta_%s_%s", (int)pT_bins[pt], (int)pT_bins[pt+1], eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
	    	h_mass[pt][eta]->Write(nome_histo);
	    	h_mass[pt][eta]->Reset("ICESM");
	    }
    }
    
    h_mass_tot->Write("Mass_Tot");
    
    outfile->Close();


}
