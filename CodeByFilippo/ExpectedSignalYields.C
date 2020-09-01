#include "TFile.h"
#include <iostream>
#include <fstream>
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
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooMyPDF_DSCB.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "ScaleUncertantiesStep2.h"

void SingleFit(TH1F *h_mass, TString canvas_name, int PM, int mP, int y, int EC, int fs);

std::vector<TString> ProdMode;
std::vector<TString> ProdMode_File;
std::vector<int> massPoint;
std::vector<TString> year;
std::vector<TString> category;
std::vector<TString> decayMode;
std::vector<float> luminosity;

void ExpectedSignalYields(){

    gROOT->Reset();
    gROOT->SetBatch();	

	Bool_t passedFullSelection;
	Int_t finalState, EventCat;
	Float_t eventWeight;
	Float_t mass4l, mass4lErr, mass4lREFIT, mass4lErrREFIT;
	Float_t mass4mu, mass4e, mass2e2mu;
	
	std::vector<int> numberEvent_ggH;
	std::vector<int> numberEvent_VBF;
	std::vector<int> numberEvent_WplusH;
	std::vector<int> numberEvent_WminusH;
	std::vector<int> numberEvent_ZH;
	std::vector<int> numberEvent_ttH;
	std::vector<std::vector<int> > numberEvent;	

	std::vector<int> numberEvent_2018;
	numberEvent_2018.push_back(1000000);
	numberEvent_2018.push_back(500000);
	numberEvent_2018.push_back(200000);
	numberEvent_2018.push_back(200000);
	numberEvent_2018.push_back(500000);
	numberEvent_2018.push_back(500000);

	
	std::vector<float> crossSection_ggH;
	crossSection_ggH.push_back(0.01212);
	crossSection_ggH.push_back(0.00792);
	crossSection_ggH.push_back(0.01122);
	crossSection_ggH.push_back(0.01308);
	crossSection_ggH.push_back(0.01706);

	std::vector<float> crossSection_VBF;
	crossSection_VBF.push_back(0.001034);
	crossSection_VBF.push_back(0.0006573);
	crossSection_VBF.push_back(0.0009606);
	crossSection_VBF.push_back(0.001133);
	crossSection_VBF.push_back(0.001506);
	
	std::vector<float> crossSection_WplusH;
	crossSection_WplusH.push_back(0.0002339);
	crossSection_WplusH.push_back(0.0001606);
	crossSection_WplusH.push_back(0.0002190);
	crossSection_WplusH.push_back(0.0002497);
	crossSection_WplusH.push_back(0.0003103);
	
	std::vector<float> crossSection_WminusH;
	crossSection_WminusH.push_back(0.0001471);
	crossSection_WminusH.push_back(0.0001015);
	crossSection_WminusH.push_back(0.0001379);
	crossSection_WminusH.push_back(0.0001570);
	crossSection_WminusH.push_back(0.0001944);
	
	std::vector<float> crossSection_ZH;
	crossSection_ZH.push_back(0.0006569);
	crossSection_ZH.push_back(0.0004467);
	crossSection_ZH.push_back(0.0006145);
	crossSection_ZH.push_back(0.0007026);
	crossSection_ZH.push_back(0.0008778);

	std::vector<float> crossSection_ttH;
	crossSection_ttH.push_back(0.00038991);
	crossSection_ttH.push_back(0.00021383);
	crossSection_ttH.push_back(0.00034867);
	crossSection_ttH.push_back(0.00041842);
	crossSection_ttH.push_back(0.00046059);
	
	std::vector<std::vector<float> > crossSection;	
	crossSection.push_back(crossSection_ggH);
	crossSection.push_back(crossSection_VBF);
	crossSection.push_back(crossSection_WplusH);
	crossSection.push_back(crossSection_WminusH);
	crossSection.push_back(crossSection_ZH);
	crossSection.push_back(crossSection_ttH);

	
	ProdMode.clear();
	ProdMode.push_back("ggH");
	ProdMode.push_back("VBF");
	ProdMode.push_back("WplusH");
	ProdMode.push_back("WminusH");
	ProdMode.push_back("ZH");
	ProdMode.push_back("ttH");
	
	ProdMode_File.clear();
	ProdMode_File.push_back("GluGluHToZZTo4L");
	ProdMode_File.push_back("VBF_HToZZTo4L");
	ProdMode_File.push_back("WplusH_HToZZTo4L");
	ProdMode_File.push_back("WminusH_HToZZTo4L");
	ProdMode_File.push_back("ZH_HToZZ_4LFilter");
	ProdMode_File.push_back("ttH_HToZZ_4LFilter");
	

	massPoint.clear();
	massPoint.push_back(125);
	massPoint.push_back(120);
	massPoint.push_back(124);
	massPoint.push_back(126);
	massPoint.push_back(130);
	
	year.clear();
	year.push_back("2016");
	year.push_back("2017");
	year.push_back("2018");	
	
	luminosity.clear();
	luminosity.push_back(35900);
	luminosity.push_back(41370);
	luminosity.push_back(60900);
	
	category.clear();
	category.push_back("");
// 	category.push_back("VBF_1jet");
// 	category.push_back("VBF_2jet");
// 	category.push_back("VH_leptonic");
// 	category.push_back("VH_hadronic");
// 	category.push_back("ttH");
// 	category.push_back("untagged");
// 	category.push_back("VH_MET");
	
	decayMode.clear();
	decayMode.push_back("4mu");
	decayMode.push_back("4e");
	decayMode.push_back("2e2mu");
	
	Double_t FitParamterValue[year.size()][ProdMode.size()][category.size()][decayMode.size()][3];
// 	Double_t FitParamterValue[year.size()][3] = {0};
	
	
	TH1F* h_yield[year.size()][ProdMode.size()][category.size()][decayMode.size()];
	TH1F* h_mass[year.size()][massPoint.size()][ProdMode.size()][category.size()][decayMode.size()];
	TH1F* h_weight[year.size()][massPoint.size()][ProdMode.size()][category.size()][decayMode.size()];
	for(int y = 0; y < year.size(); y++){
		for(int mP = 0; mP < massPoint.size(); mP++){
			for(int PM = 0; PM < ProdMode.size(); PM++){
				for(int c = 0; c < category.size(); c++){
					for(int fs = 0; fs < decayMode.size(); fs++){
						TString histo_name = Form("HistoYield_%s_%s_%s_%s", year[y].Data(), ProdMode[PM].Data(), category[c].Data(), decayMode[fs].Data());
						h_yield[y][PM][c][fs] = new TH1F(histo_name, histo_name, 20, 115.5, 135.5);
						histo_name = Form("HistoMass_%s_M%i_%s_%s_%s", year[y].Data(), massPoint.at(mP), ProdMode[PM].Data(), category[c].Data(), decayMode[fs].Data());
						h_mass[y][mP][PM][c][fs] = new TH1F(histo_name, histo_name, 100, 105, 140);
						histo_name = Form("HistoWeight_%s_M%i_%s_%s_%s", year[y].Data(), massPoint.at(mP), ProdMode[PM].Data(), category[c].Data(), decayMode[fs].Data());
						h_weight[y][mP][PM][c][fs] = new TH1F(histo_name, histo_name, 100, 0, 0.01);
// 						std::cout<<histo_name<<std::endl;
					}
				}
			}
		}
	}
	
	TLegend *legend_fs = new TLegend(0.8,0.75,0.9,0.9);

	TTree* tree;

	Color_t color;
	TString nome_file;
	
	int ciao = 0;
	
	for(int PM = 0; PM < ProdMode.size(); PM++){
		for(int mP = 0; mP <  massPoint.size(); mP++){
			for(int y = 0; y < year.size(); y++){	
			
				if(PM == 0 && mP == 0){			
					if(y == 0){
						numberEvent_ggH.push_back(1000000);
						numberEvent_ggH.push_back(500000);
						numberEvent_ggH.push_back(1000000);
						numberEvent_ggH.push_back(1000000);
						numberEvent_ggH.push_back(500000);

						numberEvent_VBF.push_back(500000);
						numberEvent_VBF.push_back(200000);
						numberEvent_VBF.push_back(500000);
						numberEvent_VBF.push_back(500000);
						numberEvent_VBF.push_back(200000);	

						numberEvent_WplusH.push_back(300000);
						numberEvent_WplusH.push_back(200000);
						numberEvent_WplusH.push_back(300000);
						numberEvent_WplusH.push_back(300000);
						numberEvent_WplusH.push_back(200000);

						numberEvent_WminusH.push_back(200000);
						numberEvent_WminusH.push_back(120000);
						numberEvent_WminusH.push_back(200000);
						numberEvent_WminusH.push_back(200000);
						numberEvent_WminusH.push_back(120000);

						numberEvent_ZH.push_back(500000);
						numberEvent_ZH.push_back(200000);
						numberEvent_ZH.push_back(500000);
						numberEvent_ZH.push_back(500000);
						numberEvent_ZH.push_back(200000);

						numberEvent_ttH.push_back(500000);
						numberEvent_ttH.push_back(200000);
						numberEvent_ttH.push_back(500000);
						numberEvent_ttH.push_back(500000);
						numberEvent_ttH.push_back(200000);
					}

					if(y == 1){
						numberEvent_ggH.push_back(1000000);
						numberEvent_ggH.push_back(500000);
						numberEvent_ggH.push_back(1000000);
						numberEvent_ggH.push_back(1000000);
						numberEvent_ggH.push_back(500000);

						numberEvent_VBF.push_back(234800);
						numberEvent_VBF.push_back(200000);
						numberEvent_VBF.push_back(500000);
						numberEvent_VBF.push_back(158100);
						numberEvent_VBF.push_back(173300);	

						numberEvent_WplusH.push_back(33400);
						numberEvent_WplusH.push_back(135400);
						numberEvent_WplusH.push_back(300000);
						numberEvent_WplusH.push_back(300000);
						numberEvent_WplusH.push_back(171100);

						numberEvent_WminusH.push_back(200000);
						numberEvent_WminusH.push_back(120000);
						numberEvent_WminusH.push_back(200000);
						numberEvent_WminusH.push_back(200000);
						numberEvent_WminusH.push_back(111200);

						numberEvent_ZH.push_back(19109);
						numberEvent_ZH.push_back(200000);
						numberEvent_ZH.push_back(500000);
						numberEvent_ZH.push_back(500000);
						numberEvent_ZH.push_back(200000);

						numberEvent_ttH.push_back(1000000);
						numberEvent_ttH.push_back(71313);
						numberEvent_ttH.push_back(155626);
						numberEvent_ttH.push_back(406086);
						numberEvent_ttH.push_back(200000);
					}
				numberEvent.push_back(numberEvent_ggH);
				numberEvent.push_back(numberEvent_VBF);
				numberEvent.push_back(numberEvent_WplusH);
				numberEvent.push_back(numberEvent_WminusH);
				numberEvent.push_back(numberEvent_ZH);
				numberEvent.push_back(numberEvent_ttH);	
				
				}

		    	nome_file = Form("/eos/user/f/ferrico/LXPLUS/ROOT_FILE_HZZ/Full_RunII/%s/%s_M%i_%s_skimmed.root", ProdMode[PM].Data(), ProdMode_File[PM].Data(), massPoint.at(mP), year[y].Data());
// 		    	std::cout<<nome_file<<std::endl;
		    	infile = new TFile(nome_file);                                                                                                                                         
				if(infile){ 
// 					infile->cd("Ana");
					tree = (TTree*)gDirectory->Get("passedEvents");		
				}
				else{ 
					std::cout<<"ERROR could not find the file"<<std::endl;
					continue;
				}
				
				Long64_t nentries = tree->GetEntries();
// 				int n_e_125;
// 				if(mP != 0) n_e_125 = nentries;
// 				else n_e_125 = 1;
// 				std::cout<<nentries<<std::endl;	
				std::cout<<nentries<<"\t"<<year.at(y)<<"\t"<<massPoint.at(mP)<<"\t"<<ProdMode.at(PM)<<std::endl; 

				TString cut_name[category.size()][decayMode.size()];
				for(int EC = 0; EC < category.size(); EC++){
					for(int fs = 0; fs < decayMode.size(); fs++){
						if(fs < 2)
// 							cut_name[EC][fs] = Form("(passedFullSelection && EventCat == %d && finalState == %d)", EC+1, fs+1);
// 							cut_name[EC][fs] = Form("eventWeight * (passedFullSelection && EventCat == %d && finalState == %d)", EC+1, fs+1);
// 							cut_name[EC][fs] = Form("eventWeight * (passedFullSelection && finalState == %d)", fs+1);

							cut_name[EC][fs] = Form("crossSection * (passedFullSelection && finalState == %d)", fs+1);
// 							cut_name[EC][fs] = Form("(passedFullSelection && finalState == %d)", fs+1);
						else
// 							cut_name[EC][fs] = Form("(passedFullSelection && EventCat == %d && (finalState == 3 || finalState == 4))", EC+1);
// 							cut_name[EC][fs] = Form("eventWeight * (passedFullSelection && EventCat == %d && (finalState == 3 || finalState == 4))", EC+1);
// 							cut_name[EC][fs] = Form("eventWeight * (passedFullSelection && (finalState == 3 || finalState == 4))");

							cut_name[EC][fs] = Form("crossSection * (passedFullSelection && (finalState == 3 || finalState == 4))");
// 							cut_name[EC][fs] = Form("(passedFullSelection && (finalState == 3 || finalState == 4))");
					}
				}

				for(int EC = 0; EC < category.size(); EC++){
// 				for(int EC = 0; EC < 2; EC++){

					for(int fs = 0; fs < decayMode.size(); fs++){
// 					for(int fs = 1; fs < 2; fs++){
						TString histo_name = Form("HistoYield_%s_%s_%s_%s_%d", year[y].Data(), ProdMode[PM].Data(), category[EC].Data(), decayMode[fs].Data(), massPoint.at(mP));
						TH1F* h_tmp;
						tree->Draw("mass4l>>h_tmp", cut_name[EC][fs],"goff");
// 						std::cout<<cut_name[EC][fs]<<std::endl;
						h_mass[y][mP][PM][EC][fs] = (TH1F*)gDirectory->Get("h_tmp");
// 						std::cout<<h_mass[y][mP][PM][EC][fs]->Integral()<<std::endl;	

						TString nome_canvas = Form("Mass %s %s %i %s %s", ProdMode_File[PM].Data(), category[EC].Data(), massPoint.at(mP), year[y].Data(), decayMode[fs].Data());
						TCanvas *c_mass = new TCanvas(canvas_name, canvas_name, 700, 500);
						c_mass->SetFrameFillColor(0);
						c_mass->cd(1)->SetBottomMargin(0.2);
						h_mass[y][mP][PM][EC][fs]->Draw();
						h_mass[y][mP][PM][EC][fs]->SetTitle(histo_name);
						h_mass[y][mP][PM][EC][fs]->GetXaxis()->SetRangeUser(105, 140);
						c_mass->Update();
						if(y == 0 && mP == 0 && PM == 0 && EC == 0 && fs == 0)
							c_mass->Print("mass.pdf[");
						c_mass->Print("mass.pdf");
						if(y == year.size()-1 && mP == massPoint.size()-1 && PM == ProdMode.size()-1 && EC == category.size()-1 && fs == decayMode.size()-1)
							c_mass->Print("mass.pdf]");

						float integral = luminosity.at(y) * h_mass[y][mP][PM][EC][fs]->Integral() * crossSection[PM][mP];
						if(y != 2) integral /= numberEvent[PM][mP];
						else integral /= numberEvent_2018[PM];

// 						h_yield[y][PM][EC][fs]->SetBinContent(massPoint.at(mP)-115, integral*luminosity.at(y)/nentries);
						h_yield[y][PM][EC][fs]->SetBinContent(massPoint.at(mP)-115, integral);
// 						h_yield[y][PM][EC][fs]->SetBinError(massPoint.at(mP)-115, pow(TMath::Sqrt(integral*luminosity.at(y)/nentries), -1));


// 						TH1F* h_tmp_2;						
// 						tree->Draw("eventWeight>>h_tmp_2", cut_name[EC][fs], "goff");
// 						h_weight[y][mP][PM][EC][fs] = (TH1F*)gDirectory->Get("h_tmp_2");
// 						nome_canvas = Form("Weight %s %s %i %s %s", ProdMode_File[PM].Data(), category[EC].Data(), massPoint.at(mP), year[y].Data(), decayMode[fs].Data());
// 						TCanvas *c_weight = new TCanvas(canvas_name, canvas_name, 700, 500);
// 						c_weight->SetFrameFillColor(0);
// 						c_weight->cd(1)->SetBottomMargin(0.2);
// 						h_weight[y][mP][PM][EC][fs]->Draw();
// 						h_weight[y][mP][PM][EC][fs]->SetTitle(histo_name);
// 						h_weight[y][mP][PM][EC][fs]->GetXaxis()->SetRangeUser(0, 0.01);
// 						c_weight->Update();
// 						if(y == 0 && mP == 0 && PM == 0 && EC == 0 && fs == 0)
// 							c_weight->Print("weight.pdf[");
// 						c_weight->Print("weight.pdf");
// 						if(y == year.size()-1 && mP == massPoint.size()-1 && PM == ProdMode.size()-1 && EC == category.size()-1 && fs == decayMode.size()-1)
// 							c_weight->Print("weight.pdf]");

// 						std::cout<<decayMode.at(fs)<<"\t"<<category.at(EC)<<"\t"<<integral<<"\t"<<h_yield[y][PM][EC][fs]->GetBinContent(massPoint.at(mP)-115)<<std::endl;
					}
				}
/*		 		for(int entry = 0; entry < nentries; entry++){
// 				for(int entry = 0; entry < (int)nentries/100; entry++){
// 				for(int entry = 0; entry < 1000; entry++){

					tree->GetEntry(entry);                                   
					if(entry % 100000 == 0)       
						std::cout<<entry<<" --- Dentro il TREE --- "<<year.at(y)<<"\t"<<massPoint.at(mP)<<"\t"<<ProdMode.at(PM)<<std::endl; 
				
						tree->SetBranchAddress("passedFullSelection", &passedFullSelection);       	
						tree->SetBranchAddress("finalState", &finalState);       	
						tree->SetBranchAddress("mass4l", &mass4l);       	
						tree->SetBranchAddress("mass4lErr", &mass4lErr);       	
						tree->SetBranchAddress("mass4lREFIT", &mass4lREFIT);       	
						tree->SetBranchAddress("mass4lErrREFIT", &mass4lErrREFIT);       	
						tree->SetBranchAddress("mass4mu", &mass4mu);       	
						tree->SetBranchAddress("mass4e", &mass4e);       	
						tree->SetBranchAddress("mass2e2mu", &mass2e2mu);       	
						tree->SetBranchAddress("EventCat", &EventCat);       	
						tree->SetBranchAddress("eventWeight", &eventWeight);       	

						if(!passedFullSelection) continue;
						if(EventCat == 0) continue;	
						
						
// 						eventWeight = eventWeight * n_e_125 / nentries;
// 						eventWeight /= nentries;
						eventWeight = 1;
						
// 						ciao++;
// 						std::cout<<n_e_125<<"\t"<<nentries<<"\t"<<eventWeight<<std::endl;
// 						if(ciao == 10){ 
// 							ciao = 0;
// 							continue;
// 						}
						
						if(finalState > 2)												
							h_mass[y][mP][PM][EventCat-1][2]->Fill(mass4l, eventWeight);
						else
							h_mass[y][mP][PM][EventCat-1][finalState-1]->Fill(mass4l, eventWeight);
				}		
*/
// 				for(int EC = 0; EC < category.size(); EC++){
// 					for(int fs = 0; fs < decayMode.size(); fs++){
						
// 						float integral = 	h_mass[y][mP][PM][EC][fs]->Integral();
// 						TString nome_canvas = Form("Mass %s %s %i %s %s", ProdMode_File[PM].Data(), category[EC].Data(), massPoint.at(mP), year[y].Data(), decayMode[fs].Data());
// // 						SingleFit(h_mass[y][mP][PM][EC][fs], nome_canvas, PM, mP, y, EC, fs);
// 						TCanvas *c_mass = new TCanvas(canvas_name, canvas_name, 700, 500);
// 						c_mass->SetFrameFillColor(0);
// 						c_mass->cd(1)->SetBottomMargin(0.2);
// 						h_mass[y][mP][PM][EC][fs]->Draw();
// 						if(y == 0 && mP == 0 && PM == 0 && EC == 0 && fs == 0)
// 							c_mass->Print("mass.pdf[");
// 						c_mass->Print("mass.pdf");
// 						if(y == year.size()-1 && mP == massPoint.size()-1 && PM == ProdMode.size()-1 && EC == category.size()-1 && fs == decayMode.size()-1)
// 							c_mass->Print("mass.pdf]");
											
// 						h_yield[y][PM][EC][fs]->SetBinContent(massPoint.at(mP)-115, integral);
// 						h_yield[y][mP][PM][EC][fs]->SetBinContent(massPoint.at(mP)-115, integral);
// 						std::cout<<decayMode.at(fs)<<"\t"<<category.at(EC)<<"\t"<<integral<<"\t"<<h_yield[y][PM][EC][fs]->GetBinContent(massPoint.at(mP)-115)<<std::endl;
// 					}
// 				}
		    }
		}
	}

	for(int y = 0; y < year.size(); y++){
		for(int PM = 0; PM < ProdMode.size(); PM++){
			for(int EC = 0; EC < category.size(); EC++){

				if(PM == 0) color = kRed + 2;
				else if(PM == 1) color = kBlue;
				else if(PM == 2) color = kGreen + 2;
				else if(PM == 3) color = kYellow + 2;
				else if(PM == 4) color = kMagenta + 2;
				else color = 1;
				TCanvas* c_yield = new TCanvas("c_yield", "c_yield", 700, 500);	
				float xmax = -999;
				std::cout<<ProdMode.at(PM)<<std::endl;
				for(int mP = 0; mP < massPoint.size(); mP++){
// 					std::cout<<massPoint.at(mP)<<"\t"<<xmax;
					if(h_yield[y][PM][EC][0]->GetBinContent(massPoint.at(mP)-115) > xmax) xmax = h_yield[y][PM][EC][0]->GetBinContent(massPoint.at(mP)-115) * 1.2;//  + h_yield[y][PM][EC][0]->GetBinError(massPoint.at(mP)-115) * 1.2;
// 					std::cout<<"\t("<<h_yield[y][PM][EC][0]->GetBinContent(massPoint.at(mP)-115)<<")\t"<<xmax;
					if(h_yield[y][PM][EC][1]->GetBinContent(massPoint.at(mP)-115) > xmax) xmax = h_yield[y][PM][EC][1]->GetBinContent(massPoint.at(mP)-115) * 1.2;// + h_yield[y][PM][EC][1]->GetBinError(massPoint.at(mP)-115) * 1.2;
// 					std::cout<<"\t("<<h_yield[y][PM][EC][1]->GetBinContent(massPoint.at(mP)-115)<<")\t"<<xmax;
					if(h_yield[y][PM][EC][2]->GetBinContent(massPoint.at(mP)-115) > xmax) xmax = h_yield[y][PM][EC][2]->GetBinContent(massPoint.at(mP)-115) * 1.2;// + h_yield[y][PM][EC][2]->GetBinError(massPoint.at(mP)-115) * 1.2;
// 					std::cout<<"\t("<<h_yield[y][PM][EC][2]->GetBinContent(massPoint.at(mP)-115)<<")\t"<<xmax<<std::endl;
				}

				TString histo_name = Form("HistoYield_%s_%s", year[y].Data(), ProdMode[PM].Data());

				for(int mP = 0; mP < massPoint.size(); mP++){														

					h_yield[y][PM][EC][0]->GetYaxis()->SetRangeUser(0, xmax);
					h_yield[y][PM][EC][1]->GetYaxis()->SetRangeUser(0, xmax);
					h_yield[y][PM][EC][2]->GetYaxis()->SetRangeUser(0, xmax);
					h_yield[y][PM][EC][0]->GetXaxis()->SetTitle("m_{H} (GeV)");
					h_yield[y][PM][EC][1]->GetXaxis()->SetTitle("m_{H} (GeV)");
					h_yield[y][PM][EC][2]->GetXaxis()->SetTitle("m_{H} (GeV)");
					h_yield[y][PM][EC][0]->GetYaxis()->SetTitle("Exp. yield");
					h_yield[y][PM][EC][1]->GetYaxis()->SetTitle("Exp. yield");
					h_yield[y][PM][EC][2]->GetYaxis()->SetTitle("Exp. yield");
					h_yield[y][PM][EC][0]->SetTitle(histo_name);
					h_yield[y][PM][EC][1]->SetTitle(histo_name);
					h_yield[y][PM][EC][2]->SetTitle(histo_name);

					if(mP == 0)
						h_yield[y][PM][EC][0]->Draw("PE");
					else
						h_yield[y][PM][EC][0]->Draw("PE same");
					h_yield[y][PM][EC][0]->SetMarkerStyle(22);
					h_yield[y][PM][EC][0]->SetMarkerSize(0.75);
					h_yield[y][PM][EC][0]->SetMarkerColor(color);
					h_yield[y][PM][EC][0]->SetLineColor(color);
					h_yield[y][PM][EC][0]->SetStats(0);
					h_yield[y][PM][EC][1]->Draw("PE same");
					h_yield[y][PM][EC][1]->SetMarkerStyle(20);
					h_yield[y][PM][EC][1]->SetMarkerSize(0.75);
					h_yield[y][PM][EC][1]->SetMarkerColor(color);
					h_yield[y][PM][EC][1]->SetLineColor(color);
					h_yield[y][PM][EC][1]->SetStats(0);
					h_yield[y][PM][EC][2]->Draw("PE same");
					h_yield[y][PM][EC][2]->SetMarkerStyle(21);
					h_yield[y][PM][EC][2]->SetMarkerSize(0.75);
					h_yield[y][PM][EC][2]->SetMarkerColor(color);
					h_yield[y][PM][EC][2]->SetLineColor(color);
					h_yield[y][PM][EC][2]->SetStats(0);
					
					c_yield->Update();

				}

				TF1* f1 = new TF1("f1","[0]*x*x + [1]*x + [2]", 115.5, 135.5);
				f1->SetLineColor(color);
				f1->SetLineWidth(1);
				

				for(int fs = 0; fs < decayMode.size(); fs++){
					h_yield[y][PM][EC][fs]->Fit("f1");
					f1->SetParNames ("P0","P1","P2");
					FitParamterValue[y][PM][EC][fs][0] = f1->GetParameter("P0");
					FitParamterValue[y][PM][EC][fs][1] = f1->GetParameter("P1");
					FitParamterValue[y][PM][EC][fs][2] = f1->GetParameter("P2");
				}

				if(y == 0 && PM == 0 && EC == 0){
					legend_fs->AddEntry(h_yield[0][0][0][0], decayMode.at(0), "lep");
					legend_fs->AddEntry(h_yield[0][0][0][1], decayMode.at(1), "lep");
					legend_fs->AddEntry(h_yield[0][0][0][2], decayMode.at(2), "lep");
				}
				legend_fs->Draw();

				if(y == 0 && PM == 0 && EC == 0)
					c_yield->Print("yield.pdf[");
				c_yield->Print("yield.pdf");
				if(y == year.size()-1 && PM == ProdMode.size()-1 && EC == category.size()-1)
					c_yield->Print("yield.pdf]");
				
			}
		}
	}
// // 			infile->ls();
// 		}
// 	
// 		Long64_t nentries_MU = tree->GetEntries();
// 		std::cout<<nentries_MU<<std::endl;	
//  		for(int entry = 0; entry < nentries_MU; entry++){
// // 		for(int entry = 0; entry < (int)nentries_MU/100; entry++){
// // 		for(int entry = 0; entry < 1000; entry++){
// 
// 			tree->GetEntry(entry);                                   
// 			if(entry % 100000 == 0)   
// 		
// 		}		
// 	}
// 	
// 	
// 
// 
// 
// 
// 
// 

// 	ofstream myfile;
// 	myfile.open ("File_Parameter_Yield.txt");
// 	for(int y = 0; y < year.size(); y++){	
// 		for(int PM = 0; PM < ProdMode.size(); PM++){
// 			for(int c = 0; c < category.size(); c++){
// 				for(int fs = 0; fs < decayMode.size(); fs++){
// 					for(int p = 0; p < 3; p++){
// // 						std::cout<<year.at(y)<<"\t"<<
// 						myfile<<year.at(y)<<"\t"<<
// 									ProdMode.at(PM)<<"\t"<<
// 									category.at(c)<<"\t"<<
// 									decayMode.at(fs)<<"\t"<<
// 									FitParamterValue[y][PM][c][fs][0]<<" - "<<
// 									FitParamterValue[y][PM][c][fs][1]<<" - "<<
// // 									FitParamterValue[y][PM][c][fs][2]<<std::endl;
// 									FitParamterValue[y][PM][c][fs][2]<<"\n";
// 					}
// 				}
// 			}
// 		}
// 	}
				

}

/*

void SingleFit(TH1F *h_mass, TString canvas_name, int PM, int mP, int y, int EC, int fs){

	RooRealVar x("x", "Mass (GeV/c^{2})", 105, 140);
	RooDataHist h_Higgs("h_Higgs", "h_Higgs", x, h_mass);
	
	RooRealVar DSCB_mean("DSCB_mean", "DSCB_mean", 105, 140);
	RooRealVar DSCB_sigma("DSCB_sigma", "DSCB_sigma", 1, 0, 30);
	RooRealVar DSCB_alphaL("DSCB_alphaL", "DSCB_alphaL", 1, 0, 30);
	RooRealVar DSCB_expL("DSCB_expL", "DSCB_expL", 1, 0, 30);
	RooRealVar DSCB_alphaR("DSCB_alphaR", "DSCB_alphaR", 1, 0, 30);
	RooRealVar DSCB_expR("DSCB_expR", "DSCB_expR", 1, 0, 30);
	RooMyPDF_DSCB DSCB("DSCB", "DSCB", x, DSCB_mean, DSCB_sigma, DSCB_alphaL, DSCB_alphaR, DSCB_expL, DSCB_expR);
	
	TCanvas *c_mass = new TCanvas(canvas_name, canvas_name, 700, 500);
	c_mass->SetFrameFillColor(0);
	c_mass->cd(1)->SetBottomMargin(0.2);
	RooPlot* xframe = x.frame(Title(canvas_name));
	h_Higgs.plotOn(xframe);
	DSCB.fitTo(h_Higgs, Range(105, 140));
	DSCB.plotOn(xframe,RooFit::LineColor(kBlue));
	DSCB.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));
	xframe->Draw();
	TLatex *tex1;
	fit_param_latex = Form("#chi^{2}/dof = %.3f", xframe->chiSquare(6));
	tex1 = new TLatex(0.65,0.85, fit_param_latex);
	tex1->SetNDC();
	tex1->Draw();
// 	RooPlot *framePull_DATA = x.frame("");
// 	framePull_DATA->addObject((TObject*)xframe->pullHist(), "p");
// 	TPad* padPull_DATA =  new  TPad("padPull","padPull",0.,0.,1.,0.2);
// 	padPull_DATA->Draw();
// 	padPull_DATA->cd(0);
// 	framePull_DATA->GetYaxis()->SetLabelSize(0.1);
// 	framePull_DATA->GetXaxis()->SetLabelSize(0.1);
// 	framePull_DATA->SetMinimum(-20.);
// 	framePull_DATA->SetMaximum(20.);
// 	framePull_DATA->Draw();
// 	lineRef_H->Draw("same");


}*/