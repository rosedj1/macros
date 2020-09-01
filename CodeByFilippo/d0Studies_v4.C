#include "ScaleUncertantiesStep2.h"
void Proiezione(TH2F* h2, TString save_name, int MC, int type, int fs, int j);
void Proiezione_deltaPt(int maximum, TH2F* h2, TString save_name, int type, int fs);
void Proiezione_deltaPt(int maximum, TH3F* h2, TString save_name, int type, int fs, int pt, int isPT);
void Draw(TH1F *h1, TH1F *h2, TString nome_canvas, TString save, TLegend *legend, TString x_name, bool mcDATA, float min = 0.1);
void DrawMCData(TH1F *hMC, TH1F *hData, TString nome_canvas, TString save, TLegend *legend, TString x_name, float weight, float min);
float FitMass(TH1D* h1, TString save_name, int fit, int save_control_1, int save_control_2);
TH2F* MakeMean(TH3F* h2, TString save_name, int fit, TH2F* h2_integral, float counting, int n_bins, float x_bounder);


TGraphErrors* Mean_vs_d0[2][2][2][3][2];
TGraphErrors* DeltaPt_vs_d0[2][2];
std::vector<TGraphErrors*> DeltaPt_vs_d0_inPt[2][2][3];
std::vector<TGraphErrors*> DeltaPt_vs_d0_inEta[2][2][3];
std::vector<TGraphErrors*> DeltaPt_vs_d0_inEta_Above20[2][2][3];
std::vector<TGraphErrors*> DeltaPt_vs_d0_inEta_Above20_Below100[2][2][3];
std::vector<TGraphErrors*> DeltaPt_vs_d0_inEta_Above20_Below200[2][2][3];
std::vector<Float_t> mean_vs_d0_tmp;
std::vector<Float_t> mean_vs_d0err_tmp;
std::vector<Float_t> mean_vs_d0_tmp2;
std::vector<Float_t> mean_vs_d0err_tmp2;


// Float_t d0_BINS[21] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
// 						10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
// 						20};
// Float_t d0_BINSerr[21] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
// 							0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
// 							0.5};
// Int_t  binnum_d0 = sizeof(d0_BINS)/sizeof(Float_t)-1;  

TLegend *legend_recoGen = new TLegend(0.75,0.75,0.9,0.9);
TLegend *legend_DataMC = new TLegend(0.75,0.75,0.9,0.9);
TLegend *legend = new TLegend(0.75,0.75,0.9,0.9);
TLegend *legend_pt = new TLegend(0.75,0.75,0.9,0.9);
TLegend *legend_eta = new TLegend(0.73,0.75,0.9,0.9);

  
void d0Studies_v4(int maximum){

	gROOT->Reset();
	gROOT->SetBatch();

	gStyle->SetOptStat(111111);
	
	TString deltapt_name = "10000 * (p_{T}^{reco} - p_{T}^{GEN}) / (p_{T}^{GEN})^{2}";

	float type_lep1, type_lep2;

	std::vector<TString> finalStates;
	finalStates.push_back("2mu");
	finalStates.push_back("2e");
	
	std::vector<TString> d0type;
	d0type.push_back("d0BS");
	d0type.push_back("d0PV");
	
	std::vector<int> number_bins;
	number_bins.push_back(20);
	number_bins.push_back(20);

	std::vector<float> extreme_BS;
	extreme_BS.push_back(0.01);
	extreme_BS.push_back(0.05);

	std::vector<float> extreme_PV;
	extreme_PV.push_back(0.01);
	extreme_PV.push_back(0.01);
	
	std::vector<Double_t> eta_bins;
	eta_bins.push_back(0);
	eta_bins.push_back(0.8);
	eta_bins.push_back(1.8);
	eta_bins.push_back(2.4);
	
	std::vector<TString> eta_bins_name;
	eta_bins_name.push_back("0p0");
	eta_bins_name.push_back("0p8");
	eta_bins_name.push_back("1p8");
	eta_bins_name.push_back("2p4");
	
	std::vector<Double_t> pT_bins;
	pT_bins.push_back(5);
	pT_bins.push_back(20);
	pT_bins.push_back(30);
	pT_bins.push_back(40);
	pT_bins.push_back(50);
	pT_bins.push_back(60);
	pT_bins.push_back(100);	
	pT_bins.push_back(250);	
	pT_bins.push_back(1000);	
	
	std::vector<TString> pT_bins_name;
	pT_bins_name.push_back("5");
	pT_bins_name.push_back("20");
	pT_bins_name.push_back("30");
	pT_bins_name.push_back("40");
	pT_bins_name.push_back("50");
	pT_bins_name.push_back("60");
	pT_bins_name.push_back("100");
	pT_bins_name.push_back("250");
	pT_bins_name.push_back("1000");

	TH2F* h_mass_vs_d0[2][2][2][3];
	TH1F* h_pt_reco[2][2][3];
	TH2F* h_pt_reco_etaBins[2][2];
	TH1F* h_pt_gen[2][3];
	TH2F* h_pt_gen_etaBins[2];
	TH1F* h_mass[2][2];
	TH1F* h_d0[2][2][2][3];
	TH2F* h_deltaPt_vs_d0[2][2];
	TH3F* h_deltaPt_vs_d0_inPt[2][2][3];
	TH3F* h_deltaPt_vs_d0_inEta[2][2][3];
	TH3F* h_deltaPt_vs_d0_inEta_Above20[2][2][3];
	TH3F* h_deltaPt_vs_d0_inEta_Above20_Below100[2][2][3];
	TH3F* h_deltaPt_vs_d0_inEta_Above20_Below200[2][2][3];
	TH1F* h_eta_reco[2][2][3];
	TH1F* h_eta_gen[2][3];
		
	TH3F* h_massZ_NEG_vs_POS[2][2][2];
	TH2F* h_integral_NEG_vs_POS[2][2][2];
	TH2F* h_mean_NEG_vs_POS[2][2][2][2];

	TH3F* h_massZ_NEG_vs_POS_gen[2][2];
	TH2F* h_integral_NEG_vs_POS_gen[2][2];
	TH2F* h_mean_NEG_vs_POS_gen[2][2][2];

	TString year;
	std::vector<TString> fitFunction;
	fitFunction.push_back("BWxCR");
	fitFunction.push_back("VOIG");
	
	std::vector<TString> monteCarlo;
	monteCarlo.push_back("MC");
	monteCarlo.push_back("Data");
	
	float charge_1, charge_2;
	float delta_pt1;
	float delta_pt2;

// 	float luminosity[3] = {35900, 41370, 60900};	
	float luminosity[3] = {33280, 40268, 55030};	
// 	float XS[3] = {5670, 7181, 7181};
	float XS[3] = {6070, 6070, 6070};
	float events[3] = {static_cast<float>(118166590), static_cast<float>(179040237), static_cast<float>(121896034)};

	float weight[2][3] = {0};
	
	float counting;
	

	for(int y = 0; y < 3; y++){
	
		if(y != 1) continue;	

		if(y == 0) year = "2016";
		else if(y == 1) year = "2017";
		else year = "2018";
		
		TString directory;
		if(maximum == 1) directory = "./prova_d0Studies_maximum_d0charge/" + year;
		else directory = "./prova_d0Studies_d0charge/" + year;
// 		directory = "./cancella/" + year;
		gSystem->Exec("mkdir " + directory);
		directory += "/";

		std::cout<<"Anno = "<<year<<std::endl;
		
		std::vector<Double_t> BS_BINS;
		std::vector<Double_t> PV_BINS;
		for(int i = 0; i < number_bins.at(0)+1; i++){
			BS_BINS.push_back(-extreme_BS.at(0) + 2 * i * extreme_BS.at(0)/number_bins.at(0));
			PV_BINS.push_back(-extreme_PV.at(0) + 2 * i * extreme_PV.at(0)/number_bins.at(0));
		}

		std::vector<Double_t> DELTAPT_BINS;
		for(int i = 0; i < 201; i++){
				DELTAPT_BINS.push_back(-500 + 2 * i * 500/200);
		}

		std::vector<Double_t> DELTAPT_BINS_MAXIMUM;
		for(int i = 0; i < 501; i++){
				DELTAPT_BINS_MAXIMUM.push_back(-50 + 2 * i * 50/500);
		}
		
		for(int MC = 0; MC < 2; MC++){
			for(int type = 0; type < 2; type++){	
				for(int fs = 0; fs < 2; fs++){	
					for(int charge = 0; charge < 3; charge++){
						TString NOME;
						if(charge == 0) NOME = Form("h_mass_vs_%s_%s_%s_%s", d0type[type].Data(), finalStates[fs].Data(), monteCarlo[MC].Data(), year.Data());
						else if(charge == 1) NOME = Form("h_mass_vs_%s_%s_NEG_%s_%s", d0type[type].Data(), finalStates[fs].Data(), monteCarlo[MC].Data(), year.Data());
						else NOME = Form("h_mass_vs_%s_%s_POS_%s_%s", d0type[type].Data(), finalStates[fs].Data(), monteCarlo[MC].Data(), year.Data());
						
						if(type == 0)
							h_mass_vs_d0[MC][type][fs][charge] = new TH2F(NOME, NOME, number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs), 50, 86, 96);//300, 60, 120);
						else
							h_mass_vs_d0[MC][type][fs][charge] = new TH2F(NOME, NOME, number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs), 50, 86, 96);//300, 60, 120);
					}
				}
			}
		}
		
		for(int type = 0; type < 2; type++){	
			for(int fs = 0; fs < 2; fs++){
					TString NOME;
					NOME = Form("h_deltaPt_vs_%s_%s_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());

					if(type == 0)
						h_deltaPt_vs_d0[type][fs] = new TH2F(NOME, NOME, number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs), 200, -500, 500);
					else
						h_deltaPt_vs_d0[type][fs] = new TH2F(NOME, NOME, number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs), 200, -500, 500);
			}
		}	
		
		for(int type = 0; type < 2; type++){	
			for(int fs = 0; fs < 2; fs++){
				for(int charge = 0; charge < 3; charge++){
					TString NOME;
					if(charge == 0)
						NOME = Form("h_deltaPt_vs_%s_%s_inPt_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());
					if(charge == 1)	
						NOME = Form("h_deltaPt_vs_%s_%s_inPt_NEG_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());
					else
						NOME = Form("h_deltaPt_vs_%s_%s_inPt_POS_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());

					if(maximum == 0){
						if(type == 0)
							h_deltaPt_vs_d0_inPt[type][fs][charge] = new TH3F(NOME, NOME, pT_bins.size()-1, &pT_bins[0], BS_BINS.size()-1, &BS_BINS[0], DELTAPT_BINS.size()-1, &DELTAPT_BINS[0]);
						else
							h_deltaPt_vs_d0_inPt[type][fs][charge] = new TH3F(NOME, NOME, pT_bins.size()-1, &pT_bins[0], PV_BINS.size()-1, &PV_BINS[0], DELTAPT_BINS.size()-1, &DELTAPT_BINS[0]);
// 						if(type == 0)
// 							h_deltaPt_vs_d0_inPt[type][fs][charge] = new TH3F(NOME, NOME, 40, -0.05, 0.05, pT_bins.size()-1, &pT_bins[0], DELTAPT_BINS.size()-1, &DELTAPT_BINS[0]);
// 						else
// 							h_deltaPt_vs_d0_inPt[type][fs][charge] = new TH3F(NOME, NOME, 40, -0.05, 0.05, pT_bins.size()-1, &pT_bins[0], DELTAPT_BINS.size()-1, &DELTAPT_BINS[0]);
					}
					else{
						if(type == 0)
							h_deltaPt_vs_d0_inPt[type][fs][charge] = new TH3F(NOME, NOME, pT_bins.size()-1, &pT_bins[0], BS_BINS.size()-1, &BS_BINS[0], DELTAPT_BINS_MAXIMUM.size()-1, &DELTAPT_BINS_MAXIMUM[0]);
						else
							h_deltaPt_vs_d0_inPt[type][fs][charge] = new TH3F(NOME, NOME, pT_bins.size()-1, &pT_bins[0], PV_BINS.size()-1, &PV_BINS[0], DELTAPT_BINS_MAXIMUM.size()-1, &DELTAPT_BINS_MAXIMUM[0]);
// 						if(type == 0)
// 							h_deltaPt_vs_d0_inPt[type][fs][charge] = new TH3F(NOME, NOME, 40, -0.05, 0.05, pT_bins.size()-1, &pT_bins[0], DELTAPT_BINS_MAXIMUM.size()-1, &DELTAPT_BINS_MAXIMUM[0]);
// 						else
// 							h_deltaPt_vs_d0_inPt[type][fs][charge] = new TH3F(NOME, NOME, 40, -0.05, 0.05, pT_bins.size()-1, &pT_bins[0], DELTAPT_BINS_MAXIMUM.size()-1, &DELTAPT_BINS_MAXIMUM[0]);
					}
				}
			}
		}
	
		for(int type = 0; type < 2; type++){	
			for(int fs = 0; fs < 2; fs++){
				for(int charge = 0; charge < 3; charge++){
					TString NOME;
					if(charge == 0)
						NOME = Form("h_deltaPt_vs_%s_%s_inEta_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());
					if(charge == 1)	
						NOME = Form("h_deltaPt_vs_%s_%s_inEta_NEG_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());
					else
						NOME = Form("h_deltaPt_vs_%s_%s_inEta_POS_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());

					if(maximum == 0){
						if(type == 0)
							h_deltaPt_vs_d0_inEta[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), BS_BINS.size()-1, &BS_BINS[0], DELTAPT_BINS.size()-1, &DELTAPT_BINS[0]);
						else
							h_deltaPt_vs_d0_inEta[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), PV_BINS.size()-1, &PV_BINS[0], DELTAPT_BINS.size()-1, &DELTAPT_BINS[0]);
					}
					else{
						if(type == 0)
							h_deltaPt_vs_d0_inEta[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), BS_BINS.size()-1, &BS_BINS[0], DELTAPT_BINS_MAXIMUM.size()-1, &DELTAPT_BINS_MAXIMUM[0]);
						else
						h_deltaPt_vs_d0_inEta[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), PV_BINS.size()-1, &PV_BINS[0], DELTAPT_BINS_MAXIMUM.size()-1, &DELTAPT_BINS_MAXIMUM[0]);
					}
				}
			}
		}

		for(int type = 0; type < 2; type++){	
			for(int fs = 0; fs < 2; fs++){
				for(int charge = 0; charge < 3; charge++){
					TString NOME;
					if(charge == 0)
						NOME = Form("h_deltaPt_vs_%s_%s_inEta_Above20_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());
					if(charge == 1)	
						NOME = Form("h_deltaPt_vs_%s_%s_inEta_Above20_NEG_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());
					else
						NOME = Form("h_deltaPt_vs_%s_%s_inEta_Above20_POS_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());

					if(maximum == 0){
						if(type == 0)
							h_deltaPt_vs_d0_inEta_Above20[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), BS_BINS.size()-1, &BS_BINS[0], DELTAPT_BINS.size()-1, &DELTAPT_BINS[0]);
						else
							h_deltaPt_vs_d0_inEta_Above20[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), PV_BINS.size()-1, &PV_BINS[0], DELTAPT_BINS.size()-1, &DELTAPT_BINS[0]);
					}
					else{
						if(type == 0)
							h_deltaPt_vs_d0_inEta_Above20[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), BS_BINS.size()-1, &BS_BINS[0], DELTAPT_BINS_MAXIMUM.size()-1, &DELTAPT_BINS_MAXIMUM[0]);
						else
							h_deltaPt_vs_d0_inEta_Above20[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), PV_BINS.size()-1, &PV_BINS[0], DELTAPT_BINS_MAXIMUM.size()-1, &DELTAPT_BINS_MAXIMUM[0]);
					}
				}
			}
		}

		for(int type = 0; type < 2; type++){	
			for(int fs = 0; fs < 2; fs++){
				for(int charge = 0; charge < 3; charge++){
					TString NOME;
					if(charge == 0)
						NOME = Form("h_deltaPt_vs_%s_%s_inEta_Above20_Below200_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());
					if(charge == 1)	
						NOME = Form("h_deltaPt_vs_%s_%s_inEta_Above20_Below200_NEG_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());
					else
						NOME = Form("h_deltaPt_vs_%s_%s_inEta_Above20_Below200_POS_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());

					if(maximum == 0){
						if(type == 0)
							h_deltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), BS_BINS.size()-1, &BS_BINS[0], DELTAPT_BINS.size()-1, &DELTAPT_BINS[0]);
						else
							h_deltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), PV_BINS.size()-1, &PV_BINS[0], DELTAPT_BINS.size()-1, &DELTAPT_BINS[0]);
					}
					else{
						if(type == 0)
							h_deltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), BS_BINS.size()-1, &BS_BINS[0], DELTAPT_BINS_MAXIMUM.size()-1, &DELTAPT_BINS_MAXIMUM[0]);
						else
							h_deltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), PV_BINS.size()-1, &PV_BINS[0], DELTAPT_BINS_MAXIMUM.size()-1, &DELTAPT_BINS_MAXIMUM[0]);
					}
				}
			}
		}
		
		for(int type = 0; type < 2; type++){	
			for(int fs = 0; fs < 2; fs++){
				for(int charge = 0; charge < 3; charge++){
					TString NOME;
					if(charge == 0)
						NOME = Form("h_deltaPt_vs_%s_%s_inEta_Above20_Below100_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());
					if(charge == 1)	
						NOME = Form("h_deltaPt_vs_%s_%s_inEta_Above20_Below100_NEG_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());
					else
						NOME = Form("h_deltaPt_vs_%s_%s_inEta_Above20_Below100_POS_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());

					if(maximum == 0){
						if(type == 0)
							h_deltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), BS_BINS.size()-1, &BS_BINS[0], DELTAPT_BINS.size()-1, &DELTAPT_BINS[0]);
						else
							h_deltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), PV_BINS.size()-1, &PV_BINS[0], DELTAPT_BINS.size()-1, &DELTAPT_BINS[0]);
					}
					else{
						if(type == 0)
							h_deltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), BS_BINS.size()-1, &BS_BINS[0], DELTAPT_BINS_MAXIMUM.size()-1, &DELTAPT_BINS_MAXIMUM[0]);
						else	
							h_deltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge] = new TH3F(NOME, NOME, eta_bins.size()-1, &(eta_bins[0]), PV_BINS.size()-1, &PV_BINS[0], DELTAPT_BINS_MAXIMUM.size()-1, &DELTAPT_BINS_MAXIMUM[0]
						);
					}
				}
			}
		}
		
		for(int MC = 0; MC < 2; MC++){	
			for(int type = 0; type < 2; type++){	
				for(int fs = 0; fs < 2; fs++){	
					for(int charge = 0; charge < 3; charge++){
						TString NOME;
						if(charge == 0) NOME = Form("h_d0_%s_%s_%s_%s", d0type[type].Data(), finalStates[fs].Data(), monteCarlo[MC].Data(), year.Data());
						else if(charge == 1) NOME = Form("h_d0_%s_%s_Leading_%s_%s", d0type[type].Data(), finalStates[fs].Data(), monteCarlo[MC].Data(), year.Data());
						else NOME = Form("h_d0_%s_%s_Subleading_%s_%s", d0type[type].Data(), finalStates[fs].Data(), monteCarlo[MC].Data(), year.Data());
					
						if(type == 0)
							h_d0[MC][type][fs][charge] = new TH1F(NOME, NOME, number_bins.at(fs)*40, -extreme_BS.at(fs)*20, extreme_BS.at(fs)*20);
						else
							h_d0[MC][type][fs][charge] = new TH1F(NOME, NOME, number_bins.at(fs)*40, -extreme_PV.at(fs)*20, extreme_PV.at(fs)*20);
					}
				}
			}
		}

		for(int MC = 0; MC < 2; MC++){		
			for(int fs = 0; fs < 2; fs++){
				TString NOME = "h_pt_reco_" + finalStates.at(fs) + "_" + monteCarlo.at(MC) + "_" + year;
				h_pt_reco[MC][fs][0] = new TH1F(NOME, NOME, 125, 0, 250);
				NOME = "h_pt_reco_" + finalStates.at(fs) + "_" + monteCarlo.at(MC) + "_Leading" + "_" + year;
				h_pt_reco[MC][fs][1] = new TH1F(NOME, NOME, 125, 0, 250);
				NOME = "h_pt_reco_" + finalStates.at(fs) + "_" + monteCarlo.at(MC) + "_Subleading" + "_" + year;
				h_pt_reco[MC][fs][2] = new TH1F(NOME, NOME, 125, 0, 250);

				NOME = "h_eta_reco_" + finalStates.at(fs) + "_" + monteCarlo.at(MC) + "_" + year;
				h_eta_reco[MC][fs][0] = new TH1F(NOME, NOME, 48, -2.4, 2.4);
				NOME = "h_eta_reco_" + finalStates.at(fs) + "_" + monteCarlo.at(MC) + "_Leading" + "_" + year;
				h_eta_reco[MC][fs][1] = new TH1F(NOME, NOME, 48, -2.4, 2.4);
				NOME = "h_eta_reco_" + finalStates.at(fs) + "_" + monteCarlo.at(MC) + "_Subleading" + "_" + year;
				h_eta_reco[MC][fs][2] = new TH1F(NOME, NOME, 48, -2.4, 2.4);

				NOME = "h_mass_" + finalStates.at(fs) + "_" + monteCarlo.at(MC) + "_" + year;
				h_mass[MC][fs] = new TH1F(NOME, NOME, 150, 60, 170);
			
				if(MC == 0){ 
					NOME = "h_pt_gen_" + finalStates.at(fs) + "_" + year;
					h_pt_gen[fs][0] = new TH1F(NOME, NOME, 125, 0, 250);
					NOME = "h_pt_gen_" + finalStates.at(fs) + "_Leading" + "_" + year;
					h_pt_gen[fs][1] = new TH1F(NOME, NOME, 125, 0, 250);
					NOME = "h_pt_gen_" + finalStates.at(fs) + "_Subleading" + "_" + year;
					h_pt_gen[fs][2] = new TH1F(NOME, NOME, 125, 0, 250);
					
					NOME = "h_eta_gen_" + finalStates.at(fs) + "_" + year;
					h_eta_gen[fs][0] = new TH1F(NOME, NOME, 48, -2.4, 2.4);
					NOME = "h_eta_gen_" + finalStates.at(fs) + "_Leading" + "_" + year;
					h_eta_gen[fs][1] = new TH1F(NOME, NOME, 48, -2.4, 2.4);
					NOME = "h_eta_gen_" + finalStates.at(fs) + "_Subleading" + "_" + year;
					h_eta_gen[fs][2] = new TH1F(NOME, NOME, 48, -2.4, 2.4);

				}
			}	
		}

		for(int MC = 0; MC < 2; MC++){
			for(int type = 0; type < 2; type++){	
				for(int fs = 0; fs < 2; fs++){	
					TString NOME = Form("h_massZ_NEG_vs_POS_%s_%s_%s_%s", d0type[type].Data(), finalStates[fs].Data(), monteCarlo[MC].Data(), year.Data());				
					if(type == 0)
						h_massZ_NEG_vs_POS[MC][type][fs] = new TH3F(NOME, NOME, number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs), number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs), 50, 86, 96);//100, 60, 120);
					else
						h_massZ_NEG_vs_POS[MC][type][fs] = new TH3F(NOME, NOME, number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs), number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs), 50, 86, 96);//100, 60, 120);

					NOME = Form("h_integral_NEG_vs_POS_%s_%s_%s_%s", d0type[type].Data(), finalStates[fs].Data(), monteCarlo[MC].Data(), year.Data());				
					if(type == 0)
						h_integral_NEG_vs_POS[MC][type][fs] = new TH2F(NOME, NOME, number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs), number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs));
					else
						h_integral_NEG_vs_POS[MC][type][fs] = new TH2F(NOME, NOME, number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs), number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs));
					for(int f = 0; f < 2; f++){
						NOME = Form("h_mean_NEG_vs_POS_%s_%s_%s_%s%s", d0type[type].Data(), finalStates[fs].Data(), monteCarlo[MC].Data(), year.Data(), fitFunction[f].Data());							
						if(type == 0)
							h_mean_NEG_vs_POS[MC][type][fs][f] = new TH2F(NOME, NOME, number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs), number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs));
						else
							h_mean_NEG_vs_POS[MC][type][fs][f] = new TH2F(NOME, NOME, number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs), number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs));
					}
				}
			}
		}
		
		for(int type = 0; type < 2; type++){	
			for(int fs = 0; fs < 2; fs++){	
				TString NOME;
				NOME = Form("h_massZgen_NEG_vs_POS_gen_%s_%s_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());				
				if(type == 0)
					h_massZ_NEG_vs_POS_gen[type][fs] = new TH3F(NOME, NOME, number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs), number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs), 50, 86, 96);//100, 60, 120);
				else
					h_massZ_NEG_vs_POS_gen[type][fs] = new TH3F(NOME, NOME, number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs), number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs), 50, 86, 96);//100, 60, 120);

				NOME = Form("h_integralGen_NEG_vs_POS_gen_%s_%s_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data());				
				if(type == 0)
					h_integral_NEG_vs_POS_gen[type][fs] = new TH2F(NOME, NOME, number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs), number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs));
				else
					h_integral_NEG_vs_POS_gen[type][fs] = new TH2F(NOME, NOME, number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs), number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs));
				for(int f = 0; f < 2; f++){
					NOME = Form("h_meanGen_NEG_vs_POS_gen_%s_%s_%s_%s", d0type[type].Data(), finalStates[fs].Data(), year.Data(), fitFunction[f].Data());							
					if(type == 0)
						h_mean_NEG_vs_POS_gen[type][fs][f] = new TH2F(NOME, NOME, number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs), number_bins.at(fs), -extreme_BS.at(fs), extreme_BS.at(fs));
					else
						h_mean_NEG_vs_POS_gen[type][fs][f] = new TH2F(NOME, NOME, number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs), number_bins.at(fs), -extreme_PV.at(fs), extreme_PV.at(fs));
				}
			}
		}

		for(int MC = 0; MC < 2; MC++){		
			for(int fs = 0; fs < 2; fs++){
// 				for(int eta = 0; eta < eta_bins.size(); eta++){
// 					TString NOME = Form("h_pt_reco_%s_%s_%s_%s_%s", finalStates[fs].Data(), monteCarlo[MC].Data(), year.Data());
					TString NOME = Form("h_pt_EtaBins_reco_%s_%s_%s", finalStates[fs].Data(), monteCarlo[MC].Data(), year.Data());
					h_pt_reco_etaBins[MC][fs] = new TH2F(NOME, NOME, 125, 0, 250, eta_bins.size()-1, &(eta_bins[0]));				

					if(MC == 0){ 
// 						NOME = Form("h_pt_gen_%s_%s_%s_%s", finalStates[fs].Data(), year.Data());
						NOME = Form("h_pt_EtaBins_gen_%s_%s", finalStates[fs].Data(), year.Data());
						h_pt_gen_etaBins[fs] = new TH2F(NOME, NOME, 125, 0, 250, eta_bins.size()-1, &(eta_bins[0]));				
					}
// 				}
			}
		}
	
// 		for(int type = 0; type < 2; type++){	
		for(int MC = 0; MC < 2; MC++){
			for(int fs = 0; fs < 1; fs++){	
			
				if(fs == 1) continue;				
				if(MC == 1) continue;
				
				TString nome_file;
				if(MC == 0)
					nome_file = "/eos/user/f/ferrico/LXPLUS/ROOT_FILE_HZZ/Full_RunII/d0studies/DY/DYJetsToLL_M-50_Full_RunII_d0studies_m" + finalStates.at(fs) + "_" + year + ".root";
				else
					nome_file = "/eos/user/f/ferrico/LXPLUS/ROOT_FILE_HZZ/Full_RunII/d0studies/Data/" + year + "/Data_m" + finalStates.at(fs) + ".root";
// 					nome_file = "/eos/user/f/ferrico/LXPLUS/ROOT_FILE_HZZ/Full_RunII/Data/" + year + "/Data_m" + finalStates.at(fs) + ".root";
					 
				infile = new TFile(nome_file);
				std::cout<<"NOME FILE = "<<nome_file<<std::endl;

				tree = (TTree*) infile->Get("passedEvents");
				tree->SetBranchAddress("massZ", &massZ);       	
				tree->SetBranchAddress("Id1", &Id1);       	
				tree->SetBranchAddress("Id2", &Id2);       	
				tree->SetBranchAddress("d0BS1", &d0BS1);       	
				tree->SetBranchAddress("d0BS2", &d0BS2);       	
				tree->SetBranchAddress("d0PV1", &d0PV1);       	
				tree->SetBranchAddress("d0PV2", &d0PV2);
				tree->SetBranchAddress("pT1", &pT1);       	
				tree->SetBranchAddress("pT2", &pT2);       	
				tree->SetBranchAddress("eta1", &eta1);       	
				tree->SetBranchAddress("eta2", &eta2);       	
				if(MC == 0){
					tree->SetBranchAddress("genLep_pt1", &genLep_pt1);
					tree->SetBranchAddress("genLep_pt2", &genLep_pt2);
					tree->SetBranchAddress("genLep_eta1", &genLep_eta1);
					tree->SetBranchAddress("genLep_eta2", &genLep_eta2);
					tree->SetBranchAddress("GENmass2l", &GENmass2l);
					
				}
			
			
		
				Long64_t nentries = tree->GetEntries();
				std::cout<<nentries<<std::endl;
				
				if(MC == 0) { 
					weight[fs][y] = luminosity[y] * XS[y] / events[y];
					std::cout<<luminosity[y]<<"\t"<<XS[y]<<"\t"<<events[y]<<"\t"<<nentries<<std::endl;
				}
				else
					weight[fs][y] = 1;

// 				std::cout<<luminosity[y]*XS[y]/nentries<<std::endl;
				std::cout<<weight[fs][y]<<std::endl;

// 				for(int entry = 0; entry < nentries; entry++){
// 				for(int entry = 0; entry < (int)nentries/10; entry++){
// 				for(int entry = 0; entry < (int)nentries/100; entry++){
				for(int entry = 0; entry < 1000; entry++){

					tree->GetEntry(entry);                                   
					if(entry % 1000000 == 0)         
						std::cout<<entry<<" --- Dentro il TREE \t"<<finalStates.at(fs)<<std::endl;
				
					if(fs == 0 && (pT1 < 5 || pT2 < 5)) continue;
					if(fs == 1 && (pT1 < 7 || pT2 < 7)) continue;
					
					if(massZ < 50) continue;
					
					counting++;
							
					h_mass[MC][fs]->Fill(massZ, weight[fs][y]);
					
					h_pt_reco[MC][fs][0]->Fill(pT1, weight[fs][y]);
					h_pt_reco[MC][fs][0]->Fill(pT2, weight[fs][y]);
					h_pt_reco[MC][fs][1]->Fill(pT1, weight[fs][y]);
					h_pt_reco[MC][fs][2]->Fill(pT2, weight[fs][y]);
					
					h_eta_reco[MC][fs][0]->Fill(eta1, weight[fs][y]);
					h_eta_reco[MC][fs][0]->Fill(eta2, weight[fs][y]);
					h_eta_reco[MC][fs][1]->Fill(eta1, weight[fs][y]);
					h_eta_reco[MC][fs][2]->Fill(eta2, weight[fs][y]);
					
					if(MC == 0){						
						h_pt_gen[fs][0]->Fill(genLep_pt1, weight[fs][y]);
						h_pt_gen[fs][0]->Fill(genLep_pt2, weight[fs][y]);
						h_pt_gen[fs][1]->Fill(genLep_pt1, weight[fs][y]);
						h_pt_gen[fs][2]->Fill(genLep_pt2, weight[fs][y]);

						h_eta_gen[fs][0]->Fill(genLep_eta1, weight[fs][y]);
						h_eta_gen[fs][0]->Fill(genLep_eta2, weight[fs][y]);
						h_eta_gen[fs][1]->Fill(genLep_eta1, weight[fs][y]);
						h_eta_gen[fs][2]->Fill(genLep_eta2, weight[fs][y]);
					}

					h_d0[MC][0][fs][0]->Fill(d0BS1, weight[fs][y]);
					h_d0[MC][0][fs][0]->Fill(d0BS2, weight[fs][y]);
					h_d0[MC][0][fs][1]->Fill(d0BS1, weight[fs][y]);
					h_d0[MC][0][fs][2]->Fill(d0BS2, weight[fs][y]);

					h_d0[MC][1][fs][0]->Fill(d0PV1, weight[fs][y]);
					h_d0[MC][1][fs][0]->Fill(d0PV2, weight[fs][y]);
					h_d0[MC][1][fs][1]->Fill(d0PV1, weight[fs][y]);
					h_d0[MC][1][fs][2]->Fill(d0PV2, weight[fs][y]);

					h_pt_reco_etaBins[MC][fs]->Fill(pT1, fabs(eta1), weight[fs][y]);
					h_pt_reco_etaBins[MC][fs]->Fill(pT2, fabs(eta2), weight[fs][y]);
					if(MC == 0){
						h_pt_gen_etaBins[fs]->Fill(genLep_pt1, fabs(eta1), weight[fs][y]);
						h_pt_gen_etaBins[fs]->Fill(genLep_pt2, fabs(eta2), weight[fs][y]);
					}
			
					if(Id1 > 0) charge_1 = -1;
					else charge_1 = 1;

					if(Id2 > 0) charge_2 = -1;
					else charge_2 = 1;
					
					if(MC == 0){
						delta_pt1 = 10000 * (pT1 - genLep_pt1) / (genLep_pt1 * genLep_pt1);
						delta_pt2 = 10000 * (pT2 - genLep_pt2) / (genLep_pt2 * genLep_pt2);
					}
					////// d0 from BS //////										
// 					std::cout<<"BS1 = "<<d0BS1<<"\t"<<Id1<<"\t"<<type_lep1<<std::endl;
// 					std::cout<<"BS2 = "<<d0BS2<<"\t"<<Id2<<"\t"<<type_lep2<<std::endl;			

// 					type_lep1 = d0BS1;
// 					type_lep2 = d0BS2;
					type_lep1 = d0BS1*charge_1;
					type_lep2 = d0BS2*charge_2;	

					h_mass_vs_d0[MC][0][fs][0]->Fill(type_lep1, massZ);
					h_mass_vs_d0[MC][0][fs][0]->Fill(type_lep2, massZ);			
					
					
					if(Id1 > 0){
						h_massZ_NEG_vs_POS[MC][0][fs]->Fill(type_lep2, type_lep1, massZ);
						h_mass_vs_d0[MC][0][fs][1]->Fill(type_lep1, massZ);
						h_mass_vs_d0[MC][0][fs][2]->Fill(type_lep2, massZ);
					}
					else{
						h_massZ_NEG_vs_POS[MC][0][fs]->Fill(type_lep1, type_lep2, massZ);
						h_mass_vs_d0[MC][0][fs][1]->Fill(type_lep2, massZ);
						h_mass_vs_d0[MC][0][fs][2]->Fill(type_lep1, massZ);
					}
					
					if(MC == 0){
					
						h_deltaPt_vs_d0[0][fs]->Fill(type_lep1, delta_pt1);
						h_deltaPt_vs_d0[0][fs]->Fill(type_lep2, delta_pt2);				

						if(Id1 > 0){
							h_massZ_NEG_vs_POS_gen[0][fs]->Fill(type_lep2, type_lep1, GENmass2l);
							h_deltaPt_vs_d0_inPt[0][fs][2]->Fill(pT2, type_lep2, delta_pt2);
							h_deltaPt_vs_d0_inPt[0][fs][0]->Fill(pT2, type_lep2, delta_pt2);							
							h_deltaPt_vs_d0_inEta[0][fs][2]->Fill(fabs(eta2), type_lep2, delta_pt2);
							h_deltaPt_vs_d0_inEta[0][fs][0]->Fill(fabs(eta2), type_lep2, delta_pt2);
							if(pT2 > 20){
								h_deltaPt_vs_d0_inEta_Above20[0][fs][2]->Fill(fabs(eta2), type_lep2, delta_pt2);
								h_deltaPt_vs_d0_inEta_Above20[0][fs][0]->Fill(fabs(eta2), type_lep2, delta_pt2);
								if(pT2 < 200){
									h_deltaPt_vs_d0_inEta_Above20_Below200[0][fs][2]->Fill(fabs(eta2), type_lep2, delta_pt2);
									h_deltaPt_vs_d0_inEta_Above20_Below200[0][fs][0]->Fill(fabs(eta2), type_lep2, delta_pt2);
								}
								if(pT2 < 100){
									h_deltaPt_vs_d0_inEta_Above20_Below100[0][fs][2]->Fill(fabs(eta2), type_lep2, delta_pt2);
									h_deltaPt_vs_d0_inEta_Above20_Below100[0][fs][0]->Fill(fabs(eta2), type_lep2, delta_pt2);
								}
							}
						}
						else{
							h_massZ_NEG_vs_POS_gen[0][fs]->Fill(type_lep1, type_lep2, GENmass2l);
							h_deltaPt_vs_d0_inPt[0][fs][1]->Fill(pT1, type_lep1, delta_pt1);
							h_deltaPt_vs_d0_inPt[0][fs][0]->Fill(pT1, type_lep1, delta_pt1);
							h_deltaPt_vs_d0_inEta[0][fs][1]->Fill(fabs(eta2), type_lep2, delta_pt2);
							h_deltaPt_vs_d0_inEta[0][fs][0]->Fill(fabs(eta2), type_lep2, delta_pt2);
							if(pT2 > 20){
								h_deltaPt_vs_d0_inEta_Above20[0][fs][1]->Fill(fabs(eta1), type_lep1, delta_pt1);
								h_deltaPt_vs_d0_inEta_Above20[0][fs][0]->Fill(fabs(eta1), type_lep1, delta_pt1);
								if(pT2 < 200){
									h_deltaPt_vs_d0_inEta_Above20_Below200[0][fs][1]->Fill(fabs(eta1), type_lep1, delta_pt1);
									h_deltaPt_vs_d0_inEta_Above20_Below200[0][fs][0]->Fill(fabs(eta1), type_lep1, delta_pt1);
								}
								if(pT2 < 100){
									h_deltaPt_vs_d0_inEta_Above20_Below100[0][fs][1]->Fill(fabs(eta1), type_lep1, delta_pt1);
									h_deltaPt_vs_d0_inEta_Above20_Below100[0][fs][0]->Fill(fabs(eta1), type_lep1, delta_pt1);
								}
							}
						}
					}
					////// d0 from BS //////					
						
						
					
					////// d0 from PV //////					
// 					std::cout<<"PV1 = "<<d0PV1<<"\t"<<Id1<<"\t"<<type_lep1<<std::endl;
// 					std::cout<<"PV2 = "<<d0PV2<<"\t"<<Id2<<"\t"<<type_lep2<<std::endl;
					
// 					type_lep1 = d0PV1;
// 					type_lep2 = d0PV2;
					type_lep1 = d0PV1*charge_1;
					type_lep2 = d0PV2*charge_2;	

					h_mass_vs_d0[MC][1][fs][0]->Fill(type_lep1, massZ);
					h_mass_vs_d0[MC][1][fs][0]->Fill(type_lep2, massZ);			

					if(Id1 > 0){
						h_massZ_NEG_vs_POS[MC][1][fs]->Fill(type_lep1, type_lep2, massZ);
						h_mass_vs_d0[MC][1][fs][1]->Fill(type_lep1, massZ);
						h_mass_vs_d0[MC][1][fs][2]->Fill(type_lep2, massZ);
					}
					else{
						h_massZ_NEG_vs_POS[MC][1][fs]->Fill(type_lep2, type_lep1, massZ);
						h_mass_vs_d0[MC][1][fs][1]->Fill(type_lep2, massZ);
						h_mass_vs_d0[MC][1][fs][2]->Fill(type_lep1, massZ);
					}
	
					if(MC == 0){					

						h_deltaPt_vs_d0[1][fs]->Fill(type_lep1, delta_pt1);
						h_deltaPt_vs_d0[1][fs]->Fill(type_lep2, delta_pt2);

						if(Id1 > 0){
							h_massZ_NEG_vs_POS_gen[1][fs]->Fill(type_lep2, type_lep1, GENmass2l);
							h_deltaPt_vs_d0_inPt[1][fs][2]->Fill(pT2, type_lep2, delta_pt2);
							h_deltaPt_vs_d0_inPt[1][fs][0]->Fill(pT2, type_lep2, delta_pt2);							
							h_deltaPt_vs_d0_inEta[1][fs][2]->Fill(fabs(eta2), type_lep2, delta_pt2);
							h_deltaPt_vs_d0_inEta[1][fs][0]->Fill(fabs(eta2), type_lep2, delta_pt2);
							if(pT2 > 20){
								h_deltaPt_vs_d0_inEta_Above20[1][fs][2]->Fill(fabs(eta2), type_lep2, delta_pt2);
								h_deltaPt_vs_d0_inEta_Above20[1][fs][0]->Fill(fabs(eta2), type_lep2, delta_pt2);
								if(pT2 < 200){
									h_deltaPt_vs_d0_inEta_Above20_Below200[1][fs][2]->Fill(fabs(eta2), type_lep2, delta_pt2);
									h_deltaPt_vs_d0_inEta_Above20_Below200[1][fs][0]->Fill(fabs(eta2), type_lep2, delta_pt2);
								}
								if(pT2 < 100){
									h_deltaPt_vs_d0_inEta_Above20_Below100[1][fs][2]->Fill(fabs(eta2), type_lep2, delta_pt2);
									h_deltaPt_vs_d0_inEta_Above20_Below100[1][fs][0]->Fill(fabs(eta2), type_lep2, delta_pt2);
								}
							}
						}
						else{
							h_massZ_NEG_vs_POS_gen[1][fs]->Fill(type_lep1, type_lep2, GENmass2l);
							h_deltaPt_vs_d0_inPt[1][fs][1]->Fill(pT1, type_lep1, delta_pt1);
							h_deltaPt_vs_d0_inPt[1][fs][0]->Fill(pT1, type_lep1, delta_pt1);
							h_deltaPt_vs_d0_inEta[1][fs][1]->Fill(fabs(eta2), type_lep2, delta_pt2);
							h_deltaPt_vs_d0_inEta[1][fs][0]->Fill(fabs(eta2), type_lep2, delta_pt2);
							if(pT2 > 20){
								h_deltaPt_vs_d0_inEta_Above20[1][fs][1]->Fill(fabs(eta1), type_lep1, delta_pt1);
								h_deltaPt_vs_d0_inEta_Above20[1][fs][0]->Fill(fabs(eta1), type_lep1, delta_pt1);
								if(pT2 < 200){
									h_deltaPt_vs_d0_inEta_Above20_Below200[1][fs][1]->Fill(fabs(eta1), type_lep1, delta_pt1);
									h_deltaPt_vs_d0_inEta_Above20_Below200[1][fs][0]->Fill(fabs(eta1), type_lep1, delta_pt1);
								}
								if(pT2 < 100){
									h_deltaPt_vs_d0_inEta_Above20_Below100[1][fs][1]->Fill(fabs(eta1), type_lep1, delta_pt1);
									h_deltaPt_vs_d0_inEta_Above20_Below100[1][fs][0]->Fill(fabs(eta1), type_lep1, delta_pt1);
								}
							}
						}

					}
					////// d0 from PV //////					


				
					
				} //for on entry

/*				
				for(int f = 0; f < 2; f++){
					for(int neg = 1; neg <= h_massZ_NEG_vs_POS[MC][0][fs]->GetNbinsY(); neg++){				
						for(int pos = 1; pos <= h_massZ_NEG_vs_POS[MC][0][fs]->GetNbinsX(); pos++){
							TString NOME;
							NOME = Form("Mass POS %.3f %.3f NEG %.3f %.3f", h_massZ_NEG_vs_POS[MC][0][fs]->GetXaxis()->GetBinLowEdge(pos), h_massZ_NEG_vs_POS[MC][0][fs]->GetXaxis()->GetBinUpEdge(pos), h_massZ_NEG_vs_POS[MC][0][fs]->GetYaxis()->GetBinLowEdge(neg), h_massZ_NEG_vs_POS[MC][0][fs]->GetYaxis()->GetBinUpEdge(neg)); 

							TH1F* proiezione = (TH1F*)h_massZ_NEG_vs_POS[MC][0][fs]->ProjectionZ(NOME,pos,pos,neg,neg);
							NOME = Form("Mass_%s_%s_d0BS", monteCarlo[MC].Data(), fitFunction[f].Data());
							float mean = FitMass(proiezione, directory + NOME, f);
							h_mean_NEG_vs_POS[MC][0][fs][f]->SetBinContent(pos, neg, mean);
							h_integral_NEG_vs_POS[MC][0][fs]->SetBinContent(neg, pos, proiezione->Integral()/counting);
							if(proiezione->Integral()/counting < 0.0001)
								h_mean_NEG_vs_POS[MC][0][fs][f]->SetBinContent(pos, neg, 0);
								
							proiezione = (TH1F*)h_massZ_NEG_vs_POS[MC][1][fs]->ProjectionZ(NOME,pos,pos,neg,neg);
							NOME = Form("Mass_%s_%s_d0PV", monteCarlo[MC].Data(), fitFunction[f].Data());
							mean = FitMass(proiezione, directory + NOME, f);
							h_mean_NEG_vs_POS[MC][0][fs][f]->SetBinContent(pos, neg, mean);
							h_integral_NEG_vs_POS[MC][1][fs]->SetBinContent(neg, pos, proiezione->Integral()/counting);
							if(proiezione->Integral()/counting < 0.0001)
								h_mean_NEG_vs_POS[MC][1][fs][f]->SetBinContent(pos, neg, 0);

	
							if(MC == 0){
								NOME = Form("Mass POS %.3f %.3f NEG %.3f %.3f", h_massZ_NEG_vs_POS[MC][0][fs]->GetXaxis()->GetBinLowEdge(pos), h_massZ_NEG_vs_POS[MC][0][fs]->GetXaxis()->GetBinUpEdge(pos), h_massZ_NEG_vs_POS[MC][0][fs]->GetYaxis()->GetBinLowEdge(neg), h_massZ_NEG_vs_POS[MC][0][fs]->GetYaxis()->GetBinUpEdge(neg)); 

								TH1F* proiezione = (TH1F*)h_massZ_NEG_vs_POS_gen[0][fs]->ProjectionZ(NOME,pos,pos,neg,neg);
								NOME = Form("Mass_GEN_%s_d0BS", fitFunction[f].Data());
								float mean = FitMass(proiezione, directory + NOME, f);
								h_mean_NEG_vs_POS_gen[0][fs][f]->SetBinContent(pos, neg, mean);	
								h_integral_NEG_vs_POS_gen[0][fs]->SetBinContent(neg, pos, proiezione->Integral()/counting);
								if(proiezione->Integral()/counting < 0.0001)
									h_mean_NEG_vs_POS_gen[0][fs][f]->SetBinContent(pos, neg, 0);						
								
								proiezione = (TH1F*)h_massZ_NEG_vs_POS_gen[1][fs]->ProjectionZ(NOME,pos,pos,neg,neg);
								NOME = Form("Mass_GEN_%s_d0PV", fitFunction[f].Data());
								mean = FitMass(proiezione, directory + NOME, f);
								h_mean_NEG_vs_POS_gen[0][fs][f]->SetBinContent(pos, neg, mean);	
								h_integral_NEG_vs_POS_gen[1][fs]->SetBinContent(neg, pos, proiezione->Integral()/counting);
								if(proiezione->Integral()/counting < 0.0001)
									h_mean_NEG_vs_POS_gen[1][fs][f]->SetBinContent(pos, neg, 0);						

							}
						}
						}
					}
*/				
// 				for(int f = 0; f < 2; f++){
// 					for(int t = 0; t < 2; t++){
				for(int f = 0; f < 2; f++){
					for(int t = 0; t < 2; t++){
						TString save_name = "Mass_" + monteCarlo.at(MC) + "_" + d0type.at(t) + "_" + fitFunction.at(f);
						if(t == 0){
							h_mean_NEG_vs_POS[MC][t][0][f] = MakeMean(h_massZ_NEG_vs_POS[MC][t][0], directory + save_name, f, h_integral_NEG_vs_POS[MC][t][0], counting, number_bins.at(fs), extreme_BS.at(fs));
							if(MC == 0){
								save_name = "MassGEN_" + d0type.at(t) + "_" + fitFunction.at(f);
								h_mean_NEG_vs_POS_gen[t][0][f] = MakeMean(h_massZ_NEG_vs_POS_gen[t][0], directory + save_name, f, h_integral_NEG_vs_POS_gen[t][0], counting, number_bins.at(fs), extreme_BS.at(fs));
							}
						}
						else{
							h_mean_NEG_vs_POS[MC][t][0][f] = MakeMean(h_massZ_NEG_vs_POS[MC][t][0], directory + save_name, f, h_integral_NEG_vs_POS[MC][t][0], counting, number_bins.at(fs), extreme_PV.at(fs));
							if(MC == 0){
								save_name = "MassGEN_" + d0type.at(t) + "_" + fitFunction.at(f);
								h_mean_NEG_vs_POS_gen[t][0][f] = MakeMean(h_massZ_NEG_vs_POS_gen[t][0], directory + save_name, f, h_integral_NEG_vs_POS_gen[t][0], counting, number_bins.at(fs), extreme_PV.at(fs));
								if(t == 1 && f == 1){
	TCanvas*c2 = new TCanvas("c2","c2", 700, 500);
	h_mean_NEG_vs_POS_gen[t][0][f]->Draw("COLZ");
	h_mean_NEG_vs_POS_gen[t][0][f]->Draw("same TEXT0");
	c2->Print("ciao2.pdf");

								}

							}
						}
					}
				}

				infile->Close();
			} // for on fs
		} // for on MC-DATA
		
	std::cout<<"QUI"<<std::endl;
	
	TCanvas*c1 = new TCanvas("c1","c1", 700, 500);
	h_mean_NEG_vs_POS_gen[1][0][0]->Draw("COLZ");
	h_mean_NEG_vs_POS_gen[1][0][0]->Draw("same TEXT0");
	c1->Print("ciao.pdf");
	
	std::cout<<"QUI 2222"<<std::endl;
/*		
		gStyle->SetPaintTextFormat("4.2f");
		for(int MC = 0; MC < 2; MC++){
			if(MC == 1) continue;
			for(int type = 0; type < 2; type++){	
// 				for(int fs = 0; fs < 2; fs++){	
				for(int f = 0; f < 2; f++){	
				
					std::cout<<"SONO DENTRO = "<<MC<<"\t"<<type<<"\t"<<f<<std::endl;

					TCanvas* c1 = new TCanvas("c1","c1", 700, 500);
					h_mean_NEG_vs_POS[MC][type][0][f]->Draw("COLZ");
					h_mean_NEG_vs_POS[MC][type][0][f]->Draw("same TEXT0");
					h_mean_NEG_vs_POS[MC][type][0][f]->SetStats(0);
					h_mean_NEG_vs_POS[MC][type][0][f]->GetZaxis()->SetRangeUser(89, 93);
					TString title = "Mean " + d0type.at(type) + " " + finalStates.at(0) + " " + monteCarlo.at(MC) + " " + fitFunction.at(f);
					h_mean_NEG_vs_POS[MC][type][0][f]->SetTitle(title);
					h_mean_NEG_vs_POS[MC][type][0][f]->GetYaxis()->SetTitleOffset(1.25);
					h_mean_NEG_vs_POS[MC][type][0][f]->GetYaxis()->SetTitle("reco d0 (#mu^{+})");
					h_mean_NEG_vs_POS[MC][type][0][f]->GetXaxis()->SetTitle("reco d0 (#mu^{-})");
					if(MC == 0 && type == 0 && f == 0)// && fs == 0)
						c1->Print(directory + "Mean_" + finalStates.at(0) + ".pdf[");
					c1->Print(directory + "Mean_" + finalStates.at(0) + ".pdf");
// 					if(MC == 1 && type == 1 && f == 1)// && fs == 1)
					if(MC == 0 && type == 1 && f == 1)// && fs == 1)
						c1->Print(directory + "Mean_" + finalStates.at(0) + ".pdf]");
					
					if(MC == 0){
						std::cout<<"GEN"<<std::endl;
						TCanvas* c2 = new TCanvas("c2","c2", 700, 500);
						h_mean_NEG_vs_POS_gen[type][0][f]->Draw("COLZ");
						h_mean_NEG_vs_POS_gen[type][0][f]->Draw("same TEXT0");
						h_mean_NEG_vs_POS_gen[type][0][f]->SetStats(0);
						h_mean_NEG_vs_POS_gen[type][0][f]->GetZaxis()->SetRangeUser(89, 93);
						title = "Mean " + d0type.at(type) + " " + finalStates.at(0) + " GENN " + fitFunction.at(f);
						h_mean_NEG_vs_POS_gen[type][0][f]->SetTitle(title);
						h_mean_NEG_vs_POS_gen[type][0][f]->GetYaxis()->SetTitleOffset(1.25);
						h_mean_NEG_vs_POS_gen[type][0][f]->GetYaxis()->SetTitle("reco d0 (#mu^{+})");
						h_mean_NEG_vs_POS_gen[type][0][f]->GetXaxis()->SetTitle("reco d0 (#mu^{-})");
						c2->Print(directory + "Mean_" + finalStates.at(0) + ".pdf");
					}
					
					std::cout<<"SONO FUORI = "<<MC<<"\t"<<type<<"\t"<<f<<std::endl;
					
// 				}
				}
			}
		}

*/
// 		gStyle->SetPaintTextFormat("0.4f");		
// 		for(int type = 0; type < 2; type++){
// 			for(int f = 0; f < 2; f++){			
// 				TCanvas* mean_MC_Data = new TCanvas("Mean (DATA - MC) / MC: " + d0type.at(type) + " " + fitFunction.at(f), "Mean (DATA - MC) / MC: " + d0type.at(type) + " " + fitFunction.at(f), 700, 500);
// 				h_mean_NEG_vs_POS[1][type][0][f]->Add(h_mean_NEG_vs_POS[0][type][0][f], -1);
// 				h_mean_NEG_vs_POS[1][type][0][f]->Divide(h_mean_NEG_vs_POS[0][type][0][f]);
// 				h_mean_NEG_vs_POS[1][type][0][f]->Draw("COLZ");
// 				h_mean_NEG_vs_POS[1][type][0][f]->Draw("same TEXT35");
// 				h_mean_NEG_vs_POS[1][type][0][f]->SetStats(0);
// 				h_mean_NEG_vs_POS[1][type][0][f]->GetZaxis()->SetRangeUser(-0.25, 0.25);
// 				h_mean_NEG_vs_POS[1][type][0][f]->SetTitle("Mean (DATA - MC) / MC: " + d0type.at(type) + " " + fitFunction.at(f));
// 				if(type == 0 && f == 0) mean_MC_Data->Print(directory + "Data_minusMC_overMC.pdf[");
// 				mean_MC_Data->Print(directory + "Data_minusMC_overMC.pdf");
// 				if(type == 1 && f == 1) mean_MC_Data->Print(directory + "Data_minusMC_overMC.pdf]");
// 			}
// 		}

// 		for(int type = 0; type < 2; type++){			
// 			for(int f = 0; f < 2; f++){	
// 				std::cout<<type<<"\t"<<f<<std::endl;		
// 				TCanvas* mean_MC_gen = new TCanvas("Mean (GEN - MC) / MC: " + d0type.at(type) + " " + fitFunction.at(f), "Mean (GEN - MC) / MC: " + d0type.at(type) + " " + fitFunction.at(f), 700, 500);
// 				h_mean_NEG_vs_POS_gen[type][0][f]->Add(h_mean_NEG_vs_POS[0][type][0][f], -1);
// 				h_mean_NEG_vs_POS_gen[type][0][f]->Divide(h_mean_NEG_vs_POS[0][type][0][f]);
// 				h_mean_NEG_vs_POS_gen[type][0][f]->Draw("COLZ");
// 				h_mean_NEG_vs_POS_gen[type][0][f]->Draw("same TEXT35");
// 				h_mean_NEG_vs_POS_gen[type][0][f]->SetStats(0);
// 				h_mean_NEG_vs_POS_gen[type][0][f]->GetZaxis()->SetRangeUser(-0.25, 0.25);
// 				h_mean_NEG_vs_POS_gen[type][0][f]->SetTitle("Mean (GEN - MC) / MC: " + d0type.at(type) + " " + fitFunction.at(f));
// 				if(type == 0 && f == 0) mean_MC_gen->Print(directory + "GEN_minusMC_overMC.pdf[");
// 				mean_MC_gen->Print(directory + "GEN_minusMC_overMC.pdf");
// 				if(type == 1 && f == 1) mean_MC_gen->Print(directory + "GEN_minusMC_overMC.pdf]");
// 			}
// 		}
										
		for(int MC = 0; MC < 2; MC++){
			if(MC == 1) continue;
			for(int type = 0; type < 2; type++){	
// 				for(int fs = 0; fs < 2; fs++){	
				TCanvas* c1 = new TCanvas("c1","c1", 700, 500);
				h_integral_NEG_vs_POS[MC][type][0]->Draw("COLZ");
				h_integral_NEG_vs_POS[MC][type][0]->Draw("same TEXT35");
				h_integral_NEG_vs_POS[MC][type][0]->SetStats(0);
				TString title = "Integral_" + d0type.at(type) + "_" + finalStates.at(0) + "_" + monteCarlo.at(MC);
				h_integral_NEG_vs_POS[MC][type][0]->SetTitle(title);
				h_integral_NEG_vs_POS[MC][type][0]->GetYaxis()->SetTitleOffset(1.25);
				h_integral_NEG_vs_POS[MC][type][0]->GetYaxis()->SetTitle("reco d0 (#mu^{+})");
				h_integral_NEG_vs_POS[MC][type][0]->GetXaxis()->SetTitle("reco d0 (#mu^{-})");

				if(MC == 0 && type == 0)// && fs == 0)
					c1->Print(directory + "Integral_" + finalStates.at(0)+ ".pdf[");
				c1->Print(directory + "Integral_" + finalStates.at(0)+ ".pdf");
				if(MC == 1 && type == 1)// && fs == 1)
					c1->Print(directory + "Integral_" + finalStates.at(0)+ ".pdf]");
// 				}
			}
		}

		for(int fs = 0; fs < 1; fs++){
			if(y == 0 && fs == 0){
				legend_recoGen->AddEntry(h_pt_reco[0][fs][0], "Reco", "f");
				legend_recoGen->AddEntry(h_pt_gen[fs][0], "GEN", "lep");

				legend_DataMC->AddEntry(h_pt_reco[0][fs][0], "MC", "f");
				legend_DataMC->AddEntry(h_pt_reco[1][fs][0], "Data", "lep");
			}
			for(int charge = 0; charge < 3; charge++){				
								
				TString NOME = 	"LeptonPt_" + finalStates.at(fs);	
				if(charge == 1) NOME = "LeptonPt_" + finalStates.at(fs) + "_POS";	
				if(charge == 2) NOME = "LeptonPt_" + finalStates.at(fs) + "_NEG";	
				if(fs == 0 && charge == 0) Draw(h_pt_gen[fs][charge], h_pt_reco[0][fs][charge], NOME,  directory + "LeptonPt.pdf[", legend_recoGen, "p_{T}", 0);
				Draw(h_pt_gen[fs][charge], h_pt_reco[0][fs][charge], NOME, directory + "LeptonPt.pdf", legend_recoGen, "p_{T}", 0);
				if(fs == 0 && charge == 2) Draw(h_pt_gen[fs][charge], h_pt_reco[0][fs][charge], NOME,  directory + "LeptonPt.pdf]", legend_recoGen, "p_{T}", 0);

				NOME = "LeptonPt_" + finalStates.at(fs) + "DataMC";
				if(charge == 1) NOME = "LeptonPt_" + finalStates.at(fs) + "DataMC" + "_POS";	
				if(charge == 2) NOME = "LeptonPt_" + finalStates.at(fs) + "DataMC" + "_NEG";	
				if(fs == 0 && charge == 0) Draw(h_pt_reco[1][fs][charge], h_pt_reco[0][fs][charge], NOME, directory + "LeptonPt_DataMC.pdf[", legend_DataMC, "p_{T}", 1, 10); 
				Draw(h_pt_reco[1][fs][charge], h_pt_reco[0][fs][charge], NOME, directory + "LeptonPt_DataMC.pdf", legend_DataMC, "p_{T}", 1, 10);
				if(fs == 0 && charge == 2) Draw(h_pt_reco[1][fs][charge], h_pt_reco[0][fs][charge], NOME, directory + "LeptonPt_DataMC.pdf]", legend_DataMC, "p_{T}", 1, 10); 

				NOME = 	"LeptonEta_" + finalStates.at(fs);	
				if(charge == 1) NOME = "LeptonEta_" + finalStates.at(fs) + "_POS";	
				if(charge == 2) NOME = "LeptonEta_" + finalStates.at(fs) + "_NEG";	
				if(fs == 0 && charge == 0) Draw(h_eta_gen[fs][charge], h_eta_reco[0][fs][charge], NOME,  directory + "LeptonEta.pdf[", legend_recoGen, "#eta", 0);
				Draw(h_eta_gen[fs][charge], h_eta_reco[0][fs][charge], NOME, directory + "LeptonEta.pdf", legend_recoGen, "#eta", 0);
				if(fs == 0 && charge == 2) Draw(h_eta_gen[fs][charge], h_eta_reco[0][fs][charge], NOME,  directory + "LeptonEta.pdf]", legend_recoGen, "#eta", 0);

				NOME = "LeptonEta_" + finalStates.at(fs) + "DataMC";
				if(charge == 1) NOME = "LeptonEta_" + finalStates.at(fs) + "DataMC" + "_POS";	
				if(charge == 2) NOME = "LeptonEta_" + finalStates.at(fs) + "DataMC" + "_NEG";	
				if(fs == 0 && charge == 0) Draw(h_eta_reco[1][fs][charge], h_eta_reco[0][fs][charge], NOME, directory + "LeptonEta_DataMC.pdf[", legend_DataMC, "#eta", 1, 10); 
				Draw(h_eta_reco[1][fs][charge], h_eta_reco[0][fs][charge], NOME, directory + "LeptonEta_DataMC.pdf", legend_DataMC, "p_{T}", 1, 10);
				if(fs == 0 && charge == 2) Draw(h_eta_reco[1][fs][charge], h_eta_reco[0][fs][charge], NOME, directory + "LeptonEta_DataMC.pdf]", legend_DataMC, "#eta", 1, 10); 


			}
	
			Draw(h_mass[1][fs], h_mass[0][fs], "Mass", directory + "Mass_" + finalStates.at(fs) + "_DataMC.pdf", legend_DataMC, "m_{ll}", 1, 10);

			for(int type = 0; type < 2; type++){	
				if(fs == 0 && type == 0) Draw(h_d0[1][type][fs][0], h_d0[0][type][fs][0], d0type.at(type) + finalStates.at(fs), directory + "d0DataMC.pdf[", legend_DataMC, "reco d0 * charge", 1, 1);
				Draw(h_d0[1][type][fs][0], h_d0[0][type][fs][0], d0type.at(type) + finalStates.at(fs), directory + "d0DataMC.pdf", legend_DataMC, "reco d0 * charge", 1, 1);
				Draw(h_d0[1][type][fs][1], h_d0[0][type][fs][1], d0type.at(type) + " " + finalStates.at(fs) + " POS", directory + "d0DataMC.pdf", legend_DataMC, "reco d0 * charge", 1, 1);
				Draw(h_d0[1][type][fs][2], h_d0[0][type][fs][2], d0type.at(type) + " " + finalStates.at(fs) + " NEG", directory + "d0DataMC.pdf", legend_DataMC, "reco d0 * charge", 1, 1);
				if(fs == 0 && type == 1) Draw(h_d0[1][type][fs][2], h_d0[0][type][fs][2], d0type.at(type) + " " + finalStates.at(fs) + " NEG", directory + "d0DataMC.pdf]", legend_DataMC, "reco d0 * charge", 1, 1);				
			}	
			
			for(int eta = 1; eta <= h_pt_reco_etaBins[0][0]->GetNbinsY(); eta++){
				TString NOME;
				NOME = Form("LepPt RECO: %s < |#eta| < %s", eta_bins_name[eta-1].Data(), eta_bins_name[eta].Data());
				TCanvas* canvas = new TCanvas(NOME, NOME, 200,10,700,500);
				TH1F* proiezione_MC = (TH1F*) h_pt_reco_etaBins[0][fs]->ProjectionX(NOME,eta,eta);
				NOME = Form("LepPt DATA: %s < |#eta| < %s", eta_bins_name[eta-1].Data(), eta_bins_name[eta].Data());
				TH1F* proiezione_DATA = (TH1F*) h_pt_reco_etaBins[1][fs]->ProjectionX(NOME,eta,eta);
				NOME = Form("LepPt GEN: %s < |#eta| < %s", eta_bins_name[eta-1].Data(), eta_bins_name[eta].Data());
				TH1F* proiezione_GEN = (TH1F*) h_pt_gen_etaBins[fs]->ProjectionX(NOME,eta,eta);
				NOME = Form("Lepton p_{T}: %s < |#eta| < %s", eta_bins_name[eta-1].Data(), eta_bins_name[eta].Data());
				if(fs == 0 && eta == 1){
					Draw(proiezione_GEN, proiezione_MC, NOME,  directory + "LeptonPt_etaBins.pdf[", legend_recoGen, "p_{T}", 0);
					Draw(proiezione_DATA, proiezione_MC, NOME,  directory + "LeptonPt_etaBins_DataMC.pdf[", legend_DataMC, "p_{T}", 1, 1);
				}
				Draw(proiezione_GEN, proiezione_MC, NOME,  directory + "LeptonPt_etaBins.pdf", legend_recoGen, "p_{T}", 0);
				Draw(proiezione_DATA, proiezione_MC, NOME,  directory + "LeptonPt_etaBins_DataMC.pdf", legend_DataMC, "p_{T}", 1, 1);
				if(fs == 0 && eta == h_pt_reco_etaBins[0][0]->GetNbinsY()){
					Draw(proiezione_GEN, proiezione_MC, NOME,  directory + "LeptonPt_etaBins.pdf]", legend_recoGen, "p_{T}", 0);
					Draw(proiezione_DATA, proiezione_MC, NOME,  directory + "LeptonPt_etaBins_DataMC.pdf]", legend_DataMC, "p_{T}", 1, 1);		
				}
			}
		}
		
		for(int type = 0; type < 2; type++){		
			for(int fs = 0; fs < 2; fs++){
// 		 		TString save_name = directory + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs);
				TString NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs);
				TCanvas* c1 = new TCanvas(NOME, NOME, 700, 500);
				h_deltaPt_vs_d0[type][fs]->Draw();
				h_deltaPt_vs_d0[type][fs]->GetYaxis()->SetTitle(deltapt_name);
				h_deltaPt_vs_d0[type][fs]->GetXaxis()->SetTitle("reco d0 * charge");
				
				if(type == 0 && fs == 0)
					c1->Print(directory + "DeltaPt_vs_d0.pdf[");
				c1->Print(directory + "DeltaPt_vs_d0.pdf");
				if(type == 1 && fs == 1)
					c1->Print(directory + "DeltaPt_vs_d0.pdf]");				
			}
		}

		for(int type = 0; type < 2; type++){		
			for(int fs = 0; fs < 2; fs++){
		 		TString save_name = directory + "DeltaPt_" + finalStates.at(fs) + "_" + d0type.at(type);
	 			Proiezione_deltaPt(maximum, h_deltaPt_vs_d0[type][fs], save_name, type, fs);
				TString NOME = "#Delta p_{T}/p_{T}^{2} in " +	finalStates.at(fs) + " " + d0type.at(type);
				TCanvas* c1 = new TCanvas(NOME, NOME, 700, 500);
				DeltaPt_vs_d0[type][fs]->Draw("AP*");
				DeltaPt_vs_d0[type][fs]->SetTitle(NOME);
				DeltaPt_vs_d0[type][fs]->GetYaxis()->SetTitle(deltapt_name);
				DeltaPt_vs_d0[type][fs]->GetXaxis()->SetTitle("reco d0 * charge");

				if(type == 0 && fs == 0)
					c1->Print(directory + "DeltaPt_Graph.pdf[");
				c1->Print(directory + "DeltaPt_Graph.pdf");
				if(type == 1 && fs == 1)
					c1->Print(directory + "DeltaPt_Graph.pdf]");
			}
		}

		TString dir_mean = "Mean/";
		gSystem->Exec("mkdir " + directory + dir_mean);
		for(int j = 10; j < 3; j ++){
			if(maximum == 1 && j < 10) continue;
			std::cout<<"SALTA = "<<maximum<<"\t"<<j<<std::endl;
			for(int type = 0; type < 2; type++){		
				for(int fs = 0; fs < 2; fs++){
					TString NOME;
					if(j == 0) NOME = directory + dir_mean + "Mass_vs_" + d0type.at(type) + "_" + finalStates.at(fs);
					else if(j == 1) NOME = directory + dir_mean + "Mass_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_POS";
					else NOME = directory + dir_mean + "Mass_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_NEG";				
// 					std::cout<<NOME<<std::endl;
					Proiezione(h_mass_vs_d0[0][type][fs][j], NOME, 0, type, fs, j);

					if(j == 0) NOME = directory + dir_mean + "Mass_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_Data";
					else if(j == 1) NOME = directory + dir_mean + "Mass_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_POS" + "_Data";
					else NOME = directory + dir_mean + "Mass_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_NEG" + "_Data";				
// 					std::cout<<NOME<<std::endl;
					Proiezione(h_mass_vs_d0[1][type][fs][j], NOME, 1, type, fs, j);

					if(j == 0) NOME = "Mean vs " + d0type.at(type) + " in " + finalStates.at(fs);
					else if(j == 1) NOME = "Mean vs " + d0type.at(type) + " in " + finalStates.at(fs) + "POS";
					else NOME = "Mean vs " + d0type.at(type) + " in " + finalStates.at(fs) + "NEG";
					if(y == 0 && type == 0 && fs == 0 and j == 0){
						legend->AddEntry(Mean_vs_d0[0][type][fs][j][0], "MC: Mean CB fit");
					legend->AddEntry(Mean_vs_d0[0][type][fs][j][1], "MC: Mean Voig fit");
	
						legend->AddEntry(Mean_vs_d0[1][type][fs][j][0], "Data: Mean CB fit");
						legend->AddEntry(Mean_vs_d0[1][type][fs][j][1], "Data: Mean Voig fit");
					}
					TCanvas* c1 = new TCanvas("ciao", "ciao", 700, 500);				
					Mean_vs_d0[0][type][fs][j][0]->Draw("AP*");
					Mean_vs_d0[0][type][fs][j][0]->SetMarkerColor(kRed+2);
					Mean_vs_d0[0][type][fs][j][0]->SetLineColor(kRed+2);
					Mean_vs_d0[0][type][fs][j][0]->SetMaximum(95);
					Mean_vs_d0[0][type][fs][j][0]->SetMinimum(87);
					Mean_vs_d0[0][type][fs][j][1]->Draw("P*");
					Mean_vs_d0[0][type][fs][j][1]->SetMarkerColor(kBlue);
					Mean_vs_d0[0][type][fs][j][1]->SetLineColor(kBlue);
					Mean_vs_d0[0][type][fs][j][0]->SetTitle(NOME);
					Mean_vs_d0[1][type][fs][j][0]->Draw("P*");
					Mean_vs_d0[1][type][fs][j][0]->SetMarkerColor(kRed+2);
					Mean_vs_d0[1][type][fs][j][0]->SetLineColor(kRed+2);
					Mean_vs_d0[1][type][fs][j][0]->SetMarkerStyle(21);
					Mean_vs_d0[1][type][fs][j][1]->Draw("P*");
					Mean_vs_d0[1][type][fs][j][1]->SetMarkerColor(kBlue);
					Mean_vs_d0[1][type][fs][j][1]->SetLineColor(kBlue);
					Mean_vs_d0[1][type][fs][j][1]->SetMarkerStyle(21);
					Mean_vs_d0[0][type][fs][j][0]->GetYaxis()->SetTitle("Mean");
					Mean_vs_d0[0][type][fs][j][0]->GetXaxis()->SetTitle("reco d0 * charge");

					legend->Draw();
				 	if(type == 0 && fs == 0 && j == 0) 
				 		c1->Print(directory + dir_mean + "Graph.pdf[");
					c1->Print(directory + dir_mean + "Graph.pdf");
				 	if(type == 1 && fs == 1 && j == 2)
				 		c1->Print(directory + dir_mean + "Graph.pdf]");
				}
			}
		} 		
		
		TString dir_Pt = "Pt/";
		gSystem->Exec("mkdir " + directory + dir_Pt);
		for(int type = 0; type < 2; type++){		
			for(int fs = 0; fs < 2; fs++){
				for(int charge = 0; charge < 3; charge++){
			 		TString save_name;
		 			save_name = directory + dir_Pt + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs);
		 			if(charge == 1)
			 			save_name = directory + dir_Pt + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_NEG";
		 			if(charge == 2)
			 			save_name = directory + dir_Pt + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_POS";
			 		Proiezione_deltaPt(maximum, h_deltaPt_vs_d0_inPt[type][fs][charge], save_name, type, fs, charge, 1);
// 					if(y == 0 && type == 0 && charge == 0){
// 						TString NOME_legend = Form("%s < p_{T} < %s", pT_bins_name[pt].Data(), pT_bins_name[pt+1].Data());
// 						legend_pt->AddEntry(DeltaPt_vs_d0_inPt[type][fs][pt], NOME_legend);
// 					}		
					TString NOME;
					NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs);
					if(charge == 1)
						NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " + finalStates.at(fs) + " with #mu^{-}";
					if(charge == 2)
						NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " + finalStates.at(fs) + " with #mu^{+}";
					TCanvas* c1 = new TCanvas(NOME, NOME, 700, 500);
					Color_t color;
					
					for(int pT = 0; pT < pT_bins.size()-1; pT++){
						
						if(pT == 1) color = kRed+2;
						else if(pT == 2) color = kGreen+2;
						else if(pT == 3) color = kMagenta+2;
						else if(pT == 4) color = kYellow+2;
						else if(pT == 5) color = kBlue+2;
						else color = kCyan+2;
						if(pT == 0){
							DeltaPt_vs_d0_inPt[type][fs][charge].at(pT)->Draw("AP*");
							if(maximum == 1){
								DeltaPt_vs_d0_inPt[type][fs][charge].at(pT)->SetMaximum(20);
								DeltaPt_vs_d0_inPt[type][fs][charge].at(pT)->SetMinimum(-20);
							}
							else{
								DeltaPt_vs_d0_inPt[type][fs][charge].at(pT)->SetMaximum(50);
								DeltaPt_vs_d0_inPt[type][fs][charge].at(pT)->SetMinimum(-100);
							}
						}
						else{
							DeltaPt_vs_d0_inPt[type][fs][charge].at(pT)->Draw("P*");
							DeltaPt_vs_d0_inPt[type][fs][charge].at(pT)->SetMarkerColor(color);
						}	
						DeltaPt_vs_d0_inPt[type][fs][charge].at(pT)->SetTitle(NOME);
						DeltaPt_vs_d0_inPt[type][fs][charge].at(pT)->GetYaxis()->SetTitle(deltapt_name);
						DeltaPt_vs_d0_inPt[type][fs][charge].at(pT)->GetXaxis()->SetTitle("reco d0 * charge");
					}

// 				legend_pt->Draw();
					if(type == 0 && fs == 0 && charge == 0)
						c1->Print(directory + dir_Pt + "DeltaPt_Graph_inPt.pdf[");
					c1->Print(directory + dir_Pt + "DeltaPt_Graph_inPt.pdf");
					if(type == 1 && fs == 1 && charge == 2)
						c1->Print(directory + dir_Pt + "DeltaPt_Graph_inPt.pdf]");
				}
			}
		}

		TString dir_Eta = "Eta/";
		gSystem->Exec("mkdir " + directory + dir_Eta);	
		for(int type = 0; type < 2; type++){		
			for(int fs = 0; fs < 2; fs++){
				for(int charge = 0; charge < 3; charge++){
			 		TString save_name;
			 		save_name = directory + dir_Eta + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs);
			 		if(charge == 1) save_name = directory + dir_Eta + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_NEG";
			 		if(charge == 2) save_name = directory + dir_Eta + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_POS";
			 		Proiezione_deltaPt(maximum, h_deltaPt_vs_d0_inEta[type][fs][charge], save_name, type, fs, charge, 0);
// 					if(y == 0 && type == 0 && fs == 0){
// 						TString NOME_legend = Form("%s < |#eta| < %s", eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
// 						legend_eta->AddEntry(DeltaPt_vs_d0_inEta[type][fs][eta], NOME_legend);
// 					}

					TString NOME;
					NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs);
					if(charge == 1)
						NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs) + " with #mu^{-}";
					if(charge == 2)
						NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs) + " with #mu^{+}";

					TCanvas* c1 = new TCanvas(NOME, NOME, 700, 500);
					Color_t color;
					for(int eta = 0; eta < eta_bins.size()-1; eta++){

						if(eta == 1) color = kRed+2;
						else if(eta == 2) color = kGreen+2;
						else if(eta == 3) color = kMagenta+2;
						else if(eta == 4) color = kYellow+2;
						else if(eta == 5) color = kBlue+2;
						else color = kCyan+2;
						if(eta == 0){
							DeltaPt_vs_d0_inEta[type][fs][charge].at(eta)->Draw("AP*");
							if(maximum == 1){
								DeltaPt_vs_d0_inEta[type][fs][charge].at(eta)->SetMaximum(20);
								DeltaPt_vs_d0_inEta[type][fs][charge].at(eta)->SetMinimum(-20);
							}
							else{
								DeltaPt_vs_d0_inEta[type][fs][charge].at(eta)->SetMaximum(30);
								DeltaPt_vs_d0_inEta[type][fs][charge].at(eta)->SetMinimum(-60);
							}
						}
						else{
							DeltaPt_vs_d0_inEta[type][fs][charge].at(eta)->Draw("P*");
							DeltaPt_vs_d0_inEta[type][fs][charge].at(eta)->SetMarkerColor(color);
						}	
						DeltaPt_vs_d0_inEta[type][fs][charge].at(eta)->SetTitle(NOME);
						DeltaPt_vs_d0_inEta[type][fs][charge].at(eta)->GetYaxis()->SetTitle(deltapt_name);
						DeltaPt_vs_d0_inEta[type][fs][charge].at(eta)->GetXaxis()->SetTitle("reco d0 * charge");
					}					

					legend_eta->Draw();
					if(type == 0 && fs == 0 && charge == 0)
						c1->Print(directory + dir_Eta + "DeltaPt_Graph_inEta.pdf[");
					c1->Print(directory + dir_Eta + "DeltaPt_Graph_inEta.pdf");
					if(type == 1 && fs == 1 && charge == 2)
						c1->Print(directory + dir_Eta + "DeltaPt_Graph_inEta.pdf]");
				}
			}
		}

		dir_Eta = "Eta_Above20GeV/";
		gSystem->Exec("mkdir " + directory + dir_Eta);	
		for(int type = 0; type < 2; type++){		
			for(int fs = 0; fs < 2; fs++){
				for(int charge = 0; charge < 3; charge++){
			 		TString save_name;
			 		save_name = directory + dir_Eta + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs);
			 		if(charge == 1) save_name = directory + dir_Eta + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_NEG";
			 		if(charge == 2) save_name = directory + dir_Eta + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_POS";
			 		Proiezione_deltaPt(maximum, h_deltaPt_vs_d0_inEta_Above20[type][fs][charge], save_name, type, fs, charge, -1);
// 					if(y == 0 && type == 0 && fs == 0){
// 						TString NOME_legend = Form("%s < |#eta| < %s", eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
// 						legend_eta->AddEntry(DeltaPt_vs_d0_inEta_Above20[type][fs][eta], NOME_legend);
// 					}

					TString NOME;
					NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs);
					if(charge == 1)
						NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs) + " with #mu^{-}";
					if(charge == 2)
						NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs) + " with #mu^{+}";

					std::cout<<NOME<<std::endl;
					TCanvas* c1 = new TCanvas(NOME, NOME, 700, 500);
					Color_t color;
					for(int eta = 0; eta < eta_bins.size()-1; eta++){

						if(eta == 1) color = kRed+2;
						else if(eta == 2) color = kGreen+2;
						else if(eta == 3) color = kMagenta+2;
						else if(eta == 4) color = kYellow+2;
						else if(eta == 5) color = kBlue+2;
						else color = kCyan+2;
						if(eta == 0){
							DeltaPt_vs_d0_inEta_Above20[type][fs][charge].at(eta)->Draw("AP*");
							if(maximum == 1){
								DeltaPt_vs_d0_inEta_Above20[type][fs][charge].at(eta)->SetMaximum(20);
								DeltaPt_vs_d0_inEta_Above20[type][fs][charge].at(eta)->SetMinimum(-20);
							}
							else{
								DeltaPt_vs_d0_inEta_Above20[type][fs][charge].at(eta)->SetMaximum(30);
								DeltaPt_vs_d0_inEta_Above20[type][fs][charge].at(eta)->SetMinimum(-60);
							}
						}
						else{
							DeltaPt_vs_d0_inEta_Above20[type][fs][charge].at(eta)->Draw("P*");
							DeltaPt_vs_d0_inEta_Above20[type][fs][charge].at(eta)->SetMarkerColor(color);
						}	
						DeltaPt_vs_d0_inEta_Above20[type][fs][charge].at(eta)->SetTitle(NOME);
						DeltaPt_vs_d0_inEta_Above20[type][fs][charge].at(eta)->GetYaxis()->SetTitle(deltapt_name);
						DeltaPt_vs_d0_inEta_Above20[type][fs][charge].at(eta)->GetXaxis()->SetTitle("reco d0 * charge");
					}					

// 					legend_eta->Draw();
					if(type == 0 && fs == 0 && charge == 0)
						c1->Print(directory + dir_Eta + "DeltaPt_Graph_inEta_Above20.pdf[");
					c1->Print(directory + dir_Eta + "DeltaPt_Graph_inEta_Above20.pdf");
					if(type == 1 && fs == 1 && charge == 2)
						c1->Print(directory + dir_Eta + "DeltaPt_Graph_inEta_Above20.pdf]");
				}
			}
		}

		dir_Eta = "Eta_Above20_Below100GeV/";
		gSystem->Exec("mkdir " + directory + dir_Eta);	
		for(int type = 0; type < 2; type++){		
			for(int fs = 0; fs < 2; fs++){
				for(int charge = 0; charge < 3; charge++){
			 		TString save_name;
			 		save_name = directory + dir_Eta + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs);
			 		if(charge == 1) save_name = directory + dir_Eta + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_NEG";
			 		if(charge == 2) save_name = directory + dir_Eta + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_POS";
			 		Proiezione_deltaPt(maximum, h_deltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge], save_name, type, fs, charge, -2);
// 					if(y == 0 && type == 0 && fs == 0){
// 						TString NOME_legend = Form("%s < |#eta| < %s", eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
// 						legend_eta->AddEntry(DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][eta], NOME_legend);
// 					}
					TString NOME;
					NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs);
					if(charge == 1)
						NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs) + " with #mu^{-}";
					if(charge == 2)
						NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs) + " with #mu^{+}";

					TCanvas* c1 = new TCanvas(NOME, NOME, 700, 500);
					Color_t color;
					for(int eta = 0; eta < eta_bins.size()-1; eta++){

						if(eta == 1) color = kRed+2;
						else if(eta == 2) color = kGreen+2;
						else if(eta == 3) color = kMagenta+2;
						else if(eta == 4) color = kYellow+2;
						else if(eta == 5) color = kBlue+2;
						else color = kCyan+2;
						if(eta == 0){
							DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge].at(eta)->Draw("AP*");
							if(maximum == 1){
								DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge].at(eta)->SetMaximum(20);
								DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge].at(eta)->SetMinimum(-20);
							}
							else{
								DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge].at(eta)->SetMaximum(30);
								DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge].at(eta)->SetMinimum(-60);
							}
						}
						else{
							DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge].at(eta)->Draw("P*");
							DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge].at(eta)->SetMarkerColor(color);
						}	
						DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge].at(eta)->SetTitle(NOME);
						DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge].at(eta)->GetYaxis()->SetTitle(deltapt_name);
						DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge].at(eta)->GetXaxis()->SetTitle("reco d0 * charge");
					}					

					legend_eta->Draw();
					if(type == 0 && fs == 0 && charge == 0)
						c1->Print(directory + dir_Eta + "DeltaPt_Graph_inEta_Above20_Below100.pdf[");
					c1->Print(directory + dir_Eta + "DeltaPt_Graph_inEta_Above20_Below100.pdf");
					if(type == 1 && fs == 1 && charge == 2)
						c1->Print(directory + dir_Eta + "DeltaPt_Graph_inEta_Above20_Below100.pdf]");
				}
			}
		}

		dir_Eta = "Eta_Above20_Below200GeV/";
		gSystem->Exec("mkdir " + directory + dir_Eta);	
		for(int type = 0; type < 2; type++){		
			for(int fs = 0; fs < 2; fs++){
				for(int charge = 0; charge < 3; charge++){
			 		TString save_name;
			 		save_name = directory + dir_Eta + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs);
			 		if(charge == 1) save_name = directory + dir_Eta + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_NEG";
			 		if(charge == 2) save_name = directory + dir_Eta + "DeltaPt_vs_" + d0type.at(type) + "_" + finalStates.at(fs) + "_POS";
			 		Proiezione_deltaPt(maximum, h_deltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge], save_name, type, fs, charge, -3);
// 					if(y == 0 && type == 0 && fs == 0){
// 						TString NOME_legend = Form("%s < |#eta| < %s", eta_bins_name[eta].Data(), eta_bins_name[eta+1].Data());
// 						legend_eta->AddEntry(DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][eta], NOME_legend);
// 					}
					TString NOME;
					NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs);
					if(charge == 1)
						NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs) + " with #mu^{-}";
					if(charge == 2)
						NOME = "#Delta p_{T}/p_{T}^{2} vs " + d0type.at(type) + " in " +	finalStates.at(fs) + " with #mu^{+}";

					TCanvas* c1 = new TCanvas(NOME, NOME, 700, 500);
					Color_t color;
					for(int eta = 0; eta < eta_bins.size()-1; eta++){

						if(eta == 1) color = kRed+2;
						else if(eta == 2) color = kGreen+2;
						else if(eta == 3) color = kMagenta+2;
						else if(eta == 4) color = kYellow+2;
						else if(eta == 5) color = kBlue+2;
						else color = kCyan+2;
						if(eta == 0){
							DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge].at(eta)->Draw("AP*");
							if(maximum == 1){
								DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge].at(eta)->SetMaximum(20);
								DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge].at(eta)->SetMinimum(-20);
							}
							else{
								DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge].at(eta)->SetMaximum(30);
								DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge].at(eta)->SetMinimum(-60);
							}
						}
						else{
							DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge].at(eta)->Draw("P*");
							DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge].at(eta)->SetMarkerColor(color);
						}	
						DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge].at(eta)->SetTitle(NOME);
						DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge].at(eta)->GetYaxis()->SetTitle(deltapt_name);
						DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge].at(eta)->GetXaxis()->SetTitle("reco d0 * charge");
					}					

					legend_eta->Draw();
					if(type == 0 && fs == 0 && charge == 0)
						c1->Print(directory + dir_Eta + "DeltaPt_Graph_inEta_Above20_Below200.pdf[");
					c1->Print(directory + dir_Eta + "DeltaPt_Graph_inEta_Above20_Below200.pdf");
					if(type == 1 && fs == 1 && charge == 2)
						c1->Print(directory + dir_Eta + "DeltaPt_Graph_inEta_Above20_Below200.pdf]");
				}
			}
		}

		counting = 0;
		
		std::cout<<"Passo allo anno successivo"<<std::endl;
		
	} //for on years

}


















void Proiezione(TH2F* h2, TString save_name, int MC, int type, int fs, int j){
	
	std::cout<<"ABCDEF = "<<save_name<<std::endl;
	TString deltapt_name = "10000 * (p_{T}^{reco} - p_{T}^{GEN}) / (p_{T}^{GEN})^{2}";	

	Int_t binnum_d0 = h2->GetNbinsX();
	float bin_min = h2->GetXaxis()->GetXmin();
	float bin_max = h2->GetXaxis()->GetXmax();
	std::vector<float> d0_BINS;
	std::vector<float> d0_BINSerr;
	std::vector<float> d0_extreme;
	for(int i = 0; i < binnum_d0; i++){
		d0_BINS.push_back(bin_min + (bin_max - bin_min) / binnum_d0 * (i + 0.5));
		d0_BINSerr.push_back((bin_max - bin_min) / (2 * binnum_d0 ) );
		
		d0_extreme.push_back(bin_min + i * ((bin_max - bin_min) / (binnum_d0 )));
	}
	d0_extreme.push_back(bin_max);

	for(int d0 = 1; d0 <= h2->GetNbinsX(); d0++){				
		TString NOME;
		NOME = Form("Mass: %.4f < d0 < %.4f", d0_extreme.at(d0-1), d0_extreme.at(d0));
		TCanvas* canvas = new TCanvas(NOME, NOME, 200,10,700,500);
		TH1D* proiezione = h2->ProjectionY(NOME,d0,d0);
	 	proiezione->Draw();
		proiezione->GetYaxis()->SetTitle("Entries");
		proiezione->GetXaxis()->SetTitle("m_{ll}");

		RooPlot* xframe;
		RooRealVar massZ = RooRealVar("massZ","massZ", massZ_min, massZ_max);
		RooDataHist histo("MC", "MC", massZ, proiezione);
		xframe = massZ.frame(Title(NOME));	 	
			
		//// BW			
		RooRealVar BW_mean("BW_mean", "BW_mean", 80, 100);
		RooRealVar BW_sigma("BW_sigma", "BW_sigma", 2.44); 
		RooBreitWigner BW("BW","BW", massZ, BW_mean, BW_sigma);
		//// BW			

		//// CB			
		RooRealVar CB_mean("CB_mean", "CB_mean", 0, -5, 5);
		RooRealVar CB_sigma("CB_sigma", "CB_sigma", 1, 0, 10);
		RooRealVar CB_alpha("CB_alpha", "CB_alpha", 1, 0, 10);
		RooRealVar CB_exp("CB_exp", "CB_exp", 5, 0, 30);
	 	RooCBShape CB("CB", "CB", massZ, CB_mean, CB_sigma, CB_alpha, CB_exp);
		//// CB	
		
		//// VOIGTIAN		
		RooRealVar VOIG_mean("VOIG_mean", "VOIG_mean", 80, 100);
		RooRealVar VOIG_sigma("VOIG_sigma", "VOIG_sigma", 2.44); 
		RooRealVar VOIG_gamma("VOIG_gamma","VOIG_gamma", 0, 10);		
		RooVoigtian VOIG("VOIG", "VOIG", massZ, VOIG_mean, VOIG_sigma, VOIG_gamma);
		//// VOIGTIAN

		//// expo
		RooRealVar tau("tau", "tau", 0, -1, 1);
		RooExponential bkg("bkg","bkg", massZ, tau);
		//// expo	
		
		RooFFTConvPdf BWxCB("BWxCB","BWxCB", massZ, BW, CB);
				
		RooRealVar fsig("fsig","signal fraction", 0.7, 0.5, 1);
		
		RooAddPdf model("model","model", BWxCB, bkg, fsig);
		RooAddPdf model_2("model_2","model_2", VOIG, bkg, fsig);
		
		histo.plotOn(xframe, MarkerSize(0.75));

		std::cout<<"MEANRANGE = "<<proiezione->GetMean()<<"\t"<<proiezione->GetRMS()<<std::endl;
		
		RooFitResult* r_model = model.fitTo(histo, Range(proiezione->GetMean() - 2 * proiezione->GetRMS(), proiezione->GetMean() + 2 * proiezione->GetRMS()), Save());
// 		RooFitResult* r_model = model.fitTo(histo, Range(86, 96), Save());
		model.plotOn(xframe, LineStyle(2), LineColor(kGreen+2));//, LineWidth(2));	
		model.paramOn(xframe, Layout(0.1, 0.4, 0.85));
		xframe->getAttText()->SetTextSize(0.03);
		xframe->getAttText()->SetTextColor(kGreen+2);
		TString t1 = Form("#chi^{2} = %.3f", xframe->chiSquare(r_model->floatParsFinal().getSize())/ r_model->floatParsFinal().getSize());
		TLatex* tex1 = new TLatex(0.2,0.9, t1);
		tex1->SetNDC();
		tex1->SetTextColor(kGreen);
		tex1->Draw();

		canvas->Update();
	 	
	 	mean_vs_d0_tmp.push_back(BW_mean.getVal());
	 	mean_vs_d0err_tmp.push_back(BW_mean.getError());		

		RooFitResult* r_model_2 = model_2.fitTo(histo, Range(proiezione->GetMean() - 2 * proiezione->GetRMS(), proiezione->GetMean() + 2 * proiezione->GetRMS()), Save());
// 		RooFitResult* r_model_2 = model_2.fitTo(histo, Range(86, 96), Save());
		model_2.plotOn(xframe, LineStyle(3), LineColor(kRed+2));//, LineWidth(2));	
		model_2.paramOn(xframe, Layout(0.6, 0.9, 0.85));
		xframe->getAttText()->SetTextSize(0.03);
		xframe->getAttText()->SetTextColor(kRed+2);
		t1 = Form("#chi^{2} = %.3f", xframe->chiSquare(r_model_2->floatParsFinal().getSize())/ r_model_2->floatParsFinal().getSize());
		tex1 = new TLatex(0.6,0.9, t1);
		tex1->SetNDC();
		tex1->SetTextColor(kRed);
		tex1->Draw();

		canvas->Update();

	 	mean_vs_d0_tmp2.push_back(VOIG_mean.getVal());
	 	mean_vs_d0err_tmp2.push_back(VOIG_mean.getError());		

		xframe->Draw();

		std::cout<<"ABCDEF = "<<save_name<<std::endl;
		if(MC == 0){
	 		if(d0 == 1) 
 				canvas->Print(save_name + ".pdf[");
			canvas->Print(save_name + ".pdf");
 			if(d0 == h2->GetNbinsX()) 
		 		canvas->Print(save_name + ".pdf]");	
		}

		if(MC == 1){
	 		if(d0 == 1) 
 				canvas->Print(save_name + "_Data.pdf[");
			canvas->Print(save_name + "_Data.pdf");
 			if(d0 == h2->GetNbinsX()) 
		 		canvas->Print(save_name + "_Data.pdf]");	
		 }
	 		
	}
	
	Mean_vs_d0[MC][type][fs][j][0] = new TGraphErrors(binnum_d0, &(d0_BINS[0]), &(mean_vs_d0_tmp[0]), &(d0_BINSerr[0]), &(mean_vs_d0err_tmp[0]));
	Mean_vs_d0[MC][type][fs][j][1] = new TGraphErrors(binnum_d0, &(d0_BINS[0]), &(mean_vs_d0_tmp2[0]), &(d0_BINSerr[0]), &(mean_vs_d0err_tmp2[0]));

	mean_vs_d0_tmp.clear();
	mean_vs_d0err_tmp.clear();
	mean_vs_d0_tmp2.clear();
	mean_vs_d0err_tmp2.clear();

}

void Proiezione_deltaPt(int maximum, TH2F* h2, TString save_name, int type, int fs){
	
	std::cout<<save_name<<std::endl;
	TString deltapt_name = "10000 * (p_{T}^{reco} - p_{T}^{GEN}) / (p_{T}^{GEN})^{2}";
		
	Int_t binnum_d0 = h2->GetNbinsX();
	float bin_min = h2->GetXaxis()->GetXmin();
	float bin_max = h2->GetXaxis()->GetXmax();
	std::vector<float> d0_BINS;
	std::vector<float> d0_BINSerr;
	std::vector<float> d0_extreme;
	for(int i = 0; i < binnum_d0; i++){
		d0_BINS.push_back(bin_min + (bin_max - bin_min) / binnum_d0 * (i + 0.5));
		d0_BINSerr.push_back((bin_max - bin_min) / (2 * binnum_d0 ) );
		
		d0_extreme.push_back(bin_min + i * ((bin_max - bin_min) / (binnum_d0 )));
	}
	d0_extreme.push_back(bin_max);

	for(int d0 = 1; d0 <= h2->GetNbinsX(); d0++){				
		TString NOME = Form("#Delta p_{T}/p_{T}^{2}: %.4f < d0 < %.4f", d0_extreme.at(d0-1), d0_extreme.at(d0));
		TCanvas* canvas = new TCanvas(NOME, NOME, 200,10,700,500);
		TH1D* proiezione = h2->ProjectionY(NOME,d0,d0);
	 	proiezione->SetStats(1111111);
	 	proiezione->Draw();
	 	proiezione->SetStats(1111111);
		proiezione->GetYaxis()->SetTitle("Entries");
		proiezione->GetXaxis()->SetTitle(deltapt_name);
	 	
// 	 	mean_vs_d0_tmp.push_back(proiezione->GetMean());
// 	 	mean_vs_d0err_tmp.push_back(proiezione->GetMeanError());	

		if(maximum == 1){
			int binmax = proiezione->GetMaximumBin(); 
			double x = proiezione->GetXaxis()->GetBinCenter(binmax);
			mean_vs_d0_tmp.push_back(x);
			mean_vs_d0err_tmp.push_back(sqrt(fabs(x)));
		}
		else{ 
		 	mean_vs_d0_tmp.push_back(proiezione->GetMaximumBin());
		 	mean_vs_d0err_tmp.push_back(sqrt(proiezione->GetMaximumBin()));	
		 }

 		if(d0 == 1) 
 			canvas->Print(save_name + ".pdf[");
		canvas->Print(save_name + ".pdf");
 		if(d0 == h2->GetNbinsX()) 
	 		canvas->Print(save_name + ".pdf]");	
	 		
	}
	
	DeltaPt_vs_d0[type][fs] = new TGraphErrors(binnum_d0, &(d0_BINS[0]), &(mean_vs_d0_tmp[0]), &(d0_BINSerr[0]), &(mean_vs_d0err_tmp[0]));

	mean_vs_d0_tmp.clear();
	mean_vs_d0err_tmp.clear();

}

void Proiezione_deltaPt(int maximum, TH3F* h2, TString save_name, int type, int fs, int charge, int isPT){
	
// 	std::cout<<save_name<<std::endl;
	TString deltapt_name = "10000 * (p_{T}^{reco} - p_{T}^{GEN}) / (p_{T}^{GEN})^{2}";
		
	Int_t binnum_d0 = h2->GetNbinsY();
	float bin_min = h2->GetYaxis()->GetXmin();
	float bin_max = h2->GetYaxis()->GetXmax();
	std::vector<float> d0_BINS;
	std::vector<float> d0_BINSerr;
	std::vector<float> d0_extreme;
	for(int i = 0; i < binnum_d0; i++){
		d0_BINS.push_back(bin_min + (bin_max - bin_min) / binnum_d0 * (i + 0.5));
		d0_BINSerr.push_back((bin_max - bin_min) / (2 * binnum_d0 ) );
		
		d0_extreme.push_back(bin_min + i * ((bin_max - bin_min) / (binnum_d0 )));
	}
	d0_extreme.push_back(bin_max);
	
	for(int pT = 1; pT <= h2->GetNbinsX(); pT++){			
		for(int d0 = 1; d0 <= h2->GetNbinsY(); d0++){
			TString NOME;
			if(isPT == 1) NOME = Form("#Delta p_{T}/p_{T}^{2}: %.4f < d0 < %.4f [%.0f, %.0f]", d0_extreme.at(d0-1), d0_extreme.at(d0), h2->GetXaxis()->GetBinLowEdge(pT), h2->GetXaxis()->GetBinUpEdge(pT));
			else NOME = Form("#Delta p_{T}/p_{T}^{2}: %.4f < d0 < %.4f [%.1f, %.1f]", d0_extreme.at(d0-1), d0_extreme.at(d0), h2->GetXaxis()->GetBinLowEdge(pT), h2->GetXaxis()->GetBinUpEdge(pT));
			std::cout<<NOME<<std::endl;
			TCanvas* canvas = new TCanvas(NOME, NOME, 200,10,700,500);
			TH1D* proiezione = h2->ProjectionZ(NOME, pT, pT, d0,d0);
		 	proiezione->SetStats(1111111);
		 	proiezione->Draw();
// 	 		proiezione->SetStats(1111111);
			proiezione->GetYaxis()->SetTitle("Entries");
			proiezione->GetXaxis()->SetTitle(deltapt_name);
	 	
			if(maximum == 1){
				int binmax = proiezione->GetMaximumBin(); 
				double x = proiezione->GetXaxis()->GetBinCenter(binmax);
			 	mean_vs_d0_tmp.push_back(x);
			 	mean_vs_d0err_tmp.push_back(sqrt(fabs(x)));	
	 		}
			else{
			 	mean_vs_d0_tmp.push_back(proiezione->GetMean());
			 	mean_vs_d0err_tmp.push_back(proiezione->GetMeanError());
			 }

	 		if(pT == 1 && d0 == 1) 
 				canvas->Print(save_name + ".pdf[");
			canvas->Print(save_name + ".pdf");
	 		if(pT == h2->GetNbinsX() && d0 == h2->GetNbinsY()) 
		 		canvas->Print(save_name + ".pdf]");		 		
		}
				
		TGraphErrors* tmp_graph = new TGraphErrors(binnum_d0, &(d0_BINS[0]), &(mean_vs_d0_tmp[0]), &(d0_BINSerr[0]), &(mean_vs_d0err_tmp[0]));
		if(isPT == 1) DeltaPt_vs_d0_inPt[type][fs][charge].push_back(tmp_graph);
		else if(isPT == 0) DeltaPt_vs_d0_inEta[type][fs][charge].push_back(tmp_graph);
		else if(isPT == -1) DeltaPt_vs_d0_inEta_Above20[type][fs][charge].push_back(tmp_graph);
		else if(isPT == -2) DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge].push_back(tmp_graph);
		else DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge].push_back(tmp_graph);
		
// 		if(isPT == 1) DeltaPt_vs_d0_inPt[type][fs][charge] = new TGraphErrors(binnum_d0, &(d0_BINS[0]), &(mean_vs_d0_tmp[0]), &(d0_BINSerr[0]), &(mean_vs_d0err_tmp[0]));
// 		else if(isPT == 0) DeltaPt_vs_d0_inEta[type][fs][charge] = new TGraphErrors(binnum_d0, &(d0_BINS[0]), &(mean_vs_d0_tmp[0]), &(d0_BINSerr[0]), &(mean_vs_d0err_tmp[0]));
// 		else if(isPT == -1) DeltaPt_vs_d0_inEta_Above20_Below100[type][fs][charge] = new TGraphErrors(binnum_d0, &(d0_BINS[0]), &(mean_vs_d0_tmp[0]), &(d0_BINSerr[0]), &(mean_vs_d0err_tmp[0]));
// 		else if(isPT == -2) DeltaPt_vs_d0_inEta_Above20_Below200[type][fs][charge] = new TGraphErrors(binnum_d0, &(d0_BINS[0]), &(mean_vs_d0_tmp[0]), &(d0_BINSerr[0]), &(mean_vs_d0err_tmp[0]));
// 		else DeltaPt_vs_d0_inEta_Above20[type][fs][charge] = new TGraphErrors(binnum_d0, &(d0_BINS[0]), &(mean_vs_d0_tmp[0]), &(d0_BINSerr[0]), &(mean_vs_d0err_tmp[0]));
		mean_vs_d0_tmp.clear();
		mean_vs_d0err_tmp.clear();
	}
	
}

void Draw(TH1F *h1, TH1F *h2, TString nome_canvas, TString save, TLegend *legend, TString x_name, bool mcDATA, float min = 0.1){

	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.23, 1, 1);
//    	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.05, 1, 1);
//    	pad11->SetGrid();
   	pad11->SetLogy();
   	pad11->SetBottomMargin(0.035);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();

	Double_t max;
	
// 	h1->Scale(1/h1->Integral());
// 	h2->Scale(1/h2->Integral());
// 	h3->Scale(1/h3->Integral());
	
	if(h1->GetMaximum() > h2->GetMaximum())
		max = h1->GetMaximum();
	else
		max = h2->GetMaximum();
	
	max = 10 * max;
	
	if(max == 0) max = 1;
		
	h1->SetTitle(nome_canvas);
	h2->SetTitle(nome_canvas);
// 	h1->SetLineColor(1);
	h1->SetMarkerColor(1);
	h1->SetMarkerStyle(20);

	h1->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h2->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	
	canvas->Update();

	h1->SetStats(0);
	h2->SetStats(0);


	h1->GetYaxis()->SetTitleSize(0.125);
	h1->GetYaxis()->SetTitle("Entries");
	h2->GetYaxis()->SetTitle("Entries");
	
	h2->SetLineColor(kBlue);
	h2->SetFillColor(kBlue);
	h2->Draw("HIST");
	h1->Draw("same PE");
	
// 	h2->Draw("Same");
	
	legend->Draw();	
	
// 	TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
// 	pt->SetBorderSize(0);
// 	pt->SetTextAlign(12);
// 	pt->SetFillStyle(0);
// 	pt->SetTextFont(42);
// 	pt->SetTextSize(0.05);
// 	TText *text;
// 	if(muon)
// 		text = pt->AddText(0.01,0.5,"Muon");
// 	else
// 		text = pt->AddText(0.01,0.5,"Electron");
// 	pt->Draw();  
	
	TH1F *ratio_2017 = (TH1F*) h1->Clone();

	legend->Draw();

	canvas->Update();
	canvas->cd();

	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.25);
// 	pad22->SetTopMargin(0);
// 	pad22->SetTopMargin(0.95);
	pad22->SetBottomMargin(0.35);
// 	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	
	ratio_2017->Divide(h2);
	ratio_2017->GetYaxis()->SetRangeUser(0.5, 1.5);
	ratio_2017->SetTitle("");
	ratio_2017->SetStats(0);
// 	ratio->GetYaxis()->SetTitleFont(43);


	ratio_2017->GetYaxis()->SetLabelSize(0.1);
	ratio_2017->GetYaxis()->SetTitleOffset(0.35);
	ratio_2017->GetYaxis()->SetTitleSize(0.11);
	if(!mcDATA) ratio_2017->GetYaxis()->SetTitle("gen / reco");
	else ratio_2017->GetYaxis()->SetTitle("data / MC");
	ratio_2017->GetYaxis()->SetNdivisions(506); 

	ratio_2017->GetXaxis()->SetLabelSize(0.12);
	ratio_2017->GetXaxis()->SetTitleSize(0.12);
	ratio_2017->GetXaxis()->SetTitle(x_name);
// 	ratio_2017->GetYaxis()->SetNdivisions(506); 


	ratio_2017->SetLineColor(kBlack);
	ratio_2017->SetMarkerColor(kBlack);

	ratio_2017->Draw();
	canvas->Update();
	pad22->Update();
	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kRed);
	line->SetLineWidth(1);
	line->Draw();

	canvas->Print(save);
	
}

void DrawMCData(TH1F *hMC, TH1F *hData, TString nome_canvas, TString save, TLegend *legend, TString x_name, float weight, float min = 0.1){

	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
//    	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.05, 1, 1);
//    	pad11->SetGrid();
   	pad11->SetLogy();
//    	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
//    	pad11->SetTicks();

	Double_t max;
	
	hMC->Scale(1/weight);
// 	hData->Scale(1/hData->Integral());
// 	h3->Scale(1/h3->Integral());
	
	if(hMC->GetMaximum() > hData->GetMaximum())
		max = hMC->GetMaximum();
	else
		max = hData->GetMaximum();
	
	max = 10 * max;
	
	if(max == 0) max = 1;
		
	hMC->SetTitle(nome_canvas);
	hData->SetTitle(nome_canvas);

	hMC->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	hData->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	
	canvas->Update();
	
	hData->SetLineColor(1);
	hData->SetMarkerColor(1);
	hData->SetMarkerStyle(20);

	hMC->SetStats(0);
	hData->SetStats(0);
	
	hMC->SetLineColor(kBlue);
	hMC->SetFillColor(kBlue);
	
	hMC->GetYaxis()->SetTitle("entries");
	hData->GetYaxis()->SetTitle("entries");


	hMC->Draw("HIST");
	hData->Draw("Same");
	
	legend->Draw();	
	
// 	TPaveText *pt = new TPaveText(0.1577181,0.9562937,0.9580537,0.9947552,"brNDC");
// 	pt->SetBorderSize(0);
// 	pt->SetTextAlign(12);
// 	pt->SetFillStyle(0);
// 	pt->SetTextFont(42);
// 	pt->SetTextSize(0.05);
// 	TText *text;
// 	if(muon)
// 		text = pt->AddText(0.01,0.5,"Muon");
// 	else
// 		text = pt->AddText(0.01,0.5,"Electron");
// 	pt->Draw();  
	
	TH1F *ratio_2017 = (TH1F*) hData->Clone();

	legend->Draw();

	canvas->Update();
	canvas->cd();

	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.25);
// 	pad22->SetTopMargin(0);
// 	pad22->SetTopMargin(0.95);
// 	pad22->SetBottomMargin(0.05);
// 	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
// 	pad22->SetTicks();
	
	ratio_2017->Divide(hMC);
	ratio_2017->GetYaxis()->SetRangeUser(0.5, 1.5);
	ratio_2017->SetTitle("");
	ratio_2017->SetStats(0);
// 	ratio->GetYaxis()->SetTitleFont(43);
	ratio_2017->GetYaxis()->SetLabelSize(0.14);
	ratio_2017->GetYaxis()->SetTitleOffset(0.50);
	ratio_2017->GetYaxis()->SetTitleSize(0.1);
	ratio_2017->GetYaxis()->SetTitle("data / MC");//Reco / Gen");
	ratio_2017->GetYaxis()->SetNdivisions(506); 

	ratio_2017->GetXaxis()->SetLabelSize(0.1);
	ratio_2017->GetXaxis()->SetTitleSize(0.1);
	ratio_2017->GetXaxis()->SetTitle(x_name);
// 	ratio_2017->GetYaxis()->SetNdivisions(506); 

	ratio_2017->SetLineColor(kRed+2);
	ratio_2017->Draw();
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

	canvas->Print(save);
	
}

float FitMass(TH1D* h1, TString save_name, int fit, int pos, int neg){

	std::cout<<"ABCDEF = "<<save_name<<std::endl;
	TString deltapt_name = "10000 * (p_{T}^{reco} - p_{T}^{GEN}) / (p_{T}^{GEN})^{2}";	

// 	Int_t binnum_d0 = h2->GetNbinsX();
// 	float bin_min = h1->GetXaxis()->GetXmin();
// 	float bin_max = h1->GetXaxis()->GetXmax();
// 	std::vector<float> d0_BINS;
// 	std::vector<float> d0_BINSerr;
// 	std::vector<float> d0_extreme;
// 	for(int i = 0; i < binnum_d0; i++){
// 		d0_BINS.push_back(bin_min + (bin_max - bin_min) / binnum_d0 * (i + 0.5));
// 		d0_BINSerr.push_back((bin_max - bin_min) / (2 * binnum_d0 ) );
// 		
// 		d0_extreme.push_back(bin_min + i * ((bin_max - bin_min) / (binnum_d0 )));
// 	}
// 	d0_extreme.push_back(bin_max);

	TString NOME;
// 	NOME = Form("Mass: %.4f < d0 < %.4f", d0_extreme.at(d0-1), d0_extreme.at(d0));
	NOME = "Mass";
	TCanvas* canvas = new TCanvas(NOME, NOME, 200,10,700,500);
 	h1->Draw();
	h1->GetYaxis()->SetTitle("Entries");
	h1->GetXaxis()->SetTitle("m_{ll}");
	

	RooPlot* xframe;
	RooRealVar massZ = RooRealVar("massZ","massZ", massZ_min, massZ_max);
	RooDataHist histo("MC", "MC", massZ, h1);
	xframe = massZ.frame(Title(NOME));	 	

			
	//// BW			
	RooRealVar BW_mean("BW_mean", "BW_mean", 80, 100);
	RooRealVar BW_sigma("BW_sigma", "BW_sigma", 2.44); 
	RooBreitWigner BW("BW","BW", massZ, BW_mean, BW_sigma);
	//// BW			

	//// CB			
	RooRealVar CB_mean("CB_mean", "CB_mean", 0, -5, 5);
	RooRealVar CB_sigma("CB_sigma", "CB_sigma", 1, 0, 10);
	RooRealVar CB_alpha("CB_alpha", "CB_alpha", 1, 0, 10);
	RooRealVar CB_exp("CB_exp", "CB_exp", 5, 0, 30);
 	RooCBShape CB("CB", "CB", massZ, CB_mean, CB_sigma, CB_alpha, CB_exp);
	//// CB	
		
	//// VOIGTIAN		
	RooRealVar VOIG_mean("VOIG_mean", "VOIG_mean", 80, 100);
	RooRealVar VOIG_sigma("VOIG_sigma", "VOIG_sigma", 2.44); 
	RooRealVar VOIG_gamma("VOIG_gamma","VOIG_gamma", 0, 10);		
	RooVoigtian VOIG("VOIG", "VOIG", massZ, VOIG_mean, VOIG_sigma, VOIG_gamma);
	//// VOIGTIAN

	//// expo
	RooRealVar tau("tau", "tau", 0, -1, 1);
	RooExponential bkg("bkg","bkg", massZ, tau);
	//// expo	
		
	RooFFTConvPdf BWxCB("BWxCB","BWxCB", massZ, BW, CB);
				
	RooRealVar fsig("fsig","signal fraction", 0.7, 0.5, 1);
		
	RooAddPdf model("model","model", BWxCB, bkg, fsig);
	RooAddPdf model_2("model_2","model_2", VOIG, bkg, fsig);
		
	histo.plotOn(xframe, MarkerSize(0.75));
/*
	std::cout<<"MEANRANGE = "<<h1->GetMean()<<"\t"<<h1->GetRMS()<<std::endl;
		
	RooFitResult* r_model = model.fitTo(histo, Range(h1->GetMean() - 2 * h1->GetRMS(), h1->GetMean() + 2 * h1->GetRMS()), Save());
// 	RooFitResult* r_model = model.fitTo(histo, Range(86, 96), Save());
	model.plotOn(xframe, LineStyle(2), LineColor(kGreen+2));//, LineWidth(2));	
	model.paramOn(xframe, Layout(0.1, 0.4, 0.85));
	xframe->getAttText()->SetTextSize(0.03);
	xframe->getAttText()->SetTextColor(kGreen+2);
	TString t1 = Form("#chi^{2} = %.3f", xframe->chiSquare(r_model->floatParsFinal().getSize())/ r_model->floatParsFinal().getSize());
	TLatex* tex1 = new TLatex(0.2,0.9, t1);
	tex1->SetNDC();
	tex1->SetTextColor(kGreen);
	tex1->Draw();

	canvas->Update();

	RooFitResult* r_model_2 = model_2.fitTo(histo, Range(h1->GetMean() - 2 * h1->GetRMS(), h1->GetMean() + 2 * h1->GetRMS()), Save());
// 	RooFitResult* r_model_2 = model_2.fitTo(histo, Range(86, 96), Save());
	model_2.plotOn(xframe, LineStyle(3), LineColor(kRed+2));//, LineWidth(2));	
	model_2.paramOn(xframe, Layout(0.6, 0.9, 0.85));
	xframe->getAttText()->SetTextSize(0.03);
	xframe->getAttText()->SetTextColor(kRed+2);
	t1 = Form("#chi^{2} = %.3f", xframe->chiSquare(r_model_2->floatParsFinal().getSize())/ r_model_2->floatParsFinal().getSize());
	tex1 = new TLatex(0.6,0.9, t1);
	tex1->SetNDC();
	tex1->SetTextColor(kRed);
	tex1->Draw();
*/
	canvas->Update();

	xframe->Draw();

	std::cout<<"ABCDEF = "<<save_name<<std::endl;
	
	if(neg == 0 && pos == 0)
		canvas->Print(save_name + ".pdf[");
	canvas->Print(save_name + ".pdf");
	if(neg == 1 && pos == 1)
		canvas->Print(save_name + ".pdf]");
	 		
// 	if(fit == 0)
// 		return BW_mean.getVal();
// 	else if(fit == 1)
// 		return VOIG_mean.getVal();	
// 	else{
// 		std::cout<<"FUNZIONE INESISTENTE"<<std::endl;
		return -999;
// 	}

}


















TH2F* MakeMean(TH3F* h2, TString save_name, int fit, TH2F* h2_integral, float counting, int n_bins, float x_bounder){

	TH2F* tmp  = new TH2F("tmp", "tmp", n_bins, -x_bounder, x_bounder, n_bins, -x_bounder, x_bounder);
	std::cout<<"INIZIO = "<<n_bins<<"\t"<<-x_bounder<<"\t"<<x_bounder<<std::endl;

	int save_control_1 = -999;
	int save_control_2 = -999;
	
	for(int neg = 1; neg <= h2->GetNbinsY(); neg++){				
		for(int pos = 1; pos <= h2->GetNbinsX(); pos++){
		
			std::cout<<"STIAMO a = "<<neg<<"\t"<<pos<<std::endl;
			std::cout<<h2->GetXaxis()->GetBinLowEdge(pos)<<"	"<<h2->GetXaxis()->GetBinUpEdge(pos)<<"	"<<h2->GetYaxis()->GetBinLowEdge(neg)<<"	"<<h2->GetYaxis()->GetBinUpEdge(neg)<<std::endl;
			
			if(neg == 1 && pos == 1){
				save_control_1 = 0;
				save_control_2 = 0;
			}
			if(neg == h2->GetNbinsY() && pos == h2->GetNbinsX()){
				save_control_1 = 1;
				save_control_2 = 1;
			}
			TString NOME;
			NOME = Form("Mass POS %.3f %.3f NEG %.3f %.3f", h2->GetXaxis()->GetBinLowEdge(pos), h2->GetXaxis()->GetBinUpEdge(pos), h2->GetYaxis()->GetBinLowEdge(neg), h2->GetYaxis()->GetBinUpEdge(neg)); 
			TH1D* proiezione = h2->ProjectionZ(NOME,pos,pos,neg,neg);
// 			NOME = Form("Mass_%s_%s_d0BS", monteCarlo[MC].Data(), fitFunction[f].Data());

			float mean = FitMass(proiezione, save_name, fit, save_control_1, save_control_2);
			tmp->SetBinContent(pos, neg, mean);
			h2_integral->SetBinContent(neg, pos, proiezione->Integral()/counting);
			if(proiezione->Integral()/counting < 0.0001)
				tmp->SetBinContent(pos, neg, 0);			
		}
	}
	
	return tmp;								
				
}
