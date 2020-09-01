#include "ScaleUncertantiesStep2.h"

std::vector<TString> ProdMode;
std::vector<TString> ProdMode_File;
std::vector<int> massPoint;
std::vector<TString> year;
std::vector<TString> category;
std::vector<TString> decayMode;

void SignalLineShape_v2(bool refit){

    gROOT->Reset();
    gROOT->SetBatch();	
    
    TString directory;

	Bool_t passedFullSelection;
	Int_t finalState, EventCat;
	Float_t eventWeight;
	Float_t mass4l, mass4lErr, mass4lREFIT, mass4lErrREFIT;
	Float_t mass4mu, mass4e, mass2e2mu;

	ProdMode.clear();
	ProdMode.push_back("ggH");
	ProdMode.push_back("VBF");
	ProdMode.push_back("WplusH");
	ProdMode.push_back("WminusH");
	ProdMode.push_back("ZH");
// 	ProdMode.push_back("ttH");
	
	ProdMode_File.clear();
	ProdMode_File.push_back("GluGluHToZZTo4L");
	ProdMode_File.push_back("VBF_HToZZTo4L");
	ProdMode_File.push_back("WplusH_HToZZTo4L");
	ProdMode_File.push_back("WminusH_HToZZTo4L");
	ProdMode_File.push_back("ZH_HToZZ_4LFilter");
// 	ProdMode_File.push_back("ttH_HToZZ_4LFilter");
	
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
	
	category.clear();
	category.push_back("VBF_1jet");
	category.push_back("VBF_2jet");
	category.push_back("VH_leptonic");
	category.push_back("VH_hadronic");
	category.push_back("ttH");
	category.push_back("untagged");
	category.push_back("VH_MET");
	
	decayMode.clear();
	decayMode.push_back("4mu");
	decayMode.push_back("4e");
	decayMode.push_back("2e2mu");
	
	Double_t FitParamterValue[year.size()][ProdMode.size()][category.size()][decayMode.size()][3];
// 	Double_t FitParamterValue[year.size()][3] = {0};
	
	
	TH1F* h_yield[year.size()][ProdMode.size()][category.size()][decayMode.size()];
// 	TH1F* h_mass[year.size()][massPoint.size()][ProdMode.size()][category.size()][decayMode.size()];
	TH1F* h_mass[year.size()][massPoint.size()][ProdMode.size()][decayMode.size()];
	for(int y = 0; y < year.size(); y++){
		for(int mP = 0; mP < massPoint.size(); mP++){
			for(int PM = 0; PM < ProdMode.size(); PM++){
// 				for(int c = 0; c < category.size(); c++){
					for(int fs = 0; fs < decayMode.size(); fs++){
						TString histo_name = Form("HistoMass_%s_M%i_%s_%s", year[y].Data(), massPoint.at(mP), ProdMode[PM].Data(), decayMode[fs].Data());
// 						h_mass[y][mP][PM][c][fs] = new TH1F(histo_name, histo_name, 100, 105, 140);
						h_mass[y][mP][PM][fs] = new TH1F(histo_name, histo_name, 100, 105, 140);
// 						std::cout<<histo_name<<std::endl;
					}
// 				}
			}
		}
	}
	
	TLegend *legend_fs = new TLegend(0.75,0.75,0.9,0.9);

	TTree* tree;

	Color_t color;
	TString nome_file;
	
	for(int y = 0; y < year.size(); y++){
		for(int PM = 0; PM < ProdMode.size(); PM++){
		    directory = "SignalLineShape/";
		    if(refit) directory = "SignalLineShape_REFIT/";
			TString execute = "mkdir " + directory;
			gSystem->Exec(execute);
			directory += year.at(y) + "/";
			execute = "mkdir " + directory;
			gSystem->Exec(execute);
			directory += ProdMode.at(PM) + "/";
			execute = "mkdir " + directory;
			gSystem->Exec(execute);

			for(int mP = 0; mP <  massPoint.size(); mP++){
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
				
/*				Long64_t nentries = tree->GetEntries();
				std::cout<<nentries<<std::endl;									
		 		for(int entry = 0; entry < nentries; entry++){
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
												
						if(finalState > 2)												
							h_mass[y][mP][PM][2]->Fill(mass4l);//, eventWeight*(mP+1));
						else
							h_mass[y][mP][PM][finalState-1]->Fill(mass4l);//, eventWeight*(mP+1));

				}*/
				
				for(int fs = 0; fs < decayMode.size(); fs++){
					TString cut_name;
					if(fs < 2)
						cut_name = Form("(passedFullSelection && finalState == %d)", fs+1);
// 						cut_name = Form("eventWeight * (passedFullSelection && finalState == %d)", fs+1);
					else
						cut_name = Form("(passedFullSelection && (finalState == 3 || finalState == 4))");
// 						cut_name = Form("eventWeight * (passedFullSelection && (finalState == 3 || finalState == 4))");
					TString histo_name = Form("HistoYield_%s_%s_%s_%d", year[y].Data(), ProdMode[PM].Data(), decayMode[fs].Data(), massPoint.at(mP));
					TH1F* h_tmp = new TH1F("h_tmp", histo_name, 100, 105, 140);
					gROOT->cd();
					if(refit)
						tree->Draw("mass4lREFIT>>h_tmp(100,105,140)", cut_name);
					else
						tree->Draw("mass4l>>h_tmp(100,105,140)", cut_name);
					TH1F* h_tmp_2 = (TH1F*)gDirectory->Get("h_tmp");
					h_mass[y][mP][PM][fs] = (TH1F*)h_tmp_2->Clone();
// 					h_mass[y][mP][PM][fs] = (TH1F*)gDirectory->Get("h_tmp");

					TString nome_canvas = Form("Mass %s %d %s %s", ProdMode_File[PM].Data(), massPoint.at(mP), year[y].Data(), decayMode[fs].Data());
					TCanvas *c_mass = new TCanvas(canvas_name, canvas_name, 700, 500);
					c_mass->SetFrameFillColor(0);
					c_mass->cd(1)->SetBottomMargin(0.2);
					h_mass[y][mP][PM][fs]->Draw();
					h_mass[y][mP][PM][fs]->SetTitle(histo_name);
// 					h_mass[y][mP][PM][fs]->SetBins(100, 0, 140);
					h_mass[y][mP][PM][fs]->GetXaxis()->SetRangeUser(105, 140);
					c_mass->Update();
					if(mP == 0 && fs == 0)
						c_mass->Print(directory + "mass.pdf[");
					c_mass->Print(directory + "mass.pdf");
					if(mP == massPoint.size()-1 && fs == decayMode.size()-1)
						c_mass->Print(directory + "mass.pdf]");
// 					std::cout<<h_mass[y][mP][PM][fs]->Integral()<<std::endl;	
				}	
		    }
		}
	}

	
	float chi_square_all[3][3];
	std::vector<float> chi_square_125;
	
	for(int y = 0; y < year.size(); y++){	
		for(int PM = 0; PM < ProdMode.size(); PM++){

		    directory = "SignalLineShape/";
		    if(refit) directory = "SignalLineShape_REFIT/";
			directory += year.at(y) + "/" + ProdMode.at(PM) + "/";;
		
			ofstream myfile;
			TString nome_file_parameter = Form("SignalLineShape_%s_%s.txt", ProdMode[PM].Data(), year[y].Data());
			nome_file_parameter = directory + "/" + nome_file_parameter;
			myfile.open(nome_file_parameter);

			for(int fs = 0; fs < decayMode.size(); fs++){
				RooRealVar x("x", "Mass (GeV/c^{2})", 105, 140);
				RooRealVar weight("weight", "weight", 0, 10);
	
				RooRealVar DSCB_mean_p0("mean_p0", "mean_p0", 125, 115, 135);
				RooRealVar DSCB_sigma_p0("sigma_p0", "sigma_p0", 1, 0, 30);
				RooRealVar DSCB_alphaL_p0("alphaL_p0", "alphaL_p0", 1, 0, 30);
				RooRealVar DSCB_expL_p0("expL_p0", "expL_p0", 1, 0, 30);
				RooRealVar DSCB_alphaR_p0("alphaR_p0", "alphaR_p0", 1, 0, 30);
				RooRealVar DSCB_expR_p0("expR_p0", "expR_p0", 1, 0, 30);

				RooRealVar DSCB_mean_p1("mean_p1", "mean_p1",-10, 10);
				RooRealVar DSCB_sigma_p1("sigma_p1", "sigma_p1", 0, -30, 30);
				RooRealVar DSCB_alphaL_p1("alphaL_p1", "alphaL_p1", 0, -30, 30);
				RooRealVar DSCB_expL_p1("expL_p1", "expL_p1", 0, -30, 30);
				RooRealVar DSCB_alphaR_p1("alphaR_p1", "alphaR_p1", 0, -30, 30);
				RooRealVar DSCB_expR_p1("expR_p1", "expR_p1", 0, -30, 30);

				RooRealVar DSCB_mean_p2("mean_p2", "mean_p2", -10, 10);
				RooRealVar DSCB_sigma_p2("sigma_p2", "sigma_p2", 0, -30, 30);
				RooRealVar DSCB_alphaL_p2("alphaL_p2", "alphaL_p2", 0, -30, 30);
				RooRealVar DSCB_expL_p2("expL_p2", "expL_p2", 0, -30, 30);
				RooRealVar DSCB_alphaR_p2("alphaR_p2", "alphaR_p2", 0, -30, 30);
				RooRealVar DSCB_expR_p2("expR_p2", "expR_p2", 0, -30, 30);
	
				RooArgList intersection_params = RooArgList(DSCB_mean_p0, DSCB_sigma_p0, DSCB_alphaL_p0, DSCB_expL_p0, DSCB_alphaR_p0, DSCB_expR_p0);
				RooArgList slope_params        = RooArgList(DSCB_mean_p1, DSCB_sigma_p1, DSCB_alphaL_p1, DSCB_expL_p1, DSCB_alphaR_p1, DSCB_expR_p1);

				RooCategory rc_signals = RooCategory("signals","signals");
		
				std::vector<RooRealVar> massPointFit;
		
				std::vector<RooFormulaVar> mean;
				std::vector<RooFormulaVar> sigma;
				std::vector<RooFormulaVar> alphaL;
				std::vector<RooFormulaVar> expL;
				std::vector<RooFormulaVar> alphaR;
				std::vector<RooFormulaVar> expR;
	
				std::vector<RooMyPDF_DSCB> DSCB;
				std::vector<RooDataHist> h_Higgs;
		
				std::vector<TString> cat_name;
				std::vector<TString> dataset_name;

				for(int mP = 0; mP <  massPoint.size(); mP++){
					cat_name.push_back(Form("cat_signal_%d", massPoint.at(mP)));
					dataset_name.push_back(Form("dataset_name_%d", massPoint.at(mP)));
					rc_signals.defineType(cat_name.at(mP));
					TString nome_var =  Form("_%d", massPoint.at(mP));		
// 					std::cout<<cat_name<<"\t"<<nome_var<<std::endl;
					massPointFit.push_back(RooRealVar("massPointFit"+nome_var, "massPointFit"+nome_var, massPoint.at(mP)));
// 					std::cout<<massPointFit.at(mP).getVal()<<std::endl;
				}
				for(int mP = 0; mP <  massPoint.size(); mP++){
					TString nome_var =  Form("_%d", massPoint.at(mP));
					mean.push_back(RooFormulaVar("mean"+nome_var,"@0 + @1 * (@3 - 125) + @2 * (@3 - 125) * (@3 - 125)", RooArgList(DSCB_mean_p0, DSCB_mean_p1, DSCB_mean_p2, massPointFit.at(mP))));
					sigma.push_back(RooFormulaVar("sigma"+nome_var,"@0 + @1 * (@3 - 125) + @2 * (@3 - 125) * (@3 - 125)", RooArgList(DSCB_sigma_p0, DSCB_sigma_p1, DSCB_sigma_p2, massPointFit.at(mP))));
					alphaL.push_back(RooFormulaVar("alphaL"+nome_var,"@0 + @1 * (@3 - 125) + @2 * (@3 - 125) * (@3 - 125)", RooArgList(DSCB_alphaL_p0, DSCB_alphaL_p1, DSCB_alphaL_p2, massPointFit.at(mP))));
					expL.push_back(RooFormulaVar("expL"+nome_var,"@0 + @1 * (@3 - 125) + @2 * (@3 - 125) * (@3 - 125)", RooArgList(DSCB_expL_p0, DSCB_expL_p1, DSCB_expL_p2, massPointFit.at(mP))));
					alphaR.push_back(RooFormulaVar("alphaR"+nome_var,"@0 + @1 * (@3 - 125) + @2 * (@3 - 125) * (@3 - 125)", RooArgList(DSCB_alphaR_p0, DSCB_alphaR_p1, DSCB_alphaR_p2, massPointFit.at(mP))));
					expR.push_back(RooFormulaVar("expR"+nome_var,"@0 + @1 * (@3 - 125) + @2 * (@3 - 125) * (@3 - 125)", RooArgList(DSCB_expR_p0, DSCB_expR_p1, DSCB_expR_p2, massPointFit.at(mP))));
// 					mean.push_back(RooFormulaVar("mean"+nome_var,"@0 + @1 * (@2 - 125)", RooArgList(DSCB_mean_p0, DSCB_mean_p1, massPointFit.at(mP))));
// 					sigma.push_back(RooFormulaVar("sigma"+nome_var,"@0 + @1 * (@2 - 125)", RooArgList(DSCB_sigma_p0, DSCB_sigma_p1, massPointFit.at(mP))));
// 					alphaL.push_back(RooFormulaVar("alphaL"+nome_var,"@0 + @1 * (@2 - 125)", RooArgList(DSCB_alphaL_p0, DSCB_alphaL_p1, massPointFit.at(mP))));
// 					expL.push_back(RooFormulaVar("expL"+nome_var,"@0 + @1 * (@2 - 125)", RooArgList(DSCB_expL_p0, DSCB_expL_p1, massPointFit.at(mP))));
// 					alphaR.push_back(RooFormulaVar("alphaR"+nome_var,"@0 + @1 * (@2 - 125)", RooArgList(DSCB_alphaR_p0, DSCB_alphaR_p1, massPointFit.at(mP))));
// 					expR.push_back(RooFormulaVar("expR"+nome_var,"@0 + @1 * (@2 - 125)", RooArgList(DSCB_expR_p0, DSCB_expR_p1, massPointFit.at(mP))));
		
// 					h_Higgs.push_back(RooDataHist("h_HiggsMass"+nome_var, "h_HiggsMass"+nome_var, RooArgSet(x, weight), h_mass[0][mP][0][0]));
					h_Higgs.push_back(RooDataHist("h_HiggsMass"+nome_var, "h_HiggsMass"+nome_var, RooArgSet(x), h_mass[y][mP][PM][fs]));
//	 				std::cout<<"3\t"<<"\t"<<nome_var<<std::endl;
				}

				for(int mP = 0; mP <  massPoint.size(); mP++){
					TString nome_var =  Form("_%d", massPoint.at(mP));
					DSCB.push_back(RooMyPDF_DSCB("DSCB"+nome_var, "DSCB"+nome_var, x, mean.at(mP), sigma.at(mP), alphaL.at(mP), expL.at(mP), alphaR.at(mP), expR.at(mP)));
//			 		std::cout<<"6\t"<<"\t"<<nome_var<<std::endl;
				}
		
				DSCB_mean_p1.setConstant();
				DSCB_sigma_p1.setConstant();
				DSCB_alphaL_p1.setConstant();
				DSCB_expL_p1.setConstant();
				DSCB_alphaR_p1.setConstant();
				DSCB_expR_p1.setConstant();
	
				DSCB_mean_p2.setConstant();
				DSCB_sigma_p2.setConstant();
				DSCB_alphaL_p2.setConstant();
				DSCB_expL_p2.setConstant();
				DSCB_alphaR_p2.setConstant();
				DSCB_expR_p2.setConstant();
								
				RooFitResult* r_125;
				TCanvas *c_mass = new TCanvas("mass", "mass", 700, 500);
// 				pad11->SetBottomMargin(0.05);
				c_mass->SetFrameFillColor(0);
				c_mass->cd(1)->SetBottomMargin(0.2);
				RooPlot* xframe = x.frame(Title("mass"));
				h_Higgs.at(0).plotOn(xframe);
				r_125 = DSCB.at(0).fitTo(h_Higgs.at(0), Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
				DSCB.at(0).plotOn(xframe,RooFit::LineColor(kBlue));
				DSCB.at(0).paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));
				xframe->Draw();
				TLatex *tex1;
				fit_param_latex = Form("#chi^{2}/dof = %.3f", xframe->chiSquare(r_125->floatParsFinal().getSize())/ r_125->floatParsFinal().getSize());
				tex1 = new TLatex(0.65,0.85, fit_param_latex);
				tex1->SetNDC();
				tex1->Draw();
				TString mass_save = directory + "Test_125_mass_" + decayMode.at(fs) + "_" + ProdMode.at(PM) + "_" + year.at(y);
				c_mass->Print(mass_save + ".png");
				c_mass->Print(mass_save + ".pdf");
				TCanvas *c_matrix = new TCanvas("matrix", "matrix", 700, 500);
				TH2F* correlation_matrix = (TH2F*)r_125->correlationHist();
				correlation_matrix->Draw("colz");
// 				SetOwnership(c,False)
				mass_save = directory + "Test_125_matrix_" + decayMode.at(fs) + "_" + ProdMode.at(PM) + "_" + year.at(y);
// 				c_matrix->Print(mass_save + ".png");
// 				c_matrix->Print(mass_save + ".pdf");
				
				chi_square_125.push_back(xframe->chiSquare(r_125->floatParsFinal().getSize())/ r_125->floatParsFinal().getSize());
		
				DSCB_mean_p1.setConstant(kFALSE);
				DSCB_sigma_p1.setConstant(kFALSE);
				DSCB_alphaL_p1.setConstant(kFALSE);
				DSCB_expL_p1.setConstant(kFALSE);
				DSCB_alphaR_p1.setConstant(kFALSE);
				DSCB_expR_p1.setConstant(kFALSE);
	
				DSCB_mean_p2.setConstant(kFALSE);
				DSCB_sigma_p2.setConstant(kFALSE);
				DSCB_alphaL_p2.setConstant(kFALSE);
				DSCB_expL_p2.setConstant(kFALSE);
				DSCB_alphaR_p2.setConstant(kFALSE);
				DSCB_expR_p2.setConstant(kFALSE);
	
				DSCB_mean_p0.setConstant();
				DSCB_sigma_p0.setConstant();
				DSCB_alphaL_p0.setConstant();
				DSCB_expL_p0.setConstant();
				DSCB_alphaR_p0.setConstant();
				DSCB_expR_p0.setConstant();		
		
				RooSimultaneous sim_pdf("sim_pdf", "Simultaneous pdf", rc_signals);
	
				for(int mP = 0; mP <  massPoint.size(); mP++)
					sim_pdf.addPdf(DSCB.at(mP), cat_name.at(mP));		
		
// 	RooDataHist rds_all_signals("combData","combined data", RooArgSet(x, weight),Index(rc_signals),
		
				RooDataHist rds_all_signals("combData","combined data", RooArgSet(x),Index(rc_signals),
// 				Import(cat_name.at(0), h_Higgs.at(0)),
							Import(cat_name.at(0), h_Higgs.at(0)),
							Import(cat_name.at(1), h_Higgs.at(1)),
							Import(cat_name.at(2), h_Higgs.at(2)),
							Import(cat_name.at(3), h_Higgs.at(3)),
							Import(cat_name.at(4), h_Higgs.at(4)));
			
				RooFitResult* r_all;
				r_all = sim_pdf.fitTo(rds_all_signals, Save(kTRUE), SumW2Error(kTRUE), Verbose(kFALSE), PrintLevel(-1), Warnings(kFALSE), NumCPU(12), Timer(kTRUE));
				r_all->Print();
				TCanvas *c_matrix_all = new TCanvas("c_matrix_all", "c_matrix_all", 700, 500);
				correlation_matrix = (TH2F*)r_all->correlationHist();
				correlation_matrix->Draw("colz");
				mass_save = directory + "Test_All_matrix_" + decayMode.at(fs) + "_" + ProdMode.at(PM) + "_" + year.at(y);
// 				c_matrix_all->Print(mass_save + ".png");
// 				c_matrix_all->Print(mass_save + ".pdf");
			
				int color_shift = 0;
		
				for(int mP = 0; mP <  massPoint.size(); mP++){
		
					TCanvas* c1 = new TCanvas(cat_name.at(mP), cat_name.at(mP), 700, 500);
// 					c1->SetFrameFillColor(0);
// 					c1->cd(1)->SetBottomMargin(0.2);
				   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
				   	pad11->Draw();
				   	pad11->cd();

					RooPlot* xframe = x.frame(Title(cat_name.at(mP)));
					TString cut = Form("signals==signals::%s", cat_name[mP].Data());
					rds_all_signals.plotOn(xframe, LineColor(kBlack),
// 									MarkerSize(0),
										Cut(cut)
// 									Cut(cat_name[mP].Data())
								);
					sim_pdf.plotOn(xframe,
							Slice(rc_signals, cat_name.at(mP)),
								ProjWData(RooArgSet(rc_signals), rds_all_signals),
								LineColor(kPink+color_shift)
								);
// 					sim_pdf.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));

					RooPlot *xframePull = x.frame("Pull");
					RooHist* hpull_CB = xframe->pullHist();
					hpull_CB->SetMarkerColor(8);
					hpull_CB->SetMarkerSize(0.75);
					xframePull->addObject(hpull_CB, "p");

					xframe->Draw();
					fit_param_latex = Form("#chi^{2}/dof = %.3f", xframe->chiSquare(r_all->floatParsFinal().getSize())/ r_all->floatParsFinal().getSize());
					tex1 = new TLatex(0.65,0.85, fit_param_latex);
					tex1->SetNDC();
					tex1->Draw();
					tex1 = new TLatex(0.15,0.85, decayMode.at(fs));
					tex1->SetNDC();
					tex1->Draw();
					color_shift++;
					
					c1->Update();
					c1->cd();

					TPad* padPull =  new  TPad("padPull","padPull", 0., 0.0, 1., 0.25);
// 	padPull->SetBottomMargin(0.05);
					padPull->Draw();
					padPull->cd();
					xframePull->SetTitle("Pull DSCB");
					xframePull->SetMarkerSize(0.05);
					xframePull->SetMarkerColor(4);
					xframePull->GetYaxis()->SetLabelSize(0.1);
					xframePull->GetXaxis()->SetLabelSize(0.1);
					xframePull->SetMinimum(-3.);
					xframePull->SetMaximum(3.);
					xframePull->Draw();

					c1->Update();
					c1->cd();
		
					mass_save = directory + "Test_All_mass_" + decayMode.at(fs) + "_" + ProdMode.at(PM) + "_" + year.at(y);
		
					if(mP == 0)
						c1->Print(mass_save + ".pdf[");
					c1->Print(mass_save + ".pdf");
					if(mP == massPoint.size()-1)
						c1->Print(mass_save + ".pdf]");		
		
					chi_square_all[fs][mP] = xframe->chiSquare(r_all->floatParsFinal().getSize())/ r_all->floatParsFinal().getSize();
				}

				myfile<<decayMode.at(fs)<<std::endl;
				myfile<<"Param CB0 = "<<DSCB_mean_p0.getVal()<<"; "<<DSCB_sigma_p0.getVal()<<"; "<<DSCB_alphaL_p0.getVal()<<"; "<<DSCB_expL_p0.getVal()<<"; "<<DSCB_alphaR_p0.getVal()<<"; "<<DSCB_expR_p0.getVal()<<std::endl;
				myfile<<"Param CB1 = "<<DSCB_mean_p1.getVal()<<"; "<<DSCB_sigma_p1.getVal()<<"; "<<DSCB_alphaL_p1.getVal()<<"; "<<DSCB_expL_p1.getVal()<<"; "<<DSCB_alphaR_p1.getVal()<<"; "<<DSCB_expR_p1.getVal()<<std::endl;
				myfile<<"Param CB2 = "<<DSCB_mean_p2.getVal()<<"; "<<DSCB_sigma_p2.getVal()<<"; "<<DSCB_alphaL_p2.getVal()<<"; "<<DSCB_expL_p2.getVal()<<"; "<<DSCB_alphaR_p2.getVal()<<"; "<<DSCB_expR_p2.getVal()<<std::endl;

			} //for on fs
		} // for on ProdMode
	} // for on year
	
		
	std::cout<<"ESCO"<<std::endl;
	
	for(int fs = 0; fs < decayMode.size(); fs++){
		std::cout<<decayMode.at(fs)<<std::endl;
		std::cout<<"First 125 = "<<chi_square_125.at(fs)<<std::endl;
		std::cout<<"Sim fit = ";
		for(int mP = 0; mP <  massPoint.size(); mP++)
			std::cout<<chi_square_all[fs][mP]<<"\t";
		std::cout<<"\n"<<std::endl;
	}
			
	
}
