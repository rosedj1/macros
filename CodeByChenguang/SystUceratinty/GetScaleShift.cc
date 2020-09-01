#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "TFile.h"
#include "TLine.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TH1F.h"
using namespace RooFit;
Double_t massZ, massZErr, eta1, eta2, pT1, pT2, pterr1, pterr2, m1, m2, phi1, phi2, RelIso1, RelIso2; Int_t lep1_ecalDriven, lep2_ecalDriven;
TString year = "2016";
TString fs = "muon";
TString plots = "";

//set up range
Double_t low = 80;
Double_t high = 100;

RooDataHist* MakeDataSet(Double_t* ptcut, Double_t* etacut, bool isData,Int_t seed);
//void DataMCComp(RooDataHist* mc, RooDataHist* data, TString name);
Double_t* GetMean(RooDataHist* dataset, TString name, bool isData, bool fixtail);
void ClosurePlot(vector<Double_t> &results);
void SetAddress(TTree* t);
void MakeLUT(Double_t* sc1, Double_t* sc2, Double_t* sc3);

void GetScaleShift(){	
	
	Double_t ptbin[7] = {7,20,30,40,50,60,100};
	Double_t etabin[4] = {0,0.9,1.4,2.4};
	Double_t* ptcut = new Double_t[2];
	Double_t* etacut = new Double_t[2];
	Double_t* tmp_mc = new Double_t[2];
	Double_t* tmp_data = new Double_t[2];
	Double_t* result = new Double_t[2];//1st scale value, 2nd scale error, 3rd pt position, 4th pt error
	vector<Double_t> results;results.clear();
	for(Int_t i=0; i<3; i++){//eta loop 
		for(Int_t j=0; j<6; j++){//pt loop
			
			ptcut[0] = ptbin[j];ptcut[1] = ptbin[j+1];
			etacut[0] = etabin[i];etacut[1] = etabin[i+1];
			
			RooDataHist* mc = MakeDataSet(ptcut, etacut, 0, 5*j+i);
			RooDataHist* data = MakeDataSet(ptcut, etacut, 1, 7*j+i);
			TString name = "pt_"+to_string(ptcut[0]).substr(0,3)+"_"+to_string(ptcut[1]).substr(0,3)+"_eta_"+to_string(etacut[0]).substr(0,3)+"_"+to_string(etacut[1]).substr(0,3);
			//DataMCComp(mc, data, name);
			tmp_mc = GetMean(mc, name, 0, 1);
			tmp_data = GetMean(data, name ,1, 1);

			result[0] = (tmp_data[0]-tmp_mc[0])/91.19;
			result[1] = sqrt(tmp_mc[1]*tmp_mc[1]+tmp_data[1]*tmp_data[1])/91.19;
			results.push_back(result[0]);
			results.push_back(result[1]);
			delete mc; delete data;
		}
	}
	ClosurePlot(results);
}

void DataMCComp(RooDataHist* mc, RooDataHist* data, TString name){
	RooRealVar massZ("massZ","massZ",low,high);
	RooHistPdf* pdf_mc = new RooHistPdf("pdf_mc","pdf_mc",massZ,*mc);
	RooHistPdf* pdf_data = new RooHistPdf("pdf_data","pdf_data",massZ,*data);
	RooPlot* frame = massZ.frame(Bins(100));
	frame->SetTitle("");
	pdf_mc->plotOn(frame,LineColor(kBlue));
	pdf_data->plotOn(frame,LineColor(kRed));
	TCanvas c("c","c",1400,1000);
	c.cd();
	frame->Draw();
	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+fs+"/"+plots+"/Comparison_"+name+".png");
	delete pdf_mc; delete pdf_data; delete frame;
}

void ClosurePlot(vector<Double_t> &results){
	gStyle->SetTitleOffset(2,"Y");	
	
	Double_t x[6] = {13.5,25,35,45,55,80};
	Double_t x_error[6] = {6.5,5,5,5,5,20};
	Double_t* sc1 = new Double_t[6];Double_t* sc1Err = new Double_t[6];
	Double_t* sc2 = new Double_t[6];Double_t* sc2Err = new Double_t[6];
	Double_t* sc3 = new Double_t[6];Double_t* sc3Err = new Double_t[6];
	for(Int_t i=0; i<6; i++){
		sc1[i] = results.at(i*2);
		sc1Err[i] = results.at(i*2+1);
		sc2[i] = results.at(i*2+12);
		sc2Err[i] = results.at(i*2+13);
		sc3[i] = results.at(i*2+24);
		sc3Err[i] = results.at(i*2+25);
	}

	MakeLUT(sc1,sc2,sc3);

	TGraph* g1 = new TGraphErrors(6,x,sc1,x_error,sc1Err);
	TGraph* g2 = new TGraphErrors(6,x,sc2,x_error,sc2Err);
	TGraph* g3 = new TGraphErrors(6,x,sc3,x_error,sc3Err);

	g1->SetMarkerStyle(25);
	g1->SetMarkerColor(kBlue);
	g1->SetTitle("");
	g1->GetXaxis()->SetTitle(fs+" p_{T} (GeV)");
	g1->GetYaxis()->SetTitle("(m_{data}-m_{MC})/m_{PDG}");
	g1->GetXaxis()->SetLimits(0,100);
	g1->GetHistogram()->SetMaximum(0.004);
	g1->GetHistogram()->SetMinimum(-0.004);

	g2->SetMarkerStyle(25);
	g2->SetMarkerColor(kBlack);
	g3->SetMarkerStyle(25);
	g3->SetMarkerColor(kRed);

	TLine* l1 = new TLine(0,0,100,0);
	l1->SetLineStyle(kDashed);
	TLine* l2 = new TLine(0,0.001,100,0.001);
	l2->SetLineStyle(kDashed);
	TLine* l3 = new TLine(0,-0.001,100,-0.001);
	l3->SetLineStyle(kDashed);
	TLine* l4 = new TLine(0,-0.0004,100,-0.0004);
	l4->SetLineStyle(kDashed);
	TLine* l5 = new TLine(0,0.0004,100,0.0004);
	l5->SetLineStyle(kDashed);

	TLegend* legend = new TLegend(0.7,0.15,0.75,0.35);
        legend->AddEntry(g1, " |#eta| 0.0-0.9", "P");
        legend->AddEntry(g2, " |#eta| 0.9-1.4", "P");
        legend->AddEntry(g3, " |#eta| 1.4-2.4", "P");
        legend->SetTextSize(0.03);
        legend->SetLineWidth(0);
        legend->SetFillColor(0);
        legend->SetBorderSize();
	
	TCanvas* c = new TCanvas("c","c",1000,1000);
	c->cd();
	c->SetLeftMargin(0.14);	
	g1->Draw("ap");
	g2->Draw("p");
	g3->Draw("p");
	legend->Draw("same");
	l1->Draw("same");
	l2->Draw("same");
	l3->Draw("same");
	l4->Draw("same");
	l5->Draw("same");

	c->SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+plots+"/"+fs+"/Closure_test_"+fs+".png");
	c->SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+plots+"/"+fs+"/Closure_test_"+fs+".pdf");
	delete c; delete g1; delete g2; delete g3; delete legend; delete l1; delete l2; delete l3;

}

Double_t* GetMean(RooDataHist* dataset, TString name, bool isData, bool fixtail){	
	Double_t Zwidth;
	if(isData==0)Zwidth=2.44;
	if(isData==1)Zwidth=2.49;
	RooRealVar massZ("massZ","massZ",low,high);
	RooRealVar PDGmZ("PDGmZ","PDGmZ",91.19);
	RooRealVar PDGwZ("PDGwZ","PDGwZ",Zwidth);
	PDGmZ.setConstant(kTRUE);
	PDGwZ.setConstant(kTRUE);
	RooBreitWigner PDGBW("PDGBW","PDGBW",massZ,PDGmZ,PDGwZ);
	
	RooRealVar mean("mean","mean",0,-5,5);
	RooRealVar sigma("sigma","sigma",1,0,10);
	RooRealVar alpha1("alpha1","alpha1",1,0,10);
	RooRealVar n1("n1","n1",5,0,60);
	RooRealVar alpha2("alpha2","alpha2",1,0,10);
	RooRealVar n2("n2","n2",5,0,60);
	RooDoubleCB DCB("DCB","DCB",massZ,mean,sigma,alpha1,n1,alpha2,n2);
	RooFFTConvPdf BW_DCB("BW_DCB","BW_DCB",massZ,PDGBW,DCB);
	RooRealVar tau("tau","tau",0,-1,1);
	RooExponential bkg("bkg","bkg",massZ,tau);
	RooRealVar fsig("fsig","fsig",0.9,0.5,1);
	RooAddPdf Model("Model","Model",BW_DCB,bkg,fsig);
	Model.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1));
	cout<<"1st mean "<<mean.getVal()<<" "<<mean.getError()<<endl;
	cout<<"1st n1 "<<n1.getVal()<<endl;	
	
	if(fixtail==1){
		alpha1.setConstant(1);
		n1.setConstant(1);
		alpha2.setConstant(1);
		n2.setConstant(1);
		tau.setConstant(1);
		fsig.setConstant(1);
		mean.setConstant(0);
		sigma.setConstant(0);
		Model.fitTo(*dataset,SumW2Error(kTRUE),PrintLevel(-1));
		cout<<"2ed mean "<<mean.getVal()<<" "<<mean.getError()<<endl;
		cout<<"2ed n1 "<<n1.getVal()<<endl;
	}

	RooPlot* frame = massZ.frame(Bins(100));
	frame->SetTitle("");
	dataset->plotOn(frame);
	Model.plotOn(frame,LineColor(kBlue),LineWidth(1));
	TString mean_s, mean_error_s,sigma_s, alpha1_s, n1_s, alpha2_s, n2_s, chi2_s;
	mean_s = to_string(mean.getVal()); sigma_s = to_string(sigma.getVal());
	alpha1_s = to_string(alpha1.getVal()); n1_s = to_string(n1.getVal());
	alpha2_s = to_string(alpha2.getVal()); n2_s = to_string(n2.getVal());
	mean_error_s = to_string(mean.getError());
	Double_t chi2;
	if(fixtail==0)chi2 = frame->chiSquare(8);
	if(fixtail==1)chi2 = frame->chiSquare(2);
	chi2_s = to_string(chi2);
	TCanvas c("c","c",1400,1000);
	frame->Draw();
	TLatex *latex=new TLatex();
	latex->SetNDC();
        latex->SetTextSize(0.05);
        latex->SetTextFont(42);
        latex->SetTextAlign(23);
	latex->DrawLatex(0.7,0.8,"#chi^{2}/DOF="+chi2_s);
	latex->DrawLatex(0.7,0.7,"mean="+mean_s+"+/-"+mean_error_s);
	latex->DrawLatex(0.7,0.6,"sigma="+sigma_s);
	latex->DrawLatex(0.7,0.5,"alpha1="+alpha1_s);
	latex->DrawLatex(0.7,0.4,"n1="+n1_s);
	latex->DrawLatex(0.7,0.3,"alpha2="+alpha2_s);
	latex->DrawLatex(0.7,0.2,"n2="+n2_s);
	TString subname;
	if(isData==0)subname = "MC_";if(isData==1)subname = "Data_";
	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+plots+"/"+fs+"/"+subname+name+".png");
	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+plots+"/"+fs+"/"+subname+name+".pdf");
	delete frame; delete latex;
	Double_t* results = new Double_t[2];
	results[0] = mean.getVal();results[1] = mean.getError();
	return results;	
}

RooDataHist* MakeDataSet(Double_t* ptcut, Double_t* etacut, bool isData, Int_t seed){
	
	RooRealVar massZ_("massZ","massZ",low,high);
	massZ_.setBins(300);
	RooArgSet argset(massZ_);
	//RooDataSet* dataset = new RooDataSet("dataset","dataset",argset);
	TFile *f=new TFile();
	
/*	//HIG-16-041
	if(isData==0&&year=="2016"&&fs=="muon")f = new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/DY_2016MC_v3_20170312_afterApproval/DYJetsToLL_M-50_kalman_v4_m2mu.root");	
	if(isData==1&&year=="2016"&&fs=="muon")f = new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/Data_2016_v1_20170312_afterApproval/DoubleLepton_m2mu.root");
	if(isData==0&&year=="2016"&&fs=="electron")f = new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/DY_2016MC_v3_20170312_afterApproval/DYJetsToLL_M-50_kalman_v4_m2e.root");
	if(isData==1&&year=="2016"&&fs=="electron")f = new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/Data_2016_v1_20170312_afterApproval/DoubleLepton_m2e.root");
	*/
	
	//Run2 Legacy
	if(isData==0&&year=="2016"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2016_m2mu.root");
	if(isData==0&&year=="2017"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2017_m2mu.root");
	if(isData==0&&year=="2018"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2018_m2mu.root");
	if(isData==0&&year=="2016"&&fs=="electron")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2016_m2e.root");
	if(isData==0&&year=="2017"&&fs=="electron")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2017_m2e.root");
	if(isData==0&&year=="2018"&&fs=="electron")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/2018_m2e.root");
	
	if(isData==1&&year=="2016"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2016/Muon/SingleMuon_m2mu.root");
	if(isData==1&&year=="2017"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2017/Muon/SingleMuon_m2mu.root");
	if(isData==1&&year=="2018"&&fs=="muon")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2018/Muon/SingleMuon_m2mu.root");
	if(isData==1&&year=="2016"&&fs=="electron")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2016/Electron/SingleElectron_m2e.root");
	if(isData==1&&year=="2017"&&fs=="electron")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2017/Electron/SingleElectron_m2e.root");
	if(isData==1&&year=="2018"&&fs=="electron")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/withRelIsoinformation/Data/2018/Electron/EGamma_m2e.root");
	

	TTree* t = (TTree*)f->Get("passedEvents");
	Int_t sum = t->GetEntries();
	SetAddress(t);
	TH1D* h = new TH1D("h","h",300,low,high);
	TRandom3 rand;
	rand.SetSeed(12345679);
	Double_t randcut;
	Double_t factor1, factor2;
	for(Int_t i=0; i<sum; i++){
		t->GetEntry(i);
		randcut=rand.Gaus(0,1);
		if(massZ<high&&massZ>low&&RelIso1<0.35&&RelIso2<0.35){
			if((pT1>ptcut[0]&&pT1<ptcut[1]&&eta1<etacut[1]&&eta1>etacut[0]&&randcut>0)||(randcut<0&&eta2<etacut[1]&&eta2>etacut[0]&&pT2>ptcut[0]&&pT2<ptcut[1])){
				h->Fill(massZ);
			}
		}		
	}
	RooDataHist* dataset = new RooDataHist("dataset","dataset",argset,Import(*h));
	f->Close();
	delete f; 
	return dataset;

}
				
void SetAddress(TTree* t){
	t->SetBranchAddress("massZ",&massZ);
	t->SetBranchAddress("massZErr",&massZErr);
	t->SetBranchAddress("eta1",&eta1);
	t->SetBranchAddress("eta2",&eta2);
	t->SetBranchAddress("pT1",&pT1);
	t->SetBranchAddress("pT2",&pT2);
	t->SetBranchAddress("pterr1",&pterr1);
	t->SetBranchAddress("pterr2",&pterr2);
	t->SetBranchAddress("m1",&m1);
	t->SetBranchAddress("m2",&m2);
	t->SetBranchAddress("phi1",&phi1);
	t->SetBranchAddress("phi2",&phi2);
	t->SetBranchAddress("lep1_ecalDriven",&lep1_ecalDriven);
	t->SetBranchAddress("lep2_ecalDriven",&lep2_ecalDriven);
	t->SetBranchAddress("RelIso1",&RelIso1);
	t->SetBranchAddress("RelIso2",&RelIso2);

}		

void MakeLUT(Double_t* sc1, Double_t* sc2, Double_t* sc3){
	
	Double_t ptbin[7] = {7,20,30,40,50,60,100};
	Double_t etabin[4] = {0,0.9,1.4,2.4};
	TH2D* LUT = new TH2D("LUT","LUT",6,ptbin,3,etabin);
	for(Int_t j=0; j<6; j++){
		LUT->SetBinContent(j+1,1,sc1[j]);
		LUT->SetBinContent(j+1,2,sc2[j]);
		LUT->SetBinContent(j+1,3,sc3[j]);
	}
	TCanvas c("c","c",1400,1000);
	c.cd();
        LUT->Draw("TEXT");
        c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+plots+"/"+fs+"/LUT.png");
	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+plots+"/"+fs+"/LUT.pdf");
        TFile* f = new TFile("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/ScaleShift/"+plots+"/"+fs+"/LUT.root","RECREATE");
        f->cd();
        LUT->Write();
        f->Close();

}
