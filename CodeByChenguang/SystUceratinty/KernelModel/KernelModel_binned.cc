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

void SetAddress(TTree* t);
RooDataHist* GetDataSet(TString year, TString fs, Double_t a1, Double_t a2, Double_t b1, Double_t b2);
Double_t* GetShift(TString year, TString fs, Double_t a1, Double_t a2, Double_t b1, Double_t b2, RooDataHist* Data);

Double_t GENmassZ, massZ, massZErr, weight, pT1, pT2, eta1, eta2, m1, m2, phi1, phi2, pterr1, pterr2;

void KernelModel_binned(){
	
	TString year = "2016";
	TString fs = "2e";
	
	RooDataHist* Data = GetDataSet(year, fs, 0.9, 1.4, 5, 20);
	Data->Print();
	
	Double_t* results = new Double_t[2];
	results = GetShift(year, fs, 0.9, 1.4, 5, 20, Data);
	cout<<results[0]<<endl;
	cout<<results[1]<<endl;	

	
	
	
	
}







RooDataHist* GetDataSet(TString year, TString fs, Double_t a1, Double_t a2, Double_t b1, Double_t b2){
	
	TFile* f = new TFile();
	
	if(year=="2016"&&fs=="2e")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2016/Electron/SingleElectron_m2e.root");
	if(year=="2017"&&fs=="2e")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2016/Electron/SingleElectron_m2e.root");
	if(year=="2018"&&fs=="2e")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2017/Electron/EGamma_m2e.root");
	if(year=="2016"&&fs=="2mu")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2017/Muon/SingleMuon_m2mu.root");
	if(year=="2017"&&fs=="2mu")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2018/Muon/SingleMuon_m2mu.root");
	if(year=="2018"&&fs=="2mu")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2018/Muon/SingleMuon_m2mu.root");
	
	RooRealVar massZ_d("massZ","massZ",80,100);
	RooArgSet argset(massZ_d);
	RooDataSet* dataset = new RooDataSet("dataset","dataset",argset);
	massZ_d.setBins(1000,"cache");

	TH1D* h = new TH1D("h","h",400,80,100);

	TTree* t = (TTree*)f->Get("passedEvents");
	SetAddress(t);
	Int_t sum = t->GetEntries();
	Int_t counter = 0;
	for(Int_t i=0; i<sum; i++){
		t->GetEntry(i);
		if(massZErr<7.2&&massZErr>0.2&&massZ<100&&massZ>80){
			if((abs(eta1)>a1&&abs(eta1)<a2&&pT1>b1&&pT1<b2)||(abs(eta2)>a1&&abs(eta2)<a2&&pT2>b1&&pT2<b2)){
				
				massZ_d.setVal(massZ);
				dataset->add(argset);
				counter = counter+1;
				if(counter==20000)break;
			}
		}
	}

	RooDataHist* dataset_binned = dataset->binnedClone();
	
	
	return dataset_binned;
	
				
	
	
	

}


Double_t* GetShift(TString year, TString fs, Double_t a1, Double_t a2, Double_t b1, Double_t b2, RooDataHist* Data){
	
		
	TFile* f = new TFile();
	if(year=="2016"&&fs=="2e")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2016_m2e.root");
	if(year=="2017"&&fs=="2e")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2017_m2e.root");
	if(year=="2018"&&fs=="2e")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2018_m2e.root");
	if(year=="2016"&&fs=="2mu")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2016_m2mu.root");
	if(year=="2017"&&fs=="2mu")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2017_m2mu.root");
	if(year=="2018"&&fs=="2mu")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2018_m2mu.root");
	
	TTree* t = (TTree*)f->Get("passedEvents");
	SetAddress(t);
	Double_t sum = t->GetEntries();
	Int_t counter = 0;
	TString index;
	
	RooArgList FuncList("FuncList");
	RooArgList CoefList("CoefList");

	RooRealVar* massZ_d = new RooRealVar("massZ","massZ",80,100);

	for(Int_t i=0; i<sum; i++){
		t->GetEntry(i);
		if(i%1000000==0)cout<<i<<endl;
		if(massZ<100&&massZ>80&&massZErr<7.2&&massZErr>0.2){
			if((abs(eta1)>a1&&abs(eta1)<a2&&b1<pT1&&pT1<b2)||(abs(eta2)>a1&&abs(eta2)<a2&&pT2<b2&&pT2>b1)){
				index = to_string(counter);	
				counter = counter+1;
				
				RooRealVar* mean = new RooRealVar("mean"+index,"mean"+index,massZ);
				RooRealVar* sigma = new RooRealVar("sigma"+index,"sigma"+index,massZErr);
				RooGaussian* kernel = new RooGaussian("kernel"+index,"kernel"+index,*massZ_d,*mean,*sigma);
				RooRealVar* coef = new RooRealVar("coef"+index,"coef"+index,0.0001);

				FuncList.add(*kernel);
				CoefList.add(*coef);
								
				
				
			}	
		}
		if(counter==20000)break;
	}
	cout<<counter<<"kernels for model building"<<endl;	
	RooRealSumPdf SumKernel("SumKernel","SumKernel",FuncList,CoefList);
	
	cout<<"convolute gauss for final model"<<endl;
	RooRealVar shift("shift","shift",0.362147);
	RooRealVar smear("smear","smear",0.012868);
	RooGaussian gauss("gauss","gauss",*massZ_d,shift,smear);
	RooFFTConvPdf Model("Model","Model",*massZ_d,SumKernel,gauss);
	cout<<"MODEL BUILDED"<<endl;
	
	//Model.fitTo(*Data,SumW2Error(kTRUE),PrintLevel(-1));
	RooPlot* frameMC = massZ_d->frame(Bins(100));
	frameMC->SetTitle("");
	Data->plotOn(frameMC);
	Model.plotOn(frameMC,LineColor(2),LineWidth(3));
	TCanvas c("c","c",1400,1000);
	c.cd();
	Double_t chi2 = frameMC->chiSquare(2);
	TString shift_s, smear_s, chi2_s;
	shift_s = to_string(shift.getVal());
	smear_s = to_string(smear.getVal());
	chi2_s = to_string(chi2);
	TLatex *latex=new TLatex();
    latex->SetNDC();
    latex->SetTextSize(0.05);
   	latex->SetTextFont(42);
   	latex->SetTextAlign(23);
	frameMC->Draw();
	latex->DrawLatex(0.7,0.8,"#chi2/DOF="+chi2_s);
	latex->DrawLatex(0.7,0.7,"shift="+shift_s);
	latex->DrawLatex(0.7,0.6,"smear="+smear_s);
	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/KernelModel/"+year+"/"+fs+"_ScaleShift00.png");
	
	
	Double_t* results = new Double_t[2];
	results[0] = shift.getVal();
	results[1] = smear.getVal();
	


	return results;
	
	
}



void SetAddress(TTree* t){
	t->SetBranchAddress("GENmass2l",&GENmassZ);
	t->SetBranchAddress("massZ",&massZ);
	t->SetBranchAddress("massZErr",&massZErr);
	t->SetBranchAddress("weight",&weight);
	t->SetBranchAddress("pT1",&pT1);
	t->SetBranchAddress("pT2",&pT2);
	t->SetBranchAddress("eta1",&eta1);
	t->SetBranchAddress("eta2",&eta2);
	t->SetBranchAddress("m1",&m1);
	t->SetBranchAddress("m2",&m2);
	t->SetBranchAddress("phi1",&phi1);
	t->SetBranchAddress("phi2",&phi2);
	t->SetBranchAddress("pterr1",&pterr1);
	t->SetBranchAddress("pterr2",&pterr2);
}
