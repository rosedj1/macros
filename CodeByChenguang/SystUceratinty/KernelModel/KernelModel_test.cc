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

Double_t GENmassZ, massZ, massZErr, pT1, pT2, eta1, eta2;

void KernelModel_test(){
	
	TString year="2017";
	TString fs = "2mu";
	Double_t a1, a2, b1, b2;
	a1 = 0;
	a2 = 0.8;
	b1 = 5;
	b2 = 20;
	
	
	TFile* f = new TFile();
	if(year=="2016"&&fs=="2e")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2016_m2e.root");
	if(year=="2017"&&fs=="2e")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2017_m2e.root");
	if(year=="2018"&&fs=="2e")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2018_m2e.root");
	if(year=="2016"&&fs=="2mu")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2016_m2mu.root");
	if(year=="2017"&&fs=="2mu")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2017_m2mu.root");
	if(year=="2018"&&fs=="2mu")f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2018_m2mu.root");
	TTree* t = (TTree*)f->Get("passedEvents");
	SetAddress(t);
	
	TString index = "";
	
	RooArgList FuncList("FuncSet");
	RooArgList CoefList("CoefSet");
	Double_t random = 0;
	TRandom3 rand;
	rand.SetSeed(123457);
	Int_t sum = t->GetEntries();
	cout<<sum<<endl;
	Int_t counter = 0;
	
	RooRealVar* massZ_d = new RooRealVar("massZ","massZ",80,100);
//	RooRealVar* sigma = new RooRealVar("sigma","sigma",0.5);
//	RooRealVar*  mu = new RooRealVar("mu","mu",1);
	
	for(Int_t i=0; i<7; i++){
//		t->GetEntry(i);
//		if(i%1000000==0)cout<<i<<endl;
//		random = rand.Gaus(0,1);
//		if(massZ<100&&massZ>80&&massZErr<7.2&&massZErr>0.2&&GENmassZ<100&&GENmassZ>80){
//			if((abs(eta1)>a1&&abs(eta1)<a2&&b1<pT1&&pT1<b2)||(abs(eta2)>a1&&abs(eta2)<a2&&pT2<b2&&pT2>b1)){
				index = to_string(i);
//				RooRealVar* massZ_d = new RooRealVar("massZ","massZ",60,120);
//				RooRealVar* sigma = new RooRealVar("sigma","sigma",1,0,10);
//				RooRealVar*  mu = new RooRealVar("mu","mu",1,-2,2);
//				RooRealVar* GENmassZ_d = new RooRealVar("GENmassZ_d"+index,"GENmassZ_d"+index,90+5*i);

				RooRealVar* mean = new RooRealVar("mean"+index,"mean"+index,90+3*i);
				RooRealVar* sigma = new RooRealVar("sigma"+index,"sigma"+index,1);
				RooGaussian* kernel = new RooGaussian("kernel"+index,"kernel"+index,*massZ_d,*mean,*sigma);
				RooRealVar* coef = new RooRealVar("coef"+index,"coef"+index,1);
				
				counter = counter+1;	
				FuncList.add(*kernel);
				CoefList.add(*coef);
		
				if(counter==6)break;
cout<<"!!!!!!!!!!!!!!!!!!!"<<massZ<<endl;
//			}	
//		}
	}
	
	cout<<counter<<"kernels"<<endl;
	RooRealSumPdf SumKernel("SumKernel","SumKernel",FuncList,CoefList);
	cout<<"MODEL BUILDED!!!!!!!!"<<endl;
		
//	sigma->setVal(2);

	RooPlot* frame = massZ_d->frame();
	frame->SetTitle("");
	SumKernel.plotOn(frame);
	TCanvas c("c","c",1000,1000);
	c.cd();
	frame->Draw();
	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/KernelModel/test_1000kernel.png");


	
}



void SetAddress(TTree* t){
	t->SetBranchAddress("GENmass2l",&GENmassZ);
	t->SetBranchAddress("massZ",&massZ);
	t->SetBranchAddress("massZErr",&massZErr);
	t->SetBranchAddress("pT1",&pT1);
	t->SetBranchAddress("pT2",&pT2);
	t->SetBranchAddress("eta1",&eta1);
	t->SetBranchAddress("eta2",&eta2);
}
