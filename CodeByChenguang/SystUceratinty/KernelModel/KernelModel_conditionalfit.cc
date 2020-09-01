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
RooDataHist* GetDataSet(TString year, TString fs, bool isData, Double_t a1, Double_t a2, Double_t b1, Double_t b2);
Double_t* GetShift(TString year, TString fs, Double_t a1, Double_t a2, Double_t b1, Double_t b2, RooDataHist* MC, RooDataHist* Data, bool isData);

Double_t GENmassZ, massZ, massZErr, weight, pT1, pT2, eta1, eta2, m1, m2, phi1, phi2, pterr1, pterr2;

void KernelModel_conditionalfit(){
	
	TString year = "2017";
	TString fs = "2mu";
	
	RooDataHist* MC = GetDataSet(year, fs, 0, 0.9, 1.4, 5, 20);
	RooDataHist* Data = GetDataSet(year, fs, 1, 0.9, 1.4, 5, 20);

	MC->Print();
	Data->Print();
	Double_t* results = new Double_t[2];
	results = GetShift(year, fs, 0.9, 1.4, 5, 20, MC, Data,0);
	cout<<results[0]<<endl;
	cout<<results[1]<<endl;	

	
	
	
	
}







RooDataHist* GetDataSet(TString year, TString fs, bool isData, Double_t a1, Double_t a2, Double_t b1, Double_t b2){
	
	TFile* f = new TFile();
	if(year=="2016"&&fs=="2e"&&isData==0)f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2016_m2e.root");
	if(year=="2017"&&fs=="2e"&&isData==0)f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2017_m2e.root");
	if(year=="2018"&&fs=="2e"&&isData==0)f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2018_m2e.root");
	if(year=="2016"&&fs=="2mu"&&isData==0)f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2016_m2mu.root");
	if(year=="2017"&&fs=="2mu"&&isData==0)f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2017_m2mu.root");
	if(year=="2018"&&fs=="2mu"&&isData==0)f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/2018_m2mu.root");
	
	if(year=="2016"&&fs=="2e"&&isData==1)f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2016/Electron/SingleElectron_m2e.root");
	if(year=="2017"&&fs=="2e"&&isData==1)f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2016/Electron/SingleElectron_m2e.root");
	if(year=="2018"&&fs=="2e"&&isData==1)f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2017/Electron/EGamma_m2e.root");
	if(year=="2016"&&fs=="2mu"&&isData==1)f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2017/Muon/SingleMuon_m2mu.root");
	if(year=="2017"&&fs=="2mu"&&isData==1)f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2018/Muon/SingleMuon_m2mu.root");
	if(year=="2018"&&fs=="2mu"&&isData==1)f = new TFile("/raid/raid9/chenguan/input/Legacy_102X_DY_2018/FromFillipo/Data/2018/Muon/SingleMuon_m2mu.root");
	
	RooRealVar massZ_d("massZ","massZ",80,100);
	RooRealVar GENmassZ_d("GENmassZ","GENmassZ",0.2,7.2);
	RooRealVar weight_d("weight","weight",0.00001,100);
	RooArgSet argset(massZ_d,GENmassZ_d,weight_d);
	RooDataSet* dataset = new RooDataSet("dataset","dataset",argset);
	massZ_d.setBins(1000,"cache");

	TH1D* h = new TH1D("h","h",400,80,100);

	TTree* t = (TTree*)f->Get("passedEvents");
	SetAddress(t);
	Int_t sum = t->GetEntries();
	Int_t counter = 0;
	Double_t random = 0;
	TRandom3 rand;
	rand.SetSeed(123457);
	for(Int_t i=0; i<sum; i++){
		t->GetEntry(i);
		random = rand.Gaus(0,1);
		if(massZErr<7.2&&massZErr>0.2&&massZ<100&&massZ>80&&GENmassZ<100&&GENmassZ>80){
			if((abs(eta1)>a1&&abs(eta1)<a2&&pT1>b1&&pT1<b2)||(abs(eta2)>a1&&abs(eta2)<a2&&pT2>b1&&pT2<b2)){
				
//				h->Fill(massZ);
				massZ_d.setVal(massZ);
				GENmassZ_d.setVal(GENmassZ);
				weight_d.setVal(weight);
				dataset->add(argset);
				counter = counter+1;
				if(counter==1000)break;
//				cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<massZ<<endl;
			}
		}
	}
	RooDataSet *dataset_w=new RooDataSet(dataset->GetName(),dataset->GetTitle(),dataset,*dataset->get(),"1","weight");
	RooDataHist* dataset_binned = dataset_w->binnedClone();
	
	
	return dataset_binned;
	
				
	
	
	

}


Double_t* GetShift(TString year, TString fs, Double_t a1, Double_t a2, Double_t b1, Double_t b2, RooDataHist* MC, RooDataHist* Data, bool isData){
	
	TString name;
	if(isData==false&&year=="2017")name="2017_MC";
	if(isData==true&&year=="2017")name="2017_Data";
	if(isData==false&&year=="2016")name="2016_MC";
	if(isData==true&&year=="2016")name="2016_Data";
	if(isData==false&&year=="2018")name="2018_MC";
	if(isData==true&&year=="2018")name="2018_Data";
	
	TFile* f_e1 = new TFile("LUT/"+name+"_e1.root");
	TFile* f_e2 = new TFile("LUT/"+name+"_e2.root");
	TFile* f_e3 = new TFile("LUT/"+name+"_e3.root");
	TFile* f_mu = new TFile("LUT/"+name+"_mu.root");
	
	TH2D* mu_corr = (TH2D*)f_mu->Get("mu");
	TH2D* e1_corr = (TH2D*)f_e1->Get("e1");
	TH2D* e2_corr = (TH2D*)f_e2->Get("e2");
	TH2D* e3_corr = (TH2D*)f_e3->Get("e3");
	
	TAxis* x_mu = mu_corr->GetXaxis();
	TAxis* x_e1 = e1_corr->GetXaxis();
	TAxis* x_e2 = e2_corr->GetXaxis();
	TAxis* x_e3 = e3_corr->GetXaxis();
	TAxis* y_mu = mu_corr->GetYaxis();
	TAxis* y_e1 = e1_corr->GetYaxis();
	TAxis* y_e2 = e2_corr->GetYaxis();
	TAxis* y_e3 = e3_corr->GetYaxis();
	
	
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

	Double_t random = 0;
	TRandom3 rand;
	rand.SetSeed(123457);
	
	TString index = "";
	Int_t counter = 0;

	RooArgList FuncList("FuncList");
	RooArgList CoefList("CoefList");

	RooRealVar* massZ_d = new RooRealVar("massZ","massZ",80,100);
	RooRealVar* sigma = new RooRealVar("sigma","sigma",0.1,0,5);
	RooRealVar* mu = new RooRealVar("mu","mu",1,0.8,1.2);
	RooRealVar* GENmassZ_d = new RooRealVar("GENmassZ","GENmassZ",80,100);
	RooFormulaVar* mean = new RooFormulaVar("mean"+index,"@0*@1",RooArgSet(*mu,*GENmassZ_d));

//	RooRealVar* mean = new RooRealVar("mean","mean",91,80,100);
//	mu->setVal(1);
//	mu->setConstant(kTRUE);	
//	RooRealVar* PDGZm = new RooRealVar("PDGZm","PDGZm",91.19);
//	RooRealVar* PDGZw = new RooRealVar("PDGZw","PDGZw",2.44);
//	RooBreitWigner* PDGZ = new RooBreitWigner("PDGZ","PDGZ",*massZ_d,*PDGZm,*PDGZw);

	for(Int_t i=0; i<sum; i++){
		t->GetEntry(i);
		if(i%1000000==0)cout<<i<<endl;
		random = rand.Gaus(0,1);
		if(massZ<100&&massZ>80&&massZErr<7.2&&massZErr>0.2&&GENmassZ>80&&GENmassZ<100){
			if((abs(eta1)>a1&&abs(eta1)<a2&&b1<pT1&&pT1<b2)||(abs(eta2)>a1&&abs(eta2)<a2&&pT2<b2&&pT2>b1)){
				index = to_string(counter);	
				counter = counter+1;
				Double_t lambda1, lambda2;
				//updata massZErr
				if(fs=="2mu"){
				Int_t xbin1 = x_mu->FindBin(pT1);
				Int_t ybin1 = y_mu->FindBin(abs(eta1));
				Int_t xbin2 = x_mu->FindBin(pT2);
				Int_t ybin2 = y_mu->FindBin(abs(eta2));
				lambda1 = mu_corr->GetBinContent(xbin1,ybin1);
				lambda2 = mu_corr->GetBinContent(xbin2,ybin2);
				}
				Double_t pt1err_corr = pterr1*lambda1;
				Double_t pt2err_corr = pterr2*lambda2;
				
				TLorentzVector lep1, lep2;
				lep1.SetPtEtaPhiM(pT1,eta1,phi1,m1);
				lep2.SetPtEtaPhiM(pT2,eta2,phi2,m2);
			
				TLorentzVector lep1p, lep2p;
				lep1p.SetPtEtaPhiM(pT1+pt1err_corr,eta1,phi1,m1);
				lep2p.SetPtEtaPhiM(pT2+pt2err_corr,eta2,phi2,m2);
				
				Double_t dm1corr = (lep1p+lep2).M()-(lep1+lep2).M();
				Double_t dm2corr = (lep1+lep2p).M()-(lep1+lep2).M();
				Double_t newmassZErr = TMath::Sqrt(dm1corr*dm1corr+dm2corr*dm2corr);
				
//				RooRealVar* sigma = new RooRealVar("sigma"+index,"sigma"+index,newmassZErr);
//				RooRealVar* GENmassZ_d = new RooRealVar("GENmassZ_d"+index,"GENmassZ_d"+index,GENmassZ);
//				RooFormulaVar* mean = new RooFormulaVar("mean"+index,"@0*@1",RooArgSet(*mu,*GENmassZ_d));
//				RooRealVar* mean = new RooRealVar("mean"+index,"mean"+index,GENmassZ,GENmassZ-5,GENmassZ+5);
				RooGaussian* kernel = new RooGaussian("kernel"+index,"kernel"+index,*massZ_d,*mean,*sigma);
//				RooFFTConvPdf* kernel = new RooFFTConvPdf("kernel"+index,"kernel"+index,*massZ_d,*PDGZ,*gauss);
				RooRealVar* coef = new RooRealVar("coef"+index,"coef"+index,0.001);

				FuncList.add(*kernel);
				if(counter<1000)CoefList.add(*coef);
				
				if(counter==1000)break;				
//			cout<<"!!!!!!!!!!!!!!!!!!!!!"<<massZ<<endl;	
				
			}	
		}
	}
	cout<<counter<<"kernels for model building"<<endl;	
	RooRealSumPdf SumKernel("SumKernel","SumKernel",FuncList,CoefList);
	cout<<"MODEL BUILDED"<<endl;		
	SumKernel.fitTo(*MC,ConditionalObservables(*GENmassZ_d),SumW2Error(kTRUE),PrintLevel(-1),Timer(kTRUE));
	RooPlot* frameMC = massZ_d->frame(Bins(100));
	frameMC->SetTitle("");
	MC->plotOn(frameMC);
	SumKernel.plotOn(frameMC,ProjWData(*GENmassZ_d,*MC),LineColor(2),LineWidth(3));
	TCanvas c("c","c",4000,1000);
	c.cd();
	Double_t chi2 = frameMC->chiSquare(2);
	TString mu_s, sigma_s, chi2_s;
	mu_s = to_string(mu->getVal());
//	sigma_s = to_string (sigma->getVal());
	chi2_s = to_string(chi2);
	TLatex *latex=new TLatex();
        latex->SetNDC();
    	latex->SetTextSize(0.05);
    	latex->SetTextFont(42);
    	latex->SetTextAlign(23);
	frameMC->Draw();
	latex->DrawLatex(0.7,0.8,"#chi2/DOF="+chi2_s);
	latex->DrawLatex(0.7,0.7,"#mu="+mu_s);
//	latex->DrawLatex(0.7,0.6,"#sigma="+sigma_s);
	TString tag = "MC";
	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/KernelModel/"+year+"/"+tag+"_"+fs+"_ScaleShift5.png");
	Double_t result1 = mu->getVal();
//	Double_t result2 = sigma->getVal();

/*	
	SumKernel.fitTo(*Data,SumW2Error(kTRUE),PrintLevel(-1),Timer(kTRUE));
	RooPlot* frameData = massZ_d.frame(Bins(100));
	frameData->SetTitle("");
	Data->plotOn(frameData);
	SumKernel.plotOn(frameData,LineColor(2),LineWidth(1));
	chi2 = frameData->chiSquare(2);
	mu_s = to_string(mu.getVal());
	sigma_s = to_string(sigma.getVal());
	chi2_s = to_string(chi2);
	c.Clear();
	c.cd();
	frameData->Draw();
	latex->DrawLatex(0.7,0.8,"#chi2/DOF="+chi2_s);
	latex->DrawLatex(0.7,0.7,"#mu="+mu_s);
	latex->DrawLatex(0.7,0.6,"#sigma="+sigma_s);
	tag = "Data";
	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/KernelModel/"+year+"/"+tag+"_"+fs+"_ScaleShift.png");
	Double_t result2 = mu.getVal();
*/	
	Double_t* results = new Double_t[2];
	results[0] = result1;
//	results[1] = result2;
	

/*	
	RooPlot* frame = massZ_d.frame();
	frame->SetTitle("");
	SumKernel.plotOn(frame);
	
	TCanvas c("c","c",1000,1000);
	c.cd();
	frame->Draw();
	c.SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/KernelModel/test.png");
*/
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
