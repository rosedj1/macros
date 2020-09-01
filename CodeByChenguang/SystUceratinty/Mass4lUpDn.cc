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

void doplot(TString fs,TString year);



void Mass4lUpDn(){
	        
			TString year="2016";

			doplot("4e", year);//alpha1 alpha2 n1 n2 sigma
//			doplot("4mu", year);
//			doplot("2e2mu", year);
			
			
			// 0.76, 1.41, 7.32, 8.62, 2.0, 0.77, 1.41, 7.25, 9.89, 2.1, 0.76, 1.39, 7.47, 10.15, 2.08 ); 
			
}

void doplot(TString fs, TString year){
	
	
	
	
	
	
	
	Int_t fs1,fs2;
	if(fs=="4e"){fs1=2;fs2=2;}
	if(fs=="4mu"){fs1=1;fs2=1;}
	if(fs=="2e2mu"){fs1=3;fs2=4;}
	


	
	bool passedFullSelection;
	Int_t finalState;
	Float_t mass4lno_v;
	Float_t mass4lup_v;
	Float_t mass4ldn_v;
	Float_t mass4lErr;	
	Int_t sum;
	Int_t index=0;
	
	
	
	
	
	RooRealVar *mass4lno=new RooRealVar("mass4lno","mass4lno",105,140); 
	RooDataSet *dataset4lno=new RooDataSet("dataset4lno","dataset4lno",RooArgSet(*mass4lno));
	
	RooRealVar *mass4lup=new RooRealVar("mass4lup","mass4lup",105,140); 
	RooDataSet *dataset4lup=new RooDataSet("dataset4lup","dataset4lup",RooArgSet(*mass4lup));
	
	RooRealVar *mass4ldn=new RooRealVar("mass4ldn","mass4ldn",105,140); 
	RooDataSet *dataset4ldn=new RooDataSet("dataset4ldn","dataset4ldn",RooArgSet(*mass4ldn));
	
	TFile* f=new TFile();
	f = new TFile("/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_2016LUT_testHIG16041/Ntuples/ggH_HIG16041.root");
//	f = new TFile("/raid/raid9/chenguan/Mass/Z1MassConstraint/liteUFHZZ4LAnalyzer_"+year+"LUT/Ntuples/GluGluHToZZTo4L_M125_"+year+".root");
	TTree *t=(TTree*)f->Get("passedEvents");
	
	sum=t->GetEntries();
	cout<<sum<<endl;
	
	
	t->SetBranchAddress("finalState",&finalState);
	t->SetBranchAddress("passedFullSelection",&passedFullSelection);
	t->SetBranchAddress("mass4l",&mass4lno_v);
	t->SetBranchAddress("mass4l_up",&mass4lup_v);
	t->SetBranchAddress("mass4l_dn",&mass4ldn_v);
	t->SetBranchAddress("mass4lErr",&mass4lErr);
	
	
	for(Int_t j=0; j<sum; j++){
		t->GetEntry(j);
		
		if(passedFullSelection==1&&105<mass4lno_v&&mass4lno_v<140&&(finalState==fs1||finalState==fs2)){
			
			
			
			*mass4lno=mass4lno_v;
			dataset4lno->add(RooArgSet(*mass4lno));
			
			
			*mass4lup=mass4lup_v;
			dataset4lup->add(RooArgSet(*mass4lup));
			
			
			*mass4ldn=mass4ldn_v;
			dataset4ldn->add(RooArgSet(*mass4ldn));
			
			index=index+1;
		}
		
	}
	cout<<index<<endl;
	//make model
	
	RooRealVar meanDCB1("meanDCB","meanDCB",125,120,130);
	RooRealVar sigmaDCB1("sigmaDCB","sigmaDCB",1,0,10);
	RooRealVar alphaDCB1("alphaDCB","alphaDCB",1,0,10);
	RooRealVar nDCB1("nDCB","nDCB",1,0,10);
    RooRealVar alpha21("alpha2","alpha2",1,0,10);
	RooRealVar n21("n2","n2",1,0,10);
	RooDoubleCB DCBno("DCB","DCB",*mass4lno,meanDCB1,sigmaDCB1,alphaDCB1,nDCB1,alpha21,n21);
	
	
	
	
	
	DCBno.fitTo(*dataset4lno,PrintLevel(-1),SumW2Error(kTRUE),Timer(kTRUE));
    RooPlot *frame1=mass4lno->frame(Bins(100));
    frame1->SetTitle(fs+"_normal");/////////////////////////////////////////////////////////////////////
    dataset4lno->plotOn(frame1);
    DCBno.plotOn(frame1);
    DCBno.paramOn(frame1,Layout(0.1,0.4,0.9));
	
	Double_t chisquare1=frame1->chiSquare(6);
	TLatex *latex1=new TLatex();
    latex1->SetNDC();
    latex1->SetTextSize(0.05);
    latex1->SetTextFont(42);
    latex1->SetTextAlign(23);
    char chi21[20];
    sprintf(chi21,"%s%1.4f","#chi^{2}/DOF=",chisquare1);
   
	TCanvas* c = new TCanvas("c","c",1400,1000);
	c->cd();
	frame1->Draw();
	latex1->DrawLatex(0.7, 0.8, chi21);
	c->SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/"+fs+"_Normal.png");
	
//	RooRealVar shift_up("shift","shift",0,-1,1);
//	RooFormulaVar meanDCB2("meanDCB2","@1+@0",RooArgList(meanDCB,shift_up);	
	RooRealVar meanDCB2("meanDCB","meanDCB",125,120,130);
//	RooFormulaVar meanDCB_up("meanDCB","@1+@0",RooArgList(meanDCB2,shift_up)); 
	RooRealVar sigmaDCB2("sigmaDCB","sigmaDCB",1,0,10);//sigmaDCB1.getVal());
	RooRealVar alphaDCB2("alphaDCB","alphaDCB",1,0,10);//alphaDCB1.getVal());
	RooRealVar nDCB2("nDCB","nDCB",1,0,10);//nDCB1.getVal());
    RooRealVar alpha22("alpha2","alpha2",1,0,10);//alpha21.getVal());
	RooRealVar n22("n2","n2",1,0,10);//n21.getVal());	
	RooDoubleCB DCBup("DCB","DCB",*mass4lup,meanDCB2,sigmaDCB2,alphaDCB2,nDCB2,alpha22,n22);
	
	DCBup.fitTo(*dataset4lup,PrintLevel(-1),SumW2Error(kTRUE),Timer(kTRUE));
    RooPlot *frame2=mass4lup->frame(Bins(100));
    frame2->SetTitle(fs+"_up");/////////////////////////////////////////////////////////////////////
    dataset4lup->plotOn(frame2);
    DCBup.plotOn(frame2);
    DCBup.paramOn(frame2,Layout(0.1,0.4,0.9));
	
	Double_t chisquare2=frame2->chiSquare(6);
	TLatex *latex2=new TLatex();
    latex2->SetNDC();
    latex2->SetTextSize(0.05);
    latex2->SetTextFont(42);
    latex2->SetTextAlign(23);
    char chi22[20];
    sprintf(chi22,"%s%1.4f","#chi^{2}/DOF=",chisquare2);
	
	c->Clear();
	c->cd();
	frame2->Draw();
	latex2->DrawLatex(0.7,0.8,chi22);
	c->SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/"+fs+"_Up.png");

	
//	RooRealVar shift_dn("shift","shift",0,1,-1);
	RooRealVar meanDCB3("meanDCB","meanDCB",125,120,130);
//	RooFormulaVar meanDCB_dn("meanDCB","@1+@0",RooArgList(meanDCB3,shift_dn));
	RooRealVar sigmaDCB3("sigmaDCB","sigmaDCB",1,0,10);//sigmaDCB1.getVal());
	RooRealVar alphaDCB3("alphaDCB","alphaDCB",1,0,10);//alphaDCB1.getVal());
	RooRealVar nDCB3("nDCB","nDCB",1,0,10);//nDCB1.getVal());
    RooRealVar alpha23("alpha2","alpha2",1,0,10);//alpha21.getVal());
	RooRealVar n23("n2","n2",1,0,10);//n21.getVal());
	RooDoubleCB DCBdn("DCB","DCB",*mass4ldn,meanDCB3,sigmaDCB3,alphaDCB3,nDCB3,alpha23,n23);
	
	DCBdn.fitTo(*dataset4ldn,PrintLevel(-1),SumW2Error(kTRUE),Timer(kTRUE));
    RooPlot *frame3=mass4ldn->frame(Bins(100));
    frame3->SetTitle(fs+"_dn");/////////////////////////////////////////////////////////////////////
    dataset4ldn->plotOn(frame3);
    DCBdn.plotOn(frame3);
    DCBdn.paramOn(frame3,Layout(0.1,0.4,0.9));
	
	Double_t chisquare3=frame3->chiSquare(6);
	TLatex *latex3=new TLatex();
    latex3->SetNDC();
    latex3->SetTextSize(0.05);
    latex3->SetTextFont(42);
    latex3->SetTextAlign(23);
    char chi23[20];
    sprintf(chi23,"%s%1.4f","#chi^{2}/DOF=",chisquare3);
	
	c->Clear();
	c->cd();
	frame3->Draw();
	latex3->DrawLatex(0.7,0.8,chi23);
	c->SaveAs("/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/"+fs+"_Down.png");
	
	
	TString inforpath="/raid/raid9/chenguan/Mass/SystUceratinty/Uncertainty/"+year+"/"+fs+"/";
	
	ofstream fout1;
	fout1.open(inforpath+"normal_noabs.txt");
	fout1<<meanDCB1.getValV();
	
	ofstream fout2;
	fout2.open(inforpath+"up_noabs.txt");
	fout2<<meanDCB2.getValV();
	
	ofstream fout3;
	fout3.open(inforpath+"down_noabs.txt");
	fout3<<meanDCB3.getValV();
	
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!result!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
cout<<"                              "<<"Scale"<<endl;
cout<<"(normal-down)/125 = "<<((meanDCB1.getVal()-meanDCB3.getVal())/125)*100<<" %"<<endl;
cout<<" "<<endl;
cout<<"(up-normal)/125 = "<<((meanDCB2.getVal()-meanDCB1.getVal())/125)*100<<" %"<<endl;
cout<<" "<<endl;
//cout<<"                             "<<"Resolution"<<endl;
//cout<<"normal-down     "<<(sigmaDCB1.getVal()-sigmaDCB2.getVal())/sigmaDCB1.getVal()<<endl;
cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;




	delete mass4lno;
	delete dataset4lno;
	delete mass4lup;
	delete dataset4lup;
	delete mass4ldn;
	delete dataset4ldn;
	
	delete t;
	delete frame1;
	delete frame2;
	delete frame3;
	delete latex1;
	delete latex2;
	delete latex3;
	delete c;
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
