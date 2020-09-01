void massZvsXXX(){
	
	TString year = "2017";
	bool isData = 0;
	


	
	Double_t massZ;
	Double_t massZErr;
	Double_t eta1;
	Double_t eta2;
	Double_t pt1;
	Double_t pt2;
	
	TString plotpath = "/home/chenguan/public_html/Syst_uncertainty/LepEnergyScale/"+year+"/2Dmap/";
	TString plotname = "";
	if(isData==0) plotname = plotpath+"MC_";
	if(isData==1) plotname = plotpath+"Data_";
	
	TH2F* h_eta1vspt1 = new TH2F("eta1vspt1","eta1vspt1",100,-2.5,2.5,100,5,100);
	TH2F* h_eta2vspt2 = new TH2F("eta2vspt2","eta2vspt2",100,-2.5,2.5,100,5,100);
	TH2F* h_massZvseta1 = new TH2F("massZvseta1","massZvseta1",100,60,120,100,-2.5,2.5);
	TH2F* h_massZvseta2 = new TH2F("massZvseta2","massZvseta2",100,60,120,100,-2.5,2.5);
	TH2F* h_massZvspt1 = new TH2F("massZvspt1","massZvspt1",100,60,120,100,5,100);
	TH2F* h_massZvspt2 = new TH2F("massZvspt2","massZvspt2",100,60,120,100,5,100);
	
	h_eta1vspt1->SetStats(0);
	h_eta2vspt2->SetStats(0);
	h_massZvseta1->SetStats(0);
	h_massZvseta2->SetStats(0);
	h_massZvspt1->SetStats(0);
	h_massZvspt2->SetStats(0);
	
	TFile* f = new TFile();
	if(year=="2016"&&isData==0)f = new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/DY_2016MC_v1_20170222/DYJetsToLL_M-50_kalman_v4_m2mu.root");
	if(year=="2016"&&isData==1)f = new TFile("/raid/raid9/mhl/HZZ4L_Run2_post2016ICHEP/outputRoot/Data_2016_v1_20170223/DoubleLepton_m2mu.root");
	if(year=="2017"&&isData==0)f = new TFile("/raid/raid9/chenguan/input/DY_"+year+"/DYJetsToLL_M50_m2mu.root");
	if(year=="2017"&&isData==1)f = new TFile("/raid/raid9/chenguan/input/2lSkim_"+year+"/DoubleMuon_Run2017-17Nov2017_m2mu.root");
	
	TTree* t = (TTree*)f->Get("passedEvents");
	
	t->SetBranchAddress("massZ",&massZ);
	t->SetBranchAddress("massZErr",&massZErr);
	t->SetBranchAddress("eta1",&eta1);
	t->SetBranchAddress("eta2",&eta2);
	t->SetBranchAddress("pT1",&pt1);
	t->SetBranchAddress("pT2",&pt2);
	
	Int_t counter=0;

	Int_t sum = (Int_t)t->GetEntries();
	for(Int_t i=0; i<sum; i++){
		t->GetEntry(i);	
		if(0.2<massZErr&&massZErr<7.2&&60<massZ&&massZ<120){
			
			h_eta1vspt1->Fill(eta1,pt2);
			h_eta2vspt2->Fill(eta2,pt2);
			h_massZvseta1->Fill(massZ,eta1);
			h_massZvseta2->Fill(massZ,eta2);
			h_massZvspt1->Fill(massZ,pt1);
			h_massZvspt2->Fill(massZ,pt2);
			counter=counter+1;

		}
	}
//	cout<<sum<<endl;
	
	TCanvas* c = new TCanvas("c","c",1000,1000);
	c->cd();
	h_eta1vspt1->Draw("COLZ");
	c->SaveAs(plotname+"eta1vspt1.png");
	c->Clear();
	c->cd();
	h_eta2vspt2->Draw("COLZ");
	c->SaveAs(plotname+"eta2vspt2.png");
	c->Clear();
	c->cd();
	h_massZvseta1->Draw("COLZ");
	c->SaveAs(plotname+"massZvseta1.png");
	c->Clear();
	c->cd();
	h_massZvseta2->Draw("COLZ");
	c->SaveAs(plotname+"massZvseta2.png");
	c->Clear();
	c->cd();
	h_massZvspt1->Draw("COLZ");
	c->SaveAs(plotname+"massZvspt1.png");
	c->Clear();
	c->cd();
	h_massZvspt2->Draw("COLZ");
	c->SaveAs(plotname+"massZvspt2.png");
	
	
	
}
