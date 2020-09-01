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

void OpenCheck(){

    TFile *infile = TFile::Open("root://cmsio5.rc.ufl.edu//store/user/t2/users/ferrico/Full_RunII/Data/2018/Electron/EGamma_RunD.root");
  TTree* tree; 
    if(infile){ 
            std::cout<<"File trovato.OK"<<std::endl;
                infile->cd("Ana");
                    tree = (TTree*)gDirectory->Get("passedEvents");
                          }
}
