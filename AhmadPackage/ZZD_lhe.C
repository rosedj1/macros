/*
  Macro to translate the (POWHEG?) H->ZZ->4l LHE events input file into the .root file
  (with coresponding tree and histograms). It is tailored to  H->ZZ->4l events.
  Example syntax:
     root -l -q ZZ.C(<input file path/name>, <number of events to read>)
*/

// C++
#include <iostream>
#include <algorithm>
#include <vector>

// ROOT
#include "Riostream.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TString.h"

// HZZ4l angles
#include "computeAngles.cc"

using namespace std;

// Globals
const int LEPSIZE = 4;
const double PI = 3.1415926;
const double mZ_PDG = 91.188;//GeV

void ZZD_lhe(TString pohegLheFile = "ZZD_lhe_events.lhe", long int Nevents = 10000, bool bApplyCuts = false, unsigned int kDebugLevel = 1) {
// define input
    ifstream in;
    in.open(pohegLheFile.Data(),ios::in);

    // define output
    string strtemp;
    strtemp.assign(pohegLheFile.Data());
    strtemp.append("syntax2");
    strtemp.append(".root");
    std::cout<<" preparing root file: "<<strtemp.c_str()<<" \n";
    std::cout<<" running for: "<<Nevents<<" events.\n";

    string str;
    int IST, moth1, moth2, col1, col2;
    double aux1, aux2, aux3, aux4, aux5, wt;
    double p3_m, p4_m, p5_m, p6_m;
    double p3_0,p3_1,p3_2,p3_3,p4_0,p4_1,p4_2,p4_3;
    double p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3;
    long int nLHEevents = 0;
    long int nLHEevents4l = 0;
    long int nLHEevents4lout = 0;
    long int nLHEeventsNon4l = 0;
    TFile *f = new TFile(strtemp.c_str(),"RECREATE");

    // define observables
    int id1,id2,id3,id4;
    double eta3,eta4,eta5,eta6;
    double pt3,pt4,pt5,pt6;
    double pto1,pto2,pto3,pto4;
    double mass4l, pT4l, eta4l, rapidity4l;
    float _costheta1, _costheta2, _costhetastar, _phi, _phi1;
    double costheta1, costheta2, costhetastar, phi, phi1;
    double massZ1, massZ2;

    // define histograms
    TH1D *mass4l_hist_tchan = new TH1D("mass4l_tchan","mass4l_tchan",200,0.0,1000.0);
    TH1D *pT4l_hist_tchan   = new TH1D("pT4l_tchan","pT4l_tchan",200,0.0,200.0);
    TH1D *rapidity4l_hist_tchan = new TH1D("rapidity4l_tchan","rapidity4l_tchan",50,0.0,5.0);
    TH1D *eta4l_hist_tchan  = new TH1D("eta4l_tchan","eta4l_tchan",50,0.0,5.0);
    TH1D *massZ1_hist_tchan  = new TH1D("massZ1_tchan","massZ1_tchan",100,20.0,120.0);
    TH1D *massZ2_hist_tchan  = new TH1D("massZ2_tchan","massZ2_tchan",100,0.0,50.0);
    TH1D *pto1_hist_tchan   = new TH1D("pto1_tchan","pto1_tchan",100,0.0,500.0);
    TH1D *pto2_hist_tchan   = new TH1D("pto2_tchan","pto2_tchan",100,0.0,500.0);
    TH1D *pto3_hist_tchan   = new TH1D("pto3_tchan","pto3_tchan",100,0.0,500.0);
    TH1D *pto4_hist_tchan   = new TH1D("pto4_tchan","pto4_tchan",100,0.0,500.0);
    TH1D *costheta1_hist_tchan  = new TH1D("costheta1_tchan","costheta1_tchan",100,-1,+1);
    TH1D *costheta2_hist_tchan  = new TH1D("costheta2_tchan","costheta2_tchan",100,-1,+1);
    TH1D *costhetastar_hist_tchan  = new TH1D("costhetastar_tchan","costhetastar_tchan",100,-1,+1);
    TH1D *phi_hist_tchan  = new TH1D("phi_tchan","phi_tchan",100,-PI,PI);
    TH1D *phi1_hist_tchan  = new TH1D("phi1_tchan","phi1_tchan",100,-PI,PI);

    // define the root tree
    TTree* TT = new TTree("lheEvents_tchan","lheEvents_tchan");
    TT->Branch("weight",&wt,"weight/D");
    TT->Branch("id3",&id1,"id3/I");
    TT->Branch("id4",&id2,"id4/I");
    TT->Branch("id5",&id3,"id5/I");
    TT->Branch("id6",&id4,"id6/I");
    TT->Branch("eta3",&eta3,"eta3/D");
    TT->Branch("eta4",&eta4,"eta4/D");
    TT->Branch("eta5",&eta5,"eta5/D");
    TT->Branch("eta6",&eta6,"eta6/D");
    TT->Branch("pt3",&pt3,"pt3/D");
    TT->Branch("pt4",&pt4,"pt4/D");
    TT->Branch("pt5",&pt5,"pt5/D");
    TT->Branch("pt6",&pt6,"pt6/D");
    TT->Branch("pto1",&pto1,"pto1/D");
    TT->Branch("pto2",&pto2,"pto2/D");
    TT->Branch("pto3",&pto3,"pto3/D");
    TT->Branch("pto4",&pto4,"pto4/D");
    TT->Branch("mass4l",&mass4l,"mass4l/D");
    TT->Branch("pT4l",&pT4l,"pT4l/D");
    TT->Branch("eta4l",&eta4l,"eta4l/D");
    TT->Branch("rapidity4l",&rapidity4l,"rapidity4l/D");
    TT->Branch("massZ1",&massZ1,"massZ1/D");
    TT->Branch("massZ2",&massZ2,"massZ2/D");
    TT->Branch("costheta1",&costheta1,"costheta1/D");
    TT->Branch("costheta2",&costheta2,"costheta2/D");
    TT->Branch("costhetastar",&costhetastar,"costhetastar/D");
    TT->Branch("phi",&phi,"phi/D");
    TT->Branch("phi1",&phi1,"phi1/D");

    // Read the lines until the first <event> line
    while (1) {
        getline (in,str);if (!in.good()) break;
        if (str=="<event>") break;
    }

    // Read events assuming H->ZZ->4l process and fill histograms

    while (1) {
        nLHEevents++;
        cout<<"nLHEevents"<<nLHEevents<<endl;
        if(nLHEevents > Nevents) break;
        // get the event lines before leptons
        // first line for the weight
        //in >> aux1 >> aux2 >> wt >> aux3 >> aux4 >> aux5; if (!in.good()) break;        
        
        // skip 7 lines for ggH, 5 lines for qqZZ, 3 lines for ggZZ
        /*
        getline (in,str);if (!in.good()) break;
        getline (in,str);if (!in.good()) break;
        getline (in,str);if (!in.good()) break;
        getline (in,str);if (!in.good()) break;
        getline (in,str);if (!in.good()) break;*/
        // get the event lines for lepton (4 lines)
        while(1) {
            getline (in,str);
            in >> id1 >> IST >> moth1 >> moth2 >> col1 >> col2 >>p3_0>>p3_1>>p3_2>>p3_3>>p3_m>>aux1>>aux2; if (!in.good()) break;
            //cout<<id1<<"\t"<<IST<<"\t"<<moth1<<"\t"<<moth2<<"\t"<<col1<<"\t"<<col2<<"\t"<<p3_0<<"\t"<<p3_1<<"\t"<<p3_2<<"\t"<<p3_3<<"\t"<<p3_m<<"\t"<<aux1<<"\t"<<aux2<<endl;
            if (fabs(id1)>=11 && fabs(id1)<=15) break;
        }
        getline (in,str);
        in >> id2 >> IST >> moth1 >> moth2 >> col1 >> col2 >>p4_0>>p4_1>>p4_2>>p4_3>>p4_m>>aux1>>aux2; if (!in.good()) break;
        getline (in,str);
        in >> id3 >> IST >> moth1 >> moth2 >> col1 >> col2 >>p5_0>>p5_1>>p5_2>>p5_3>>p5_m>>aux1>>aux2; if (!in.good()) break;
        getline (in,str);
        in >> id4 >> IST >> moth1 >> moth2 >> col1 >> col2 >>p6_0>>p6_1>>p6_2>>p6_3>>p6_m>>aux1>>aux2;if (!in.good()) break; 

        //wt = 1; // weight = 1 at the momen
        // print some progress info
        if (kDebugLevel>0){
            if ((nLHEevents%100) == 0)printf("Analysed %lu events\n",nLHEevents);
        }

        // check if it is 4e, 4mu or 2e2mu event (move to the next one if it is not)
        int sumAbsId = fabs(id1) + fabs(id2) + fabs(id3) + fabs(id4);
        if ( (sumAbsId!=44 && sumAbsId!=48 && sumAbsId!=52) || (sumAbsId==52 && (fabs(id1)==15 || fabs(id2)==15 || fabs(id3)==15 || fabs(id4)==15)) ) {
            nLHEeventsNon4l++;
            if (kDebugLevel>1){
                printf("Event %lu is not a 4e, 4mu or 2e2mu event\n",nLHEevents);
            }
        } else {
            if (kDebugLevel>1){
                printf("Event %lu :\n %e, %e, %e, %e \n %e, %e, %e, %e \n %e, %e, %e, %e \n %e, %e, %e, %e \n",nLHEevents,p3_0,p3_1,p3_2,p3_3,p4_0,p4_1,p4_2,p4_3,p5_0,p5_1,p5_2,p5_3,p6_0,p6_1,p6_2,p6_3);
            }
            //std::cout<<"weight: "<<wt<<std::endl;
            // prepare TLorentzVector(s)
            int idL1,idL2,idL3,idL4;
            TLorentzVector L11P4,L12P4,L21P4,L22P4;
            // convention for leptons - example for 2mu2e event:    pp -> H -> ZZ -> mu-(p3)+mu+(p4)+e-(p5)+e+(p6)
            L11P4.SetPxPyPzE(p3_0, p3_1, p3_2, p3_3); idL1 = id1;
            L12P4.SetPxPyPzE(p4_0, p4_1, p4_2, p4_3); idL2 = id2;
            L21P4.SetPxPyPzE(p5_0, p5_1, p5_2, p5_3); idL3 = id3;
            L22P4.SetPxPyPzE(p6_0, p6_1, p6_2, p6_3); idL4 = id4;
            cout<<"pT,eta,phi,m L1: "<<L11P4.Pt()<<" "<<L11P4.Eta()<<" "<<L11P4.Phi()<<" "<<L11P4.M()<<endl;
            cout<<"pT,eta,phi,m L2: "<<L12P4.Pt()<<" "<<L12P4.Eta()<<" "<<L12P4.Phi()<<" "<<L12P4.M()<<endl;
            cout<<"pT,eta,phi,m L3: "<<L21P4.Pt()<<" "<<L21P4.Eta()<<" "<<L21P4.Phi()<<" "<<L21P4.M()<<endl;
            cout<<"pT,eta,phi,m L4: "<<L22P4.Pt()<<" "<<L22P4.Eta()<<" "<<L22P4.Phi()<<" "<<L22P4.M()<<endl;
            // prepare vectors of TLorentzVector(s

            // production-related jet-inclusive observables            
            TLorentzVector ZZ = (L11P4 + L12P4 + L21P4 + L22P4);
            mass4l     = ZZ.M();
            pT4l       = ZZ.Pt();
            rapidity4l = ZZ.Rapidity();
            eta4l      = ZZ.Eta();

            // form the Z1 and Z2 candidates
            TLorentzVector Z1, Z2;

            TLorentzVector tmpZ1 = (L11P4 + L12P4);
            TLorentzVector tmpZ2 = (L21P4 + L22P4);

            if (sumAbsId==48){// for 2e2mu
              Z1 = tmpZ1; Z2 = tmpZ2;
              mela::computeAngles(L11P4, idL1, L12P4, idL2, L21P4, idL3, L22P4, idL4, _costhetastar,_costheta1,_costheta2,_phi,_phi1);
              if (fabs(tmpZ1.M()-mZ_PDG)>fabs(tmpZ2.M()-mZ_PDG)){ //swap Z1 and Z2 if needed
                  Z2 = tmpZ1; Z1 = tmpZ2;
                  mela::computeAngles(L21P4, idL3, L22P4, idL4, L11P4, idL1, L12P4, idL2, _costhetastar,_costheta1,_costheta2,_phi,_phi1);
              }
            } else if (sumAbsId==44 || sumAbsId==52){// for 4e
              TLorentzVector tmpZ3 = (L11P4 + L22P4);
              TLorentzVector tmpZ4 = (L21P4 + L12P4);
              double min_dM = min(min(fabs(tmpZ1.M()-mZ_PDG),fabs(tmpZ2.M()-mZ_PDG)),min(fabs(tmpZ3.M()-mZ_PDG),fabs(tmpZ4.M()-mZ_PDG)));
              Z1 = tmpZ1; Z2 = tmpZ2;
              mela::computeAngles(L11P4, idL1, L12P4, idL2, L21P4, idL3, L22P4, idL4, _costhetastar,_costheta1,_costheta2,_phi,_phi1);
              if (fabs(tmpZ2.M()-mZ_PDG)<=min_dM){
                  Z2 = tmpZ1; Z1 = tmpZ2;
                  mela::computeAngles(L21P4, idL3, L22P4, idL4, L11P4, idL1, L12P4, idL2, _costhetastar,_costheta1,_costheta2,_phi,_phi1);
              }
              if (fabs(tmpZ3.M()-mZ_PDG)<=min_dM){
                  Z1 = tmpZ3; Z2 = tmpZ4;
                  mela::computeAngles(L11P4, idL1, L22P4, idL4, L21P4, idL3, L12P4, idL2, _costhetastar,_costheta1,_costheta2,_phi,_phi1);
              }
              if (fabs(tmpZ4.M()-mZ_PDG)<=min_dM){
                  Z2 = tmpZ3; Z1 = tmpZ4;
                  mela::computeAngles(L21P4, idL3, L12P4, idL2, L11P4, idL1, L22P4, idL4, _costhetastar,_costheta1,_costheta2,_phi,_phi1);
              }
            } else {
              cout << "WARNING: problem in the type of the final state!" << endl;
              return;
            }
            // decay-related jet-inclusive observables
            massZ1=Z1.M();
            massZ2=Z2.M();
            cout<<"massZ1"<<massZ1<<endl;
            cout<<"massZ2"<<massZ2<<endl;
            costheta1 = (double) _costheta1;
            costheta2 = (double) _costheta2;
            costhetastar = (double) _costhetastar;
            phi = (double) _phi;
            phi1 = (double) _phi1;

            // other observables
            pt3 = L11P4.Pt(); eta3 = L11P4.Eta();
            pt4 = L12P4.Pt(); eta4 = L12P4.Eta();
            pt5 = L21P4.Pt(); eta5 = L21P4.Eta();
            pt6 = L22P4.Pt(); eta6 = L22P4.Eta();

        
            // sort pts

            double ptArray[4];
            ptArray[0]=pt3; ptArray[1]=pt4; ptArray[2]=pt5; ptArray[3]=pt6;
            sort(ptArray, ptArray + LEPSIZE);
            pto1 = ptArray[3]; pto2 = ptArray[2]; pto3 = ptArray[1]; pto4 = ptArray[0];
            // check if event satifies the fiducial cuts
            if ( (!bApplyCuts) || ( (pto1>20.0 && pto2>10.0 && pto3>5.0 && pto4>5.0) 
                                    && (fabs(eta3)<2.5 && fabs(eta4)<2.5 && fabs(eta5)<2.5 && fabs(eta6)<2.5) 
                                    && (40.0<massZ1 && 12.0<massZ2)) ) {

                nLHEevents4l++;
                
                // Fill the histograms
                mass4l_hist_tchan->Fill(mass4l,wt);
                pT4l_hist_tchan->Fill(pT4l,wt);
                rapidity4l_hist_tchan->Fill(rapidity4l,wt);
                eta4l_hist_tchan->Fill(eta4l,wt);
                massZ1_hist_tchan->Fill(massZ1,wt);
                massZ2_hist_tchan->Fill(massZ2,wt);
                pto1_hist_tchan->Fill(pto1,wt);
                pto2_hist_tchan->Fill(pto2,wt);
                pto3_hist_tchan->Fill(pto3,wt);
                pto4_hist_tchan->Fill(pto4,wt);
                costheta1_hist_tchan->Fill(costheta1,wt);
                costheta2_hist_tchan->Fill(costheta2,wt);
                costhetastar_hist_tchan->Fill(costhetastar,wt);
                phi_hist_tchan->Fill(phi,wt);
                phi1_hist_tchan->Fill(phi1,wt);
                // Fill the tree
                TT->Fill();

                // print-outs
                //      cout << "[nLHEevents: " << nLHEevents << "]: mass4l: " << mass4l  << ", pT4l: " << pT4l  << ", rapidity4l: " << rapidity4l  << ", eta4l: " << eta4l << ", mZ1: " << mZ1 << ", mZ2: " << mZ2 << endl;
                //      cout << "pto1: " << pto1 << ", pto2: " << pto2 << "pto3: " << pto3 << ", pto4: " << pto4 << endl;
                //      cout << "costhetastar: " << costhetastar << ", costheta1: " << costheta1 << "costheta2: " << costheta2 << ", phi: " << phi << ", phi1: " << phi1 << endl << endl;
            } else {
                //nLHEevents4lout++;
            }

        }

        // move to the next event
        while (1) {
            getline (in,str); if (!in.good()) break;
            if (str=="<event>") break;
        }
    }

    // Normalise the histogram with the 1/binsize
    mass4l_hist_tchan->Scale(1./(mass4l_hist_tchan->GetBinWidth(1)));
    pT4l_hist_tchan->Scale(1./(pT4l_hist_tchan->GetBinWidth(1)));
    rapidity4l_hist_tchan->Scale(1./(rapidity4l_hist_tchan->GetBinWidth(1)));
    eta4l_hist_tchan->Scale(1./(eta4l_hist_tchan->GetBinWidth(1)));
    massZ1_hist_tchan->Scale(1./(massZ1_hist_tchan->GetBinWidth(1)));
    massZ2_hist_tchan->Scale(1./(massZ2_hist_tchan->GetBinWidth(1)));
    pto1_hist_tchan->Scale(1./(pto1_hist_tchan->GetBinWidth(1)));
    pto2_hist_tchan->Scale(1./(pto2_hist_tchan->GetBinWidth(1)));
    pto3_hist_tchan->Scale(1./(pto3_hist_tchan->GetBinWidth(1)));
    pto4_hist_tchan->Scale(1./(pto4_hist_tchan->GetBinWidth(1)));
    costheta1_hist_tchan->Scale(1./(costheta1_hist_tchan->GetBinWidth(1)));
    costheta2_hist_tchan->Scale(1./(costheta2_hist_tchan->GetBinWidth(1)));
    costhetastar_hist_tchan->Scale(1./(costhetastar_hist_tchan->GetBinWidth(1)));
    phi_hist_tchan->Scale(1./(phi_hist_tchan->GetBinWidth(1)));
    phi1_hist_tchan->Scale(1./(phi1_hist_tchan->GetBinWidth(1)));

    // Write tree and histograms, close & delete
    mass4l_hist_tchan->Write();
    pT4l_hist_tchan->Write();
    rapidity4l_hist_tchan->Write();
    eta4l_hist_tchan->Write();
    massZ1_hist_tchan->Write();
    massZ2_hist_tchan->Write();
    costheta1_hist_tchan->Write();
    costheta2_hist_tchan->Write();
    costhetastar_hist_tchan->Write();
    phi_hist_tchan->Write();
    phi1_hist_tchan->Write();
    pto1_hist_tchan->Write();
    pto2_hist_tchan->Write();
    pto3_hist_tchan->Write();
    pto4_hist_tchan->Write();

    TT->Write();

    printf("Found %lu events\n",--nLHEevents);
    if (kDebugLevel>0){
        printf("Found %lu non-4l events\n",nLHEeventsNon4l);
        printf("Found %lu 4l events out-of-fiducial cuts\n",nLHEevents4lout);
        printf("Found %lu 4l events within selection cuts\n",nLHEevents4l);
    }
    in.close();
    delete f;
}
