#ifndef ANALYZETPROXYTBSM_H
#define ANALYZETPROXYTBSM_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "NtupleVarsTProxy.h"
#include "TH1F.h"
#include "TTree.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

//#pragma link C++ class std::vector< std::vector >+; 
//#pragma link C++ class std::vector< TLorentzVector >+;
//#pragma link C++ class NtupleVarsTProxy+;

//Define Root LorentzVector
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > myLV;

class AnalyzeTProxytBSM : public NtupleVarsTProxy{

 public:
  AnalyzeTProxytBSM(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data", const char *sample="sample", const char* LostlepFlag ="Flag", const char* phoID="phoID");
  ~AnalyzeTProxytBSM();
  //void     EventLoop(const char *,const char *,const char *,const char *, const char*, const char*);
  void   EventLoop(std::string buffer);
  void   BookHistogram(const char *);//, const char *);
  Bool_t Process(Long64_t entry);
  myLV   getBestPhoton(int);
  //TLorentzVector  getBestPhoton(int);
  int    bestPhotonIndxAmongPhotons=-100;

  TFile *oFile;
  TH1F *h_MET;
  TH1F *h_MET2;
  TH1F *h_EvtBrk;
  TH1F *h_hadJets_Pt;
  TH1F *h_hadJets_Pt1;
  TH1F *h_hadJets_Eta;
  TH1F *h_hadJets_Eta1;
  TH1F *h_NJets0;   //jets with no cuts
  TH1F *h_NJets1;   // jets with cut on Pt
  TH1F *h_NJets2;   // along with eta cut
  TH1F *h_NJets3;   // along with MET cut
  TH1F *h_Jet_pT0, *h_Jet_Eta0;   //0 stands for no cut 
  TH1F *h_Jet_pT1, *h_Jet_Eta1;   //1 stands for cut in the pT 
  TH1F *h_Jet_pT2, *h_Jet_Eta2;   //2 stands for cut in pT and Eta
  TH1F *h_Jet_pT3, *h_Jet_Eta3;   //3 stands for cut in pT and Eta and MET
  
  TH2F *h_NJet_MET, *h_NJet_PhoPt, *h_NJet_genPhoPt; 

  TH1F *h_ele_pT0, *h_ele_eta0, *h_ele_phi0; 
  TH1F *h_ele_pT3, *h_ele_eta3, *h_ele_phi3;  
  TH1F *h_pho_pT0, *h_pho_eta0, *h_pho_phi0;
  TH1F *h_pho_pT3, *h_pho_eta3, *h_pho_phi3;
  TH1F *h_gen_pT0, *h_gen_eta0, *h_gen_phi0;
  TH1F *h_gen_pT3, *h_gen_eta3, *h_gen_phi3;
};
#endif


#ifdef ANALYZETPROXYTBSM_cxx

//void AnalyzeLightBSM::BookHistogram(const char *outFileName, const char *N2_mass) {
void AnalyzeTProxytBSM::BookHistogram(const char *outFileName) {
  std::cout << "AnalyzeLightBSM::BookHistogram " << std::endl;

  oFile = new TFile(outFileName, "recreate");
  oFile->cd();

  // Book your histograms & summary counters here 
  h_MET = new TH1F("h_MET", "h_MET", 100, 0.0, 1000.0);
  //h_MET2 = new TH1F("h_MET2", "hMET2", 100, 0.0, 1000.0);
  h_hadJets_Pt = new TH1F("h_hadJets_Pt", "h_hadJets_Pt", 100, 0.0, 1000.0);
  h_hadJets_Pt1 = new TH1F("h_hadJets_Pt1", "h_hadJets_Pt1", 100, 0.0, 1000.0);
  h_hadJets_Eta = new TH1F("h_hadJets_Eta", "h_hadJets_Eta", 100, 0.0, 10.0);
  h_hadJets_Eta1 = new TH1F("h_hadJets_Eta1", "h_hadJets_Eta1", 100, 0.0, 10.0);
  
  h_Jet_pT0  = new TH1F("h_Jet_pT0",  "h_Jet_pT0", 100, 0.0, 1000.0);
  h_Jet_pT1  = new TH1F("h_Jet_pT1",  "h_Jet_pT1", 100, 0.0, 1000.0);
  h_Jet_pT2  = new TH1F("h_Jet_pT2",  "h_Jet_pT2", 100, 0.0, 1000.0);
  h_Jet_pT3  = new TH1F("h_Jet_pT3",  "h_Jet_pT3", 100, 0.0, 1000.0);
  h_Jet_Eta0 = new TH1F("h_Jet_Eta0", "h_Jet_Eta0", 100, -10.0, 10.0);
  h_Jet_Eta1 = new TH1F("h_Jet_Eta1", "h_Jet_Eta1", 100, -10.0, 10.0);
  h_Jet_Eta2 = new TH1F("h_Jet_Eta2", "h_Jet_Eta2", 100, -10.0, 10.0);
  h_Jet_Eta3 = new TH1F("h_Jet_Eta3", "h_Jet_Eta3", 100, -10.0, 10.0);

  h_NJets0 = new TH1F("h_NJets0", "h_NJets0",50, 0, 50);
  h_NJets1 = new TH1F("h_NJets1", "h_NJets1",50, 0, 50);
  h_NJets2 = new TH1F("h_NJets2", "h_NJets2",50, 0, 50);
  h_NJets3 = new TH1F("h_NJets3", "h_NJets3",50, 0, 50);
   
  h_ele_pT0  = new TH1F("h_ele_pT0",  "h_ele_pT0", 100, 0.0, 1000.0);
  h_ele_eta0 = new TH1F("h_ele_eta0", "h_ele_eta0",100, -3.0, 3.0);
  h_ele_phi0 = new TH1F("h_ele_phi0", "h_ele_phi0",100, -3.2, 3.2);
  h_ele_pT3  = new TH1F("h_ele_pT3",  "h_ele_pT3", 100, 0.0, 1000.0);
  h_ele_eta3 = new TH1F("h_ele_et3", "h_ele_eta3",100, -3.0, 3.0);
  h_ele_phi3 = new TH1F("h_ele_phi3", "h_ele_phi3",100, -3.2, 3.2);

  h_pho_pT0  = new TH1F("h_pho_pT0",  "h_pho_pT0", 100, 0.0, 1000.0);
  //h_pho_pT1  = new TH1F("h_pho_pT1",  "h_pho_pT1", 100, 0.0, 1000.0);
  h_pho_eta0 = new TH1F("h_pho_eta0", "h_pho_eta0",100, -3.0, 3.0);
  h_pho_phi0 = new TH1F("h_pho_phi0", "h_pho_phi0",100, -3.2, 3.2);
  h_pho_pT3 = new TH1F("h_pho_pT3",  "h_pho_pT3", 100, 0.0, 1000.0);
  //h_pho_pT1  = new TH1F("h_pho_pT1",  "h_pho_pT1", 100, 0.0, 1000.0);
  h_pho_eta3 = new TH1F("h_pho_eta3", "h_pho_eta3",100, -3.0, 3.0);
  h_pho_phi3= new TH1F("h_pho_phi3", "h_pho_phi3",100, -3.2, 3.2);

  h_NJet_MET = new TH2F("h_NJet_MET", "h_NJet_MET", 100, 0.0, 1000.0, 50, 0, 50);
  h_NJet_MET->GetXaxis()->SetTitle("MET");
  h_NJet_MET->GetYaxis()->SetTitle("NJet");

  h_NJet_PhoPt = new TH2F("h_NJet_PhoPt", "h_PhoPt", 100, 0.0, 1000.0, 50, 0, 50);
  h_NJet_PhoPt->GetXaxis()->SetTitle("Pho_Pt");
  h_NJet_PhoPt->GetYaxis()->SetTitle("NJet");
  
  h_NJet_genPhoPt = new TH2F("h_genNJet_PhoPt", "h_genPhoPt", 100, 0.0, 1000.0, 50, 0, 50);
  h_NJet_genPhoPt->GetXaxis()->SetTitle("genPho_Pt");
  h_NJet_genPhoPt->GetYaxis()->SetTitle("NJet");
  
  h_gen_pT0  = new TH1F("h_gen_pT0",  "h_gen_pT0", 100, 0.0, 1000.0);
  h_gen_eta0 = new TH1F("h_gen_eta0", "h_gen_eta0",100, -6.0, 6.0);
  h_gen_phi0 = new TH1F("h_gen_phi0", "h_gen_phi0",100, -3.2, 3.2);
  h_gen_pT3  = new TH1F("h_gen_pT3",  "h_gen_pT3", 100, 0.0, 1000.0);
  h_gen_eta3 = new TH1F("h_gen_eta3", "h_gen_eta3",100, -6.0, 6.0);
  h_gen_phi3 = new TH1F("h_gen_phi3", "h_gen_phi3",100, -3.2, 3.2);

  h_EvtBrk = new TH1F("h_EvtBrk", "h_EvtBrk", 11, 0, 11);
  h_EvtBrk->GetXaxis()->SetBinLabel(2,"e");
  h_EvtBrk->GetXaxis()->SetBinLabel(3,"mu");
  h_EvtBrk->GetXaxis()->SetBinLabel(4,"tau");
  h_EvtBrk->GetXaxis()->SetBinLabel(5,"e-e");
  h_EvtBrk->GetXaxis()->SetBinLabel(6,"mu-mu");
  h_EvtBrk->GetXaxis()->SetBinLabel(7,"tau-tau");
  h_EvtBrk->GetXaxis()->SetBinLabel(8,"e-mu");
  h_EvtBrk->GetXaxis()->SetBinLabel(9,"e-tau");
  h_EvtBrk->GetXaxis()->SetBinLabel(10,"mu-tau");
  h_EvtBrk->GetXaxis()->SetBinLabel(11,"q-q");

    
}

AnalyzeTProxytBSM::AnalyzeTProxytBSM(const TString &inputFileList, const char *outFileName,const char *dataset, const char *sample, const char* LostlepFlag, const char* phoID) {
  
  std::cout << outFileName << std::endl;

  string nameData=dataset;//vvv

  if(nameData!="signalH") nameData="BG";
  if(nameData=="signalH") nameData="signal";
  cout<<"Treating the input files as "<<nameData<<" for setting tree branches"<<endl;
  BookHistogram(outFileName); //, N2_mass);
  //CrossSection_Map_Init();
}

AnalyzeTProxytBSM::~AnalyzeTProxytBSM() { 
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();
}

#endif // AnalyzeTProxytBSM_cxx
