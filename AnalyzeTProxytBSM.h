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

  TH1F *h_ele_pT, *h_ele_eta, *h_ele_phi;
  TH1F *h_pho_pT, *h_pho_pT1, *h_pho_eta, *h_pho_phi;
  TH1F *h_gen_pT, *h_gen_eta, *h_gen_phi;
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

  h_ele_pT  = new TH1F("h_ele_pT",  "h_ele_pT", 100, 0.0, 1000.0);
  h_ele_eta = new TH1F("h_ele_eta", "h_ele_eta",100, -3.0, 3.0);
  h_ele_phi = new TH1F("h_ele_phi", "h_ele_phi",100, -3.2, 3.2);
  h_pho_pT  = new TH1F("h_pho_pT",  "h_pho_pT", 100, 0.0, 1000.0);
  h_pho_pT1  = new TH1F("h_pho_pT1",  "h_pho_pT1", 100, 0.0, 1000.0);
  h_pho_eta = new TH1F("h_pho_eta", "h_pho_eta",100, -3.0, 3.0);
  h_pho_phi = new TH1F("h_pho_phi", "h_pho_phi",100, -3.2, 3.2);
  h_gen_pT  = new TH1F("h_gen_pT",  "h_gen_pT", 100, 0.0, 1000.0);
  h_gen_eta = new TH1F("h_gen_eta", "h_gen_eta",100, -6.0, 6.0);
  h_gen_phi = new TH1F("h_gen_phi", "h_gen_phi",100, -3.2, 3.2);
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
