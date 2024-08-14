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
  //void   EventLoop(std::string buffer);
  void   EventLoop(std::string buffer,const char *,const char *);
  void   BookHistogram(const char *);//, const char *);
  Bool_t Process(Long64_t entry);
  myLV   getBestPhoton(int);
  //TLorentzVector  getBestPhoton(int);
  int    bestPhotonIndxAmongPhotons=-100;
  vector<string> selection = {"no_cut", "MET", "Pho_pT", "NHadjets", "ST", "Lep_veto", "Iso_Lep_Trk_veto","TrigEff","EvtCln","JetMetPhi"};
  vector<string> genparticle = {"Electron", "Muon", "Tau"};
  vector<string> recoparticle = {"Electron", "Muon", "Tau"};
  void CrossSection_Map_Init();
  TFile *oFile;
  // TH1D *h_MET[10];
  // TH1D *h_NHadJets[10];  
  // TH1D *h_Pho_pT[10];

  TH1D *h_MET[10];
  TH1D *h_NHadJets[10];  
  TH1D *h_Pho_pT[10];
  
  TH1F *h_Jet_pT[10], *h_Jet_eta[10], *h_Jet_phi[10];
  TH1F *h_Pho_eta[10], *h_Pho_phi[10];
  TH1F *h_Gen_pT[5][10], *h_Gen_eta[5][10], *h_Gen_phi[5][10];
  TH1F *h_Reco_pT[5][10], *h_Reco_eta[5][10], *h_Reco_phi[5][10];
  TH2F *h_NHadJets_pTSum, *h_GenRecoE;
  TH1F *h_Gen_MET, *h_EFakePho_MET, *h_LostElectron_MET, *h_LostMuon_MET, *h_HadTau_MET, *h_Rest_MET, *h_EFakePho_eta, *h_LostElectron_eta, *h_LostMuon_eta, *h_HadTau_eta, *h_Rest_eta;
};
#endif


#ifdef ANALYZETPROXYTBSM_cxx
 
//void AnalyzeLightBSM::BookHistogram(const char *outFileName, const char *N2_mass) {
void AnalyzeTProxytBSM::BookHistogram(const char *outFileName) {
  std::cout << "AnalyzeLightBSM::BookHistogram " << std::endl;

  oFile = new TFile(outFileName, "recreate");
  oFile->cd();
  //char hname_NHadJets[100], hname_Jet_Pt[100], hname_Jet_Eta[100], hname_Jet_Phi[100], hname_Met[100], hname_PhoPt[100], hname_PhoEta[100], hname_PhoPhi[100], hname_ElectronPt[100], hname_Muon_Pt[100], hname_tau_Pt[100], hname_ElectronEta[100], hname_Muon_Eta[100], hname_tau_Eta[100], hname_ElectronPhi[100], hname_Muon_Phi[100], hname_tau_Phi[100];
  char hname_NHadJets[100], hname_Jet_Pt[100], hname_Jet_Eta[100], hname_Jet_Phi[100], hname_Met[100], hname_PhoPt[100], hname_PhoEta[100], hname_PhoPhi[100], hname_GenPt[100], hname_GenEta[100], hname_GenPhi[100], hname_RecoPt[100], hname_RecoEta[100], hname_RecoPhi[100]; 
  // vector<char> hname_GenPtcl_ID = {hname_ElectronPt[100], hname_Muon_Pt[100], hname_tau_Pt[100], hname_ElectronEta[100], hname_Muon_Eta[100], hname_tau_Eta[100], hname_ElectronPhi[100], hname_Muon_Phi[100], hname_tau_Phi[100]} 
  // Book your histograms & summary counters here 
  for (int i=0; i<selection.size();i++)
    {
      sprintf(hname_NHadJets,"h_NHadJets_%s",selection[i].c_str());
      sprintf(hname_Jet_Pt,"h_Jet_Pt_%s",selection[i].c_str());
      sprintf(hname_Jet_Eta,"h_Jet_Eta_%s",selection[i].c_str());
      sprintf(hname_Jet_Phi,"h_Jet_Phi_%s",selection[i].c_str());
      sprintf(hname_Met,"h_MET_%s",selection[i].c_str());
      sprintf(hname_PhoPt,"h_Pho_Pt_%s",selection[i].c_str());
      sprintf(hname_PhoEta,"h_Pho_Eta_%s",selection[i].c_str());
      sprintf(hname_PhoPhi,"h_Pho_Phi_%s",selection[i].c_str());

      // h_NHadJets[i]= new TH1D(hname_NHadJets, hname_NHadJets, 50,0,50);
      // h_MET[i] = new TH1D(hname_Met,hname_Met,100,0,5000);
      // h_Pho_pT[i]= new TH1D(hname_PhoPt,hname_PhoPt,100,0,1000);

      h_NHadJets[i]= new TH1D(hname_NHadJets, hname_NHadJets, 50,0,50);
      h_MET[i] = new TH1D(hname_Met,hname_Met,100,0,5000);
      h_Pho_pT[i]= new TH1D(hname_PhoPt,hname_PhoPt,100,0,1000);
      
      h_Jet_pT[i]  = new TH1F(hname_Jet_Pt,hname_Jet_Pt, 100,0.0, 1000.0);
      h_Jet_eta[i] = new TH1F(hname_Jet_Eta,hname_Jet_Eta, 100, -10.0, 10.0);
      h_Jet_phi[i] = new TH1F(hname_Jet_Phi,hname_Jet_Phi,100, -3.2, 3.2);      
      h_Pho_eta[i] = new TH1F(hname_PhoEta,hname_PhoEta,100, -3.0, 3.0);
      h_Pho_phi[i] = new TH1F(hname_PhoPhi,hname_PhoPhi,100, -3.2, 3.2);

      //Defing histogram for different gen particles for different cuts
      for (int j=0; j<genparticle.size(); j++)
      {
         // sprintf(hname_GenPt[i][4],"h_Gen%s_Pt_ST",genparticle[i]);
	 // sprintf(hname_GenEta[i][4],"h_Gen%s_eta_ST",genparticle[i]);
	 // sprintf(hname_GenParticlePhi[i][4],"h_Gen%s_phi_ST",genparticle[i]);
	if (i==0 || i>3){      //i>3 conditon is to take histograms after sT cut                              
	   
	   sprintf(hname_GenPt,"h_Gen%s_Pt_%s",genparticle[j].c_str(),selection[i].c_str());
	   sprintf(hname_GenEta,"h_Gen%s_Eta_%s",genparticle[j].c_str(),selection[i].c_str());
	   sprintf(hname_GenPhi,"h_Gen%s_Phi_%s",genparticle[j].c_str(),selection[i].c_str());

	   sprintf(hname_RecoPt,"h_Reco%s_Pt_%s",recoparticle[j].c_str(),selection[i].c_str());
	   sprintf(hname_RecoEta,"h_Reco%s_Eta_%s",recoparticle[j].c_str(),selection[i].c_str());
	   sprintf(hname_RecoPhi,"h_Reco%s_Phi_%s",recoparticle[j].c_str(),selection[i].c_str());

	   
	   h_Gen_pT[j][i] = new TH1F(hname_GenPt, hname_GenPt, 100,0.0,1000.0);
	   //h_Gen_eta[j][i] = new TH1F(hname_GenEta, hname_GenEta, 500,-10,10);
	   h_Gen_eta[j][i] = new TH1F(hname_GenEta, hname_GenEta, 50,-10,10);
	   h_Gen_phi[j][i] = new TH1F(hname_GenPhi, hname_GenPhi, 100,-5.0,5.0);
       

	   if (j==2) continue;
	   h_Reco_pT[j][i] = new TH1F(hname_RecoPt, hname_RecoPt, 100,0.0,1000.0);
	   h_Reco_eta[j][i] = new TH1F(hname_RecoEta, hname_RecoEta, 500,-10.0,10.0);
	   h_Reco_phi[j][i] = new TH1F(hname_RecoPhi, hname_RecoPhi, 100,-5.0,5.0);
	   
	}
      }      
    }
  
     

  h_NHadJets_pTSum = new TH2F("h_NHadJets_pTSum","h_NHadJets_pTSum",50,0,50,10000,0,3500);
  h_GenRecoE = new TH2F("h_GenRecoE", "h_GenRecoE",10,0,10,10,0,10);
  h_GenRecoE->SetXTitle("Gen");
  h_GenRecoE->SetYTitle("Reco");

  
  h_EFakePho_eta = new TH1F("h_EFakePho_Eta","hname_EFakePho_Eta",50, -10.0, 10.0);
  h_LostElectron_eta = new TH1F("h_LostElectron_Eta","hname_LostElectron_Eta",50, -10.0, 10.0);
  h_LostMuon_eta = new TH1F("h_LostMuon_Eta","hname_LostMuon_Eta",50, -10.0, 10.0);
  h_HadTau_eta = new TH1F("h_HadTau_Eta", "h_HadTau_Eta",50,-10.0,10.0);
  h_Rest_eta = new TH1F ("h_Rest_Eta", "h_Rest_Eta", 50,-10.0,10.0);

  h_Gen_MET = new TH1F("h_Gen_MET", "h_Gen_MET", 100,0,5000);
  h_EFakePho_MET = new TH1F("h_EFakePho_MET","hname_EFakePho_MET",100, 0.0, 5000.0);
  h_LostElectron_MET = new TH1F("h_LostElectron_MET","hname_LostElectron_MET",100, 0.0, 5000.0);
  h_LostMuon_MET = new TH1F("h_LostMuon_MET","hname_LostMuon_MET",100, 0.0, 5000.0);
  h_HadTau_MET = new TH1F("h_HadTau_MET", "h_HadTau_MET",100, 0.0, 5000.0);
  h_Rest_MET = new TH1F ("h_Rest_MET", "h_Rest_MET", 100, 0.0, 5000.0);

}

AnalyzeTProxytBSM::AnalyzeTProxytBSM(const TString &inputFileList, const char *outFileName,const char *dataset, const char *sample, const char* LostlepFlag, const char* phoID) {
  
  std::cout << outFileName << std::endl;

  string nameData=dataset;//vvv

  if(nameData!="signalH") nameData="BG";
  if(nameData=="signalH") nameData="signal";
  cout<<"Treating the input files as "<<nameData<<" for setting tree branches"<<endl;
  BookHistogram(outFileName); //, N2_mass);
  CrossSection_Map_Init();
}

void AnalyzeTProxytBSM::CrossSection_Map_Init()
{
  char *f_name_EH = new char[2000];
  sprintf(f_name_EH,"./map_crosssection_SMprocess_v1.txt");//,chi2_method);
  std::ifstream in_EH(f_name_EH);
  if(!in_EH) {
    cout<<"ERROR => "<<f_name_EH<<" Not found"<<endl;
    //return;                                                                                                                                
    exit(0);
  }
  string process_name;
  float value, entries;
  cout<<"File name = "<<f_name_EH<<endl;
  while(in_EH>>process_name>>value>>entries){
    std::pair<std::string, float> temp_pair;    
    float weight =value/entries;
    cout << "values: " << value << " entries: " << entries << endl;
    temp_pair = std::make_pair(process_name,weight);
    cross_sectionValues.insert(temp_pair);
  }
}


AnalyzeTProxytBSM::~AnalyzeTProxytBSM() { 
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();
}

#endif // AnalyzeTProxytBSM_cxx






//backups

/* TH1F *h_Jet_pT0, *h_Jet_eta0, *h_Jet_phi0;  //0 stands for no cut  */
  /* TH1F *h_Jet_pT1, *h_Jet_eta1;   //1 stands for cut in the pT  */
  /* TH1F *h_Jet_pT2, *h_Jet_eta2;   //2 stands for cut in pT and Eta */
  /* TH1F *h_Jet_pT3, *h_Jet_eta3, *h_Jet_phi3;   //3 stands for cut in pT and Eta and MET */
  /* TH1F *h_Jet_pT4, *h_Jet_eta4, *h_Jet_phi4; */
  /* TH1F *h_Jet_pT5, *h_Jet_eta5, *h_Jet_phi5; */
  
  
  //TH2F *h_NJet_MET, *h_NJet_PhoPt, *h_NJet_genPhoPt; 

  //TH1F *h_ele_pT0, *h_ele_eta0, *h_ele_phi0; 
  // TH1F *h_ele_pT3, *h_ele_eta3, *h_ele_phi3;  
  
  /* TH1F *h_pho_eta0,*h_pho_eta3, *h_pho_eta4, *h_pho_eta5; */
  /* TH1F *h_pho_pT0, *h_pho_pT3, *h_pho_pT4, *h_pho_pT5; */
  /* TH1F *h_pho_phi0, *h_pho_phi3, *h_pho_phi4, *h_pho_phi5; */
  //TH1F *h_gen_pT0, *h_gen_eta0, *h_gen_phi0;
  //TH1F *h_gen_pT3, *h_gen_eta3, *h_gen_phi3;

 //h_MET0 = new TH1F("h_MET0", "h_MET0", 100, 0.0, 1000.0);
  /* h_MET3 = new TH1F("h_MET3", "hMET3", 100, 0.0, 1000.0); */
  /* h_MET4 = new TH1F("h_MET4", "hMET4", 100, 0.0, 1000.0); */
  /* h_MET5 = new TH1F("h_MET5", "hMET5", 100, 0.0, 1000.0); */
  // h_hadJets_Pt = new TH1F("h_hadJets_Pt", "h_hadJets_Pt", 100, 0.0, 1000.0);
  // h_hadJets_Pt1 = new TH1F("h_hadJets_Pt1", "h_hadJets_Pt1", 100, 0.0, 1000.0);
  // h_hadJets_Eta = new TH1F("h_hadJets_Eta", "h_hadJets_Eta", 100, 0.0, 10.0);
  // h_hadJets_Eta1 = new TH1F("h_hadJets_Eta1", "h_hadJets_Eta1", 100, 0.0, 10.0);
  
  /* h_Jet_pT0  = new TH1F("h_Jet_pT0",  "h_Jet_pT0", 100, 0.0, 1000.0); */
  /* h_Jet_pT1  = new TH1F("h_Jet_pT1",  "h_Jet_pT1", 100, 0.0, 1000.0); */
  /* h_Jet_pT2  = new TH1F("h_Jet_pT2",  "h_Jet_pT2", 100, 0.0, 1000.0); */
  /* h_Jet_pT3  = new TH1F("h_Jet_pT3",  "h_Jet_pT3", 100, 0.0, 1000.0); */
  /* h_Jet_pT4  = new TH1F("h_Jet_pT4",  "h_Jet_pT4", 100, 0.0, 1000.0); */
  /* h_Jet_pT5  = new TH1F("h_Jet_pT5",  "h_Jet_pT5", 100, 0.0, 1000.0); */
  /* h_Jet_eta0 = new TH1F("h_Jet_eta0", "h_Jet_eta0", 100, -10.0, 10.0); */
  /* h_Jet_eta1 = new TH1F("h_Jet_eta1", "h_Jet_eta1", 100, -10.0, 10.0); */
  /* h_Jet_eta2 = new TH1F("h_Jet_eta2", "h_Jet_eta2", 100, -10.0, 10.0); */
  /* h_Jet_eta3 = new TH1F("h_Jet_eta3", "h_Jet_eta3", 100, -10.0, 10.0); */
  /* h_Jet_eta4 = new TH1F("h_Jet_eta4", "h_Jet_eta4", 100, -10.0, 10.0); */
  /* h_Jet_eta5 = new TH1F("h_Jet_eta5", "h_Jet_eta5", 100, -10.0, 10.0); */
  /* h_Jet_phi0 = new TH1F("h_Jet_phi0", "h_Jet_phi0",100, -3.2, 3.2); */
  /* h_Jet_phi3 = new TH1F("h_Jet_phi3", "h_Jet_phi3",100, -3.2, 3.2); */
  /* h_Jet_phi4 = new TH1F("h_Jet_phi4", "h_Jet_phi4",100, -3.2, 3.2); */
  /* h_Jet_phi5 = new TH1F("h_Jet_phi5", "h_Jet_phi5",100, -3.2, 3.2); */
  

  /* h_NHadJets0 = new TH1F("h_NHadJets0", "h_NHadJets0",50, 0, 50); */
  /* h_NHadJets1 = new TH1F("h_NHadJets1", "h_NHadJets1",50, 0, 50); */
  /* h_NHadJets = new TH1F("h_NHadJets", "h_NHadJets",50, 0, 50); */
  /* h_NHadJets3 = new TH1F("h_NHadJets3", "h_NHadJets3",50, 0, 50); */
  /* h_NHadJets4 = new TH1F("h_NHadJets4", "h_NHadJets4",50, 0, 50); */
  /* h_NHadJets5 = new TH1F("h_NHadJets5", "h_NHadJets5",50, 0, 50); */

  // h_NHadJets0 = new TH1F("h_NHadJets0", "h_NHadJets0",50, 0, 50);
  // h_NHadJets3 = new TH1F("h_NHadJets3", "h_NHadJets3",50, 0, 50);
  // h_NHadJets4 = new TH1F("h_NHadJets4", "h_NHadJets4",50, 0, 50);
  // h_NHadJets5 = new TH1F("h_NHadJets5", "h_NHadJets5",50, 0, 50);
   
  // h_ele_pT0  = new TH1F("h_ele_pT0",  "h_ele_pT0", 100, 0.0, 1000.0);
  // h_ele_eta0 = new TH1F("h_ele_eta0", "h_ele_eta0",100, -3.0, 3.0);
  // h_ele_phi0 = new TH1F("h_ele_phi0", "h_ele_phi0",100, -3.2, 3.2);
  // h_ele_pT3  = new TH1F("h_ele_pT3",  "h_ele_pT3", 100, 0.0, 1000.0);
  // h_ele_eta3 = new TH1F("h_ele_et3", "h_ele_eta3",100, -3.0, 3.0);
  // h_ele_phi3 = new TH1F("h_ele_phi3", "h_ele_phi3",100, -3.2, 3.2);

  /* h_pho_pT0  = new TH1F("h_pho_pT0",  "h_pho_pT0", 100, 0.0, 1000.0); */
  /* //h_pho_pT1  = new TH1F("h_pho_pT1",  "h_pho_pT1", 100, 0.0, 1000.0); */
  /* h_pho_eta0 = new TH1F("h_pho_eta0", "h_pho_eta0",100, -3.0, 3.0); */
  /* h_pho_phi0 = new TH1F("h_pho_phi0", "h_pho_phi0",100, -3.2, 3.2); */
  /* h_pho_pT3 = new TH1F("h_pho_pT3",  "h_pho_pT3", 100, 0.0, 1000.0); */
  /* h_pho_pT4 = new TH1F("h_pho_pT4",  "h_pho_pT4", 100, 0.0, 1000.0); */
  /* h_pho_pT5 = new TH1F("h_pho_pT5",  "h_pho_pT5", 100, 0.0, 1000.0); */
  /* //h_pho_pT1  = new TH1F("h_pho_pT1",  "h_pho_pT1", 100, 0.0, 1000.0); */
  /* h_pho_eta3 = new TH1F("h_pho_eta3", "h_pho_eta3",100, -3.0, 3.0); */
  /* h_pho_eta4 = new TH1F("h_pho_eta4", "h_pho_eta4",100, -3.0, 3.0); */
  /* h_pho_eta5 = new TH1F("h_pho_eta5", "h_pho_eta5",100, -3.0, 3.0); */
  
  /* h_pho_phi3= new TH1F("h_pho_phi3", "h_pho_phi3",100, -3.2, 3.2); */
  /* h_pho_phi4= new TH1F("h_pho_phi4", "h_pho_phi4",100, -3.2, 3.2); */
  /* h_pho_phi5= new TH1F("h_pho_phi5", "h_pho_phi5",100, -3.2, 3.2); */

  // h_NHadJet_MET = new TH2F("h_NHadJet_MET", "h_NHadJet_MET", 100, 0.0, 1000.0, 50, 0, 50);
  // h_NHadJet_MET->GetXaxis()->SetTitle("MET");
  // h_NHadJet_MET->GetYaxis()->SetTitle("NJet");

  // h_NJet_PhoPt = new TH2F("h_NJet_PhoPt", "h_PhoPt", 100, 0.0, 1000.0, 50, 0, 50);
  // h_NJet_PhoPt->GetXaxis()->SetTitle("Pho_Pt");
  // h_NJet_PhoPt->GetYaxis()->SetTitle("NJet");
  
  // h_NJet_genPhoPt = new TH2F("h_genNJet_PhoPt", "h_genPhoPt", 100, 0.0, 1000.0, 50, 0, 50);
  // h_NJet_genPhoPt->GetXaxis()->SetTitle("genPho_Pt");
  // h_NJet_genPhoPt->GetYaxis()->SetTitle("NJet");
  
  // h_gen_pT0  = new TH1F("h_gen_pT0",  "h_gen_pT0", 100, 0.0, 1000.0);
  // h_gen_eta0 = new TH1F("h_gen_eta0", "h_gen_eta0",100, -6.0, 6.0);
  // h_gen_phi0 = new TH1F("h_gen_phi0", "h_gen_phi0",100, -3.2, 3.2);
  // h_gen_pT3  = new TH1F("h_gen_pT3",  "h_gen_pT3", 100, 0.0, 1000.0);
  // h_gen_eta3 = new TH1F("h_gen_eta3", "h_gen_eta3",100, -6.0, 6.0);
  // h_gen_phi3 = new TH1F("h_gen_phi3", "h_gen_phi3",100, -3.2, 3.2);

  // h_EvtBrk = new TH1F("h_EvtBrk", "h_EvtBrk", 11, 0, 11);
  // h_EvtBrk->GetXaxis()->SetBinLabel(2,"e");
  // h_EvtBrk->GetXaxis()->SetBinLabel(3,"mu");
  // h_EvtBrk->GetXaxis()->SetBinLabel(4,"tau");
  // h_EvtBrk->GetXaxis()->SetBinLabel(5,"e-e");
  // h_EvtBrk->GetXaxis()->SetBinLabel(6,"mu-mu");
  // h_EvtBrk->GetXaxis()->SetBinLabel(7,"tau-tau");
  // h_EvtBrk->GetXaxis()->SetBinLabel(8,"e-mu");
  // h_EvtBrk->GetXaxis()->SetBinLabel(9,"e-tau");
  // h_EvtBrk->GetXaxis()->SetBinLabel(10,"mu-tau");
  // h_EvtBrk->GetXaxis()->SetBinLabel(11,"q-q");

  /* //TH1F *h_MET3; */
  /* TH1F *h_MET4; */
  /* TH1F *h_MET5; */
  /* TH1F *h_EvtBrk; */
  /* TH1F *h_hadJets_Pt; */
  /* TH1F *h_hadJets_Pt1; */
  /* TH1F *h_hadJets_Eta; */
  /* TH1F *h_hadJets_Eta1; */


/* TH1F *h_NHadJets0; */
  /* TH1F *h_NHadJets3; */
  /* TH1F *h_NHadJets4; */
  /* TH1F *h_NHadJets5; */

 //jets with no cuts
  /* TH1F *h_NHadJets1;   // jets with cut on Pt */
  /* TH1F *h_NHadJets;   // along with eta cut */
  /* TH1F *h_NHadJets3;   // along with MET cut */
  /* TH1F *h_NHadJets4;   // cut on photon pt */
  /* TH1F *h_NHadJets5;   // cut on NHadJets */
