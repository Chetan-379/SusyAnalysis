#define ANALYZETPROXYTBSM_cxx
#include "AnalyzeTProxytBSM.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <fstream>
#include"TGraphErrors.h"
#include"TGraphAsymmErrors.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
//#pragma link C++ class std::vector< std::vector >+; 
//#pragma link C++ class std::vector< TLorentzVector >+;
//#ifdef __MAKECINT__
//#pragma link C++ class NtupleVarsTProxy+;
//#endif
using namespace TMVA;
int main(int argc, char* argv[])
{ 
  if (argc < 6) {
    cerr << "Please give 6 arguments " << "runList " << " " << "outputFileName" << " " << "which year dataset" <<" "<<"which Process"<< " "<<"which Lostlep bkg"<< " "<<"Which pho_ID"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *sample=argv[4];
  const char *elec = argv[5];
  const char *phoID = argv[6];

  //TString pho_ID = phoID;
  AnalyzeTProxytBSM ana(inputFileList, outFileName, data,sample, elec,phoID);

  //=== === Loop over input files === ====
  int iFile = 0; 
  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;
  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    iFile++;
    std::cout << "===========================================================================" << std::endl;
    std::cout << "iFile " << iFile << "  Analyzing tree from " << buffer.c_str() << std::endl;
    std::cout << "===========================================================================" << std::endl;

    // for skimmed tree
    TFile *fin = new TFile(buffer.c_str());
    TTree *chain = (TTree*) fin->FindObjectAny("PreSelection");
    std::cout << "main(): chain->GetEntries() "<<  chain->GetEntries() <<std::endl;
    ana.Init(chain);
    ana.EventLoop(buffer.c_str(),data,sample);
    delete chain; 
    delete fin;
  }

  // === === some random summary === ===
  cout << "dataset " << data << " " << endl;
  cout<<"If analyzing the lost electron estimation ? "<<"  "<<elec<<endl;
  cout<<"Which pho_ID: "<<"\t"<<phoID<<endl;
  Tools::Instance();  
  return 0;
}

void AnalyzeTProxytBSM::EventLoop(std::string buffer, const char *data, const char *sample) {
  
  std::cout << "AnalyzeTProxytBSM::EventLoop() " << std::endl;
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  cout << "Analyzing " << buffer.c_str() << " nentries " << nentries << std::endl;  
  char* s_cross = new char[100];
  sprintf(s_cross,"%s.%s",data,sample);
  std::string s_process = s_cross;
  TString s_Process = s_process;
  TString s_sample= sample;
  TString s_data= data;
  double cross_section = getCrossSection(s_process);
  double wt = 0;
  double wt1 = 0;

  std::cout << cross_section << "\t" <<"analyzed process"<<"\t"<<s_cross<<endl;
 
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
  double sumwt = 0;
  double lumiInfb = 0;
  
  int nEvents=0, NGenL=0, NLostElectrons=0, NLostMuons=0, NEFakePho=0;

  // string src_file1, src_file2, src_file3;
  // //if (s_data.Contains("20")) cout << "half detected!" << endl; 
  // // if(s_data.Contains("2016preVFP")){
  // //  src_file1 = "2016APV_lost_e_TF_NBJet.root";
  // //  src_file2 = "2016APV_lost_mu_TF_NBJet.root";
  // //  src_file3 = "2016APV_LL_TF_NBJet.root";
  // // }
  
  // // if(s_data.Contains("2016postVFP")){
  // //   src_file1 = "2016_lost_e_TF_NBJet.root";
  // //   src_file2 = "2016_lost_mu_TF_NBJet.root";
  // //   src_file3 = "2016_LL_TF_NBJet.root";
  // // }

  // // if(s_data.Contains("2017")){
  // //   src_file1 = "2017_lost_e_TF_NBJet.root";
  // //   src_file2 = "2017_lost_mu_TF_NBJet.root";
  // //   src_file3 = "2017_LL_TF_NBJet.root";
  // // }
  
  // // if(s_data.Contains("2018")){
  // //   src_file1 = "final_2018_lost_e_TF_NBJet.root";
  // //   src_file2 = "final_2018_lost_mu_TF_NBJet.root";
  // //   src_file3 = "final_2018_LL_TF_NBJet.root";
  // // }
  




  // src_file1 = "year_avg_lost_e_TF_NBJet.root";
  // src_file2 = "year_avg_lost_mu_TF_NBJet.root";
  // src_file3 = "year_avg_LL_TF_NBJet.root";
  
  
  // TFile *root_file1 = new TFile(src_file1.c_str());
  // TFile *root_file2 = new TFile(src_file2.c_str());
  // TFile *root_file3 = new TFile(src_file3.c_str());
  
  // TH1F *h_TF1, *h_TF2, *h_TF3;	   
  // h_TF1 = (TH1F*) root_file1 ->Get("combined_TF");
  // h_TF2 = (TH1F*) root_file2 ->Get("combined_TF");
  // h_TF3 = (TH1F*) root_file3 ->Get("combined_TF");
  // float TF1 = 0, TF2 = 0, TF3 = 0;
  // double wt2 = 0, wt3 = 0, wt4 = 0;
  
  
  for (Long64_t jentry=0; jentry<fChain->GetEntries(); jentry++){
    fDirector.SetReadEntry(jentry);
    
    // == == print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);      
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" << endl;
    decade = k;

    vector<myLV> hadJets, bjets;
    int BTags; 
    bool Debug=false;
    vector<int> jetMatchindx;
    int bJet1Idx=-100; 
    double deepCSVvalue = 0.4168,p0=1.787e+02,p1=6.657e+01,p2=9.47e-01;
    double minDR=99999;
    int phoMatchingJetIndx = -100, photonMatchingJetIndx = -100; 
    int minDRindx=-100;
    double qmulti = 0.0;
    int pho_ID=0;  //for simplicity taking only soft one
    myLV bestPhoton=getBestPhoton(pho_ID);
    int hadJetID=-999;
    int NJets0=Jets->size();
    int NHadJets = 0;
    float Jets_pT_Sum=0;
    float ST=0;
    int NEMu;
    double dPhi_METjet1, dPhi_METjet2;
    bool applyHEMveto = true, applyL1TrigFire_prob = true;
    bool apply_pixelveto = true;

    if(s_data.Contains("2016preVFP")){lumiInfb=19.5;deepCSVvalue = 0.6001; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01;}// APV
    if(s_data.Contains("2016postVFP")) {lumiInfb=16.5; deepCSVvalue = 0.5847; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01;} //2016
  
    if(s_data.Contains("2017")) {lumiInfb=41.48; p0=1.82e+02; p1=6.336e+01; p2=9.171e-01; deepCSVvalue = 0.4506;}
    if(s_data.Contains("2018")) {lumiInfb=59.83; p0=1.787e+02; p1=6.657e+01; p2=9.47e-01; deepCSVvalue = 0.4168;}

    wt = getEventWeight(s_process, lumiInfb, cross_section);
  
    wt = wt*puWeight;
  
    //applying hemveto                                                                  
    bool HEMaffected=false;
    if(s_data.Contains("2018") && applyHEMveto){
      for(int i=0; i<Electrons->size();i++){
	if(Electrons[i].Pt() >30 && Electrons[i].Eta() > -3.0 && Electrons[i].Eta() < -1.4 && Electrons[i].Phi() > -1.57 && Electrons[i].Phi() < -0.87) {HEMaffected = true; break;}
      }
      for(int i=0; i<Jets->size();i++){
	if(Jets[i].Pt() > 30 && Jets[i].Eta() > -3.2 && Jets[i].Eta() < -1.2 && Jets[i].Phi() > -1.77 && Jets[i].Phi() < -0.67 && DeltaPhi(Jets[i].Pt(),METPhi)<0.5) {HEMaffected = true; break;}
      }
      if(HEMaffected == true) continue;           
    }
    
    //adding l1trigger prefire issue probability
    if(applyL1TrigFire_prob && (s_data.Contains("2016postVFP") ||  s_data.Contains("2017") ))
      {
	wt =wt*NonPrefiringProb;
      }

    //=======================makeing some gen level and reco level collectons==============================================
    int nGenEle1=0, nGenMu1=0, nGenTau1=0, nGentau_had1=0;
    TLorentzVector genEle1, genMu1, genTau1;
    vector<TLorentzVector> v_genTau2,v_genEle2,v_genMu2;
    v_genEle2.clear();
    v_genMu2.clear();
    v_genTau2.clear();
    
    for(int i=0 ; i < GenElectrons->size(); i++){
      if(GenElectrons[i].Pt()!=0){
	nGenEle1++;
	genEle1.SetPtEtaPhiE(GenElectrons[i].Pt(), GenElectrons[i].Eta(), GenElectrons[i].Phi(), GenElectrons[i].E());
	v_genEle2.push_back(genEle1);
      }
    }
    
    for(int i=0 ; i < GenMuons->size(); i++){
      if(GenMuons[i].Pt()!=0){
	nGenMu1++;
	genMu1.SetPtEtaPhiE(GenMuons[i].Pt(), GenMuons[i].Eta(), GenMuons[i].Phi(), GenMuons[i].E());
	v_genMu2.push_back(genMu1);
      }
    }
    
    for(int i=0 ; i < GenTaus->size(); i++){
      if(GenTaus[i].Pt()!=0){
	nGenTau1++;
	genTau1.SetPtEtaPhiE(GenTaus[i].Pt(), GenTaus[i].Eta(), GenTaus[i].Phi(), GenTaus[i].E());
	v_genTau2.push_back(genTau1);
	if((*GenTaus_had)[i])
	  nGentau_had1++;	
      }
    }
    
    sortTLorVec(&v_genTau2);
    sortTLorVec(&v_genMu2);
    sortTLorVec(&v_genEle2);


    //Reco:
    TLorentzVector RecoE_LV, bestPhoton_v1;
    vector<TLorentzVector> Electrons_v1;
    for(int i=0; i< Electrons ->size(); i++){
      RecoE_LV.SetPtEtaPhiE(Electrons[i].Pt(), Electrons[i].Eta(), Electrons[i].Phi(), Electrons[i].E());
      Electrons_v1.push_back(RecoE_LV);
    }

    sortTLorVec(&Electrons_v1);
    //sortTLorVec(&Electrons);
    //sortmyLV(&Electrons);

    bestPhoton_v1.SetPtEtaPhiE(bestPhoton.Pt(), bestPhoton.Eta(), bestPhoton.Phi(), bestPhoton.E());
    
    //selecting good reco electrons for Electrons vector
    bool Elec_passAccep =false, Elec_passEtacut =false, Elec_passpTcut =false, Elec_passIso =false;
    int e_index = -1, nelec_reco =0, nlep=0;
    TLorentzVector emobj_elec;
    for(int i=0; i < Electrons_v1.size(); i++){
	// if(nelec_reco>0) continue;
	if((Electrons_v1[i].Pt() >10) && abs(Electrons_v1[i].Eta()) < 2.5){
	  if( (*Electrons_passIso)[i]==1)
    	    {	      
    	      //pass_isoElec++; Elec_passIso = true; nelec_reco++; nlep++; e_index=i; recElec=Electrons_v1[i]; v_recEle.push_back(Electrons_v1[i]);
	      Elec_passIso = true; e_index=i; nelec_reco++; nlep++;
	      emobj_elec=Electrons_v1[i];
    	    }	    
	}
    }

    //=====================================================================================================================

    //==================================Some new variables for FR=============================
    bool bestEMObjIsEle = false, bestEMObjIsEle_px=false, bestEMobj=false;
    TLorentzVector bestEMObj;
    //myLV bestEMObj;
    myLV bestEMObj_myLV;
    int  hasPho =0, hasPho_px =0, hasEle =0, lep =0;

   if(nlep==1 && emobj_elec.Pt()>40){
     hasEle =1;
     lep++;    
   }

   bool bestPhoHasPxlSeed=true, noPhoHasPxlSeed=true;
   if(bestPhotonIndxAmongPhotons>=0){
     if((*Photons_hasPixelSeed)[bestPhotonIndxAmongPhotons]<0.001) bestPhoHasPxlSeed=false;

     if(!bestPhoHasPxlSeed && bestPhoton.Pt()>40) {hasPho=1; hasPho_px=0;}
     else if(bestPhoHasPxlSeed && bestPhoton.Pt()>40) {hasPho_px=true; hasPho=0;}
     else {hasPho_px=0; hasPho=0;}
   }
  

   //if (hasEle == 1 && hasPho == 0) {bestEMObjIsEle =true; bestEMObj = Electrons_v1[e_index]; bestEMobj=true; }
   if (hasEle == 1 && hasPho == 0) {bestEMObjIsEle =true; bestEMObj = Electrons_v1[e_index]; bestEMobj=true; }
   //else if(hasEle==0 &&  hasPho==1) {bestEMObjIsEle=false;  bestEMObjIsEle_px=false;bestEMObj = bestPhoton_v1; bestEMobj=true;}
   else if(hasEle==0 &&  hasPho==1) {bestEMObjIsEle=false;  bestEMObjIsEle_px=false; bestEMObj = bestPhoton_v1; bestEMobj=true;}
   else continue;
   if(bestEMobj==false) continue;
   //h_selectBaselineYields_v1->Fill("no reco #gamma or e",wt);
   if(bestEMObj.Pt()<40) continue;
   bestEMObj_myLV.SetPxPyPzE(bestEMObj.Px(), bestEMObj.Py(), bestEMObj.Pz(), bestEMObj.E());

  
  
//===============================================================================    
   if(apply_pixelveto){
      if(s_data.Contains("2017") && ((bestEMObj.Eta() > 0.0 && bestEMObj.Eta() < 1.5 && bestEMObj.Phi() > 2.7))) continue;
      if(s_data.Contains("2018") && ((bestEMObj.Eta() > 0.0 && bestEMObj.Eta() < 1.5 && bestEMObj.Phi() > 0.4 && bestEMObj.Phi() < 0.8))) continue;
    }
    
    
    //selecting Hadronic and b jets
    for(int i=0;i<Jets->size();i++){
      if( (Jets[i].Pt() > 30.0) && (abs(Jets[i].Eta()) <= 2.4) ){	
	  double dR = DeltaR(bestEMObj.Eta(), bestEMObj.Phi(),Jets[i].Eta(),Jets[i].Phi());
	  if(dR<minDR) {minDR=dR;minDRindx=i;}	
      } 
    }
    
    if(Debug)
      cout<<"===load tree entry  ==="<<"\t"<<jentry<<"\t"<<"Jets check == "<<minDR<<endl;
    
    for(int i=0;i<Jets->size();i++){
      if( (Jets[i].Pt() > 30.0) && (abs(Jets[i].Eta()) <= 2.4) ){	  
	if( !(minDR < 0.3 && i==minDRindx) )
	  {		
	    hadJetID= (*Jets_ID)[i];
	    if(hadJetID)
	      {
		hadJets.push_back(Jets[i]);
		if((*Jets_bJetTagDeepCSVBvsAll)[i] > deepCSVvalue){
		  bjets.push_back(Jets[i]); bJet1Idx = i;}		  
		jetMatchindx.push_back(i);
	      }
	  }	
      }
    }
    
    BTags = bjets.size();
    
    if( minDR<0.3){
      photonMatchingJetIndx = minDRindx;
      qmulti=(*Jets_chargedMultiplicity)[photonMatchingJetIndx];      
    }

    if( photonMatchingJetIndx>=0 ){
      if( (Jets[photonMatchingJetIndx].Pt()) > 1.1*(bestEMObj.Pt()) ){
	if( ((Jets[photonMatchingJetIndx] - bestEMObj_myLV).Pt())>30){
	  hadJets.push_back(Jets[photonMatchingJetIndx] - bestEMObj_myLV);
	}
      }
    }

    
    if(hadJets.size()==0) continue;
    
    for(int i=0;i<hadJets.size();i++){
      
      if( (abs(hadJets[i].Eta()) < 2.4) ){
	NHadJets++;
	  ST=ST+(hadJets[i].Pt());
      }
    }

    //sortTLorVec(&hadJets);
    
    if( minDR<0.3) {
      ST=ST+bestEMObj.Pt();
      phoMatchingJetIndx = minDRindx;
    }
    
    if(NHadJets>=1)
      dPhi_METjet1 = abs(DeltaPhi(METPhi,hadJets[0].Phi()));
    if(NHadJets>=2)
      dPhi_METjet2 = abs(DeltaPhi(METPhi,hadJets[1].Phi()));
    
    NEMu = NElectrons + NMuons;

    bool EvtCln_18 =false, EvtCln_17 =false, EvtCln_16 = false;
    if ((s_data.Contains("2017") || s_data.Contains("2018")) && PrimaryVertexFilter==1 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 &&HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter==1 && BadPFMuonDzFilter==1 && eeBadScFilter==1 && ecalBadCalibFilter==1 && NVtx>0 && PFCaloMETRatio < 5){
      if((!(phoMatchingJetIndx>=0 && (Jets[phoMatchingJetIndx].Pt())/(bestPhoton.Pt()) < 1.0)) && phoMatchingJetIndx >= 0) {
	EvtCln_17 = true;
	EvtCln_18 = true;}
    }

    if (s_data.Contains("2016") && PrimaryVertexFilter==1 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 &&HBHEIsoNoiseFilter==1 &&EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter==1 && BadPFMuonDzFilter==1 && eeBadScFilter==1 && PFCaloMETRatio < 5){
      if((!(phoMatchingJetIndx>=0 && (Jets[phoMatchingJetIndx].Pt())/(bestPhoton.Pt()) < 1.0)) && phoMatchingJetIndx >= 0) EvtCln_16 = true;}

      
    
    //defining flags for applying baseline selections
    bool Pass_EMu_veto=false, Pass_Iso_trk_veto=false, Pass_Pho_pT=false, Pass_MET=false, Pass_NHadJets=false, Pass_ST=false, applyTrgEff = false, EvtCln=false,JetMetPhi=false,rmOvrlp=false, Pass_MET2=false;
    if (bestEMObj.Pt() >40) Pass_Pho_pT = true;
    if (Pass_Pho_pT && MET > 100) Pass_MET = true;
    if (Pass_MET && NHadJets >=2) Pass_NHadJets = true;		
    if (Pass_NHadJets && ST > 300) {Pass_ST = true; applyTrgEff = true;}
    if (applyTrgEff && (EvtCln_18 || EvtCln_17 || EvtCln_16)) EvtCln = true;
    if (EvtCln && dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3) JetMetPhi = true;		      
    if (JetMetPhi && RemoveSampleOverlap(s_sample, bestPhoton)) rmOvrlp = true;
    if (rmOvrlp && MET>200) Pass_MET2=true;
    if (Pass_MET2 && NEMu == 0) Pass_EMu_veto = true;	
    if (Pass_EMu_veto && !(isoElectronTracks || isoMuonTracks || isoPionTracks)) Pass_Iso_trk_veto = true;									
    
    //if(jentry<20) cout << wt << endl;      
    //h_NHadJets[0]->Fill(0.0,wt);      
    h_NHadJets[0]->Fill(NHadJets,wt);
    h_NbJets[0]->Fill(BTags,wt);
    h_MET[0]->Fill(MET,wt);
    h_Pho_pT[0] -> Fill(bestPhoton.Pt(),wt);
    h_Pho_eta[0] -> Fill(bestPhoton.Eta(),wt);
    h_Pho_phi[0] -> Fill(bestPhoton.Phi(),wt);
    
    if (Pass_Pho_pT){
      h_MET[2] -> Fill(MET,wt);
      h_Pho_pT[2] -> Fill(bestPhoton.Pt(),wt);
      h_Pho_eta[2] -> Fill(bestPhoton.Eta(),wt);
      h_Pho_phi[2] -> Fill(bestPhoton.Phi(),wt);
      h_NHadJets[2]-> Fill(NHadJets,wt);
      h_NbJets[2]->Fill(BTags,wt);
    }
    
    if (Pass_MET){ 
      h_MET[1]-> Fill(MET,wt);
      h_Pho_pT[1]  ->Fill(bestPhoton.Pt(),wt);
      h_Pho_eta[1]  ->Fill(bestPhoton.Eta(),wt);
      h_Pho_phi[1]  ->Fill(bestPhoton.Phi(),wt);
      h_NHadJets[1]-> Fill(NHadJets,wt);
      h_NbJets[1]->Fill(BTags,wt);
    }
    
    if (Pass_NHadJets){
      h_MET[3] ->Fill(MET,wt);
      h_Pho_pT[3] ->Fill(bestPhoton.Pt(),wt);
      h_Pho_eta[3] -> Fill(bestPhoton.Eta(),wt);
      h_Pho_phi[3] -> Fill(bestPhoton.Phi(),wt);
      h_NHadJets[3]-> Fill(NHadJets,wt);
      h_NbJets[3]->Fill(BTags,wt);
    }
    
    if (Pass_ST){
      h_MET[4] ->Fill(MET,wt);
      h_Pho_pT[4] ->Fill(bestPhoton.Pt(),wt);
      h_Pho_eta[4] -> Fill(bestPhoton.Eta(),wt);
      h_Pho_phi[4] -> Fill(bestPhoton.Phi(),wt);
      h_NHadJets[4]-> Fill(NHadJets,wt);
      h_NbJets[4]->Fill(BTags,wt);
    }
    
    if (applyTrgEff)
      {
	wt1 = wt * (((TMath::Erf((MET - p0)/p1)+1)/2.0)*p2);
	h_MET[7] ->Fill(MET,wt1);
	h_Pho_pT[7] ->Fill(bestPhoton.Pt(),wt1);
	h_Pho_eta[7] -> Fill(bestPhoton.Eta(),wt1);
	h_Pho_phi[7] -> Fill(bestPhoton.Phi(),wt1);
	h_NHadJets[7]-> Fill(NHadJets,wt1);
	h_NbJets[7]->Fill(BTags,wt1);
      }
    
    if (EvtCln)
      {
	h_MET[8] ->Fill(MET,wt1);
	h_Pho_pT[8] ->Fill(bestPhoton.Pt(),wt1);
	h_Pho_eta[8] -> Fill(bestPhoton.Eta(),wt1);
	h_Pho_phi[8] -> Fill(bestPhoton.Phi(),wt1);
	h_NHadJets[8]-> Fill(NHadJets,wt1);
	h_NbJets[8]->Fill(BTags,wt1);
      }
    
    if (JetMetPhi)
      {
	h_MET[9] ->Fill(MET,wt1);
	h_Pho_pT[9] ->Fill(bestPhoton.Pt(),wt1);
	h_Pho_eta[9] -> Fill(bestPhoton.Eta(),wt1);
	h_Pho_phi[9] -> Fill(bestPhoton.Phi(),wt1);
	h_NHadJets[9]-> Fill(NHadJets,wt1);
	h_NbJets[9]->Fill(BTags,wt1);
      }
    
    if (rmOvrlp)
      {
	h_MET[10] ->Fill(MET,wt1);
	h_Pho_pT[10] ->Fill(bestPhoton.Pt(),wt1);
	h_Pho_eta[10] -> Fill(bestPhoton.Eta(),wt1);
	h_Pho_phi[10] -> Fill(bestPhoton.Phi(),wt1);
	h_NHadJets[10]-> Fill(NHadJets,wt1);
	h_NbJets[10]->Fill(BTags,wt1);
      }
	  
    if (Pass_MET2)
      {
	h_MET[11] ->Fill(MET,wt1);
	h_Pho_pT[11] ->Fill(bestPhoton.Pt(),wt1);
	h_Pho_eta[11] -> Fill(bestPhoton.Eta(),wt1);
	h_Pho_phi[11] -> Fill(bestPhoton.Phi(),wt1);
	h_NHadJets[11]-> Fill(NHadJets,wt1);
	h_NbJets[11]->Fill(BTags,wt1);
      }
	  
    if (Pass_EMu_veto) {
      h_MET[5] ->Fill(MET,wt1);
      h_Pho_pT[5] ->Fill(bestPhoton.Pt(),wt1);
      h_Pho_eta[5] -> Fill(bestPhoton.Eta(),wt1);
      h_Pho_phi[5] -> Fill(bestPhoton.Phi(),wt1);
      h_NHadJets[5]-> Fill(NHadJets,wt1);
      h_NbJets[5]->Fill(BTags,wt1);
    }
	  
    if (Pass_Iso_trk_veto){
      h_MET[6] ->Fill(MET,wt1);
      h_Pho_pT[6] ->Fill(bestPhoton.Pt(),wt1);
      h_Pho_eta[6] -> Fill(bestPhoton.Eta(),wt1);
      h_Pho_phi[6] -> Fill(bestPhoton.Phi(),wt1);
      h_NHadJets[6] -> Fill(NHadJets,wt1);
      h_NbJets[6]->Fill(BTags,wt1);
    }

//     //h_Electrons_size_NElectrons->Fill(Electrons->size(),NElectrons);        
	 

    //for stitching part
    double minDR_pho_gen_lep, minDR_pho_qg;
    minDR_pho_gen_lep = getGenLep(bestPhoton);
    minDR_pho_qg = madMinPhotonDeltaR;

    if (JetMetPhi){
      h_mindR_pho_gen_lep_Ovrlp->Fill(minDR_pho_gen_lep,wt1);
      h_mindR_pho_qg_Ovrlp->Fill(minDR_pho_qg,wt1);	  
    }
      
    if (rmOvrlp){
      h_mindR_pho_gen_lep_rmOvrlp->Fill(minDR_pho_gen_lep,wt1);
      h_mindR_pho_qg_rmOvrlp->Fill(minDR_pho_qg,wt1);
    }
      
    //putting conditions for hasgenPromptPhotons
    if(hasGenPromptPhoton){
      if (JetMetPhi){
	h_mindR_pho_gen_lep_Ovrlp_genPromptPho->Fill(minDR_pho_gen_lep,wt1);
	h_mindR_pho_qg_Ovrlp_genPromptPho->Fill(minDR_pho_qg,wt1);	  
      }
	
      if (rmOvrlp){
	h_mindR_pho_gen_lep_rmOvrlp_genPromptPho->Fill(minDR_pho_gen_lep,wt1);
	h_mindR_pho_qg_rmOvrlp_genPromptPho->Fill(minDR_pho_qg,wt1);
      }
    }
      
//     // //sorting Gen_e vector
// //     // TLorentzVector GenE_LV;
// //     // vector<TLorentzVector> GenElectrons_v1;
// //     // for (int i = 0; i < GenElectrons->size(); i++){
// //     //   GenE_LV.SetPtEtaPhiE(GenElectrons[i].Pt(), GenElectrons[i].Eta(), GenElectrons[i].Phi(), GenElectrons[i].E());
// //     //   GenElectrons_v1.push_back(GenE_LV);
// //     // }
// //     // sortTLorVec(&GenElectrons_v1);


//==================================================================FR SR=============================================================================
    bool elec_CR= false; 
    bool pho_SR= false;
    double mt_ele = 0.0;
    
    bool fakeElectron=false, proc=false, proc1=false;
    if(Pass_MET2 && hasPho && bestEMobj && !bestEMObjIsEle){
      //h_selectBaselineYields_v1->Fill("Enter Pho SR",wt);
      if(isoMuonTracks !=0 || isoPionTracks!=0 || isoElectronTracks!=0) continue; // veto muon/pions/electrons 
      if(NMuons!=0 || NElectrons!=0) continue;
      if(nGenMu1==0 && nGenEle1==0 && v_genTau2.size()==0) continue; //reject w->qq events or all hadronic decay of TTbar                                                    
      if(nGenMu1>0) continue;
      if(nGenEle1==0) continue;
      if(((s_sample.Contains("WGJets") || s_sample.Contains("WJets"))&& nGentau_had1>0) ||( (s_sample.Contains("TTGJets") || s_sample.Contains("TTJets")) && nGentau_had1>1)) continue;
      if(nGenEle1==0) continue;
      
      proc = v_genEle2.size()==1;
      proc1 = v_genEle2.size()>1;
      fakeElectron =(proc && (bestPhoton_v1.DeltaR(v_genEle2[0])<0.2 && bestPhoton_v1.Pt()/v_genEle2[0].Pt()<1.2 && bestPhoton_v1.Pt()/v_genEle2[0].Pt()>0.8))||(proc1 && ((bestPhoton_v1.DeltaR(v_genEle2[0])<0.2 && bestPhoton_v1.Pt()/v_genEle2[0].Pt()<1.2 && bestPhoton_v1.Pt()/v_genEle2[0].Pt()>0.8)|| (v_genEle2.size()>1 && (bestPhoton_v1.DeltaR(v_genEle2[1])<0.2 && bestPhoton_v1.Pt()/v_genEle2[1].Pt()<1.2 && bestPhoton_v1.Pt()/v_genEle2[1].Pt()>0.8))));
      
      if(fakeElectron) pho_SR= true;
      else continue;
      //h_selectBaselineYields_v1->Fill("After Pho SR",wt);

      if(bestEMObj.Pt()<=40) continue;
      if(pho_SR){
	 h_FR_SR_Pho_Pt->Fill(bestEMObj.Pt(),wt1);
	 h_FR_SR_MET->Fill(MET,wt1);
	 h_FR_SR_NHadJets->Fill(NHadJets,wt1);
	 h_FR_SR_NbJets->Fill(BTags,wt1);
	 h_FR_SR_QMult->Fill(qmulti,wt1);
      }
    }
//====================================================================================================================================================

    
  

    

//===========================================================FR CR====================================================================================
    
    if(Pass_MET2 && bestEMobj && bestEMObjIsEle && !hasPho){
      //h_selectBaselineYields_v1->Fill("Enter Elec CR",wt);
      if(isoMuonTracks !=0 || isoPionTracks !=0) continue; // veto muon/pions from 1 electron CR
      if(NMuons!=0 || NElectrons>1) continue;     
      mt_ele = sqrt(2*bestEMObj.Pt()*MET*(1-cos(DeltaPhi(METPhi,bestEMObj.Phi()))));
      if(mt_ele>100) { continue;}//h_selectBaselineYields_CR->Fill("e-CR: mT<100",wt);continue;} // remove signal contamination 
      elec_CR = true;
      // h_selectBaselineYields_v1->Fill("After Elec CR",wt);
      // int searchBin = getBinNoV6_WithOnlyBLSelec(nHadJets,BTags);
      // if(nHadJets>=7 && Debug)
      // 	cout<<"Entry: "<<jentry<<"\t"<<searchBin<<"\t"<<nHadJets<<"\t"<<BTags<<"\t"<<MET<<"\t"<<bestPhoton.Pt()<<"\t"<<wt<<endl;
      //if (bestEMObj.Pt()<=40) continue;

      if (elec_CR){
	 h_FR_CR_Pho_Pt->Fill(bestEMObj.Pt(),wt1);
	 h_FR_CR_MET->Fill(MET,wt1);
	 h_FR_CR_NHadJets->Fill(NHadJets,wt1);
	 h_FR_CR_NbJets->Fill(BTags,wt1);
	 h_FR_CR_QMult->Fill(qmulti,wt1);
      }
    }
//====================================================================================================================================================    

    

        

//     if (lost_elec_SR) {
//       h_Lost_e_SR_Pho_Pt->Fill(bestPhoton.Pt(),wt1);
//       h_Lost_e_SR_MET->Fill(MET,wt1);
//       h_Lost_e_SR_NHadJets->Fill(NHadJets,wt1);
//       h_Lost_e_SR_NbJets->Fill(BTags,wt1);
	
//       h_Lost_e_SR_binned ->Fill(TFbins,wt1);
//       h_Lost_e_SR_srch_binned -> Fill(SrchBins,wt1);
//     }
      
//     if (lost_mu_SR){
//       h_Lost_mu_SR_Pho_Pt->Fill(bestPhoton.Pt(),wt1);
//       h_Lost_mu_SR_MET->Fill(MET,wt1);
//       h_Lost_mu_SR_NHadJets->Fill(NHadJets,wt1);
//       h_Lost_mu_SR_NbJets->Fill(BTags,wt1);
	
//       h_Lost_mu_SR_binned ->Fill(TFbins,wt1);
//       h_Lost_mu_SR_srch_binned -> Fill(SrchBins,wt1);
//     }

//     if (lost_elec_SR || lost_mu_SR) {
//       h_LL_SR_Pho_Pt ->Fill(bestPhoton.Pt(),wt1);
//       h_LL_SR_MET->Fill(MET,wt1);
//       h_LL_SR_NHadJets->Fill(NHadJets,wt1);
//       h_LL_SR_NbJets->Fill(BTags,wt1);
	 
//       h_LL_SR_binned ->Fill(TFbins,wt1);
//       h_LL_SR_srch_binned -> Fill(SrchBins,wt1);
//     }
                  
//     if (FR_SR) {
//       h_FR_SR_binned -> Fill(TFbins,wt1);
//       h_FR_SR_NHadJets -> Fill(NHadJets,wt1);
//     }
      
//     if (rest_SR) {
	
//       // cout << "Event: " << jentry << endl;
//       // for(Long64_t ii=0; ii<GenParticles->size(); ii++) {
//       //   cout << "PdgId: " << GenParticles_PdgId[(int)ii] << ", ParenId: " << GenParticles_ParentId[(int)ii] << ", status: " << GenParticles_Status[(int)ii] << endl;
//       //   if (abs(GenParticles_ParentId[(int)ii] == 24)) h_ptcl_W_rest->Fill(GenParticles_PdgId[(int)ii]);    
//       // }
//       // cout << endl;

//       h_Pho_Pt_rest ->Fill(bestPhoton.Pt(),wt1);
//       h_MET_rest ->Fill(MET,wt1);
//       h_rest_NHadJets ->Fill(NHadJets,wt1);
//       h_NbJets_rest ->Fill(BTags,wt1);
//       h_dhi_JetMET1_rest ->Fill(dPhi_METjet1,wt1);
//       h_dhi_JetMET2_rest ->Fill(dPhi_METjet2,wt1);		
//     }

     
//     if (lost_elec_CR) {
//       h_Lost_e_CR_Pho_Pt ->Fill(bestPhoton.Pt(),wt1);
//       h_Lost_e_CR_MET->Fill(MET,wt1);
//       h_Lost_e_CR_NHadJets->Fill(NHadJets,wt1);
//       h_Lost_e_CR_NbJets->Fill(BTags,wt1);
	
//       h_Lost_e_CR_binned ->Fill(TFbins,wt1);
//       h_Lost_e_CR_srch_binned -> Fill(SrchBins,wt1);
//     }

//     if (lost_mu_CR) {
//       h_Lost_mu_CR_Pho_Pt ->Fill(bestPhoton.Pt(),wt1);
//       h_Lost_mu_CR_MET->Fill(MET,wt1);
//       h_Lost_mu_CR_NHadJets->Fill(NHadJets,wt1);
//       h_Lost_mu_CR_NbJets->Fill(BTags,wt1);
	
//       h_Lost_mu_CR_binned ->Fill(TFbins,wt1);
//       h_Lost_mu_CR_srch_binned -> Fill(SrchBins,wt1);
//     }

//     if (lost_elec_CR || lost_mu_CR) {	
//       h_LL_CR_binned ->Fill(TFbins,wt1);
//       h_LL_CR_srch_binned -> Fill(SrchBins,wt1);
//     }
      
//     //validating TF Method
//     if (lost_elec_CR){
// 	TF1 = h_TF1 ->GetBinContent(TFbins+1);
// 	wt2 = TF1*wt1;
// 	h_Lost_e_SR_srch_binned_pred->Fill(SrchBins,wt2);
//     }
//     //cout << "\n\n";
//     if (lost_mu_CR){
// 	TF2 = h_TF2 ->GetBinContent(TFbins+1);
// 	wt3 = TF2*wt1;
// 	h_Lost_mu_SR_srch_binned_pred->Fill(SrchBins,wt3);
//     }

//     if (lost_elec_CR || lost_mu_CR){
// 	TF3 = h_TF3 ->GetBinContent(TFbins+1);
// 	wt4 = TF3*wt1;
// 	h_LL_SR_srch_binned_pred->Fill(SrchBins,wt4);
//     }
	
    //     //================================================================================================================
    //     //testing
    //     //if (GenElectrons->size() == 0 && GenMuons->size() ==0 && GenTaus->size() > 1) h_Gen_taus_size->Fill(GenTaus_had->size());
      
      
      
      
    //     // if (Pass_Iso_trk_veto && GenElectrons->size()>0){        
    //     // 	double dR_gen_e_gamma = DeltaR(GenElectrons[0].Eta(),GenElectrons[0].Phi(),bestPhoton.Eta(),bestPhoton.Phi());
    //     // 	h_dR_gen_e_reco_pho ->Fill(dR_gen_e_gamma);
    //     // 	h_gen_e_reco_pho_ratio ->Fill(GenElectrons[0].Pt()/bestPhoton.Pt());
    //     // 	h_dRvsRatio->Fill(dR_gen_e_gamma,GenElectrons[0].Pt()/bestPhoton.Pt());
    //     // }
      
    //     //==============================================================================================================
      
      
    //     // for(int i=0;i<Jets->size();i++){
    //     // 	if( (Jets[i].Pt() > 30.0) && (abs(Jets[i].Eta()) <= 2.4) ){
    //     // 	  if (!(minDR < 0.3 && i==minDRindx)){
    //     // 	    if (Jets_ID[i]){
    //     // 	      Jets_pT_Sum += Jets[i].Pt();
    //     // 	      h_Jet_pT[0]->Fill(Jets[i].Pt());
    //     // 	      h_Jet_eta[0]->Fill(Jets[i].Eta());
    //     // 	      h_Jet_phi[0]->Fill(Jets[i].Phi());
	      
    //     // 	      if (MET>200){
    //     // 		h_Jet_pT[1]->Fill(Jets[i].Pt());
    //     // 		h_Jet_eta[1]->Fill(Jets[i].Eta());
    //     // 		h_Jet_phi[1]->Fill(Jets[i].Phi());
		
    //     // 		if (bestPhoton.Pt()<20) continue;
    //     // 		h_Jet_pT[2]->Fill(Jets[i].Pt());
    //     // 		h_Jet_eta[2]->Fill(Jets[i].Eta());
    //     // 		h_Jet_phi[2]->Fill(Jets[i].Phi());
		
    //     // 		if (NHadJets < 2) continue;
    //     // 		h_Jet_pT[3]->Fill(Jets[i].Pt());
    //     // 		h_Jet_eta[3]->Fill(Jets[i].Eta());
    //     // 		h_Jet_phi[3]->Fill(Jets[i].Phi());
		
    //     // 		if (!(NEMu == 0)) continue;
    //     // 		h_Jet_pT[5]->Fill(Jets[i].Pt());
    //     // 		h_Jet_eta[5]->Fill(Jets[i].Eta());
    //     // 		h_Jet_phi[5]->Fill(Jets[i].Phi());
		
    //     // 		if (isoElectronTracks || isoMuonTracks || isoPionTracks) continue;
    //     // 		h_Jet_pT[6]->Fill(Jets[i].Pt());
    //     // 		h_Jet_eta[6]->Fill(Jets[i].Eta());
    //     // 		h_Jet_phi[6]->Fill(Jets[i].Phi());
		
    //     // 		}
    //     // 	      }
    //     // 	    }
    //     // 	  }
    //     // }
      
    //     // 	//if(hadJets.size()==0) continue;

    //     // 	if(Debug)
    //     // 	  cout<<"===load tree entry ===  "<<"\t"<<jentry<<"\t"<<"No of B-Jets ===  "<<bjets.size()<<endl;

    //     // 	//defining flags for different categories of detector inefficiencies
    //     // 	bool LostMu_flag=false, LostE_flag=false, EfakePho_flag=false, hadTau_flag=false, Rest_flag=false;
    //     // 	if (GenMuons->size() > 0 && NMuons == 0) LostMu_flag = true;	    
    //     // 	else if(GenElectrons -> size() > 0 && NElectrons == 0 && bestPhotonIndxAmongPhotons > 0)
    //     // 	  { double dR = DeltaR(bestPhoton.Eta(),bestPhoton.Phi(),GenElectrons[0].Eta(),GenElectrons[0].Phi());
    //     // 	    if (dR > 0.1) LostE_flag = true;	      
    //     // 	    else EfakePho_flag = true;	     
    //     // 	  }
    //     // 	else if(GenElectrons -> size() > 0 && NElectrons == 0 && bestPhotonIndxAmongPhotons < 0) LostE_flag = true;	
    //     // 	else if (GenTaus->size() > 0 && GenTaus_had[0]) hadTau_flag = true;
    //     // 	else if (GenElectrons->size()>0 || GenMuons->size()>0 || GenTaus->size()>0) Rest_flag = true;
	
    //     // 	//filling histo to check stack plot:
    //     // 	if (MET>200 && bestPhoton.Pt()>40 && NHadJets>=2 && ST>300 && NEMu==0)
    //     // 	  {
    //     // 	    //filling the histogram for total
    //     // 	    if(GenElectrons->size()>0 && GenMuons->size()==0 && GenTaus->size()==0) h_Gen_eta[0][6]->Fill(GenElectrons[0].Eta());
    //     // 	    if(GenElectrons->size()==0 && GenMuons->size()>0 && GenTaus->size()==0) h_Gen_eta[0][6]->Fill(GenMuons[0].Eta());
    //     // 	    if(GenElectrons->size()==0 && GenMuons->size()==0 && GenTaus->size()>0)h_Gen_eta[0][6]->Fill(GenTaus[0].Eta());
    //     // 	    if(GenElectrons->size()>0 && GenMuons->size()==0 && GenTaus->size()>0)h_Gen_eta[0][6]->Fill(GenElectrons[0].Eta());
    //     // 	    if(GenElectrons->size()==0 && GenMuons->size()>0 && GenTaus->size()>0)h_Gen_eta[0][6]->Fill(GenMuons[0].Eta());

	    
    //     // 	    //filling the histogrms for the different categories
    //     // 	    if (LostE_flag) h_LostElectron_eta -> Fill(GenElectrons[0].Eta());
    //     // 	    if (EfakePho_flag) h_EFakePho_eta-> Fill(GenElectrons[0].Eta());
    //     // 	    if (LostMu_flag) h_LostMuon_eta -> Fill(GenMuons[0].Eta());
    //     // 	    if (hadTau_flag) h_HadTau_eta -> Fill(GenTaus[0].Eta());  
    //     // 	    if (Rest_flag) h_Rest_eta-> FillRandom("gaus",1);
    //     // 	    }

	  
	     
    //     // 	// Loop for plotting the Gen information
    //     // 	for(Long64_t ii=0; ii<GenParticles->size(); ii++) {
    //     // 	  if (abs(GenParticles_PdgId[(int)ii])==24 && abs(GenParticles_ParentId[(int)ii])==13) cout << "W from mu" << endl;
    //     // 	}   //end genparticle loop
	
    //     // 	for (Long64_t ii =0; ii<GenElectrons->size(); ii++){
    //     // 	  h_Gen_pT[0][0]->Fill(GenElectrons[(int)ii].Pt(),wt1);
    //     // 	  h_Gen_eta[0][0]->Fill(GenElectrons[(int)ii].Eta(),wt1);
    //     // 	  h_Gen_phi[0][0]->Fill(GenElectrons[(int)ii].Phi(),wt1);
	  
	  
	  
    //     // 	  if (ST<300) continue;
    //     // 	  h_Gen_pT[0][4]->Fill(GenElectrons[(int)ii].Pt(),wt1);  
    //     // 	  h_Gen_eta[0][4]->Fill(GenElectrons[(int)ii].Eta(),wt1); 	    
    //     // 	  h_Gen_phi[0][4]->Fill(GenElectrons[(int)ii].Phi(),wt1); 
	  
	  
    //     // 	  if (NEMu!=0) continue;
    //     // 	  h_Gen_pT[0][5]->Fill(GenElectrons[(int)ii].Pt(),wt1); 	    
    //     // 	  h_Gen_eta[0][5]->Fill(GenElectrons[(int)ii].Eta(),wt1); 
    //     // 	  h_Gen_phi[0][5]->Fill(GenElectrons[(int)ii].Phi(),wt1); 
	  

    //     // 	  if (isoElectronTracks != 0 || isoMuonTracks != 0 || isoPionTracks != 0) continue;
    //     // 	  h_Gen_pT[0][6]->Fill(GenElectrons[(int)ii].Pt(),wt1);  
    //     // 	 //h_Gen_eta[0][6]->Fill(GenElectrons[(int)ii].Eta(),wt1); 	    
    //     // 	  h_Gen_phi[0][6]->Fill(GenElectrons[(int)ii].Phi(),wt1);
    //     // 	}  //end Gen electron loop

	
    //     // 	for (Long64_t ii =0; ii<GenMuons->size(); ii++){
    //     // 	  h_Gen_pT[1][0]->Fill(GenMuons[(int)ii].Pt(),wt1);
    //     // 	  h_Gen_eta[1][0]->Fill(GenMuons[(int)ii].Eta(),wt1);
    //     // 	  h_Gen_phi[1][0]->Fill(GenMuons[(int)ii].Phi(),wt1);
	  
	  
    //     // 	  if (ST<300) continue;
    //     // 	  h_Gen_pT[1][4]->Fill(GenMuons[(int)ii].Pt(),wt1);  
    //     // 	  h_Gen_eta[1][4]->Fill(GenMuons[(int)ii].Eta(),wt1); 	    
    //     // 	  h_Gen_phi[1][4]->Fill(GenMuons[(int)ii].Phi(),wt1); 
	  
	  
    //     // 	  if (NEMu!=0) continue;
    //     // 	  h_Gen_pT[1][5]->Fill(GenMuons[(int)ii].Pt(),wt1); 	    
    //     // 	  h_Gen_eta[1][5]->Fill(GenMuons[(int)ii].Eta(),wt1); 
    //     // 	  h_Gen_phi[1][5]->Fill(GenMuons[(int)ii].Phi(),wt1); 
	  
	  
    //     // 	  if (isoElectronTracks || isoMuonTracks || isoPionTracks) continue;
    //     // 	  h_Gen_pT[1][6]->Fill(GenMuons[(int)ii].Pt(),wt1);  
    //     // 	  h_Gen_eta[1][6]->Fill(GenMuons[(int)ii].Eta(),wt1); 	    
    //     // 	  h_Gen_phi[1][6]->Fill(GenMuons[(int)ii].Phi(),wt1);
    //     // 	} //end Gen Muon loop
	
    //     // 	for (Long64_t ii =0; ii<GenTaus->size(); ii++){
    //     // 	  h_Gen_pT[2][0]->Fill(GenTaus[(int)ii].Pt(),wt1);
    //     // 	  h_Gen_eta[2][0]->Fill(GenTaus[(int)ii].Eta(),wt1);
    //     // 	  h_Gen_phi[2][0]->Fill(GenTaus[(int)ii].Phi(),wt1);
	  
    //     // 	  if (ST<300) continue;
    //     // 	  h_Gen_pT[2][4]->Fill(GenTaus[(int)ii].Pt(),wt1);  
    //     // 	  h_Gen_eta[2][4]->Fill(GenTaus[(int)ii].Eta(),wt1); 	    
    //     // 	  h_Gen_phi[2][4]->Fill(GenTaus[(int)ii].Phi(),wt1); 
	  
    //     // 	  if (NEMu!=0) continue;
    //     // 	  h_Gen_pT[2][5]->Fill(GenTaus[(int)ii].Pt(),wt1); 	    
    //     // 	  h_Gen_eta[2][5]->Fill(GenTaus[(int)ii].Eta(),wt1); 
    //     // 	  h_Gen_phi[2][5]->Fill(GenTaus[(int)ii].Phi(),wt1); 
	  
    //     // 	  if (isoElectronTracks != 0 || isoMuonTracks != 0 || isoPionTracks != 0) continue;
    //     // 	  h_Gen_pT[2][6]->Fill(GenTaus[(int)ii].Pt(),wt1);  
    //     // 	  h_Gen_eta[2][6]->Fill(GenTaus[(int)ii].Eta(),wt1); 	    
    //     // 	  h_Gen_phi[2][6]->Fill(GenTaus[(int)ii].Phi(),wt1);	  
    //     // 	} //end Gen Tau loop
	
    //     // 	for(Long64_t ii=0; ii<Electrons->size(); ii++){
    //     // 	  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > myele = Electrons[(int)ii];
    //     // 	  h_Reco_pT[0][0]->Fill(Electrons[(int)ii].Pt(),wt1);
    //     // 	  h_Reco_eta[0][0]->Fill(Electrons[(int)ii].Eta(),wt1);
    //     // 	  h_Reco_phi[0][0]->Fill(Electrons[(int)ii].Phi(),wt1);
	  	  
    //     // 	  if (ST<300) continue;
    //     // 	  h_Reco_pT[0][4]->Fill(Electrons[(int)ii].Pt(),wt1);
    //     // 	  h_Reco_eta[0][4]->Fill(Electrons[(int)ii].Eta(),wt1);
    //     // 	  h_Reco_phi[0][4]->Fill(Electrons[(int)ii].Phi(),wt1);
	  
    //     // 	  vector<myLV> v_electron;
    //     // 	  if (Electrons_passIso) v_electron.push_back(Electrons[(int)ii]);
    //     // 	  if (v_electron.size()!=0) continue;
    //     // 	  h_Reco_pT[0][5]->Fill(v_electron[0].Pt(),wt1);
    //     // 	  h_Reco_eta[0][5]->Fill(v_electron[0].Eta(),wt1);
    //     // 	  h_Reco_phi[0][5]->Fill(v_electron[0].Phi(),wt1);
	    
    //     // 	  if (isoElectronTracks != 0 || isoMuonTracks != 0 || isoPionTracks != 0) continue;
    //     // 	  h_Reco_pT[0][6]->Fill(Electrons[(int)ii].Pt(),wt1);
    //     // 	  //h_Reco_eta[0][6]->Fill(Electrons[(int)ii].Eta(),wt1);
    //     // 	  h_Reco_phi[0][6]->Fill(Electrons[(int)ii].Phi(),wt1);	    
    //     // 	} //end electron loop
    //     // 	//if (jentry<1000)cout << endl;

    //     // 	for(Long64_t ii=0; ii<Muons->size(); ii++){
    //     // 	    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > mymu = Muons[(int)ii];
    //     // 	    h_Reco_pT[1][0]->Fill(Muons[(int)ii].Pt(),wt1);
    //     // 	    h_Reco_eta[1][0]->Fill(Muons[(int)ii].Eta(),wt1);
    //     // 	    h_Reco_phi[1][0]->Fill(Muons[(int)ii].Phi(),wt1);

    //     // 	    if (ST<300) continue;
    //     // 	    h_Reco_pT[1][4]->Fill(Muons[(int)ii].Pt(),wt1);
    //     // 	    h_Reco_eta[1][4]->Fill(Muons[(int)ii].Eta(),wt1);
    //     // 	    h_Reco_phi[1][4]->Fill(Muons[(int)ii].Phi(),wt1);

    //     // 	    vector<myLV> v_muon;
    //     // 	    if (Muons_passIso) v_muon.push_back(Muons[(int)ii]);
    //     // 	    if (v_muon.size()!=0) continue;
    //     // 	    h_Reco_pT[1][5]->Fill(v_muon[0].Pt(),wt1);
    //     // 	    h_Reco_eta[1][5]->Fill(v_muon[0].Eta(),wt1);
    //     // 	    h_Reco_phi[1][5]->Fill(v_muon[0].Phi(),wt1);

    //     // 	    if (isoElectronTracks != 0 || isoMuonTracks != 0 || isoPionTracks != 0) continue;
    //     // 	    h_Reco_pT[1][6]->Fill(Muons[(int)ii].Pt(),wt1);
    //     // 	    h_Reco_eta[1][6]->Fill(Muons[(int)ii].Eta(),wt1);
    //     // 	    h_Reco_phi[1][6]->Fill(Muons[(int)ii].Phi(),wt1);
    //     // 	} //end Muon loop
	

    //     // 	for(Long64_t ii=0; ii<Photons->size(); ii++){
    //     // 	  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > mypho = Photons[(int)ii];
    //     // 	} //end photon loop
	

    } // end jentry loop
  
  // //combining lost e and lost_mu tau_had
  // //hists for SR
  // h_LL_SR_Pho_Pt ->Add(h_Lost_e_SR_Pho_Pt, h_Lost_mu_SR_Pho_Pt);
  // h_LL_SR_MET ->Add(h_Lost_e_SR_MET, h_Lost_mu_SR_MET);
  // h_LL_SR_NHadJets ->Add(h_Lost_e_SR_NHadJets, h_Lost_mu_SR_NHadJets);
  // h_LL_SR_NbJets ->Add(h_Lost_e_SR_NbJets, h_Lost_mu_SR_NbJets);

  // //hists for CR
  // h_LL_CR_Pho_Pt ->Add(h_Lost_e_CR_Pho_Pt, h_Lost_mu_CR_Pho_Pt);
  // h_LL_CR_MET ->Add(h_Lost_e_CR_MET, h_Lost_mu_CR_MET);
  // h_LL_CR_NHadJets ->Add(h_Lost_e_CR_NHadJets, h_Lost_mu_CR_NHadJets);
  // h_LL_CR_NbJets ->Add(h_Lost_e_CR_NbJets, h_Lost_mu_CR_NbJets);  

  cout << "Cutflows for the year: " << s_data << endl;
  cout << h_NHadJets[0]->Integral() << endl;  
  cout << h_NHadJets[2]->Integral() << endl;
  cout << h_NHadJets[1]->Integral() << endl;
  cout << h_NHadJets[3]->Integral() << endl; 
  cout << h_NHadJets[4]->Integral() << endl;
  cout << h_NHadJets[7]->Integral()<< endl;
  cout << h_NHadJets[8]->Integral()<< endl;
  cout << h_NHadJets[9]->Integral()<< endl;
  cout << h_NHadJets[10]->Integral()<< endl;
  cout << h_NHadJets[11]->Integral()<< endl;
  cout << h_NHadJets[5]->Integral() << endl;
  cout << h_NHadJets[6]->Integral() << "\n\n";
  
  // cout << "Lost electron: " << h_Lost_e_SR_NHadJets->Integral() << endl;
  // cout << "Lost muon: " << h_Lost_mu_SR_NHadJets->Integral() << endl;
  // cout << "Hadronic tau: " << h_Had_tau_SR_NHadJets->Integral() << endl;
  // cout << "e fake photon: " << h_FR_SR_NHadJets->Integral() << endl;
  // cout << "Rest: " << h_rest_NHadJets->Integral() << endl;

  
}// End Eventloop



myLV AnalyzeTProxytBSM::getBestPhoton(int pho_ID){
  vector<myLV> goodPho;
  vector<int> goodPhoIndx;
  for(int iPho=0;iPho<Photons->size();iPho++){
    if(((*Photons_hasPixelSeed)[iPho]<0.001) && ( (*Photons_fullID)[iPho]))
      //if(abs(Photons[iPho].Eta())<2.4 && (*Photons_hasPixelSeed)[iPho]<0.001 && ( ((*Photons_fullID)[iPho] && pho_ID==0) || (pho_ID==1 &&(((*Photons_cutBasedID)[iPho]==1 || (*Photons_cutBasedID)[iPho]==2))) || (pho_ID==2 && (*Photons_cutBasedID)[iPho]==2) || (pho_ID==3 && (*Photons_mvaValuesID)[iPho]>-0.02) || (pho_ID==4 && (*Photons_mvaValuesID)[iPho]>0.42)))
      {
	goodPho.push_back(Photons[iPho]);
	goodPhoIndx.push_back(iPho);
      }
  }
  
  int highPtIndx=-100;
  for(int i=0;i<goodPho.size();i++){
    if(i==0) highPtIndx=0;
    else if( (goodPho[highPtIndx].Pt()) < (goodPho[i].Pt()) ){highPtIndx=i;}
  }
   
  if(highPtIndx>=0){
    bestPhotonIndxAmongPhotons = goodPhoIndx[highPtIndx];
  }
  else bestPhotonIndxAmongPhotons = -100;
  if(highPtIndx==-100){myLV v0;return v0;}
  else return goodPho[highPtIndx];
   
}

//double AnalyzeLightBSM::getGenLep(TLorentzVector bestPhoton){
double AnalyzeTProxytBSM::getGenLep(myLV bestPhoton){
  //vector<TLorentzVector> v_genLep2;
  vector<myLV> v_genLep2;
  //TLorentzVector genMu1, genEle1;
  myLV genMu1, genEle1;
  // if(flag)
  //   {
  for(int i=0 ; i < GenElectrons->size(); i++)
    {
      if(GenElectrons[i].Pt()!=0)
	{
	  genEle1 = (GenElectrons[i]);
	  v_genLep2.push_back(genEle1);
	}
      
    }
  //   }
  // else
  //   {
  for(int i=0 ; i < GenMuons->size(); i++)
    {
      if(GenMuons[i].Pt()!=0)
	{
	  genMu1 = (GenMuons[i]);
	  v_genLep2.push_back(genMu1);
	}
    }
  //  }
  return MinDr_myLV(bestPhoton,v_genLep2);
}

bool AnalyzeTProxytBSM::RemoveSampleOverlap(TString s_sample, myLV bestPhoton) {

  //script to define conditions to remove ovrlap
  bool genphocheck = false;
  int genphomatch_before=0, genphomatch_after=0; 
  bool cont1=true, cont2=true, cont3=true, cont4=true, cont5=true,
    cont6=true, cont7=true, cont8=true, cont9=true, cont10=true; 

  if((s_sample.Contains("TTJets_HT") || s_sample.Contains("TTJets-HT")) && madHT<600)
    cont1=false;

  if((s_sample.Contains("TTJets_inc")|| s_sample.Contains("TTJets_SingleLept") ||
      s_sample.Contains("TTJets_DiLept") || s_sample.Contains("TTJets_Leptons") ||
      s_sample.Contains("TTJets_Leptons")) && madHT>600)
    cont2=false;

  if(!genphocheck) {
    genphomatch_before++;
    double mindr_Pho_genlep=getGenLep(bestPhoton);
    
    if( s_sample.Contains("TTG") ) {
      if(!hasGenPromptPhoton) {
	//h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt1);
	//if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
      } else if(hasGenPromptPhoton) {
	//h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt1);
	if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )) {
	  //h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt1);
	  //if(madMinPhotonDeltaR<0.5)h_selectBaselineYields_v1->Fill("madMinPhotonDeltaR <0.5",wt1);
	  //if(mindr_Pho_genlep<0.5) h_selectBaselineYields_v1->Fill("mindr_Pho_genlep<0.5",wt1);
	  cont3=false;
	} else {
	  //if(madMinPhotonDeltaR >= 0.5) h_selectBaselineYields_v1->Fill("mindR(q/g, #gamma)",wt1);
	  //if(mindr_Pho_genlep >=0.5)    h_selectBaselineYields_v1->Fill("mindR(l, #gamma)",wt1);
	}
      }
    }

    if(s_sample.Contains("WGJets_MonoPhoton_PtG-40to130UL") ||
       s_sample.Contains("WGJets_MonoPhoton_PtG-130UL")) {
      //if(s_sample.Contains("WGJets_MonoPhoton_PtG-40to130UL"||"WGJets_MonoPhoton_PtG-130UL"))
      if(!hasGenPromptPhoton) {
	//h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt1);
	//if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
      } else if(hasGenPromptPhoton) {
	//h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt1);
	if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 ))
	  {//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt1);
	    //if(madMinPhotonDeltaR<0.5) h_selectBaselineYields_v1->Fill("madMinPhotonDeltaR <0.5",wt1);
	    //if(mindr_Pho_genlep<0.5)   h_selectBaselineYields_v1->Fill("mindr_Pho_genlep<0.5",wt1);
	    cont4=false;
	  } else {
	  //if(madMinPhotonDeltaR >= 0.5) h_selectBaselineYields_v1->Fill("mindR(q/g, #gamma)",wt1);
	  //if(mindr_Pho_genlep >=0.5)    h_selectBaselineYields_v1->Fill("mindR(l, #gamma)",wt1);
	}
      }
    }
    
    if(s_sample.Contains("WJets")) {
      if(!hasGenPromptPhoton) {
	//h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt1);
	//if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl; 
      } else if(hasGenPromptPhoton) {
	//h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt1);
	if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)) {
	  //h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt1);
	  cont5=false;
	} else {
	  //if(madMinPhotonDeltaR >= 0.5) h_selectBaselineYields_v1->Fill("pass_mindR(q/g, #gamma)",wt1); 
	  //if(mindr_Pho_genlep >=0.5)    h_selectBaselineYields_v1->Fill("pass_mindR(l, #gamma)",wt1); 
	}
      }
    }
    
    if(s_sample.Contains("TTJets_HT") || s_sample.Contains("TTJets-HT")||
       s_sample.Contains("TTJets-inc")|| s_sample.Contains("TTJets_inc") ||
       s_sample.Contains("TTJets2_v17")||s_sample.Contains("TTJets")  ||
       s_sample.Contains("TTJets_Leptons")) {
      if(hasGenPromptPhoton) {	
	if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)) {
	  cont6=false;
	}
      }
    }
	
    if(hasGenPromptPhoton && (s_sample.Contains("GJets"))) {
      if(!(madMinPhotonDeltaR>0.4)) cont7=false;
    }
	
    if(hasGenPromptPhoton && (s_sample.Contains("QCD"))) {
      if((madMinPhotonDeltaR>0.4 && hasGenPromptPhoton))
	cont8=false;
    }
	
    if(hasGenPromptPhoton && ((s_sample.Contains("ZG"))|| (s_sample.Contains("ZNuNuG"))
			      || s_sample.Contains("ZNuNuGJets"))) {
      if(!(madMinPhotonDeltaR>0.5))
	cont9=false;
    }
	
    if(hasGenPromptPhoton && ((s_sample.Contains("ZJets"))|| (s_sample.Contains("ZNuNuJets")))) {
      if(!(madMinPhotonDeltaR<=0.5))
	cont10=false;
    }
    genphomatch_after++;
  }

  bool rmOvrlp = cont1 && cont2 && cont3 && cont4 && cont5 && cont6 && cont7 && cont8 && cont9 && cont10;

  return rmOvrlp;
}

// void AnalyzeLightBSM::FillHistogram_Kinematics(int i, int Njets, int btags, double pho_Pt, double mt_phoMET, double dPhi, double ST,float Eta , float Phi, double E, float METPhi,double qmulti,double leadJets_qmulti, double leadJet_Pt, double leadbjet_tag, double mindR, TLorentzVector Jet_matched, vector<TLorentzVector>hadJets, TLorentzVector leadJet,int nvrtx,double mindr_Pho_genElec, double wt){
//   //cout<<"Alps "<<i<<"\t"<<Njets<<"\t"<<btags<<"\t"<<pho_Pt<<"\t"<<mt_phoMET<<"\t"<<dPhi<<"\t"<<ST<<endl;
  
//   h_Njets[i]->Fill(Njets,wt);
//   h_Nbjets[i]->Fill(btags,wt);
//   h_MET_[i]->Fill(MET,wt);
//   h_PhotonPt[i]->Fill(pho_Pt,wt);
//   h_Mt_PhoMET[i]->Fill(mt_phoMET,wt);
//   h_dPhi_PhoMET[i]->Fill(dPhi,wt);
//   h_St[i]->Fill(ST,wt);
//   h_HT[i]->Fill(HT,wt);
//   h_Photon_Eta[i]->Fill(Eta,wt);
//   h_Photon_Phi[i]->Fill(Phi,wt);
  
//   h_Photon_E[i]->Fill(E,wt);
//   h_MET_Phi[i]->Fill(METPhi,wt);
//   h_qmulti_1[i]->Fill(qmulti,wt);
//   h_qmultiVsEmobjPT[i]->Fill(pho_Pt,qmulti,wt);
//   h_qmultiVsnJets[i]->Fill(Njets,qmulti,wt);
//   h_ST_vs_EMObjPt[i]->Fill( pho_Pt,ST,wt);
//   h_Emobj_PtvsEta[i]->Fill(pho_Pt,Eta,wt);
//   h_Emobj_PtvsPhi[i]->Fill(pho_Pt,Phi,wt);
//   h_Emobj_EtavsPhi[i]->Fill(Eta,Phi,wt);
//   //  cout<<"just after eta vs phi  "<<endl;
//   h_nvrtx[i]->Fill(nvrtx,wt);
//   h_minDR_Jets_EMObject[i]->Fill(mindR,wt);
//   h_Phi_matchedJet[i]->Fill(Jet_matched.Phi(),wt);
//   h_Pt_matchedJet[i]->Fill(Jet_matched.Pt(),wt);
//   h_Eta_matchedJet[i]->Fill(Jet_matched.Eta(),wt);
//   h_EtavsPhi_matchedJet[i]->Fill(Jet_matched.Eta(),Jet_matched.Phi(),wt);
//   h_PtvsEta_matchedJet[i]->Fill(Jet_matched.Pt(),Jet_matched.Eta(),wt);
//   h_PtvsPhi_matchedJet[i]->Fill(Jet_matched.Pt(),Jet_matched.Phi(),wt);
//   h_minDR_Jets_vs_Em_Pt[i]->Fill(Jet_matched.Pt(),mindR,wt);
//   h_btaggervalue_vs_qmulti[i]->Fill(qmulti,leadbjet_tag,wt);
//   h_btaggervalue_vs_minDR_Jets[i]->Fill(mindR,leadbjet_tag,wt);
//   h_minDR_Jets_vsqMulti[i]->Fill(qmulti,mindR,wt);
//   h_HT5HT[i]->Fill(HT5/HT,wt);
//   h_Emobje_pt_vs_Jet_Pt[i]->Fill(pho_Pt,Jet_matched.Pt(),wt);
//   for(int j=0; j<hadJets.size();j++){
//     if (j>=4) continue;
//     h_Phi_leadJet[j][i]->Fill(hadJets[j].Phi(),wt);
//     h_Pt_leadJet[j][i]->Fill(hadJets[j].Pt(),wt);
//     h_Eta_leadJet[j][i]->Fill(hadJets[j].Eta(),wt);
//     h_EtavsPhi_leadJet[j][i]->Fill(hadJets[j].Eta(),hadJets[j].Phi(),wt);
//     h_PtvsPhi_leadJet[j][i]->Fill(hadJets[j].Pt(),hadJets[j].Phi(),wt);
//     h_PtvsEta_leadJet[j][i]->Fill(hadJets[j].Pt(),hadJets[j].Eta(),wt);
//     TLorentzVector Met;
//     Met.SetPtEtaPhiE(MET,0,METPhi,0);
//     double dPhi_METjet = abs(Met.DeltaPhi(hadJets[j])); 
//     h_HT5HT_vsdPhi_METJet[j][i]->Fill(dPhi_METjet,HT5/HT,wt);
//     h_dPhi_METJet[j][i]->Fill(dPhi_METjet,wt);
//     h_dPhi_METJet_vsMET[j][i]->Fill(MET,dPhi_METjet,wt);
//   }

    
//   h_nBjets_vs_qmulti[i]->Fill(btags,qmulti,wt);
//   h_qmultiVs_MET[i]->Fill(MET,qmulti,wt);
//   h_qmultiVs_ST[i]->Fill(ST,qmulti,wt);
//   h_leadJets_qmulti[i]->Fill(leadJets_qmulti,wt);
//   h_leadJet_Pt[i]->Fill(leadJet_Pt,wt);
//   h_leadbjet_tag[i]->Fill(leadbjet_tag,wt);
//   h_leadbjet_tag_vs_leadQmulti[i]->Fill(leadbjet_tag,leadJets_qmulti,wt);
//   h_leadbjet_tag_vsQmulti[i]->Fill(leadbjet_tag,qmulti,wt);
//   h_leadjet_ptvsqmulti[i]->Fill(leadJet_Pt,qmulti,wt);

//   int TFbins_v2 = getBinNoV1_le(pho_Pt, qmulti);
//   int TFbins_v3 = getBinNoV16_le(Njets,btags,pho_Pt);
//   int searchBin = getBinNoV6_WithOnlyBLSelec(Njets,btags);
//   int TFbins_v4 = getBinNoV7_le(Njets,btags);
//   int TFbins_v5 = getBinNo_v0FR(pho_Pt, qmulti, 0.1);
//   int TFbins_v6 = getBinNo_v1FR(pho_Pt,Njets);
//   int TFbins_v7 = getBinNo_v2FR(qmulti,Njets,btags);

//   h_TFbins_LL_v4[i]->Fill(TFbins_v4,wt);
//   h_TFbins_LL_v2[i]->Fill(TFbins_v2,wt);
  
//   h_Sbins_LL[i]->Fill(searchBin,wt);
//   h_TFbins_LL_v3[i]->Fill(TFbins_v3,wt);
//   h_TFbins_LL_v5[i]->Fill(TFbins_v5,wt);
//   h_TFbins_LL_v6[i]->Fill(TFbins_v6,wt);
//   h_TFbins_LL_v7[i]->Fill(TFbins_v7,wt);
  
  
//   searchBin = getBinNoV3_ST_MET_bjets(ST,MET,btags);
//   h_Sbins_LL_newSbins_v3[i]->Fill(searchBin,wt);

//   searchBin = getBinNoV7_ST_MET_bjets_phopT(ST,MET,btags,pho_Pt);
//   h_Sbins_LL_newSbins_v7[i]->Fill(searchBin,wt);
//   searchBin = getBin_ST_MET_bjets_phopT_merge(MET,ST,btags,pho_Pt);
//   h_Sbins_LL_newSbins_v7_merge[i]->Fill(searchBin,wt);
  
//   //cout<<"Alps "<<i<<"\t"<<"TFbins_v1  "<<TFbins_v1<<" TFbins_v2 "<< TFbins_v2<<"\t MET "<<MET<<"\t NJets  "<<Njets<<"  btags "<<btags<<"\t"<<h_TFbins_LL_v1[i]->GetBinContent(TFbins_v1)<<"\t"<<wt<<"\t"<<searchBin<<endl;
//   //cout<<h_Njets[i]->GetMean()<<endl;

// }

int AnalyzeTProxytBSM::getBinNoV1_le( int nHadJets, int nbjets){
  int sBin=-100,m_i=0;
  if(nbjets==0){
    if(nHadJets==2)     { sBin=1;}
    else if(nHadJets==3)     { sBin=2;}
    else if(nHadJets==4)     { sBin=3;}
    else if((nHadJets==5 || nHadJets==6)){ sBin=4;}
    else if(nHadJets>=7)   { sBin=5;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)      { sBin=6;}
    else if((nHadJets==5 || nHadJets==6)){ sBin=7;}
    else if(nHadJets>=7)   { sBin=8;}
  }
  return sBin;
}

int AnalyzeTProxytBSM::getBinNoV6_WithOnlyBLSelec(int nHadJets,int nbjets)
{
  
  int sBin=-100,m_i=0;
  if(nbjets==0 ){
    if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
    else if(nHadJets==5 || nHadJets==6){ sBin=7;}
    else if(nHadJets>=7)               { sBin=13;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)     { sBin=18;}
    else if(nHadJets==5 || nHadJets==6){ sBin=23;}
    else if(nHadJets>=7)               { sBin=28;}
  }
  if(sBin==0){
    for(int i=0;i<METLowEdge.size()-1;i++){
      if(METLowEdge[i]<199.99) continue;
      int sBin1=sBin;
      m_i++;
      if(MET >= METLowEdge[i] && MET < METLowEdge[i+1]){ sBin = sBin+m_i;
	break; }
      else if(MET >= METLowEdge[METLowEdge.size()-1])  { sBin = 7         ;
        break; }
    }
  }
  else if(sBin==7 || sBin==33 || sBin==39){
    int sBin1=sBin;
    for(int i=0;i<METLowEdge_1.size()-1;i++){
      if(METLowEdge_1[i]<199.99) continue;
      m_i++;
      if(MET >= METLowEdge_1[i] && MET < METLowEdge_1[i+1]){ sBin = sBin+m_i;break;}
      else if(MET >= METLowEdge_1[METLowEdge_1.size()-1])  { sBin = sBin+6; break; }
    }
  }

  else 
    {
      for(int i=0;i<METLowEdge_2.size()-1;i++){
	if(METLowEdge_2[i]<199.99) continue;
	m_i++;
	if(MET >= METLowEdge_2[i] && MET < METLowEdge_2[i+1]){ sBin = sBin+m_i;break; }
	else if(MET >= METLowEdge_2[METLowEdge_2.size()-1])  { sBin = sBin+5; break; }
      }
    }
  return sBin;
}

// void AnalyzeLightBSM::FillHistogram_Kinematics(int i, int Njets, int btags, double pho_Pt, float MET, double ST, float METPhi, double qmulti, double wt){
//   //cout<<"Alps "<<i<<"\t"<<Njets<<"\t"<<btags<<"\t"<<pho_Pt<<"\t"<<mt_phoMET<<"\t"<<dPhi<<"\t"<<ST<<endl;
  
//   h_Njets[i]->Fill(Njets,wt);
//   h_Nbjets[i]->Fill(btags,wt);
//   h_MET[i]->Fill(MET,wt);
//   h_PhotonPt[i]->Fill(pho_Pt,wt);
//   h_Mt_PhoMET[i]->Fill(mt_phoMET,wt);
//   h_dPhi_PhoMET[i]->Fill(dPhi,wt);
//   h_St[i]->Fill(ST,wt);
//   h_HT[i]->Fill(HT,wt);
//   h_Photon_Eta[i]->Fill(Eta,wt);
//   h_Photon_Phi[i]->Fill(Phi,wt);
  
//   h_Photon_E[i]->Fill(E,wt);
//   h_MET_Phi[i]->Fill(METPhi,wt);
//   h_qmulti_1[i]->Fill(qmulti,wt);
//   h_qmultiVsEmobjPT[i]->Fill(pho_Pt,qmulti,wt);
//   h_qmultiVsnJets[i]->Fill(Njets,qmulti,wt);
//   h_ST_vs_EMObjPt[i]->Fill( pho_Pt,ST,wt);
//   h_Emobj_PtvsEta[i]->Fill(pho_Pt,Eta,wt);
//   h_Emobj_PtvsPhi[i]->Fill(pho_Pt,Phi,wt);
//   h_Emobj_EtavsPhi[i]->Fill(Eta,Phi,wt);
//   //  cout<<"just after eta vs phi  "<<endl;
//   h_nvrtx[i]->Fill(nvrtx,wt);
//   h_minDR_Jets_EMObject[i]->Fill(mindR,wt);
//   h_Phi_matchedJet[i]->Fill(Jet_matched.Phi(),wt);
//   h_Pt_matchedJet[i]->Fill(Jet_matched.Pt(),wt);
//   h_Eta_matchedJet[i]->Fill(Jet_matched.Eta(),wt);
//   h_EtavsPhi_matchedJet[i]->Fill(Jet_matched.Eta(),Jet_matched.Phi(),wt);
//   h_PtvsEta_matchedJet[i]->Fill(Jet_matched.Pt(),Jet_matched.Eta(),wt);
//   h_PtvsPhi_matchedJet[i]->Fill(Jet_matched.Pt(),Jet_matched.Phi(),wt);
//   h_minDR_Jets_vs_Em_Pt[i]->Fill(Jet_matched.Pt(),mindR,wt);
//   h_btaggervalue_vs_qmulti[i]->Fill(qmulti,leadbjet_tag,wt);
//   h_btaggervalue_vs_minDR_Jets[i]->Fill(mindR,leadbjet_tag,wt);
//   h_minDR_Jets_vsqMulti[i]->Fill(qmulti,mindR,wt);
//   h_HT5HT[i]->Fill(HT5/HT,wt);
//   h_Emobje_pt_vs_Jet_Pt[i]->Fill(pho_Pt,Jet_matched.Pt(),wt);
//   for(int j=0; j<hadJets.size();j++){
//     if (j>=4) continue;
//     h_Phi_leadJet[j][i]->Fill(hadJets[j].Phi(),wt);
//     h_Pt_leadJet[j][i]->Fill(hadJets[j].Pt(),wt);
//     h_Eta_leadJet[j][i]->Fill(hadJets[j].Eta(),wt);
//     h_EtavsPhi_leadJet[j][i]->Fill(hadJets[j].Eta(),hadJets[j].Phi(),wt);
//     h_PtvsPhi_leadJet[j][i]->Fill(hadJets[j].Pt(),hadJets[j].Phi(),wt);
//     h_PtvsEta_leadJet[j][i]->Fill(hadJets[j].Pt(),hadJets[j].Eta(),wt);
//     TLorentzVector Met;
//     Met.SetPtEtaPhiE(MET,0,METPhi,0);
//     double dPhi_METjet = abs(Met.DeltaPhi(hadJets[j])); 
//     h_HT5HT_vsdPhi_METJet[j][i]->Fill(dPhi_METjet,HT5/HT,wt);
//     h_dPhi_METJet[j][i]->Fill(dPhi_METjet,wt);
//     h_dPhi_METJet_vsMET[j][i]->Fill(MET,dPhi_METjet,wt);
//   }

    
//   h_nBjets_vs_qmulti[i]->Fill(btags,qmulti,wt);
//   h_qmultiVs_MET[i]->Fill(MET,qmulti,wt);
//   h_qmultiVs_ST[i]->Fill(ST,qmulti,wt);
//   h_leadJets_qmulti[i]->Fill(leadJets_qmulti,wt);
//   h_leadJet_Pt[i]->Fill(leadJet_Pt,wt);
//   h_leadbjet_tag[i]->Fill(leadbjet_tag,wt);
//   h_leadbjet_tag_vs_leadQmulti[i]->Fill(leadbjet_tag,leadJets_qmulti,wt);
//   h_leadbjet_tag_vsQmulti[i]->Fill(leadbjet_tag,qmulti,wt);
//   h_leadjet_ptvsqmulti[i]->Fill(leadJet_Pt,qmulti,wt);

//   int TFbins_v2 = getBinNoV1_le(pho_Pt, qmulti);
//   int TFbins_v3 = getBinNoV16_le(Njets,btags,pho_Pt);
//   int searchBin = getBinNoV6_WithOnlyBLSelec(Njets,btags);
//   int TFbins_v4 = getBinNoV7_le(Njets,btags);
//   int TFbins_v5 = getBinNo_v0FR(pho_Pt, qmulti, 0.1);
//   int TFbins_v6 = getBinNo_v1FR(pho_Pt,Njets);
//   int TFbins_v7 = getBinNo_v2FR(qmulti,Njets,btags);

//   h_TFbins_LL_v4[i]->Fill(TFbins_v4,wt);
//   h_TFbins_LL_v2[i]->Fill(TFbins_v2,wt);
  
//   h_Sbins_LL[i]->Fill(searchBin,wt);
//   h_TFbins_LL_v3[i]->Fill(TFbins_v3,wt);
//   h_TFbins_LL_v5[i]->Fill(TFbins_v5,wt);
//   h_TFbins_LL_v6[i]->Fill(TFbins_v6,wt);
//   h_TFbins_LL_v7[i]->Fill(TFbins_v7,wt);
  
  
//   searchBin = getBinNoV3_ST_MET_bjets(ST,MET,btags);
//   h_Sbins_LL_newSbins_v3[i]->Fill(searchBin,wt);

//   searchBin = getBinNoV7_ST_MET_bjets_phopT(ST,MET,btags,pho_Pt);
//   h_Sbins_LL_newSbins_v7[i]->Fill(searchBin,wt);
//   searchBin = getBin_ST_MET_bjets_phopT_merge(MET,ST,btags,pho_Pt);
//   h_Sbins_LL_newSbins_v7_merge[i]->Fill(searchBin,wt);
  
 
// }



//SS//== not using this functions == 
Bool_t AnalyzeTProxytBSM::Process(Long64_t entry) {

  std::cout << entry << std::endl;
  fDirector.SetReadEntry(entry);
  std::cout<< "entry " << entry << " RunNum " << RunNum << std::endl;
  std::cout << "GenParticles->size() "<< GenParticles->size() << std::endl;
  return 0;
}



























//===============================================================================================================================================================================================================
//                                                                                                BACKUP
//===============================================================================================================================================================================================================
	  // if (Ch_e == 1 && Ch_mu == 0 && Ch_tau == 0) {
	  //   int identifier = 1;
	  //   h_EvtBrk->Fill(identifier);
	  // }
    // if (Ch_mu == 1 && Ch_e == 0 && Ch_tau == 0) {
    //   int identifier = 2;
    //   h_EvtBrk->Fill(identifier);
    // }
    // if (Ch_tau == 1 && Ch_e == 0 && Ch_mu == 0) {
    //   int identifier = 3;
    //   h_EvtBrk->Fill(identifier);
    // }
    
    // if (Ch_e == 2) {
    //   int identifier = 4;
    //   h_EvtBrk->Fill(identifier);
    // }
    // if (Ch_mu == 2) {
    //   int identifier = 5;
    //   h_EvtBrk->Fill(identifier);
    // }
    // if (Ch_tau == 2) {
    //   int identifier = 6;
    //   h_EvtBrk->Fill(identifier);
      //    }
	  
    // if (Ch_e == 1 && Ch_mu == 1 && Ch_tau == 0) {
    //   int identifier = 7;
    //   h_EvtBrk->Fill(identifier);
    // }
    // if (Ch_e == 1 && Ch_tau == 1 && Ch_mu == 0) {
    //   int identifier = 8;
    //   h_EvtBrk->Fill(identifier);
    // }
    // if (Ch_mu == 1 && Ch_tau == 1 && Ch_e == 0) {
    //   int identifier = 9;
    //   h_EvtBrk->Fill(identifier);
    // }
    // if (Ch_e == 0 && Ch_mu == 0 && Ch_tau == 0) {
    //   int identifier = 10;
    //   h_EvtBrk->Fill(identifier);
    // }
	             
    //std::cout << std::endl; 
    //std::cout << "Electrons->size() "<< Electrons->size() << std::endl;




  // std::cout << "No. of events with 0 Leptons " << NEvtlep0 << std::endl;
  // std::cout << "No. of events with 1 Leptons " << NEvtlep1 << std::endl;
  // std::cout << "No. of events with 2 Leptons " << NEvtlep2 << std::endl;
  // std::cout << "No. of events with 3 Leptons " << NEvtlep3 << std::endl;
  // std::cout << "No. of events with 4 Leptons " << NEvtlep4 << std::endl;
  // std::cout <<"total Events (sum of all cases) " << NEvtlep0 + NEvtlep1 + NEvtlep2 + NEvtlep3 + NEvtlep4 << std::endl;

      // h_NHadJets->Fill(NHadJets,Jets_pT_Sum);
      //if (jentry < 100 && IsoTracks==0) cout << NHadJets << endl;
      //if (jentry <100) cout << Weight << endl;
      // if (jentry < 20000 && GenElectrons->size()>1){
      // 	for (int ii = 0; ii < GenElectrons->size(); ii++){
      // 	  cout << "genelec pt: " << GenTaus[(int)ii].Pt(),wt << ", ";
      // 	  //cout << "Event No.: " << jentry << "\n\n";
	
      // 	}
      // 	cout << "\n\n";
      // }
      //if (GenElectrons->size()>1) cout << GenElectrons->size() << endl;
      //if (NElectrons!=Electrons->size()) cout << "Alert!!" << endl;
      // if (PdgId == 22) h_NJet_genPhoPt->Fill(GenTaus[(int)ii].Pt(),wt,NHadJets);
      //if (jentry < 20000 && abs(PdgId) == 11) cout << "Genparticle_electron  pt:" << GenTaus[(int)ii].Pt(),wt << endl ; 
      //if (abs(PdgId)==11) {
	//Ch_e++;
	//h_Gen_pT[2][4]->Fill(GenTaus[(int)ii].Pt(),wt1);
	
      //}
      //   if (abs(PdgId)==13) Ch_mu++;
	    //   if (abs(PdgId)==15) Ch_tau++;
	    //if (jentry < 100)  cout << "particle id: " << GenParticles_PdgId[(int)ii] << ", Mother Id: " << GenParticles_ParentId[(int)ii] << ", particle status: " << GenParticles_Status[(int)ii] << endl;

//cout << "No. of events with MET>100: " << nEvents << endl;
  //cout << "nentries: " << nentries << endl;


// h_ele_pT0  ->Fill(Electrons[(int)ii].Pt(),wt1);
	    // h_ele_eta0 ->Fill(Electrons[(int)ii].Eta(),wt1);
	    // h_ele_phi0 ->Fill(Electrons[(int)ii].Phi(),wt1);
	    // if (MET>200){
	    //   h_ele_pT3  ->Fill(Electrons[(int)ii].Pt(),wt1);
	    //   h_ele_eta3 ->Fill(Electrons[(int)ii].Eta(),wt1);
	    //   h_ele_phi3 ->Fill(Electrons[(int)ii].Phi(),wt1);
	    // }

	    //if (Electrons->size()>1 && jentry <10000) cout << Electrons[(int)ii].Pt(),wt << ", ";

// 	if(GenTaus[(int)ii].Pt(),wt>1.0){
	    // 	h_gen_pT  ->Fill(GenTaus[(int)ii].Pt(),wt1);
	    // 	h_gen_eta ->Fill(GenTaus[(int)ii].Eta(),wt1);
	    // 	h_gen_phi ->Fill(GenTaus[(int)ii].Phi(),wt1);
	    // 	}
	  


	  //   //if (jentry < 1000)     std::cout << "pdgid: " <<GenParticles_PdgId[(int)ii] << ", parentid: " << GenParticles_ParentId[(int)ii] << endl; 
	  //   //for filling the the leptons in histogram
	  //   int PdgId = GenParticles_PdgId[(int)ii];
	  //     double dR=DeltaR(bestPhoton.Eta(),bestPhoton.Phi(),GenParticles[(int)ii].Eta(),wt,GenParticles[(int)ii].Phi(),wt1);
	  //   // if (GenParticles[(int)ii].Pt(),wt>10 && abs(GenParticles[(int)ii].Eta(),wt)<2.5){
	  //   if (abs(PdgId) == 11) {
	  //     //if (jentry<100 &&  GenElectrons->size()>0) cout << GenElectrons[(int)ii].Pt(),wt << endl;
	  //     h_Gen_pT[0][0]->Fill(GenParticles[(int)ii].Pt(),wt1);
	  //     NGenE++;
	  //   }
	      
	  //   if (abs(PdgId) == 13){
	  //     h_Gen_pT[1][0]->Fill(GenParticles[(int)ii].Pt(),wt1);
	  //     NGenM++;
	  //   }
	  //   if (abs(PdgId) == 15) {
	  //     h_Gen_pT[2][0]->Fill(GenParticles[(int)ii].Pt(),wt1);
	  //     NGenT++;
	  //   }

	  //   if (abs(PdgId) == 11) h_Gen_eta[0][0]->Fill(GenParticles[(int)ii].Eta(),wt1); 
	  //   if (abs(PdgId) == 13) h_Gen_eta[1][0]->Fill(GenParticles[(int)ii].Eta(),wt1);
	  //   if (abs(PdgId) == 15) h_Gen_eta[2][0]->Fill(GenParticles[(int)ii].Eta(),wt1);

	  //   if (abs(PdgId) == 11) h_Gen_phi[0][0]->Fill(GenParticles[(int)ii].Phi(),wt1); 
	  //   if (abs(PdgId) == 13) h_Gen_phi[1][0]->Fill(GenParticles[(int)ii].Phi(),wt1);
	  //   if (abs(PdgId) == 15) h_Gen_phi[2][0]->Fill(GenParticles[(int)ii].Phi(),wt1);


	  //   if (ST<300) continue;
	  //   if (abs(PdgId) == 11) h_Gen_pT[0][4]->Fill(GenParticles[(int)ii].Pt(),wt1); 
	  //   if (abs(PdgId) == 13) h_Gen_pT[1][4]->Fill(GenParticles[(int)ii].Pt(),wt1);
	  //   if (abs(PdgId) == 15) h_Gen_pT[2][4]->Fill(GenParticles[(int)ii].Pt(),wt1);

	  //   if (abs(PdgId) == 11) h_Gen_eta[0][4]->Fill(GenParticles[(int)ii].Eta(),wt1); 
	  //   if (abs(PdgId) == 13) h_Gen_eta[1][4]->Fill(GenParticles[(int)ii].Eta(),wt1);
	  //   if (abs(PdgId) == 15) h_Gen_eta[2][4]->Fill(GenParticles[(int)ii].Eta(),wt1);

	  //   if (abs(PdgId) == 11) h_Gen_phi[0][4]->Fill(GenParticles[(int)ii].Phi(),wt1); 
	  //   if (abs(PdgId) == 13) h_Gen_phi[1][4]->Fill(GenParticles[(int)ii].Phi(),wt1);
	  //   if (abs(PdgId) == 15) h_Gen_phi[2][4]->Fill(GenParticles[(int)ii].Phi(),wt1);


	  //   if (NEMu!=0) continue;
	  //   if (abs(PdgId) == 11) h_Gen_pT[0][5]->Fill(GenParticles[(int)ii].Pt(),wt1); 
	  //   if (abs(PdgId) == 13) h_Gen_pT[1][5]->Fill(GenParticles[(int)ii].Pt(),wt1);
	  //   if (abs(PdgId) == 15) h_Gen_pT[2][5]->Fill(GenParticles[(int)ii].Pt(),wt1);

	  //   if (abs(PdgId) == 11) h_Gen_eta[0][5]->Fill(GenParticles[(int)ii].Eta(),wt1); 
	  //   if (abs(PdgId) == 13) h_Gen_eta[1][5]->Fill(GenParticles[(int)ii].Eta(),wt1);
	  //   if (abs(PdgId) == 15) h_Gen_eta[2][5]->Fill(GenParticles[(int)ii].Eta(),wt1);

	  //   if (abs(PdgId) == 11) h_Gen_phi[0][5]->Fill(GenParticles[(int)ii].Phi(),wt1); 
	  //   if (abs(PdgId) == 13) h_Gen_phi[1][5]->Fill(GenParticles[(int)ii].Phi(),wt1);
	  //   if (abs(PdgId) == 15) h_Gen_phi[2][5]->Fill(GenParticles[(int)ii].Phi(),wt1);


	  //   if (!(isoElectronTracks == 0 && isoMuonTracks == 0 && isoPionTracks == 0)) continue;
	  //   if (abs(PdgId) == 11) h_Gen_pT[0][6]->Fill(GenParticles[(int)ii].Pt(),wt1); 
	  //   if (abs(PdgId) == 13) h_Gen_pT[1][6]->Fill(GenParticles[(int)ii].Pt(),wt1);
	  //   if (abs(PdgId) == 15) h_Gen_pT[2][6]->Fill(GenParticles[(int)ii].Pt(),wt1);

	  //   if (abs(PdgId) == 11) h_Gen_eta[0][6]->Fill(GenParticles[(int)ii].Eta(),wt1); 
	  //   if (abs(PdgId) == 13) h_Gen_eta[1][6]->Fill(GenParticles[(int)ii].Eta(),wt1);
	  //   if (abs(PdgId) == 15) h_Gen_eta[2][6]->Fill(GenParticles[(int)ii].Eta(),wt1);

	  //   if (abs(PdgId) == 11) h_Gen_phi[0][6]->Fill(GenParticles[(int)ii].Phi(),wt1); 
	  //   if (abs(PdgId) == 13) h_Gen_phi[1][6]->Fill(GenParticles[(int)ii].Phi(),wt1);
	  //   if (abs(PdgId) == 15) h_Gen_phi[2][6]->Fill(GenParticles[(int)ii].Phi(),wt1);


	//for lost e, lost mu, and e faking photon
	    // if (abs(PdgId) == 13 && NMuons == 0) h_LostMuon_eta->Fill(GenParticles[(int)ii].Eta(),wt1);
	   
	    // else if (abs(PdgId) == 11 && dR > 0.1 &&  NElectrons == 0) h_LostElectron_eta->Fill(GenParticles[(int)ii].Eta(),wt1);
	     
	    // else if (abs(PdgId) == 11 && bestPhotonIndxAmongPhotons >= 0 && dR < 0.1 && NElectrons == 0) h_EFakePho_eta->Fill(GenParticles[(int)ii].Eta(),wt1);
	
	    // else if (abs(PdgId) == 15 && GenTaus_had[0]) h_HadTau_eta -> Fill[GenParticles(int)ii].Eta(),wt;

	    // else continue;

	    //  else cout << "nothing found" << endl;
   
	   




	  //if (NGenE+NGenM+NGenT == 0) NGenL++;
	  
	  //if (jentry < 100) cout << ("\n\n");
	  // if (jentry < 1000) cout << "Event: " << jentry << endl;
	  // for(Long64_t ii=0; ii<GenParticles->size(); ii++){
	  //   if (NGenE != Electrons->size() && abs(GenParticles_PdgId[(int)ii]) == 16 && jentry < 1000) cout << "tau found"<< endl; 
					     
	  // }

// cout << "no. of events with no e,mu,tau: " << NGenL << endl;
  // cout << "no. of events electron faking photon= " << NEFakePho << endl;
  // cout << "no. of events with lost electrons= " << NLostElectrons << endl;
  // cout << "no. of events with lost electrons= " << NLostMuons << endl;

  //if (GenElectrons->size()>1) cout << GenElectrons->size() << endl;

	//bool LostMu_flag=false, LostE_flag=false, EfakePho_flag=false, hadTau_flag=false, Rest_flag=false;
	
	// if(GenElectrons -> size() > 0 && NElectrons == 0 && bestPhotonIndxAmongPhotons > 0)
	//   { double dR = DeltaR(bestPhoton.Eta(),bestPhoton.Phi(),GenElectrons[0].Eta(),GenElectrons[0].Phi());
	//     if (dR > 0.1) LostE_flag = true;	      
	//     else EfakePho_flag = true;	     
	//   }
	
	//  if (!((ST<300) && (NEMu!=0) && (isoElectronTracks || isoMuonTracks || isoPionTracks)))
	//    {h_Gen_eta[0][6] -> Fill(GenElectrons[0].Eta());
	//      //h_Reco_eta[0][6] -> Fill(GenElectrons[0].Eta());
	//      if (LostE_flag) h_LostElectron_eta -> Fill(GenElectrons[0].Eta());
	//      if (EfakePho_flag) h_EFakePho_eta-> Fill(GenElectrons[0].Eta());
	//    }


 //filling the histogram for total
	    //   if(GenElectrons->size()>0 && GenMuons->size()==0 && GenTaus->size()==0) h_Gen_MET->Fill(MET);
	    //   if(GenElectrons->size()==0 && GenMuons->size()>0 && GenTaus->size()==0) h_Gen_MET->Fill(MET);
	    //   if(GenElectrons->size()==0 && GenMuons->size()==0 && GenTaus->size()>0)h_Gen_MET->Fill(MET);
	    //   if(GenElectrons->size()>0 && GenMuons->size()==0 && GenTaus->size()>0)h_Gen_MET->Fill(MET);
	    //   if(GenElectrons->size()==0 && GenMuons->size()>0 && GenTaus->size()>0)h_Gen_MET->Fill(MET);

	    //   //filling the histogrms for the different categories
	    //   if (LostE_flag) h_LostElectron_MET -> Fill(MET);
	    //   if (EfakePho_flag) h_EFakePho_MET-> Fill(MET);
	    //   if (LostMu_flag) h_LostMuon_MET -> Fill(MET);
	    //   if (hadTau_flag) h_HadTau_MET -> Fill(MET);  
	    //   if (Rest_flag) h_Rest_MET-> Fill(MET);
	    // }
		
// if (LostE_flag) h_LostElectron_eta -> Fill(GenElectrons[0].Eta(),wt1);
	  // if (EfakePho_flag) h_EFakePho_eta->Fill(GenElectrons[0].Eta(),wt1);
	  
	  
	    
	  // // for lost e, mu and e fake photon
	  // //for (Long64_t ii=0; ii<GenElectrons->size(); ii++){
	  // bool LostMu_flag, LostE_flag, EfakePho_flag, hadTau_flag, Rest_flag;
	  // // Double_t rndm;
	  // // rndm.Uniform(-10,10);		
	  
	  // //if (GenMuons->size() > 0 && NMuons == 0) h_LostMuon_eta -> Fill(GenMuons[0].Eta(),wt1);
	  // if (GenMuons->size() > 0 && NMuons == 0); 
	  //    { h_LostMuon_eta -> Fill(GenElectrons[0].Eta(),wt1);
	  //      cout << "jfskla" << endl;
	      
	  //   }
	  // //LostMu_flag = true;
	    
	  // else if(GenElectrons -> size() > 0 && NElectrons == 0 && bestPhotonIndxAmongPhotons > 0)
	  //   { double dR = DeltaR(bestPhoton.Eta(),bestPhoton.Phi(),GenElectrons[0].Eta(),GenElectrons[0].Phi());
	  //     if (dR > 0.1) h_LostElectron_eta -> Fill(GenElectrons[0].Eta(),wt1);
	  //     //LostE_flag = true;	      
	  //     else h_EFakePho_eta->Fill(GenElectrons[0].Eta(),wt1);
	     
	  //   }
	  
	  // //else if (GenTaus->size() > 0 && GenTaus_had[0]) h_HadTau_eta -> Fill(GenParticles[0].Eta(),wt1);
	  // else if (GenTaus->size() > 0 && GenTaus_had[0]) h_HadTau_eta -> Fill(GenElectrons[0].Eta(),wt1);
	  // // else h_Rest_eta->FillRandom("gaus",1);
	  // else h_Rest_eta->Fill(GenElectrons[0].Eta(),wt1);


// if (hadTau_flag) h_HadTau_eta -> Fill(GenTaus[0].Eta(),wt1);  
	  // if (Rest_flag) h_Rest_eta->FillRandom("gaus",1);

	  
	 
	    //AnalyzeTProxytBSM ana(inputFileList, outFileName, data,sampl;
  //ana.EventLoop(data,inputFileList,sample,outFileName,phoID);
  //ana.EventLoop(inputFileList,data,sample);


  // int NEvtlep0 = 0;
  // int NEvtlep1 = 0;
  // int NEvtlep2 = 0;
  // int NEvtlep3 = 0;
  // int NEvtlep4 = 0;



	// h_NHadJets[0]->Fill(0.0000,wt1);
       // 	//h_NHadJets[0]->Fill(NHadJets,wt1);
       // sumwt += wt;
       // if (k > decade) {
       // 	 cout << "sum weight" << sumwt << endl;
       // 	 cout << "NHadjets Integral: " << h_NHadJets[0]->Integral() << endl;
       // 	 cout << "overflow NhadJets: " << h_NHadJets[0]->h_TF ->GetBinContent(h_NHadJets[0]->GetNbinsX() + 1) << endl;
       // 	 cout << "underflow NHadJets: " << h_NHadJets[0]->h_TF ->GetBinContent(-1) << "\n\n";
       // }

       // //h_MET[0]->Fill(0.0000,wt1);
       // h_MET[0]->Fill(MET,wt1);
       // if (k > decade) {
       // 	 cout << "MET Integral: " << h_MET[0]->Integral() << endl;
       // 	 cout << "overflow MET: " << h_MET[0]->h_TF ->GetBinContent(h_MET[0]->GetNbinsX() + 1) << endl;
       // 	 cout << "underflow MET: " << h_MET[0]->h_TF ->GetBinContent(-1) << "\n\n";  
       // }

       // //h_Pho_pT[0] -> Fill(0.0000,wt1);
       // h_Pho_pT[0] -> Fill(bestPhoton.Pt(),wt1);
       // if (k > decade) {
       // 	 cout << "Photon pt Integral: " << h_Pho_pT[0]->Integral() << endl;
       // 	 cout << "overflow Pho Pt: " << h_Pho_pT[0]->h_TF ->GetBinContent(h_Pho_pT[0]->GetNbinsX() + 1) << endl;
       // 	 cout << "underflow Photon pt: " << h_Pho_pT[0]->h_TF ->GetBinContent(-1) << "\n\n";
       // }



// h_pho_eta0  ->Fill(bestPhoton.Eta());
      // h_pho_phi0 ->Fill(bestPhoton.Phi());
      
      // h_NJet_PhoPt->Fill(bestPhoton.Pt(),NHadJets);
      

      
      // h_NHadJets[0]->Fill(NHadJets,wt1);
      // h_MET[0]->Fill(MET,wt1);
      // h_Pho_pT[0]  ->Fill(bestPhoton.Pt(),wt1);
      // h_Pho_eta[0]  ->Fill(bestPhoton.Eta(),wt1);
      // h_Pho_phi[0]  ->Fill(bestPhoton.Phi(),wt1);

      // if (NEMu == 0) {
      // 	h_MET[5] ->Fill(MET);
      // 	h_Pho_pT[5] ->Fill(bestPhoton.Pt());
      // 	h_Pho_eta[5] -> Fill(bestPhoton.Eta());
      // 	h_Pho_phi[5] -> Fill(bestPhoton.Phi());
      // 	h_NHadJets[5]-> Fill(NHadJets);
	
      // 	if (MET  < 200) continue;
      // 	nEvents++;
      // 	h_MET[1]-> Fill(MET),wt;
      // 	h_Pho_pT[1]  ->Fill(bestPhoton.Pt(),wt1);
      // 	h_Pho_eta[1]  ->Fill(bestPhoton.Eta(),wt1);
      // 	h_Pho_phi[1]  ->Fill(bestPhoton.Phi(),wt1);
      // 	h_NHadJets[1]-> Fill(NHadJets,wt1);

      // 	if (bestPhoton.Pt()<20) continue;
      // 	h_MET[2] -> Fill(MET);
      // 	h_Pho_pT[2] -> Fill(bestPhoton.Pt());
      // 	h_Pho_eta[2] -> Fill(bestPhoton.Eta());
      // 	h_Pho_phi[2] -> Fill(bestPhoton.Phi());
      // 	h_NHadJets[2]-> Fill(NHadJets);

      // 	if (NHadJets < 2) continue;
      // 	h_MET[3] ->Fill(MET);
      // 	h_Pho_pT[3] ->Fill(bestPhoton.Pt());
      // 	h_Pho_eta[3] -> Fill(bestPhoton.Eta());
      // 	h_Pho_phi[3] -> Fill(bestPhoton.Phi());
      // 	h_NHadJets[3]-> Fill(NHadJets);

      // 	if (ST < 300) continue;
      // 	h_MET[4] ->Fill(MET);
      // 	h_Pho_pT[4] ->Fill(bestPhoton.Pt());
      // 	h_Pho_eta[4] -> Fill(bestPhoton.Eta());
      // 	h_Pho_phi[4] -> Fill(bestPhoton.Phi());
      // 	h_NHadJets[4]-> Fill(NHadJets);

	
	
      // 	if (isoElectronTracks != 0 || isoMuonTracks != 0 || isoPionTracks != 0) continue;
      // 	h_MET[6] ->Fill(MET);
      // 	h_Pho_pT[6] ->Fill(bestPhoton.Pt());
      // 	h_Pho_eta[6] -> Fill(bestPhoton.Eta());
      // 	h_Pho_phi[6] -> Fill(bestPhoton.Phi());
      // 	h_NHadJets[6] -> Fill(NHadJets);
      // }


	
      
