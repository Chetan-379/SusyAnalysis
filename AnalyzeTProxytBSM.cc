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
    ana.EventLoop(buffer.c_str());
    delete chain; 
    delete fin;
  }

  // === === some random summary === ===
  cout << "dataset " << data << " " << endl;
  cout<<"If analyzing the lost electron estimation ? "<<"  "<<elec<<endl;
  cout<<"Which pho_ID: "<<"\t"<<phoID<<endl;
  return 0;
}

//void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList, const char *sample , const char *outFileName, const char *elec, const char* phoID) {
void AnalyzeTProxytBSM::EventLoop(std::string buffer) {
  
  std::cout << "AnalyzeTProxytBSM::EventLoop() " << std::endl;

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  std::cout << "Analyzing " << buffer.c_str() << " nentries " << nentries << std::endl;  

  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
  // int NEvtlep0 = 0;
  // int NEvtlep1 = 0;
  // int NEvtlep2 = 0;
  // int NEvtlep3 = 0;
  // int NEvtlep4 = 0;
  int nEvents=0, NGenL=0, NLostElectrons=0, NLostMuons=0, NEFakePho=0;;
  for (Long64_t jentry=0; jentry<fChain->GetEntries(); jentry++){
  //for (Long64_t jentry=0; jentry<10000; jentry++){
    
      fDirector.SetReadEntry(jentry);
      
      // == == print number of events done == == == == == == == =
      double progress = 10.0 * jentry / (1.0 * nentries);
      int k = int (progress);
      if (k > decade)
	cout << 10 * k << " %" << endl;
      decade = k;
      
      //std::cout << jentry << " " << MET << std::endl;
      
   
      vector<myLV> hadJets, bjets;
      int BTags = bjets.size();
      bool Debug=false;
      vector<int> jetMatchindx;
      int bJet1Idx=-100;
      double deepCSVvalue = 0.4148;
      double minDR=99999; 
      int minDRindx=-100;
      //TLorentzVector AnalyzeTProxytBSM::getBestPhoton(int pho_ID);
      int pho_ID=0;  //for simplicity taking only soft one
      //TLorentzVector bestPhoton=getBestPhoton(pho_ID);
      myLV bestPhoton=getBestPhoton(pho_ID);
      int hadJetID=-999;
      int NJets0=Jets->size();
      int NHadJets = 0;
      float Jets_pT_Sum=0;
      float ST=0;
      int Iso_Lep_Tracks, NEMu;
      // double mT; 

      //======================
      //ss//      h_MET[0]->Fill(MET);
      //======================

      //bool passST = (bool)(ST<300);

      
      for(int i=0;i<Jets->size();i++){
	if( (Jets[i].Pt() > 30.0) && (abs(Jets[i].Eta()) <= 2.4) ){
	  double dR=DeltaR(bestPhoton.Eta(),bestPhoton.Phi(),Jets[i].Eta(),Jets[i].Phi());
	  if(dR<minDR) minDR=dR;minDRindx=i;
	  
	} 
      }
      
      

      if(Debug)
	cout<<"===load tree entry  ==="<<"\t"<<jentry<<"\t"<<"Jets check == "<<minDR<<endl;
     
      Iso_Lep_Tracks = isoElectronTracks + isoMuonTracks + isoPionTracks;
      NEMu = NElectrons + NMuons;
      for(int i=0;i<Jets->size();i++){
      	//zif(Debug)
      	  //cout<<"  = Jets.Pt()  ==  "<<Jets_v1[i].Pt()<<"\t"<< " = Jets.Eta() == "<<Jets_v1[i].Eta()<<endl;
	  
      	  if( (Jets[i].Pt() > 30.0) && (abs(Jets[i].Eta()) <= 2.4) ){
      	    //if(Debug)
      	    //cout<< "==== loadjets ==="<<"\t"<<i<<"\t"<<minDR<<endl;
	    
      	    if( !(minDR < 0.3 && i==minDRindx) )
      	      {
		
      		hadJetID= (*Jets_ID)[i];
		//if (jentry < 100) cout << "hadjetId: " << hadJetID << endl;
      		if(hadJetID)
      		  {
      		    hadJets.push_back(Jets[i]);
		    NHadJets++;
		    
		    
	      
      		    // hadJets_hadronFlavor.push_back((*Jets_hadronFlavor)[i]);
      		    // hadJets_HTMask.push_back((*Jets_HTMask)[i]);
      		    // hadJets_bJetTagDeepCSVBvsAll.push_back((*Jets_bJetTagDeepCSVBvsAll)[i]);
      		    // if(q==1) leadjet_qmulti=(*Jets_chargedMultiplicity)[q];
      		    // if(q==1) leadjet_Pt=(*Jets)[q].Pt();
		      if((*Jets_bJetTagDeepCSVBvsAll)[i] > deepCSVvalue){
			bjets.push_back(Jets[i]); bJet1Idx = i;}
		      //hadJets.push_back((*Jets)[i]);
		      jetMatchindx.push_back(i);
		      
		  }
	      }
	  }	  
      }
      
	for(int i=0;i<Jets->size();i++){
	  if( (Jets[i].Pt() > 30.0) && (abs(Jets[i].Eta()) <= 2.4) ){
	    if (!(minDR < 0.3 && i==minDRindx)){
	      if (Jets_ID[i]){
		Jets_pT_Sum = Jets_pT_Sum + Jets[i].Pt();
		ST = Jets_pT_Sum + bestPhoton.Pt();
		
		h_Jet_pT[0]->Fill(Jets[i].Pt());
		h_Jet_eta[0]->Fill(Jets[i].Eta());
		h_Jet_phi[0]->Fill(Jets[i].Phi());
		
		if (MET>200){
		  h_Jet_pT[1]->Fill(Jets[i].Pt());
		  h_Jet_eta[1]->Fill(Jets[i].Eta());
		  h_Jet_phi[1]->Fill(Jets[i].Phi());
		
		  if (bestPhoton.Pt()<20) continue;
		  h_Jet_pT[2]->Fill(Jets[i].Pt());
		  h_Jet_eta[2]->Fill(Jets[i].Eta());
		  h_Jet_phi[2]->Fill(Jets[i].Phi());
		  
		  if (NHadJets < 2) continue;
		  h_Jet_pT[3]->Fill(Jets[i].Pt());
		  h_Jet_eta[3]->Fill(Jets[i].Eta());
		  h_Jet_phi[3]->Fill(Jets[i].Phi());
		  
		  if (!(NEMu == 0)) continue;
		  h_Jet_pT[5]->Fill(Jets[i].Pt());
		  h_Jet_eta[5]->Fill(Jets[i].Eta());
		  h_Jet_phi[5]->Fill(Jets[i].Phi());

		  if (!(Iso_Lep_Tracks == 0)) continue;
		  h_Jet_pT[6]->Fill(Jets[i].Pt());
		  h_Jet_eta[6]->Fill(Jets[i].Eta());
		  h_Jet_phi[6]->Fill(Jets[i].Phi());
		  
		}
	      }
	    }
	  }
	}
	
	if(hadJets.size()==0) continue;
	if(Debug)
	  cout<<"===load tree entry ===  "<<"\t"<<jentry<<"\t"<<"No of B-Jets ===  "<<bjets.size()<<endl;
	
	  //if(jentry<10 ) {
	  // std::cout<< "jentry " << jentry << " RunNum " << RunNum << std::endl;
    // std::cout << "GenParticles->size() "<< GenParticles->size() << std::endl;
    // int Nlep = 0;
    // int Ch_e = 0;
    // int Ch_mu = 0;
    // int Ch_tau = 0;
	  
          	if (jentry < 1000) cout << "Event: " << jentry << endl;
    // begin genparticle loop
	//if (GenTaus_had->size() > 1) cout << GenTaus_had->size() << endl;
	//cout << GenTaus_had << endl;
	int NGenE=0, NGenM=0, NGenT=0;
	  for(Long64_t ii=0; ii<GenParticles->size(); ii++) {
	    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > mygen = GenParticles[(int)ii];
	    //std::cout <<" ii, Pt, Eta, Phi, E " << ii << " " << GenParticles[(int)ii].Pt() << " " << GenParticles[(int)ii].Eta() << " " << GenParticles[(int)ii].Phi() << " " << GenParticles[(int)ii].E() << " pdgid, parentid, status " << GenParticles_PdgId[(int)ii] << " " << GenParticles_ParentId[(int)ii] << " " << GenParticles_Status[(int)ii]  << std::endl;
	    if (jentry < 1000)     std::cout << "pdgid: " <<GenParticles_PdgId[(int)ii] << ", parentid: " << GenParticles_ParentId[(int)ii] << endl; 
	    //for filling the the leptons in histogram
	    int PdgId = GenParticles_PdgId[(int)ii];
	    double dR=DeltaR(bestPhoton.Eta(),bestPhoton.Phi(),GenParticles[(int)ii].Eta(),GenParticles[(int)ii].Phi());
	    // if (GenParticles[(int)ii].Pt()>10 && abs(GenParticles[(int)ii].Eta())<2.5){
	    if (abs(PdgId) == 11) {
	      //if (jentry<100 &&  GenElectrons->size()>0) cout << GenElectrons[(int)ii].Pt() << endl;
	      h_Gen_pT[0][0]->Fill(GenParticles[(int)ii].Pt());
	      NGenE++;
	    }
	      
	    if (abs(PdgId) == 13){
	      h_Gen_pT[1][0]->Fill(GenParticles[(int)ii].Pt());
	      NGenM++;
	    }
	    if (abs(PdgId) == 15) {
	      h_Gen_pT[2][0]->Fill(GenParticles[(int)ii].Pt());
	      NGenT++;
	    }

	    if (abs(PdgId) == 11) h_Gen_eta[0][0]->Fill(GenParticles[(int)ii].Eta()); 
	    if (abs(PdgId) == 13) h_Gen_eta[1][0]->Fill(GenParticles[(int)ii].Eta());
	    if (abs(PdgId) == 15) h_Gen_eta[2][0]->Fill(GenParticles[(int)ii].Eta());

	    if (abs(PdgId) == 11) h_Gen_phi[0][0]->Fill(GenParticles[(int)ii].Phi()); 
	    if (abs(PdgId) == 13) h_Gen_phi[1][0]->Fill(GenParticles[(int)ii].Phi());
	    if (abs(PdgId) == 15) h_Gen_phi[2][0]->Fill(GenParticles[(int)ii].Phi());


	    if (ST<300) continue;
	    if (abs(PdgId) == 11) h_Gen_pT[0][4]->Fill(GenParticles[(int)ii].Pt()); 
	    if (abs(PdgId) == 13) h_Gen_pT[1][4]->Fill(GenParticles[(int)ii].Pt());
	    if (abs(PdgId) == 15) h_Gen_pT[2][4]->Fill(GenParticles[(int)ii].Pt());

	    if (abs(PdgId) == 11) h_Gen_eta[0][4]->Fill(GenParticles[(int)ii].Eta()); 
	    if (abs(PdgId) == 13) h_Gen_eta[1][4]->Fill(GenParticles[(int)ii].Eta());
	    if (abs(PdgId) == 15) h_Gen_eta[2][4]->Fill(GenParticles[(int)ii].Eta());

	    if (abs(PdgId) == 11) h_Gen_phi[0][4]->Fill(GenParticles[(int)ii].Phi()); 
	    if (abs(PdgId) == 13) h_Gen_phi[1][4]->Fill(GenParticles[(int)ii].Phi());
	    if (abs(PdgId) == 15) h_Gen_phi[2][4]->Fill(GenParticles[(int)ii].Phi());


	    if (NEMu!=0) continue;
	    if (abs(PdgId) == 11) h_Gen_pT[0][5]->Fill(GenParticles[(int)ii].Pt()); 
	    if (abs(PdgId) == 13) h_Gen_pT[1][5]->Fill(GenParticles[(int)ii].Pt());
	    if (abs(PdgId) == 15) h_Gen_pT[2][5]->Fill(GenParticles[(int)ii].Pt());

	    if (abs(PdgId) == 11) h_Gen_eta[0][5]->Fill(GenParticles[(int)ii].Eta()); 
	    if (abs(PdgId) == 13) h_Gen_eta[1][5]->Fill(GenParticles[(int)ii].Eta());
	    if (abs(PdgId) == 15) h_Gen_eta[2][5]->Fill(GenParticles[(int)ii].Eta());

	    if (abs(PdgId) == 11) h_Gen_phi[0][5]->Fill(GenParticles[(int)ii].Phi()); 
	    if (abs(PdgId) == 13) h_Gen_phi[1][5]->Fill(GenParticles[(int)ii].Phi());
	    if (abs(PdgId) == 15) h_Gen_phi[2][5]->Fill(GenParticles[(int)ii].Phi());


	    if (Iso_Lep_Tracks!=0) continue;
	    if (abs(PdgId) == 11) h_Gen_pT[0][6]->Fill(GenParticles[(int)ii].Pt()); 
	    if (abs(PdgId) == 13) h_Gen_pT[1][6]->Fill(GenParticles[(int)ii].Pt());
	    if (abs(PdgId) == 15) h_Gen_pT[2][6]->Fill(GenParticles[(int)ii].Pt());

	    if (abs(PdgId) == 11) h_Gen_eta[0][6]->Fill(GenParticles[(int)ii].Eta()); 
	    if (abs(PdgId) == 13) h_Gen_eta[1][6]->Fill(GenParticles[(int)ii].Eta());
	    if (abs(PdgId) == 15) h_Gen_eta[2][6]->Fill(GenParticles[(int)ii].Eta());

	    if (abs(PdgId) == 11) h_Gen_phi[0][6]->Fill(GenParticles[(int)ii].Phi()); 
	    if (abs(PdgId) == 13) h_Gen_phi[1][6]->Fill(GenParticles[(int)ii].Phi());
	    if (abs(PdgId) == 15) h_Gen_phi[2][6]->Fill(GenParticles[(int)ii].Phi());

	    if (abs(PdgId) == 13 && bestPhotonIndxAmongPhotons < 0 &&  Muons->size() == 0)
	      {
		NLostMuons++;	    
		h_LostMuon_eta->Fill(GenParticles[(int)ii].Eta());
	      }

	    else if (abs(PdgId) == 11 && bestPhotonIndxAmongPhotons < 0 &&  Electrons->size() == 0)
	      {
		NLostElectrons++;
		h_LostElectron_eta->Fill(GenParticles[(int)ii].Eta());
	      }
	    
	    // else if (abs(PdgId) == 11 && bestPhotonIndxAmongPhotons >= 0 && dR < 0.1 && Electrons->size() == 0)
	    //    {
	    // 	NEFakePho++;
	    // 	h_EFakePho_eta->Fill(GenParticles[(int)ii].Eta());
	    //    }

	    // else if 
	    //   {for (Long64_t jj=0; jj<GenTaus_had->size(); jj++) {

		 
	    // 	}
	    //   }
	    //else if (GenTaus_had[(int)ii]){
	    // else if (ii<GenTaus_had->size())
	    //   {
	    //   if (GenTaus_had[(int)ii])
	    // 	 h_HadTau_eta->Fill(GenParticles[(int)ii].Eta());
		
	    //   }
	  }   //end genparticle loop
	  //	  cout << endl;
	  //h_GenRecoE->Fill(NGenE,Electrons->size());
	
	  if (NGenE+NGenM+NGenT == 0) NGenL++;
	  
	  //if (jentry < 100) cout << ("\n\n");
	  // if (jentry < 1000) cout << "Event: " << jentry << endl;
	  // for(Long64_t ii=0; ii<GenParticles->size(); ii++){
	  //   if (NGenE != Electrons->size() && abs(GenParticles_PdgId[(int)ii]) == 16 && jentry < 1000) cout << "tau found"<< endl; 
					     
	  // }

	  for(Long64_t ii=0; ii<Electrons->size(); ii++){
	    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > myele = Electrons[(int)ii];
	    //std::cout <<" ii, Pt, Eta, Phi, E " << ii << " " << Electrons[(int)ii].Pt() 		  << " " << Electrons[(int)ii].Eta() << " " << Electrons[(int)ii].Phi() 		  << " " << Electrons[(int)ii].E()  		  << " iso, mediumID " << Electrons_iso[(int)ii] << " " << Electrons_mediumID[(int)ii]		  << std::endl;
	    
	    
	    h_Reco_pT[0][0]->Fill(Electrons[(int)ii].Pt());
	    h_Reco_eta[0][0]->Fill(Electrons[(int)ii].Eta());
	    h_Reco_phi[0][0]->Fill(Electrons[(int)ii].Phi());
	    
	    if (ST<300) continue;
	    h_Reco_pT[0][4]->Fill(Electrons[(int)ii].Pt());
	    h_Reco_eta[0][4]->Fill(Electrons[(int)ii].Eta());
	    h_Reco_phi[0][4]->Fill(Electrons[(int)ii].Phi());
	    
	    
	    // if (NElectrons!=0) continue;
	    // h_Reco_pT[0][5]->Fill(Electrons[(int)ii].Pt());
	    // h_Reco_eta[0][5]->Fill(Electrons[(int)ii].Eta());
	    // h_Reco_phi[0][5]->Fill(Electrons[(int)ii].Phi());
	    
	    vector<myLV> v_electron;
	    if (Electrons_passIso) v_electron.push_back(Electrons[(int)ii]);
	    //sortTLorVec(&v_lepton);
	    //else v_lepton.push_back(0);
	    //if (v_electron.size()!=0) continue;
	    // h_Reco_pT[0][5]->Fill(v_electron[0].Pt());
	    // h_Reco_eta[0][5]->Fill(v_electron[0].Eta());
	    // h_Reco_phi[0][5]->Fill(v_electron[0].Phi());

	    if (NEMu!=0) continue;
	    h_Reco_pT[0][5]->Fill(Electrons[(int)ii].Pt());
	    h_Reco_eta[0][5]->Fill(Electrons[(int)ii].Eta());
	    h_Reco_phi[0][5]->Fill(Electrons[(int)ii].Phi());
	    
	    
	    if (Iso_Lep_Tracks!=0) continue;
	    h_Reco_pT[0][6]->Fill(Electrons[(int)ii].Pt());
	    h_Reco_eta[0][6]->Fill(Electrons[(int)ii].Eta());
	    h_Reco_phi[0][6]->Fill(Electrons[(int)ii].Phi());
	   
	    
	    
	  } //end electron loop
	  //if (Electrons->size()>1 && jentry <10000) cout << endl;
	  for(Long64_t ii=0; ii<Muons->size(); ii++){
	    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > mymu = Muons[(int)ii];

	    h_Reco_pT[1][0]->Fill(Muons[(int)ii].Pt());
	    h_Reco_eta[1][0]->Fill(Muons[(int)ii].Eta());
	    h_Reco_phi[1][0]->Fill(Muons[(int)ii].Phi());

	    if (ST<300) continue;

	    h_Reco_pT[1][4]->Fill(Muons[(int)ii].Pt());
	    h_Reco_eta[1][4]->Fill(Muons[(int)ii].Eta());
	    h_Reco_phi[1][4]->Fill(Muons[(int)ii].Phi());

	    vector<myLV> v_muon;
	    if (Muons_passIso) v_muon.push_back(Muons[(int)ii]);
	    if (v_muon.size()!=0) continue;
	    h_Reco_pT[1][5]->Fill(v_muon[0].Pt());
	    h_Reco_eta[1][5]->Fill(v_muon[0].Eta());
	    h_Reco_phi[1][5]->Fill(v_muon[0].Phi());

	    if (Iso_Lep_Tracks!=0) continue;
	    h_Reco_pT[1][6]->Fill(Muons[(int)ii].Pt());
	    h_Reco_eta[1][6]->Fill(Muons[(int)ii].Eta());
	    h_Reco_phi[1][6]->Fill(Muons[(int)ii].Phi());
	  } //end Muon loop

      
	     //std::cout << std::endl; 
      //std::cout << "Photons->size() "<< Photons->size() << std::endl;
      for(Long64_t ii=0; ii<Photons->size(); ii++){
	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > mypho = Photons[(int)ii];
	//std::cout <<" ii, Pt, Eta, Phi, E " << ii << " " << Photons[(int)ii].Pt() 		  << " " << Photons[(int)ii].Eta() << " " << Photons[(int)ii].Phi() 		  << " " << Photons[(int)ii].E()  		  << " mvavalueID, pfGammaIso " << Photons_mvaValuesID[(int)ii] << " " << Photons_pfGammaIso[(int)ii]		  << std::endl;
	
	
	//h_pho_eta ->Fill(Photons[(int)ii].Eta());
	//h_pho_phi ->Fill(Photons[(int)ii].Phi());
      } //end photon loop
      // h_pho_eta0  ->Fill(bestPhoton.Eta());
      // h_pho_phi0 ->Fill(bestPhoton.Phi());
      
      // h_NJet_PhoPt->Fill(bestPhoton.Pt(),NHadJets);


      
      h_NHadJets[0]->Fill(NHadJets);
      h_MET[0]->Fill(MET);
      h_Pho_pT[0]  ->Fill(bestPhoton.Pt());
      h_Pho_eta[0]  ->Fill(bestPhoton.Eta());
      h_Pho_phi[0]  ->Fill(bestPhoton.Phi());
      if (MET  > 200) {nEvents++;
	h_MET[1]-> Fill(MET);
	h_Pho_pT[1]  ->Fill(bestPhoton.Pt());
	h_Pho_eta[1]  ->Fill(bestPhoton.Eta());
	h_Pho_phi[1]  ->Fill(bestPhoton.Phi());
	h_NHadJets[1]-> Fill(NHadJets);
	if (bestPhoton.Pt()<20) continue;
	h_MET[2] -> Fill(MET);
	h_Pho_pT[2] -> Fill(bestPhoton.Pt());
	h_Pho_eta[2] -> Fill(bestPhoton.Eta());
	h_Pho_phi[2] -> Fill(bestPhoton.Phi());
	h_NHadJets[2]-> Fill(NHadJets);
	if (NHadJets < 2) continue;
	h_MET[3] ->Fill(MET);
	h_Pho_pT[3] ->Fill(bestPhoton.Pt());
	h_Pho_eta[3] -> Fill(bestPhoton.Eta());
	h_Pho_phi[3] -> Fill(bestPhoton.Phi());
	h_NHadJets[3]-> Fill(NHadJets);
	if (ST < 300) continue;
	h_MET[4] ->Fill(MET);
	h_Pho_pT[4] ->Fill(bestPhoton.Pt());
	h_Pho_eta[4] -> Fill(bestPhoton.Eta());
	h_Pho_phi[4] -> Fill(bestPhoton.Phi());
	h_NHadJets[4]-> Fill(NHadJets);
	if (!(NEMu == 0)) continue;
	h_MET[5] ->Fill(MET);
	h_Pho_pT[5] ->Fill(bestPhoton.Pt());
	h_Pho_eta[5] -> Fill(bestPhoton.Eta());
	h_Pho_phi[5] -> Fill(bestPhoton.Phi());
	h_NHadJets[5]-> Fill(NHadJets);
	
	if (!(Iso_Lep_Tracks == 0)) continue;
	h_MET[6] ->Fill(MET);
	h_Pho_pT[6] ->Fill(bestPhoton.Pt());
	h_Pho_eta[6] -> Fill(bestPhoton.Eta());
	h_Pho_phi[6] -> Fill(bestPhoton.Phi());
	h_NHadJets[6] -> Fill(NHadJets);
      }
  } // end jentry loop 
  cout << "no. of events with no e,mu,tau: " << NGenL << endl;
  cout << "no. of events electron faking photon= " << NEFakePho << endl;
  cout << "no. of events with lost electrons= " << NLostElectrons << endl;
  cout << "no. of events with lost electrons= " << NLostMuons << endl;
}

  myLV AnalyzeTProxytBSM::getBestPhoton(int pho_ID){
  //TLorentzVector AnalyzeTProxytBSM::getBestPhoton(int pho_ID){
  //vector<TLorentzVector> goodPho;
  vector<myLV> goodPho;
  vector<int> goodPhoIndx;
  for(int iPho=0;iPho<Photons->size();iPho++){
    //if(((*Photons_hasPixelSeed)[iPho]<0.001) && ( (*Photons_fullID)[iPho]))
    if(((*Photons_hasPixelSeed)[iPho]<0.001)  && ( (*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001) &&( pho_ID==0 || (pho_ID==1 &&(((*Photons_cutBasedID)[iPho]==1 || (*Photons_cutBasedID)[iPho]==2))) || (pho_ID==2 && (*Photons_cutBasedID)[iPho]==2) || (pho_ID==3 && (*Photons_mvaValuesID)[iPho]>-0.02) || (pho_ID==4 && (*Photons_mvaValuesID)[iPho]>0.42))) ) 
      {
	goodPho.push_back(Photons[iPho] );
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
   //if(highPtIndx==-100){TLorentzVector v0;return v0;}
   else return goodPho[highPtIndx];
   
  }
  
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
      // 	  cout << "genelec pt: " << GenElectrons[(int)ii].Pt() << ", ";
      // 	  //cout << "Event No.: " << jentry << "\n\n";
	
      // 	}
      // 	cout << "\n\n";
      // }
      //if (GenElectrons->size()>1) cout << GenElectrons->size() << endl;
      //if (NElectrons!=Electrons->size()) cout << "Alert!!" << endl;
      // if (PdgId == 22) h_NJet_genPhoPt->Fill(GenParticles[(int)ii].Pt(),NHadJets);
      //if (jentry < 20000 && abs(PdgId) == 11) cout << "Genparticle_electron  pt:" << GenParticles[(int)ii].Pt() << endl ; 
      //if (abs(PdgId)==11) {
	//Ch_e++;
	//h_Gen_pT[0][4]->Fill(GenParticles[(int)ii].Pt());
	
      //}
      //   if (abs(PdgId)==13) Ch_mu++;
	    //   if (abs(PdgId)==15) Ch_tau++;
	    //if (jentry < 100)  cout << "particle id: " << GenParticles_PdgId[(int)ii] << ", Mother Id: " << GenParticles_ParentId[(int)ii] << ", particle status: " << GenParticles_Status[(int)ii] << endl;

//cout << "No. of events with MET>100: " << nEvents << endl;
  //cout << "nentries: " << nentries << endl;


// h_ele_pT0  ->Fill(Electrons[(int)ii].Pt());
	    // h_ele_eta0 ->Fill(Electrons[(int)ii].Eta());
	    // h_ele_phi0 ->Fill(Electrons[(int)ii].Phi());
	    // if (MET>200){
	    //   h_ele_pT3  ->Fill(Electrons[(int)ii].Pt());
	    //   h_ele_eta3 ->Fill(Electrons[(int)ii].Eta());
	    //   h_ele_phi3 ->Fill(Electrons[(int)ii].Phi());
	    // }

	    //if (Electrons->size()>1 && jentry <10000) cout << Electrons[(int)ii].Pt() << ", ";

// 	if(GenParticles[(int)ii].Pt()>1.0){
	    // 	h_gen_pT  ->Fill(GenParticles[(int)ii].Pt());
	    // 	h_gen_eta ->Fill(GenParticles[(int)ii].Eta());
	    // 	h_gen_phi ->Fill(GenParticles[(int)ii].Phi());
	    // 	}
	  
	    
