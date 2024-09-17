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
  double cross_section = getCrossSection(s_process);
  //double wt = getEventWeight(s_process, cross_section);
  double wt = 0;
  //cout << "weight is: "<< wt << endl;
  double wt1 = 0;
  //bool genphocheck=false;
  //int genphomatch_before=0, genphomatch_after=0;

  std::cout << cross_section << "\t" <<"analyzed process"<<"\t"<<s_cross<<endl;
  
  Long64_t nbytes = 0, nb = 0;
  int decade = 0;
  double sumwt = 0;
  
  int nEvents=0, NGenL=0, NLostElectrons=0, NLostMuons=0, NEFakePho=0;
  string src_file="TF_NBJet.root";
  TFile *root_file = new TFile(src_file.c_str());
  TH1F* h_TF;	   
  if(s_sample.Contains("TTG")) h_TF = (TH1F*) root_file ->Get("h_TF_TTGJets");
  if(s_sample.Contains("TTJets") || s_sample.Contains("TTJets_Leptons")) h_TF = (TH1F*) root_file ->Get("h_TF_TTJets");
  if(s_sample.Contains("WJets")) h_TF = (TH1F*) root_file ->Get("h_TF_WJets");
  if(s_sample.Contains("WGJets_MonoPhoton_PtG-40to130UL")||s_sample.Contains("WGJets_MonoPhoton_PtG-130UL")) h_TF = (TH1F*) root_file ->Get("h_TF_WGJets");

	 

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
      double deepCSVvalue = 0.4148,p0=1.787e+02,p1=6.657e+01,p2=9.47e-01;
      double minDR=99999;
      int phoMatchingJetIndx = -100; 
      int minDRindx=-100;
      int pho_ID=0;  //for simplicity taking only soft one
      myLV bestPhoton=getBestPhoton(pho_ID);
      int hadJetID=-999;
      int NJets0=Jets->size();
      int NHadJets = 0;
      float Jets_pT_Sum=0;
      float ST=0;
      int NEMu;
      double dPhi_METjet1, dPhi_METjet2;
      

      NEMu = NElectrons + NMuons;

      
      //if (jentry < 20) cout << wt << endl; 	    
      // if(s_Process.Contains("2018.WGJets_MonoPhoton_PtG-40to130UL")|| s_Process.Contains("2018.WGJets_MonoPhoton_PtG-130UL")|| s_Process.Contains("2016preVFP.WGJets_MonoPhoton_PtG-40to130UL") ||s_Process.Contains("2016preVFP.WGJets_MonoPhoton_PtG-130UL") || s_Process.Contains("2017.WGJets_MonoPhoton_PtG-40to130UL")||s_Process.Contains("2017.WGJets_MonoPhoton_PtG-130UL")|| s_Process.Contains("2016postVFP.WGJets_MonoPhoton_PtG-130UL")||s_Process.Contains("2016postVFP.WGJets_MonoPhoton_PtG-40to130UL")) wt = cross_section*59.83*1000.0;
      // else wt = Weight*59.83*1000.0;
      wt = getEventWeight(s_process, cross_section);
      //selecting Hadronic and b jets
      for(int i=0;i<Jets->size();i++){
	if( (Jets[i].Pt() > 30.0) && (abs(Jets[i].Eta()) <= 2.4) ){
	  if (Photons->size()!=0) {
	    double dR=DeltaR(bestPhoton.Eta(),bestPhoton.Phi(),Jets[i].Eta(),Jets[i].Phi());
	    if(dR<minDR) {minDR=dR;minDRindx=i;}
	  }
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
      
      for(int i=0;i<hadJets.size();i++){
	
	if( (abs(hadJets[i].Eta()) < 2.4) ){
	  NHadJets++;
	  ST=ST+(hadJets[i].Pt());
	}
      }

      if( minDR<0.3) {
	ST=ST+bestPhoton.Pt();
	phoMatchingJetIndx = minDRindx;
      }

      if(NHadJets>=1)
	dPhi_METjet1 = abs(DeltaPhi(METPhi,hadJets[0].Phi()));
      if(NHadJets>=2)
	dPhi_METjet2 = abs(DeltaPhi(METPhi,hadJets[1].Phi()));

      
      //defining flags for applying baseline selections
      bool Pass_EMu_veto=false, Pass_Iso_trk_veto=false, Pass_Pho_pT=false, Pass_MET=false, Pass_NHadJets=false, Pass_ST=false, applyTrgEff = false, EvtCln=false,JetMetPhi=false,rmOvrlp=false, Pass_MET2=false;
      	  
	  if (bestPhoton.Pt() >40){
	    Pass_Pho_pT = true;
	    if (MET > 100){
	      Pass_MET = true;
	      if (NHadJets >=2){
		Pass_NHadJets = true;		
		if (ST > 300){
		  Pass_ST = true;		      	
		  applyTrgEff = true;
		  if(PrimaryVertexFilter==1 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 &&HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter==1 && BadPFMuonDzFilter==1 && eeBadScFilter==1 && ecalBadCalibFilter==1 && NVtx>0 && PFCaloMETRatio < 5){
		    if((!(phoMatchingJetIndx>=0 && (Jets[phoMatchingJetIndx].Pt())/(bestPhoton.Pt()) < 1.0)) && phoMatchingJetIndx >= 0){
		      EvtCln = true;			  
		      if(dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3){
			JetMetPhi = true;						  
			if(RemoveSampleOverlap(s_sample, bestPhoton)){
			  rmOvrlp = true;
			  if(MET>200) {
			    Pass_MET2=true;
			    if (NEMu == 0) {
			      Pass_EMu_veto = true;	
			      if (!(isoElectronTracks || isoMuonTracks || isoPionTracks)){
				Pass_Iso_trk_veto = true;			  	    							    
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  
	  //if(jentry<20) cout << wt << endl;      
	  //h_NHadJets[0]->Fill(0.0,wt);      
	  h_NHadJets[0]->Fill(NHadJets,wt);      
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
	  }
	  
	  if (Pass_MET){ 
	    h_MET[1]-> Fill(MET,wt);
	    h_Pho_pT[1]  ->Fill(bestPhoton.Pt(),wt);
	    h_Pho_eta[1]  ->Fill(bestPhoton.Eta(),wt);
	    h_Pho_phi[1]  ->Fill(bestPhoton.Phi(),wt);
	    h_NHadJets[1]-> Fill(NHadJets,wt);
	  }
	  
	  if (Pass_NHadJets){
	    h_MET[3] ->Fill(MET,wt);
	    h_Pho_pT[3] ->Fill(bestPhoton.Pt(),wt);
	    h_Pho_eta[3] -> Fill(bestPhoton.Eta(),wt);
	    h_Pho_phi[3] -> Fill(bestPhoton.Phi(),wt);
	    h_NHadJets[3]-> Fill(NHadJets,wt);
	  }
	  
	  if (Pass_ST){
	    h_MET[4] ->Fill(MET,wt);
	    h_Pho_pT[4] ->Fill(bestPhoton.Pt(),wt);
	    h_Pho_eta[4] -> Fill(bestPhoton.Eta(),wt);
	    h_Pho_phi[4] -> Fill(bestPhoton.Phi(),wt);
	    h_NHadJets[4]-> Fill(NHadJets,wt);
	  }
	  
	  if (applyTrgEff)
	    {
	      wt1 = wt * (((TMath::Erf((MET - p0)/p1)+1)/2.0)*p2);
	      h_MET[7] ->Fill(MET,wt1);
	      h_Pho_pT[7] ->Fill(bestPhoton.Pt(),wt1);
	      h_Pho_eta[7] -> Fill(bestPhoton.Eta(),wt1);
	      h_Pho_phi[7] -> Fill(bestPhoton.Phi(),wt1);
	      h_NHadJets[7]-> Fill(NHadJets,wt1);
	    }
      
	  if (EvtCln)
	    {
	      h_MET[8] ->Fill(MET,wt1);
	      h_Pho_pT[8] ->Fill(bestPhoton.Pt(),wt1);
	      h_Pho_eta[8] -> Fill(bestPhoton.Eta(),wt1);
	      h_Pho_phi[8] -> Fill(bestPhoton.Phi(),wt1);
	      h_NHadJets[8]-> Fill(NHadJets,wt1);   
	    }
      
	  if (JetMetPhi)
	    {
	      h_MET[9] ->Fill(MET,wt1);
	      h_Pho_pT[9] ->Fill(bestPhoton.Pt(),wt1);
	      h_Pho_eta[9] -> Fill(bestPhoton.Eta(),wt1);
	      h_Pho_phi[9] -> Fill(bestPhoton.Phi(),wt1);
	      h_NHadJets[9]-> Fill(NHadJets,wt1);
	    }
	  
	  if (rmOvrlp)
	    {
	      h_MET[10] ->Fill(MET,wt1);
	      h_Pho_pT[10] ->Fill(bestPhoton.Pt(),wt1);
	      h_Pho_eta[10] -> Fill(bestPhoton.Eta(),wt1);
	      h_Pho_phi[10] -> Fill(bestPhoton.Phi(),wt1);
	      h_NHadJets[10]-> Fill(NHadJets,wt1);
	    }
	  
	  if (Pass_MET2)
	    {
	      h_MET[11] ->Fill(MET,wt1);
	      h_Pho_pT[11] ->Fill(bestPhoton.Pt(),wt1);
	      h_Pho_eta[11] -> Fill(bestPhoton.Eta(),wt1);
	      h_Pho_phi[11] -> Fill(bestPhoton.Phi(),wt1);
	      h_NHadJets[11]-> Fill(NHadJets,wt1);
	    }
	  
	  if (Pass_EMu_veto) {
	    h_MET[5] ->Fill(MET,wt1);
	    h_Pho_pT[5] ->Fill(bestPhoton.Pt(),wt1);
	    h_Pho_eta[5] -> Fill(bestPhoton.Eta(),wt1);
	    h_Pho_phi[5] -> Fill(bestPhoton.Phi(),wt1);
	    h_NHadJets[5]-> Fill(NHadJets,wt1);
	  }
	  
	  if (Pass_Iso_trk_veto){
	    h_MET[6] ->Fill(MET,wt1);
	    h_Pho_pT[6] ->Fill(bestPhoton.Pt(),wt1);
	    h_Pho_eta[6] -> Fill(bestPhoton.Eta(),wt1);
	    h_Pho_phi[6] -> Fill(bestPhoton.Phi(),wt1);
	    h_NHadJets[6] -> Fill(NHadJets,wt1);
	    
	  }
	 

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

      int TFbins = getBinNoV1_le(NHadJets,BTags);
      int SrchBins = getBinNoV6_WithOnlyBLSelec(NHadJets,BTags);

          //selecting SR for lost e
      bool lost_elec_CR = false, lost_elec_SR = false, FR_SR = false, FR_CR = false;      
      if (Pass_Iso_trk_veto && GenElectrons->size()>0 && ((GenElectrons[0].Pt()/bestPhoton.Pt())<=0.8 || (GenElectrons[0].Pt()/bestPhoton.Pt())>=1.2)){        
	double dR_gen_e_gamma = DeltaR(GenElectrons[0].Eta(),GenElectrons[0].Phi(),bestPhoton.Eta(),bestPhoton.Phi());
	if (dR_gen_e_gamma > 0.2) lost_elec_SR = true;
      }

      //Fake rate SR
      else if (Pass_Iso_trk_veto && GenElectrons->size()>0 && !((GenElectrons[0].Pt()/bestPhoton.Pt())<=0.8 || (GenElectrons[0].Pt()/bestPhoton.Pt())>=1.2)){
	double dR_gen_e_gamma = DeltaR(GenElectrons[0].Eta(),GenElectrons[0].Phi(),bestPhoton.Eta(),bestPhoton.Phi());
	if (dR_gen_e_gamma < 0.2) FR_SR = true;	    
	  }
      
      //selecting CR for lost e  
      if (Pass_MET2 && NElectrons==1 && NMuons==0 && isoElectronTracks && (!(isoMuonTracks || isoPionTracks))){
	double dR_reco_e_gamma = DeltaR(Electrons[0].Eta(),Electrons[0].Phi(),bestPhoton.Eta(),bestPhoton.Phi());
	    float m_T=sqrt(2*((Electrons[0].Pt()*MET)-(Electrons[0].Pt()*MET*cos(DeltaPhi(Electrons[0].Phi(),METPhi)))));
	    h_mT_reco_e_G->Fill(m_T);
	    if (dR_reco_e_gamma >= 0.2 && m_T<100) lost_elec_CR = true;
      }

     
      //selecting lost muon SR
      bool lost_mu_CR = false, lost_mu_SR = false;
      if (Pass_Iso_trk_veto){        
	lost_mu_SR = true;
	  }
      
      //selecting lost muon CR
      float m_T = 0;
	  if (Pass_MET2 && NElectrons==0 && NMuons==1 && isoMuonTracks && (!( isoElectronTracks || isoPionTracks))){	    
	    m_T = sqrt(2*((Muons[0].Pt()*MET)-(Muons[0].Pt()*MET*cos(DeltaPhi(Muons[0].Phi(),METPhi)))));
	    if (m_T<100) lost_mu_CR = true;
	  }
	  
	  
	  if (lost_elec_SR) {
	    h_Lost_e_SR_Pho_Pt->Fill(bestPhoton.Pt(),wt1);
	    h_Lost_e_SR_binned ->Fill(TFbins,wt1);
	    h_Lost_e_SR_srch_binned -> Fill(SrchBins,wt1);
	  }

	  if (lost_elec_CR) {
	    h_Lost_e_CR_Pho_Pt ->Fill(bestPhoton.Pt(),wt1);
	    h_Lost_e_CR_binned ->Fill(TFbins,wt1);
	    h_Lost_e_CR_srch_binned -> Fill(SrchBins,wt1);
	    h_mT_reco_e_G -> Fill(m_T);
	  }

	  if (FR_SR) {
	    h_FR_SR_binned -> Fill(TFbins,wt1);
	  }
	  
	  if (lost_mu_SR){
	    h_Lost_mu_SR_binned -> Fill(TFbins,wt1);
	  }
	
	  
	  //cout << h_TF->GetName() << endl;
	  
	  if (lost_elec_CR){		    
	    //calculating the SR yield in diff bins
	    float TF=0;	           
	   
	   
	    for(int ibin = 0;ibin < 7; ibin++){
	      double wt2=0; 

	      if (ibin == 1){
		if (NHadJets==2) TF = h_TF ->GetBinContent(1);
		if (NHadJets==3) TF = h_TF ->GetBinContent(2);
		if (NHadJets==4) TF = h_TF ->GetBinContent(3);		
		wt2 = TF*wt1;		  
		for(int jbin =1; jbin <=7; jbin++){
		  h_Lost_e_SR_srch_binned_pred->Fill(SrchBins,wt2);
		  //SR_pred[jbin] = TF*h_Lost_e_CR_binned->h_TF ->GetBinContent(jbin);
		}
	      }
	      
	      if (ibin == 2){
		TF = h_TF ->GetBinContent(4);
		wt2 = TF*wt1;		 
		for(int jbin =8; jbin <=13; jbin++){
		  h_Lost_e_SR_srch_binned_pred->Fill(SrchBins,wt2);
		  //SR_pred[jbin] = TF*h_Lost_e_CR_binned->h_TF ->GetBinContent(jbin);
		}
	      }
	      
	      if (ibin == 3){
		TF = h_TF ->GetBinContent(5);
		wt2 = TF*wt1;	
		for(int jbin =14; jbin <=18; jbin++){
		  h_Lost_e_SR_srch_binned_pred->Fill(SrchBins,wt2);
		  //SR_pred[jbin] = TF*h_Lost_e_CR_binned->h_TF ->GetBinContent(jbin);
		}
	      }
	      
	      if (ibin == 4){
		TF = h_TF ->GetBinContent(6);
		wt2 = TF*wt1;	
		for(int jbin =19; jbin <=23; jbin++){
		  h_Lost_e_SR_srch_binned_pred->Fill(SrchBins,wt2);
		  //SR_pred[jbin] = TF*h_Lost_e_CR_binned->h_TF ->GetBinContent(jbin);
		}
	      }
	      
	      if (ibin == 5){
		TF = h_TF ->GetBinContent(7);
		wt2 = TF*wt1;
		for(int jbin = 24; jbin <= 28; jbin++){
		  h_Lost_e_SR_srch_binned_pred->Fill(SrchBins,wt2);
		  //SR_pred[jbin] = TF*h_Lost_e_CR_binned->h_TF ->GetBinContent(jbin);
		}
	      }
	      
	      if (ibin == 6){
		TF = h_TF ->GetBinContent(7);
		wt2 = TF*wt1;
		for(int jbin = 24; jbin <= 28; jbin++){
		  h_Lost_e_SR_srch_binned_pred->Fill(SrchBins,wt2);
		  //SR_pred[jbin] = TF*h_Lost_e_CR_binned->h_TF ->GetBinContent(jbin);
		}
	      }	    
	      //root_file->Close();
	    }
	  }


	  //================================================================================================================
	  //testing
	  if (Pass_Iso_trk_veto && GenElectrons->size()>0){        
	    double dR_gen_e_gamma = DeltaR(GenElectrons[0].Eta(),GenElectrons[0].Phi(),bestPhoton.Eta(),bestPhoton.Phi());
	    h_dR_gen_e_reco_pho ->Fill(dR_gen_e_gamma);
	    h_gen_e_reco_pho_ratio ->Fill(GenElectrons[0].Pt()/bestPhoton.Pt());
	    h_dRvsRatio->Fill(dR_gen_e_gamma,GenElectrons[0].Pt()/bestPhoton.Pt());
	  }

	  //==============================================================================================================

            	                             
      for(int i=0;i<Jets->size();i++){
	if( (Jets[i].Pt() > 30.0) && (abs(Jets[i].Eta()) <= 2.4) ){
	  if (!(minDR < 0.3 && i==minDRindx)){
	    if (Jets_ID[i]){
	      Jets_pT_Sum += Jets[i].Pt();
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
		
		if (isoElectronTracks || isoMuonTracks || isoPionTracks) continue;
		h_Jet_pT[6]->Fill(Jets[i].Pt());
		h_Jet_eta[6]->Fill(Jets[i].Eta());
		h_Jet_phi[6]->Fill(Jets[i].Phi());
		
		}
	      }
	    }
	  }
      }
      
	//if(hadJets.size()==0) continue;

	if(Debug)
	  cout<<"===load tree entry ===  "<<"\t"<<jentry<<"\t"<<"No of B-Jets ===  "<<bjets.size()<<endl;

	//defining flags for different categories of detector inefficiencies
	bool LostMu_flag=false, LostE_flag=false, EfakePho_flag=false, hadTau_flag=false, Rest_flag=false;
	if (GenMuons->size() > 0 && NMuons == 0) LostMu_flag = true;	    
	else if(GenElectrons -> size() > 0 && NElectrons == 0 && bestPhotonIndxAmongPhotons > 0)
	  { double dR = DeltaR(bestPhoton.Eta(),bestPhoton.Phi(),GenElectrons[0].Eta(),GenElectrons[0].Phi());
	    if (dR > 0.1) LostE_flag = true;	      
	    else EfakePho_flag = true;	     
	  }
	else if(GenElectrons -> size() > 0 && NElectrons == 0 && bestPhotonIndxAmongPhotons < 0) LostE_flag = true;	
	else if (GenTaus->size() > 0 && GenTaus_had[0]) hadTau_flag = true;
	else if (GenElectrons->size()>0 || GenMuons->size()>0 || GenTaus->size()>0) Rest_flag = true;
	
	//filling histo to check stack plot:
	if (MET>200 && bestPhoton.Pt()>40 && NHadJets>=2 && ST>300 && NEMu==0)
	  {
	    //filling the histogram for total
	    if(GenElectrons->size()>0 && GenMuons->size()==0 && GenTaus->size()==0) h_Gen_eta[0][6]->Fill(GenElectrons[0].Eta());
	    if(GenElectrons->size()==0 && GenMuons->size()>0 && GenTaus->size()==0) h_Gen_eta[0][6]->Fill(GenMuons[0].Eta());
	    if(GenElectrons->size()==0 && GenMuons->size()==0 && GenTaus->size()>0)h_Gen_eta[0][6]->Fill(GenTaus[0].Eta());
	    if(GenElectrons->size()>0 && GenMuons->size()==0 && GenTaus->size()>0)h_Gen_eta[0][6]->Fill(GenElectrons[0].Eta());
	    if(GenElectrons->size()==0 && GenMuons->size()>0 && GenTaus->size()>0)h_Gen_eta[0][6]->Fill(GenMuons[0].Eta());

	    
	    //filling the histogrms for the different categories
	    if (LostE_flag) h_LostElectron_eta -> Fill(GenElectrons[0].Eta());
	    if (EfakePho_flag) h_EFakePho_eta-> Fill(GenElectrons[0].Eta());
	    if (LostMu_flag) h_LostMuon_eta -> Fill(GenMuons[0].Eta());
	    if (hadTau_flag) h_HadTau_eta -> Fill(GenTaus[0].Eta());  
	    if (Rest_flag) h_Rest_eta-> FillRandom("gaus",1);
	    }

	  
	     
	// Loop for plotting the Gen information
	for(Long64_t ii=0; ii<GenParticles->size(); ii++) {
	}   //end genparticle loop
	
	for (Long64_t ii =0; ii<GenElectrons->size(); ii++){
	  h_Gen_pT[0][0]->Fill(GenElectrons[(int)ii].Pt(),wt1);
	  h_Gen_eta[0][0]->Fill(GenElectrons[(int)ii].Eta(),wt1);
	  h_Gen_phi[0][0]->Fill(GenElectrons[(int)ii].Phi(),wt1);
	  
	  
	  if (ST<300) continue;
	  h_Gen_pT[0][4]->Fill(GenElectrons[(int)ii].Pt(),wt1);  
	  h_Gen_eta[0][4]->Fill(GenElectrons[(int)ii].Eta(),wt1); 	    
	  h_Gen_phi[0][4]->Fill(GenElectrons[(int)ii].Phi(),wt1); 
	  
	  
	  if (NEMu!=0) continue;
	  h_Gen_pT[0][5]->Fill(GenElectrons[(int)ii].Pt(),wt1); 	    
	  h_Gen_eta[0][5]->Fill(GenElectrons[(int)ii].Eta(),wt1); 
	  h_Gen_phi[0][5]->Fill(GenElectrons[(int)ii].Phi(),wt1); 
	  

	  if (isoElectronTracks != 0 || isoMuonTracks != 0 || isoPionTracks != 0) continue;
	  h_Gen_pT[0][6]->Fill(GenElectrons[(int)ii].Pt(),wt1);  
	 //h_Gen_eta[0][6]->Fill(GenElectrons[(int)ii].Eta(),wt1); 	    
	  h_Gen_phi[0][6]->Fill(GenElectrons[(int)ii].Phi(),wt1);
	}  //end Gen electron loop

	
	for (Long64_t ii =0; ii<GenMuons->size(); ii++){
	  h_Gen_pT[1][0]->Fill(GenMuons[(int)ii].Pt(),wt1);
	  h_Gen_eta[1][0]->Fill(GenMuons[(int)ii].Eta(),wt1);
	  h_Gen_phi[1][0]->Fill(GenMuons[(int)ii].Phi(),wt1);
	  
	  
	  if (ST<300) continue;
	  h_Gen_pT[1][4]->Fill(GenMuons[(int)ii].Pt(),wt1);  
	  h_Gen_eta[1][4]->Fill(GenMuons[(int)ii].Eta(),wt1); 	    
	  h_Gen_phi[1][4]->Fill(GenMuons[(int)ii].Phi(),wt1); 
	  
	  
	  if (NEMu!=0) continue;
	  h_Gen_pT[1][5]->Fill(GenMuons[(int)ii].Pt(),wt1); 	    
	  h_Gen_eta[1][5]->Fill(GenMuons[(int)ii].Eta(),wt1); 
	  h_Gen_phi[1][5]->Fill(GenMuons[(int)ii].Phi(),wt1); 
	  
	  
	  if (isoElectronTracks || isoMuonTracks || isoPionTracks) continue;
	  h_Gen_pT[1][6]->Fill(GenMuons[(int)ii].Pt(),wt1);  
	  h_Gen_eta[1][6]->Fill(GenMuons[(int)ii].Eta(),wt1); 	    
	  h_Gen_phi[1][6]->Fill(GenMuons[(int)ii].Phi(),wt1);
	} //end Gen Muon loop
	
	for (Long64_t ii =0; ii<GenTaus->size(); ii++){
	  h_Gen_pT[2][0]->Fill(GenTaus[(int)ii].Pt(),wt1);
	  h_Gen_eta[2][0]->Fill(GenTaus[(int)ii].Eta(),wt1);
	  h_Gen_phi[2][0]->Fill(GenTaus[(int)ii].Phi(),wt1);
	  
	  if (ST<300) continue;
	  h_Gen_pT[2][4]->Fill(GenTaus[(int)ii].Pt(),wt1);  
	  h_Gen_eta[2][4]->Fill(GenTaus[(int)ii].Eta(),wt1); 	    
	  h_Gen_phi[2][4]->Fill(GenTaus[(int)ii].Phi(),wt1); 
	  
	  if (NEMu!=0) continue;
	  h_Gen_pT[2][5]->Fill(GenTaus[(int)ii].Pt(),wt1); 	    
	  h_Gen_eta[2][5]->Fill(GenTaus[(int)ii].Eta(),wt1); 
	  h_Gen_phi[2][5]->Fill(GenTaus[(int)ii].Phi(),wt1); 
	  
	  if (isoElectronTracks != 0 || isoMuonTracks != 0 || isoPionTracks != 0) continue;
	  h_Gen_pT[2][6]->Fill(GenTaus[(int)ii].Pt(),wt1);  
	  h_Gen_eta[2][6]->Fill(GenTaus[(int)ii].Eta(),wt1); 	    
	  h_Gen_phi[2][6]->Fill(GenTaus[(int)ii].Phi(),wt1);	  
	} //end Gen Tau loop
	
	for(Long64_t ii=0; ii<Electrons->size(); ii++){
	  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > myele = Electrons[(int)ii];
	  h_Reco_pT[0][0]->Fill(Electrons[(int)ii].Pt(),wt1);
	  h_Reco_eta[0][0]->Fill(Electrons[(int)ii].Eta(),wt1);
	  h_Reco_phi[0][0]->Fill(Electrons[(int)ii].Phi(),wt1);
	  	  
	  if (ST<300) continue;
	  h_Reco_pT[0][4]->Fill(Electrons[(int)ii].Pt(),wt1);
	  h_Reco_eta[0][4]->Fill(Electrons[(int)ii].Eta(),wt1);
	  h_Reco_phi[0][4]->Fill(Electrons[(int)ii].Phi(),wt1);
	  
	  vector<myLV> v_electron;
	  if (Electrons_passIso) v_electron.push_back(Electrons[(int)ii]);
	  if (v_electron.size()!=0) continue;
	  h_Reco_pT[0][5]->Fill(v_electron[0].Pt(),wt1);
	  h_Reco_eta[0][5]->Fill(v_electron[0].Eta(),wt1);
	  h_Reco_phi[0][5]->Fill(v_electron[0].Phi(),wt1);
	    
	  if (isoElectronTracks != 0 || isoMuonTracks != 0 || isoPionTracks != 0) continue;
	  h_Reco_pT[0][6]->Fill(Electrons[(int)ii].Pt(),wt1);
	  //h_Reco_eta[0][6]->Fill(Electrons[(int)ii].Eta(),wt1);
	  h_Reco_phi[0][6]->Fill(Electrons[(int)ii].Phi(),wt1);	    
	} //end electron loop
	//if (jentry<1000)cout << endl;

	for(Long64_t ii=0; ii<Muons->size(); ii++){
	    ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > mymu = Muons[(int)ii];
	    h_Reco_pT[1][0]->Fill(Muons[(int)ii].Pt(),wt1);
	    h_Reco_eta[1][0]->Fill(Muons[(int)ii].Eta(),wt1);
	    h_Reco_phi[1][0]->Fill(Muons[(int)ii].Phi(),wt1);

	    if (ST<300) continue;
	    h_Reco_pT[1][4]->Fill(Muons[(int)ii].Pt(),wt1);
	    h_Reco_eta[1][4]->Fill(Muons[(int)ii].Eta(),wt1);
	    h_Reco_phi[1][4]->Fill(Muons[(int)ii].Phi(),wt1);

	    vector<myLV> v_muon;
	    if (Muons_passIso) v_muon.push_back(Muons[(int)ii]);
	    if (v_muon.size()!=0) continue;
	    h_Reco_pT[1][5]->Fill(v_muon[0].Pt(),wt1);
	    h_Reco_eta[1][5]->Fill(v_muon[0].Eta(),wt1);
	    h_Reco_phi[1][5]->Fill(v_muon[0].Phi(),wt1);

	    if (isoElectronTracks != 0 || isoMuonTracks != 0 || isoPionTracks != 0) continue;
	    h_Reco_pT[1][6]->Fill(Muons[(int)ii].Pt(),wt1);
	    h_Reco_eta[1][6]->Fill(Muons[(int)ii].Eta(),wt1);
	    h_Reco_phi[1][6]->Fill(Muons[(int)ii].Phi(),wt1);
	} //end Muon loop
	

	for(Long64_t ii=0; ii<Photons->size(); ii++){
	  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > mypho = Photons[(int)ii];
	} //end photon loop
	

  } // end jentry loop 

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
  cout << h_NHadJets[6]->Integral() << endl;


  
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
  // if(sBin==0){
  //   for(int i=0;i<METLowEdge.size()-1;i++){
  //     if(METLowEdge[i]<199.99) continue;
  //     int sBin1=sBin;
  //     m_i++;
  //     if(MET >= METLowEdge[i] && MET < METLowEdge[i+1]){ sBin = sBin+m_i;
  //       break; }
  //     else if(MET >= METLowEdge[METLowEdge.size()])  { sBin = 6;
  //       break; }
  //   }
  //   else if (sBin==7){
  //     for(int i=0;i<METLowEdge.size()-1;i++){
  // 	if(METLowEdge[i]<199.99) continue;
  // 	int sBin1=sBin;
  // 	m_i++;
  // 	if(MET >= METLowEdge[i] && MET < METLowEdge[i+1]){ sBin = sBin+m_i;
  // 	  break; }
  // 	else if(MET >= METLowEdge[METLowEdge.size()-1])  { sBin = 6;
  // 	  break; }

  //   }
    
  // }
  // else if(sBin==7 || sBin==13 || sBin==19 || sBin==25 || sBin==31){
  //   int sBin1=sBin;
  //   for(int i=0;i<METLowEdge.size()-1;i++){
  //     if(METLowEdge_1[i]<199.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_1[i] && MET < METLowEdge_1[i+1]){ sBin = sBin+m_i;
  //       break;}
  //     else if(MET >= METLowEdge_1[METLowEdge_1.size()-1])  { sBin = sBin+6;
  //       break; }
  //   }
  // }

  // else if(sBin==37){
  //   for(int i=0;i<METLowEdge.size()-1;i++){
  //     if(METLowEdge[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge[i] && MET < METLowEdge[i+1]){ sBin = sBin+m_i;break; }
  //     else if(MET >= METLowEdge[METLowEdge.size()-1])  { sBin = 44   ;break; }
  //     // else if(MET >= METLowEdge[METLowEdge.size()-1])  { sBin = 44   ;break; }                                                                                                                    

  //   }
  // }

  // else if(sBin==44){
  //   for(int i=0;i<METLowEdge.size()-1;i++){
  //     if(METLowEdge[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge[i] && MET < METLowEdge[i+1]){ sBin = sBin+m_i;break; }
  //     else if(MET >= METLowEdge[METLowEdge.size()-1])  { sBin = 52   ;break; }
  //     // else if(MET >= METLowEdge[METLowEdge.size()-1])  { sBin = 51   ;break; }                                                                                                                    

//     }
//   }
// -
  // int sBin=-100,m_i=0;
  // if(nbjets==0 ){
  //   if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
  //   else if(nHadJets==5 || nHadJets==6){ sBin=8;}
  //   else if(nHadJets>=7)               { sBin=15;}
  // }
  // else{
  //   if(nHadJets>=2 && nHadJets<=4)     { sBin=22;}
  //   else if(nHadJets==5 || nHadJets==6){ sBin=29;}
  //   else if(nHadJets>=7)               { sBin=36;}
  // }
  // if(sBin==0){
  //   for(int i=0;i<METLowEdge.size()-1;i++){
  //     if(METLowEdge[i]<99.99) continue;
  //     int sBin1=sBin;
  //     m_i++;
  //     if(MET >= METLowEdge[i] && MET < METLowEdge[i+1]){ sBin = sBin+m_i;
  // 	break; }
  //     else if(MET >= METLowEdge[METLowEdge.size()-1])  { sBin = 7;
  // 	break; }
  //   }
  // }
  // else if(sBin==8 || sBin==15 || sBin==22 || sBin==29 || sBin==36){
  //   int sBin1=sBin;
  //   for(int i=0;i<METLowEdge_1.size()-1;i++){
  //     if(METLowEdge_1[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_1[i] && MET < METLowEdge_1[i+1]){ sBin = sBin+m_i;
  // 	break;}
  //     else if(MET >= METLowEdge_1[METLowEdge_1.size()-1])  { sBin = sBin+6;
  // 	break; }
  //   }
  // }

  // else if(sBin==37){
  //   for(int i=0;i<METLowEdge.size()-1;i++){
  //     if(METLowEdge[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge[i] && MET < METLowEdge[i+1]){ sBin = sBin+m_i;break; }
  //     else if(MET >= METLowEdge[METLowEdge.size()-1])  { sBin = 44   ;break; }
  //     else if(MET >= METLowEdge[METLowEdge.size()-1])  { sBin = 44   ;break; }

  // }
  // }

  // else if(sBin==44){
  //   for(int i=0;i<METLowEdge.size()-1;i++){
  //     if(METLowEdge[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge[i] && MET < METLowEdge[i+1]){ sBin = sBin+m_i;break; }
  //     else if(MET >= METLowEdge[METLowEdge.size()-1])  { sBin = 52   ;break; }
  //     // else if(MET >= METLowEdge[METLowEdge.size()-1])  { sBin = 51   ;break; }

  //   }
  //}

  return sBin;
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


	
      
