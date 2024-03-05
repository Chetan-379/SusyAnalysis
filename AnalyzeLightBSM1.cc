#define AnalyzeLightBSM_cxx
#include "AnalyzeLightBSM.h"
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
#pragma link C++ class std::vector< std::vector >+; 
#pragma link C++ class std::vector< TLorentzVector >+;

using namespace TMVA;
int main(int argc, char* argv[])
{

  if (argc < 6) {
    cerr << "Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "which year dataset" <<" "<<"which Process"<< " "<<"which Lostlep bkg"<< " "<<"Which pho_ID"<<endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *sample=argv[4];
  const char *elec = argv[5];
  const char *phoID = argv[6];
  //TString pho_ID = phoID;

  AnalyzeLightBSM ana(inputFileList, outFileName, data,sample, elec,phoID);
  cout << "dataset " << data << " " << endl;
  cout<<"If analyzing the lost electron estimation ? "<<"  "<<elec<<endl;
  cout<<"Which pho_ID: "<<"\t"<<phoID<<endl;
  ana.EventLoop(data,inputFileList,sample,outFileName,elec,phoID);
  Tools::Instance();
  return 0;
}

void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList, const char *sample , const char *outFileName, const char *elec, const char* phoID) {
  cout<<"inside event loop"<<endl;
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;
  TString Elec_flag = elec;
  TString pho_ID_str = phoID;
  int pho_ID=-1;
  if(pho_ID_str.Contains("loose"))
    pho_ID=0;
  else if (pho_ID_str.Contains("medium"))
    pho_ID=1;
  else if (pho_ID_str.Contains("tight"))    
    pho_ID=2;
  else if (pho_ID_str.Contains("mva_wp90"))
    pho_ID=3;

  else if (pho_ID_str.Contains("mva_wp80"))
    pho_ID=4;
  //  TString pho_ID = phoID;
  bool lost_elec_flag=false;
  if(Elec_flag.Contains("Electron"))
    lost_elec_flag = true;
  else
    lost_elec_flag = false;
  
  float counter=0.0;
  TString s_data=data;

  TString s_sample= sample;
  
  Long64_t nbytes = 0, nb = 0; 
  int decade = 0;
  int count=0;
  double lumiInfb=137.19;
  int count_QCD=0;
  Long64_t nSurvived = 0,bkg_comp=0,MET_rej=0,nocut=0;
  int genphomatch_after=0,genphomatch_before=0;
  bool genphocheck=false;
  bool v17=true, v12=false;
  bool EWselection=true;
  bool Debug=false;
  bool higMET=false, highdphi=false;
  double deepCSVvalue=0,p0=0,p1=0,p2=0;
  if(s_data.Contains("2016")) lumiInfb=35.922;
  if(s_data.Contains("2017")) lumiInfb=41.529;
  if(s_data.Contains("2018")) lumiInfb=59.74;
  if(s_data.Contains("signal"))lumiInfb= 137.19; //since we only have 2018 dataset
  if(!s_sample.Contains("data"))
    {

      if(s_data.Contains("2016") && highdphi) {lumiInfb=35.922; p0=1.66539e+02; p1=8.13254e+01; p2=9.71152e-01; deepCSVvalue = 0.6321;}
      if(s_data.Contains("2017") && highdphi) {lumiInfb=41.529; p0=1.86744e+02; p1=6.74978e+01; p2=9.65333e-01; deepCSVvalue = 0.4941;}
      if(s_data.Contains("2018") && highdphi) {lumiInfb=59.74; p0=1.89868e+02; p1=6.60434e+01; p2=9.79618e-01; deepCSVvalue = 0.4184;}

      if(s_data.Contains("2016") && !highdphi) {lumiInfb=35.922; p0=1.67229e+02; p1=8.52729e+01; p2=8.29784e-01; deepCSVvalue = 0.6321;}
      if(s_data.Contains("2017") && !highdphi) {lumiInfb=41.529; p0=1.67641e+02; p1=1.21487e+02; p2=9.23864e-01; deepCSVvalue = 0.4941;}
      if(s_data.Contains("2018") && !highdphi) {lumiInfb=59.74; p0=1.45529e+02; p1=1.08431e+02; p2=9.27220e-01; deepCSVvalue = 0.4184;}

      if(s_data.Contains("FastSim") && s_data.Contains("2016")) lumiInfb=137.19;
    }
  if(s_data.Contains("data"))
    {
      if(s_data.Contains("2016")) {deepCSVvalue = 0.6321;}
      if(s_data.Contains("2017")) {deepCSVvalue = 0.4941;}
      if(s_data.Contains("2018")) {deepCSVvalue = 0.4184;}
    }

    
  std::string s_process = sample;
  double cross_section = getCrossSection(s_process);
  cout<<cross_section<<"\t"<<"analyzed process"<<endl;
 //  cout<<"Event"<<"\t"<<"par-ID "<<"\t"<<"parentID"<<"\t"<<"GenMET"<<"\t"<<"MET"<<"\t"<<"pT"<<"\t"<<"Eta"<<"\t"<<"Phi"<<"\t"<<"E"<<endl;
 //  float met=0.0,st=0.0, njets=0, btags=0,mTPhoMET=0.0,dPhi_PhoMET=0.0,dPhi_MetJet=0.0;
  
 //  // TMVA::Reader *reader1 = new TMVA::Reader();
 //  // reader1->AddVariable( "MET", &met );
 //  // reader1->AddVariable( "NhadJets", &njets );
 //  // reader1->AddVariable( "BTags", &btags );
 //  // reader1->AddVariable("mTPhoMET_",&mTPhoMET);
 //  // reader1->AddVariable("dPhi_PhoMET_",&dPhi_PhoMET);
 //  // reader1->AddVariable("dPhi_Met_Jet",&dPhi_MetJet);
 //  // reader1->AddVariable( "ST", &st );
 //  // reader1->BookMVA( "BDT_100trees_2maxdepth method", "TMVAClassification_BDT_100trees_2maxdepth.weights.xml");
  float nsurVived=0.0;
  int searchBin=0, Tfbins=0;
  float nCR_elec =0,nCR_mu=0,nCR_Tau=0,nSR_elec =0,nSR_mu=0,nSR_Tau=0, FailIso_Elec=0,FailIso_Mu=0, FailAccept_Elec=0,FailAccept_Mu=0, FailId_Elec=0,FailId_Mu=0, PassIso_Elec=0,PassIso_Mu=0, PassAccept_Elec=0,PassAccept_Mu=0, PassId_Elec=0,PassId_Mu=0, nfakeRatePho=0,wt_LL=0.0;

 //  // counters for events yields after each selection/rejection
  float nEvents_Selec[100]={};
  const char* out_nEventsTags[100] ={};
  const char* filetag[8]={"TTGJets_2018","TTGJets_2017","TTGJets_2016","Run2_TTGJets","WGJets_2018","WGJets_2017","WGJets_2016","Run2_WGJets"};
 // //  TFile* f_muon= new TFile("TF_allin1_LLEstimation_muon_newTFbins.root");
 // // //reading TF for electrons
 // //  TFile* f_elec =new TFile("TF_allin1_LLEstimation_electron_newTFbins.root");
 // //  char* histname = new char[2000];
 // //  TH1F* h_TF;
 // // auto sample1="";
 // // if(s_sample.Contains("Autumn18.WGJets_MonoPhoton_PtG-130")|| s_sample.Contains("Autumn18.WGJets_MonoPhoton_PtG-40to130"))
 // //   sample1 = "WGJets";
 // // else
 // //   sample1 = sample;
 // // if(lost_elec_flag)
 // //   {
 // //     sprintf(histname,"h_LL_TFbins_MEtNhadJetsBjetsBins_Elec__WGJets_2018");//%s_%s",sample1,data);
 // //     h_TF = (TH1F*)f_elec->Get(histname);
 // //     cout<<"Reading TF for lost electron estimation:  "<<"\t"<<histname<<endl;
 // //     for(int i=0; i<h_TF->GetNbinsX();i++)
 // //       {cout<<i<<"\t"<<h_TF->GetBinContent(i)<<endl;}
 // //   }
 // // else
 // //   {
 // //     sprintf(histname,"h_LL_TFbins_MEtNhadJetsBjetsBins_FailAcep_TauHadronic_WGJets_2018");//%s_%s",sample1,data);
 // //     cout<<"Reading TF for lost muon estimation:  "<<"\t"<<histname<<endl;
 // //     h_TF = (TH1F*)f_muon->Get(histname);
 // //     for(int i=0; i<h_TF->GetNbinsX();i++)
 // //       {cout<<i<<"\t"<<h_TF->GetBinContent(i)<<endl;}
 // //   }


 //  //  nentries=100;
  for (Long64_t jentry=0; jentry<nentries;jentry++)
     {
      double progress = 10.0 * jentry / (1.0 * nentries);
      int k = int (progress);
      if (k > decade)
        cout << 10 * k << " %" << endl;
      decade = k;
      
      // ===============read this entry == == == == == == == == == == ==                                                                        
      Long64_t ientry = LoadTree(jentry);
      if(Debug)
	cout<<"===load tree entry ==="<<endl;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(Debug)
        cout<<"===load tree entry ==="<<jentry<<endl;

      // ==========================================================================================\\
      // ====================================== Calculating weight per event ======================\\
      // ==========================================================================================\\
      cout<<"weight setting"<<"\t"<<jentry<<endl;
      if(s_sample.Contains("MCMC_86_7257_1750"))
        wt = (0.002991*lumiInfb*1000)/nentries; //gluino mass = 1750 GeV                                                                                                      
      else if(s_sample.Contains("MCMC_86_7257_1400"))
        wt = (0.0284*lumiInfb*1000)/nentries; //need to check luminosity                                                                                                      
      else if(s_sample.Contains("MCMC_86_7257_1500"))
        wt = (0.0157*lumiInfb*1000)/nentries;
      else if(s_sample.Contains("MCMC_86_7257_1600"))
    	wt = (0.00887*lumiInfb*1000)/nentries;
      else if(s_sample.Contains("MCMC_86_7257_1700"))
        wt = (0.00507*lumiInfb*1000)/nentries;
      else  if(s_sample.Contains("MCMC_86_7257_1800"))
        wt = (0.00293*lumiInfb*1000)/nentries;
      else if(s_sample.Contains("MCMC_86_7257_1900"))
        wt = (0.00171*lumiInfb*1000)/nentries;
      else if(s_sample.Contains("MCMC_86_7257_2000"))
        wt = (0.00101*lumiInfb*1000)/nentries;
      // else if(s_sample.Contains("WGJets")|| s_sample.Contains("ZNuNuGJets")|| s_sample.Contains("ZJets") ||s_sample.Contains("QCD") || s_sample.Contains("GJets")) //(s_data.Contains("2016") || s_data.Contains("2017") || s_data.Contains("2018") || s_sample.Contains("WGJets"))
      // 	{
	  
      // 	  std::string s_process = sample;
      // 	  double cross_section = getCrossSection(s_process);
      // 	  if(jentry==10)
      // 	    cout<<cross_section<<"\t"<<"analyzed process"<<"\t"<<s_process<<endl;
      // 	  wt = cross_section*lumiInfb*1000.0;//)/nentries; // Weight*lumiInfb*1000.0; //(cross_section*lumiInfb*1000.0)/nentries;
      // 	}
      else if (s_data.Contains("2016") || s_data.Contains("2017") || s_data.Contains("2018"))// || s_sample.Contains("WGJets"))
    	{
    	  wt = Weight*lumiInfb*1000.0;
    	}
      if(Debug)
      cout<<"alpa"<<"\t"<<jentry<<endl;
      //      cout<<"===load tree entry ==="<<"\t"<<jentry<<endl;
      // if(jentry%100000==0)
      // 		cout<<"event wieght=  "<<"\t"<<wt<<endl;
      //cout<<"Event-wt  "<<"\t"<<Weight<<"  calculated wt  "<<wt<<"\t cross section "<<CrossSection<<endl;
      if(jentry==10)
    	cout<<"event weight"<<"\t"<<wt<<"\t"<<CrossSection<<endl;

      h_selectBaselineYields_->Fill("No cuts, evt in 1/fb",wt);
      nocut++;
      if(Debug)
    	cout<<"Before filling lorentz vectors"<<endl;
      //filling the lorentz vector variables 

      if(s_sample.Contains("UL")){
    	//set objetc pointer

    	Electrons_v1.clear();
    	GenElectrons_v1.clear();
    	GenParticles_v1.clear();
    	GenMuons_v1.clear();
    	GenTaus_v1.clear();
    	GenJets_v1.clear();//_v1.clear();
    	Electrons_v1.clear();
    	Photons_v1.clear();
    	Muons_v1.clear();
    	Taus_v1.clear();
    	Jets_v1.clear();
    	if(Debug)
    	  cout<<jentry<<"\t"<<"After setting pointer  lorentz vectors"<<endl;

      for(int i=0; i<Electrons_;i++)
    	{
    	  TLorentzVector P1;
    	  P1.SetPtEtaPhiE(Electrons_fCoordinates_fPt[i],Electrons_fCoordinates_fEta[i],Electrons_fCoordinates_fPhi[i],Electrons_fCoordinates_fE[i]);
	  
    	  Electrons_v1.push_back(P1);
    	  //cout<<"Electrons_v1[0].Pt()  "<<"\t"<<Electrons_v1[0].Pt()<<endl;
    	  //	  cout<<(*Electrons).size()<<endl;
    	  //SetPtEtaPhiE(Electrons_fCoordinates_fPt[i],Electrons_fCoordinates_fEta[i],Electrons_fCoordinates_fPhi[i],Electrons_fCoordinates_fE[i]);
    	}
      if(Debug)
    	cout<<jentry<<"\t"<<"After filling reco nelectrons lorentz vectors"<<endl;

      for(int i=0; i<GenElectrons_;i++)
    	{
    	  TLorentzVector P1;
    	  P1.SetPtEtaPhiE(GenElectrons_fCoordinates_fPt[i],GenElectrons_fCoordinates_fEta[i],GenElectrons_fCoordinates_fPhi[i],GenElectrons_fCoordinates_fE[i]);
    	  GenElectrons_v1.push_back(P1);
    	}
      if(Debug)
        cout<<jentry<<"\t"<<"After filling gen electrons lorentz vectors"<<"\t"<<GenJets_<<endl;
      //cout<<"GenJets size - "<<"\t"<<GenJets_<<endl;
      for(int i=0; i<GenJets_;i++)
    	{
    	  // cout<<i<<"\t"<<GenJets_fCoordinates_fPt[i]<<"\t"<<GenJets_fCoordinates_fEta[i]<<"\t"<<GenJets_fCoordinates_fPhi[i]<<"\t"<<GenJets_fCoordinates_fE[i]<<endl;
    	  TLorentzVector P1;
    	  P1.SetPtEtaPhiE(GenJets_fCoordinates_fPt[i],GenJets_fCoordinates_fEta[i],GenJets_fCoordinates_fPhi[i],GenJets_fCoordinates_fE[i]); //Electrons_fCoordinates_fPt[i],Electrons_fCoordinates_fEta[i],Electrons_fCoordinates_fPhi[i],Electrons_fCoordinates_fE[i]);
    	  //cout<<"P1.Pt() = "<<P1.Pt()<<"\t"<<endl;
    	  GenJets_v1.push_back(P1);
    	  //cout<<"GenJets_v1[0].Pt()  "<<"\t"<<GenJets_v1[0].Pt()<<endl;

    	}
      if(Debug)
        cout<<jentry<<"\t"<<"After filling gen Jets lorentz vectors"<<endl;

      for(int i=0;i<GenMuons_;i++)
      	{
    	  TLorentzVector P1;
          P1.SetPtEtaPhiE(GenMuons_fCoordinates_fPt[i],GenMuons_fCoordinates_fEta[i],GenMuons_fCoordinates_fPhi[i],GenMuons_fCoordinates_fE[i]);
    	  GenMuons_v1.push_back(P1);
    	  //cout<<"GenMuons_v1[0].Pt() "<<GenMuons_v1[0].Pt()<<endl;
      	}
      for(int i=0;i<GenTaus_;i++){
    	TLorentzVector P1;
          P1.SetPtEtaPhiE(GenTaus_fCoordinates_fPt[i],GenTaus_fCoordinates_fEta[i],GenTaus_fCoordinates_fPhi[i],GenTaus_fCoordinates_fE[i]);
    	  GenTaus_v1.push_back(P1);
      }
      for(int i=0;i<GenParticles_;i++)
      	{
    	  TLorentzVector P1;
          P1.SetPtEtaPhiE(GenParticles_fCoordinates_fPt[i],GenParticles_fCoordinates_fEta[i],GenParticles_fCoordinates_fPhi[i],GenParticles_fCoordinates_fE[i]);
    	  //if(Debug)
    	  //cout<<"Genparticles P1.Pt()  "<<P1.Pt()<<"\t"<<GenParticles_fCoordinates_fPt[i]<<endl;
    	  GenParticles_v1.push_back(P1);
    	  //if(Debug)
    	  //cout<<"Genparticles_v1.Pt()  "<<GenParticles_v1[int(GenParticles_v1.size())-1].Pt()<<endl;

      	}
      for(int i=0;i<Photons_;i++)
    	{
    	  TLorentzVector P1;
          P1.SetPtEtaPhiE(Photons_fCoordinates_fPt[i],Photons_fCoordinates_fEta[i],Photons_fCoordinates_fPhi[i],Photons_fCoordinates_fE[i]);       Photons_v1.push_back(P1);
    	}
      if(Debug)
        cout<<jentry<<"\t"<<"After filling photons lorentz vectors"<<endl;

      for(int i=0;i<Muons_;i++)
      	{
    	  TLorentzVector P1;
          P1.SetPtEtaPhiE(Muons_fCoordinates_fPt[i],Muons_fCoordinates_fEta[i],Muons_fCoordinates_fPhi[i],Muons_fCoordinates_fE[i]);
      	  Muons_v1.push_back(P1);      	}

      for(int i=0;i<Jets_;i++)	
      	{ TLorentzVector P1;//(0,0,0,0);
    	  P1.SetPtEtaPhiE(Jets_fCoordinates_fPt[i],Jets_fCoordinates_fEta[i],Jets_fCoordinates_fPhi[i],Jets_fCoordinates_fE[i]);
    	  Jets_v1.push_back(P1);
    	  //	  Jets->push_back(P1);
    	}
      if(Debug)
    	cout<<jentry<<"\t"<<"After filling  reco Jets"<<endl;

      }
      int count=0;
      int branch_size = GenParticles_v1.size();
      if(Debug)
        cout<<"Gen particles  "<<"\t"<<branch_size<<"\t"<<GenParticles_v1[0].Pt()<<endl;

      int ele_branch=Electrons_;
      if(Debug)
        cout<<"Electrons "<<"\t"<<Electrons_<<endl;

      int pho_branch = Photons_;
      if(Debug)
    	        cout<<"Photons "<<"\t"<<Photons_<<endl;
      
    //   //variables to be used in getting dR
      float iGenEle_eta =99999.0, iGenEle_phi=99999.0, iRecEle_eta =99999.0, iRecEle_phi=99999.0, iRecphot_eta =99999.0, iRecphot_phi=99999.0,iGen_Wpt=99999.0,iGen_Wp=99999.0,iGen_WEta=99999.0,iGen_WPhi=99999.0;
      float dR_Ele=0, min_dR=9999, min_dR_pho=9999;
      int pdgID=0, parentId=0, status;
      int Ids[4]={11,13,15,22};
      const char* ids[4]={"e","mu","tau","pho"};
      // loop over gen particles
      int count_genElec=0,count_genMu=0,count_genTau=0, igenPho=0, count_genEl_tau=0,count_genMu_tau=0,count_haddecay=0;
      TLorentzVector genPho_W,genElec_W,genMu_W,genTau_W,genElec_Tau,genMu_Tau,genW,genNuElec_W,genNuMu_W,genNuTau_W,genElecNu_Tau,genMuNu_Tau;
      int elec_reco=0,elec_reco0_before=0,elec_reco1_before=0,muon_reco=0,elec_gen3=0,elec_gen2=0, elec_gen=0, muon_gen=0,elec_reco0=0,elec_reco1=0,evtSurvived_preselec=0,elec_reco2=0,elec2_reco=0,survived_vetohad=0,elec_reco1_CR=0,survived_elecge1=0,events_cr=0,events_sr=0,total=0,remain=0,elec_reco0_genel=0,efakepho=0,ele=0,genphomatch_after=0,genphomatch_before=0,elec_gen4=0,gentauhad2=0,lep2=0,lep=0;
      TLorentzVector leadGenPho,genPho,genLep,v_lep1,v_lep2,v_genMu,v_genNuMu,v_genElec,v_genNuElec,v_genTau,v_genNuTau,v_genW,v_genNu;
      //      vector<TLorentzVector> v_genPho, v_genLep,v_genElec,v_genMu,v_genTau, v_genAccepLep,v_genAccepElec,v_genAccepMu,v_genAccepTau;
    // //   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // //   ///////////////////////////////////   Making gen level collections of all particles  ////////////////////////////////////////////////
    // //   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if(Debug)
    	cout<<"cp - before gen particles"<<endl;
      // for (int igen=0; igen<branch_size; igen++)
      // 	{	  
      // 	  pdgID= (*GenParticles_PdgId)[igen];
      // 	  parentId=(*GenParticles_ParentId)[igen];
      // 	  // if(Debug)
      // 	  //   {	    cout<<"printing Gen particles parent ID and pdgID"<<endl;
      // 	  //     cout<<parentId<<"\t"<<pdgID<<endl;
      // 	  //     cout<<GenParticles_fCoordinates_fPt[igen]<<"\t"<<GenParticles_fCoordinates_fEta[igen]<<endl;}
	  
      // 	  if(abs(pdgID)==24)
      // 	    {h_W_pT->Fill(GenParticles_v1[igen].Pt(),wt); h_W_Eta->Fill(GenParticles_v1[igen].Eta(),wt);h_W_Phi->Fill(GenParticles_v1[igen].Phi(),wt); h_W_E->Fill(GenParticles_v1[igen].E(),wt); v_genW = GenParticles_v1[igen];}
      // 	  if(abs(pdgID)==11 && abs(parentId)==24){h_elec_pT->Fill(GenParticles_v1[igen].Pt(),wt); h_elec_Eta->Fill(GenParticles_v1[igen].Eta(),wt);h_elec_Phi->Fill(GenParticles_v1[igen].Phi(),wt); h_elec_E->Fill(GenParticles_v1[igen].E(),wt);v_genElec = GenParticles_v1[igen];}
      // 	  if(abs(pdgID)==12 && abs(parentId)==24){h_Nu_pT->Fill(GenParticles_v1[igen].Pt(),wt); h_Nu_Eta->Fill(GenParticles_v1[igen].Eta(),wt);h_Nu_Phi->Fill(GenParticles_v1[igen].Phi(),wt); h_Nu_E->Fill(GenParticles_v1[igen].E(),wt);v_genNuElec = GenParticles_v1[igen];}
      // 	  if(abs(pdgID)==13 && abs(parentId)==24){h_mu_pT->Fill(GenParticles_v1[igen].Pt(),wt); h_mu_Eta->Fill(GenParticles_v1[igen].Eta(),wt);h_mu_Phi->Fill(GenParticles_v1[igen].Phi(),wt); h_mu_E->Fill(GenParticles_v1[igen].E(),wt); v_genMu=GenParticles_v1[igen];}
      // 	  if(abs(pdgID)==15 && abs(parentId)==24){h_tau_pT->Fill(GenParticles_v1[igen].Pt(),wt); h_tau_Eta->Fill(GenParticles_v1[igen].Eta(),wt);h_tau_Phi->Fill(GenParticles_v1[igen].Phi(),wt); h_tau_E->Fill(GenParticles_v1[igen].E(),wt); v_genTau=GenParticles_v1[igen];}
      // 	  if(abs(pdgID)==14 && abs(parentId)==24) v_genNuMu = GenParticles_v1[igen];
      // 	  if(abs(pdgID)==16 && abs(parentId)==24) v_genNuTau = GenParticles_v1[igen];
      // 	  if((abs(pdgID)==14 || abs(pdgID)==16 || abs(pdgID)==12) && abs(parentId)==24) v_genNu= GenParticles_v1[igen];
      // 	}

      // h_2d_elec_nu_pT->Fill(v_genElec.Pt(),v_genNuElec.Pt(),wt);
      // h_2d_mu_nu_pT->Fill(v_genMu.Pt(),v_genNuMu.Pt(),wt);
      // h_2d_Tau_nu_pT->Fill(v_genTau.Pt(),v_genNuTau.Pt(),wt);
      // h_2d_WvsElec_pT->Fill(v_genW.Pt(),v_genElec.Pt(),wt);
      // h_2d_WvsElecNu_pT->Fill(v_genW.Pt(),v_genNuElec.Pt(),wt);
      // h_2d_WvsMu_pT->Fill(v_genW.Pt(),v_genMu.Pt(),wt);
      // h_2d_WvsMuNu_pT->Fill(v_genW.Pt(),v_genNuMu.Pt(),wt);
      // h_2d_WvsTau_pT->Fill(v_genW.Pt(),v_genTau.Pt(),wt);
      // h_2d_WvsTauNu_pT->Fill(v_genW.Pt(),v_genNuTau.Pt(),wt);
      
      // h_2d_elec_nu_Eta->Fill(v_genElec.Eta(),v_genNuElec.Eta(),wt);
      // h_2d_mu_nu_Eta->Fill(v_genMu.Eta(),v_genNuMu.Eta(),wt);
      // h_2d_Tau_nu_Eta->Fill(v_genTau.Eta(),v_genNuTau.Eta(),wt);
      // h_2d_WvsElec_Eta->Fill(v_genW.Eta(),v_genElec.Eta(),wt);
      // h_2d_WvsElecNu_Eta->Fill(v_genW.Eta(),v_genNuElec.Eta(),wt);
      // h_2d_WvsMu_Eta->Fill(v_genW.Eta(),v_genMu.Eta(),wt);
      // h_2d_WvsMuNu_Eta->Fill(v_genW.Eta(),v_genNuMu.Eta(),wt);
      // h_2d_WvsTau_Eta->Fill(v_genW.Eta(),v_genTau.Eta(),wt);
      // h_2d_WvsTauNu_Eta->Fill(v_genW.Eta(),v_genNuTau.Eta(),wt);
      // h_2d_nu_MET_pT->Fill(v_genNu.Pt(),MET,wt);
      h_MET->Fill(MET,wt);
      h_GenMET->Fill(GenMET,wt);
      h_GenJets->Fill(GenJets_,wt);
      if(Debug)
      cout<<"===load tree entry check1 - after filling gen level at entry ==="<<"\t"<<jentry<<endl;           

    //   //selecting electron and muons
      int nGenMu=0,nGenEle=0,nGenTau=0,nGenHad=0,nGenLep=0,nGenEle_tau=0,nGenMu_tau=0,nGentau_lep=0,nGentau_had=0,nGentau_had1=0,nGenTau1=0,nGenEle1=0,nGenEle_tau1=0,nGenMu1=0,nGenMu_tau1=0,o=0,nelec_reco=0,nmu_reco=0;
      TLorentzVector genPho1,genEle1,genMu1,genTau1,recElec, recMu;//,genHad1,genLep1,gentau_lep1,gentau_had1,recElec, recMu;
      vector<TLorentzVector> v_genEle1, v_genPho1, v_genMu1,v_genTau1,v_genHad1,v_gentau_lep1,v_gentau_had1,v_gentau_had2,v_genTau2,v_genEle2,v_genMu2,v_genLep2,v_genLep1,v_recEle,v_recMu;
      v_genEle2.clear();
      v_genEle1.clear();
      v_genMu1.clear();
      v_genMu2.clear();
      v_genTau1.clear();
      v_genTau2.clear();

      int leadGenPhoIdx=-100, mu_index=-100;
      int pass_accep_elec=0,fail_accep_elec=0,fail_isoEle=0,pass_isoElec=0,fail_IdElec=0,pass_IdElec=0;
      int pass_accep_mu=0,fail_accep_mu=0,fail_isoMu=0,pass_isoMu=0,fail_IdMu=0,pass_IdMu=0;
      if(!s_sample.Contains("data")){
    	for(int i=0;i<GenParticles_;i++){
    	  if(GenParticles_v1[i].Pt()!=0){
    	    if((abs((*GenParticles_PdgId)[i])==22) && ((abs((*GenParticles_ParentId)[i])<=100) || ((*GenParticles_ParentId)[i]==2212) ) && (*GenParticles_Status)[i]==1 )
    	      {
    		leadGenPhoIdx = i;
    		genPho1 = (GenParticles_v1[i]);
    		v_genPho1.push_back(genPho1);
    	      }
    	  }
    	}
    	bool hadtau = false;
    	for(int i=0 ; i < GenElectrons_; i++)
    	  {
    	    if(GenElectrons_v1[i].Pt()!=0)// && (*GenElectrons)[i].Pt()!=inf)                                                                                                                                          
    	      {
    		nGenEle1++;
    		genEle1 = (GenElectrons_v1[i]);
    		v_genEle2.push_back(genEle1);
    		v_genLep2.push_back(genEle1);
    		//h_GenpT[3]->Fill(genEle1.Pt()); h_GenEta[3]->Fill(genEle1.Eta());
    	      }
    	  }
    	for(int i=0 ; i < GenMuons_; i++)
    	  {
    	    if(GenMuons_v1[i].Pt()!=0)// && (*GenMuons)[i].Pt()!=inf)                                                                                                                                                  
    	      {
    		nGenMu1++;
    		genMu1 = (GenMuons_v1[i]);
    		v_genMu2.push_back(genMu1);
    		v_genLep2.push_back(genMu1);
    		//h_GenpT[4]->Fill(genMu1.Pt()); h_GenEta[4]->Fill(genMu1.Eta());
    	      }
    	  }
    	//cout<<"entry:"<<"\t"<<jentry<<"\t"<<nGenMu1<<"\t"<<"genMuon size"<<"\t"<<v_genMu2.size()<<endl;                                                                                                              

    	for(int i=0 ; i < GenTaus_; i++)
    	  {
    	    if(GenTaus_v1[i].Pt()!=0) // && (*GenTaus)[i].Pt()!=inf)                                                                                                                                                   
    	      {
    		nGenTau1++;
    		genTau1 = (GenTaus_v1[i]);
    		v_genTau2.push_back(genTau1);
    		//v_genLep2.push_back(genTau1);                                                                                                                                                                        
    		//h_GenpT[5]->Fill(genTau1.Pt()); h_GenEta[5]->Fill(genTau1.Eta());
    		if((*GenTaus_had)[i])
    		  nGentau_had1++;

    	      }
    	  }
      }
      sortTLorVec(&v_genTau2);
      sortTLorVec(&v_genMu2);
      sortTLorVec(&v_genEle2);
      if(Debug)
    	cout<<"=== load entry "<<jentry<<"\t"<<"GenTaus_  "<<GenTaus_<<endl;
      //      int nlep=0;
      int total_lost_el = 0,cr_el=0,sr_el,e_index=-1,nlep=0, NgenElec=0,leadGenPhoIdx1=0;
      bool Elec_passEtacut=false,Elec_passpTcut=false,Elec_passAccep = false,Elec_failAccep= false, Elec_passId= false, Elec_failId= false,Elec_passIso = false,Elec_failIso = false;
      bool Mu_passEtacut=false,Mu_passpTcut=false, Mu_passAccep= false,Mu_failAccep= false, Mu_passId= false, Mu_failId= false,Mu_passIso = false,Mu_failIso = false;
      int fail_isoElec=0;
      nlep=0; nelec_reco=0;
      for(int i=0;i<Electrons_;i++)
    	{
    	  if(nelec_reco>0) continue;
    	  if((Electrons_v1[i].Pt()>10) && abs(Electrons_v1[i].Eta()) < 2.5 && (*Electrons_passIso)[i]==1)
    	    {
    	      pass_isoElec++; Elec_passIso = true;nelec_reco++; nlep++; e_index=i; recElec=Electrons_v1[i]; v_recEle.push_back(Electrons_v1[i]) ;
    	      h_recoElec_pT->Fill(recElec.Pt(),wt);
    	      h_recoElec_Eta->Fill(recElec.Eta(),wt);
    	      h_recoElec_Eta->Fill(recElec.Eta(),wt);
    	      h_recoElec_Phi->Fill(recElec.Phi(),wt);

    	    }
    	  else if((*Electrons_passIso)[i]!=1)
    	    {
    	      fail_isoElec++; Elec_passIso = false;}

    	}
      sortTLorVec(&v_recEle);
      //      if(Debug)      cout<<"entry: "<<"\t"<<jentry<<" "<<nlep<<"\t"<<"reco e size"<<" "<<v_recEle.size()<<" "<<Electrons->size()<<"\t"<<NElectrons<<endl;

      for(int i=0;i<Muons_;i++)
    	{
    	  if((Muons_v1[i].Pt()>10) && abs(Muons_v1[i].Eta()) < 2.5 && nmu_reco<1 && (*Muons_passIso)[i]==1)
    	    {
    	      pass_isoMu++; Mu_passIso = true;nmu_reco++; nlep++; mu_index=i; recMu=Muons_v1[i]; v_recMu.push_back(Muons_v1[i]);
    	      h_recoMu_pT->Fill(recMu.Pt(),wt);
    	      h_recoMu_Eta->Fill(recMu.Eta(),wt);
    	      h_recoMu_Eta->Fill(recMu.Eta(),wt);
    	      h_recoMu_Phi->Fill(recMu.Phi(),wt);

    	    }
    	  else if((*Muons_passIso)[i]!=1)
    	    {
    	      fail_isoMu++; Mu_passIso = false;}
    	}

      //      if(Debug && v_recMu.size()==1)      cout<<"entry: "<<"\t"<<jentry<<" "<<nlep<<"\t"<<"reco e size"<<" "<<v_recEle.size()<<" "<<Muons->size()<<"\t"<<NMuons<<endl;
      if(Electrons_==0 && nelec_reco!=0) //checking if the contribution to SR with no reco e- is coming up or not?                              
    	cout<<nelec_reco<<endl;
      float ratio =0.0, ratio1=0.0, mindr_genElecPho=-999, mindr_=-9999;


      int count_genEle=0,count_recEle=0;
      // TLorentzVector genPho1,genLep1;
      // int leadGenPhoIdx=-100;
      vector<TLorentzVector> goodPho_n;
      vector<int> goodPhoIndx_n;
      for(int iPho=0;iPho<Photons_;iPho++){
    	if((*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001)) {
    	  goodPho_n.push_back( Photons_v1[iPho] );
    	  goodPhoIndx_n.push_back(iPho);
    	}
      }
      int nPhotons = goodPho_n.size();
      TLorentzVector bestPhoton=getBestPhoton(pho_ID);      
      // *******************  Selecting Jet objects ********************************//
       int minDRindx=-100,phoMatchingJetIndx=-100,hadJetID=-999,bJet1Idx=-100,nHadJets=0;
      double minDR=99999,ST=0,Ht=0;
      vector<TLorentzVector> hadJets,bjets;
      vector<int> jetMatchindx;
      bool recoJetMatch_recoPho=false, genJetMatch_recoPho=false;
      if(Debug)
        cout<<"===load tree entry check2 at entry ==="<<"\t"<<jentry<<endl;
      for(int i=0;i<Jets_;i++)
    	{
    	  if( (Jets_v1[i].Pt() > 30.0) && (abs(Jets_v1[i].Eta()) <= 2.4) ){
    	    double dR=bestPhoton.DeltaR(Jets_v1[i]);
    	    if(dR<minDR){minDR=dR;minDRindx=i;}
    	  }
      }
      if(Debug)
    	cout<<"===load tree entry  ==="<<"\t"<<jentry<<"\t"<<"Jets check == "<<minDR<<endl;
      
    for(int i=0;i<Jets_;i++)
      { //if(Debug)
    	  //cout<<"  = Jets.Pt()  ==  "<<Jets_v1[i].Pt()<<"\t"<< " = Jets.Eta() == "<<Jets_v1[i].Eta()<<endl;
    	if( (Jets_v1[i].Pt() > 30.0) && (abs(Jets_v1[i].Eta()) <= 2.4) ){
    	  //if(Debug)
    	  //cout<< "==== loadjets ==="<<"\t"<<i<<"\t"<<minDR<<endl;
          if( !(minDR < 0.3 && i==minDRindx) )
    	    {
	      
    	      hadJetID= (*Jets_ID)[i];
    	      if(hadJetID)
    		{		 
    		  hadJets.push_back(Jets_v1[i]);
    		  // hadJets_hadronFlavor.push_back((*Jets_hadronFlavor)[i]);
    		  // hadJets_HTMask.push_back((*Jets_HTMask)[i]);
    		  // hadJets_bJetTagDeepCSVBvsAll.push_back((*Jets_bJetTagDeepCSVBvsAll)[i]);
    		  // if(q==1) leadjet_qmulti=(*Jets_chargedMultiplicity)[q];
    		  // if(q==1) leadjet_Pt=(*Jets)[q].Pt();
    		  if((*Jets_bJetTagDeepCSVBvsAll)[i] > deepCSVvalue){
    		    bjets.push_back(Jets_v1[i]); bJet1Idx = i;}
    	      // hadJets.push_back((*Jets)[i]);
    		  jetMatchindx.push_back(i);
    		}	  
    	    }
    	}
      }
    if(hadJets.size()==0) continue;
    if(Debug)
      cout<<"===load tree entry ===  "<<"\t"<<jentry<<"\t"<<"No of B-Jets ===  "<<bjets.size()<<endl;

    BTags = bjets.size();
    if( minDR<0.3 ) {phoMatchingJetIndx=minDRindx; recoJetMatch_recoPho=true;}
    double genmindr=99999, recojetmindr=99999;
    if(Debug)
      cout<<"====load entry === "<<jentry<<"\t"<<"phoMatchingJetIndx ===  "<<phoMatchingJetIndx<<endl;
    genmindr = MinDr(bestPhoton,GenJets_v1);
    recojetmindr =  minDR;// MinDr(bestPhoton,hadJets);
    if(genmindr<0.3)
      genJetMatch_recoPho = true;
    if(Debug)
      cout<<"===load tree entry === "<<"\t"<<jentry<<" hadJets.size()  "<< hadJets.size()<<endl;    
    for(int i=0;i<hadJets.size();i++){
      if( (abs(hadJets[i].Eta()) < 2.4) ){ST=ST+(hadJets[i].Pt());Ht=Ht+(hadJets[i].Pt());}
      if( (abs(hadJets[i].Eta()) < 2.4) ){nHadJets++;}
    }
    sortTLorVec(&hadJets);
    // // ********* This is to account into all visible energy: Adding photon matching with Jet ******************//
    if( minDR<0.3 ) ST=ST+bestPhoton.Pt(); 
    // //transverse mass of best photon+MET                                                                                  
    double mTPhoMET=sqrt(2*(bestPhoton.Pt())*MET*(1-cos(DeltaPhi(METPhi,bestPhoton.Phi()))));
    // // ******* Dphi between MET and Best Photon     ******************** //                                                                 
    double dPhi_PhoMET= abs(DeltaPhi(METPhi,bestPhoton.Phi()));
    if(Debug)
      cout<<"===load tree entry ==="<<"\t"<<jentry<<" mTPhoMET "<<mTPhoMET<<" METPhi  "<<bestPhoton.Eta()<<" Photons_ "<<Photons_<<" hadJets.size() " <<"\t"<<hadJets.size()<<" NElectrons "<<NElectrons<< "  Jets  "<<Jets_<<"\t"<< "Jets size  "<<Jets_v1.size()<<endl;
    h_leadJet_pT->Fill(hadJets[0].Pt(),wt);
    h_subleadJet_pT->Fill(hadJets[1].Pt(),wt);    
    h_leadJet_Eta->Fill(hadJets[0].Eta(),wt);
    h_subleadJet_Eta->Fill(hadJets[1].Eta(),wt);
    h_Njets[0]->Fill(hadJets.size(),wt);
    h_Nbjets[0]->Fill(BTags,wt);
    h_MET_[0]->Fill(MET,wt);
    //h_PhotonPt[0]->Fill(bestPhoton.Phi(),wt);
    h_Mt_PhoMET[0]->Fill(mTPhoMET,wt);
    h_dPhi_PhoMET[0]->Fill(dPhi_PhoMET,wt);
    h_St[0]->Fill(ST,wt);
    h_HT[0]->Fill(HT,wt);

    TLorentzVector Met;
    Met.SetPtEtaPhiE(MET,0,METPhi,0);
    Double_t deta_jet_pho= 0.0,deta_jet_met=0.0,deta_met_pho=0.0;
    if(    goodPho_n.size()==0 || Photons_==0) continue;
    h_recoPho_pT->Fill(bestPhoton.Pt(),wt); h_recoPho_Eta->Fill(bestPhoton.Eta(),wt);
    h_recoPho_Phi->Fill(bestPhoton.Phi(),wt);
    h_PhotonPt[0]->Fill(bestPhoton.Pt(),wt);
    deta_jet_pho = abs(bestPhoton.Eta()-hadJets[0].Eta());
    double mT= 0.0, dPhi_METjet1=5, dPhi_METjet2=5, dPhi_phojet1=5, dPhi_phojet2=5, dPhi_phoMET=5;
    if(hadJets.size() > 0) dPhi_phojet1 = abs(bestPhoton.DeltaPhi(hadJets[0]));
    if(hadJets.size() > 1) dPhi_phojet2 = abs(bestPhoton.DeltaPhi(hadJets[1]));
    dPhi_phojet2 = abs(bestPhoton.DeltaPhi(hadJets[1]));
    dPhi_phoMET = abs(bestPhoton.DeltaPhi(Met));

    dPhi_METjet1 = abs(Met.DeltaPhi(hadJets[0]));
    dPhi_METjet2 = abs(Met.DeltaPhi(hadJets[1]));
    dPhi_phoMET = abs(bestPhoton.DeltaPhi(Met));
    dPhi_METjet1=3.8,dPhi_METjet2=3.8;                                                                                                         
    Met.SetPtEtaPhiE(MET,0,METPhi,0);
    dPhi_phojet1 = abs(bestPhoton.DeltaPhi(hadJets[0]));
    dPhi_phojet2 = abs(bestPhoton.DeltaPhi(hadJets[1]));
    dPhi_phoMET = abs(bestPhoton.DeltaPhi(Met));
    dPhi_METjet1 = abs(DeltaPhi(METPhi,(hadJets)[0].Phi()));
    dPhi_METjet2 = abs(DeltaPhi(METPhi,(hadJets)[1].Phi()));
    if(Debug)
      cout<<"===load tree entry ==="<<"\t"<<jentry<<" hadJets.size() "<<hadJets.size()<<endl;



    h_selectBaselineYields_->Fill("all",wt);
    if(bestPhotonIndxAmongPhotons<0) continue;
    bool bestPhoHasPxlSeed=true;
    if((*Photons_hasPixelSeed)[bestPhotonIndxAmongPhotons]<0.001) bestPhoHasPxlSeed=false;
    if( bestPhoHasPxlSeed ) continue;
    if(Debug)
      cout<<"===load tree entry ==="<<"\t"<<jentry<<" bestPhoHasPxlSeed "<<bestPhoHasPxlSeed<<endl;
    if(hadJets.size() == 0) continue;

    h_selectBaselineYields_->Fill("hadsize!=0",wt);
    if(bestPhoton.Pt()>40)    h_selectBaselineYields_->Fill("#gamma_{pt}>40",wt);
    else continue;
    if (nHadJets>=2)h_selectBaselineYields_->Fill("njets>=2",wt);
    else continue;
    if(MET>100) 
      h_selectBaselineYields_->Fill("MET>100",wt);
    else continue;
    if(ST>300)
      h_selectBaselineYields_->Fill("ST>300",wt);
    else continue;
    if(s_data.Contains("2016") || s_data.Contains("2017") || s_data.Contains("2018"))
      {
        if(PFCaloMETRatio >=  5) continue;
        if(s_sample.Contains("data"))
          {
            double ecalBadCalibReducedExtraFilter=1;
    	    if(!(PrimaryVertexFilter==1 && globalSuperTightHalo2016Filter == 1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter && BadChargedCandidateFilter && eeBadScFilter==1 && NVtx>0 && ecalBadCalibReducedExtraFilter==1)) continue;
          }
        if(!(s_sample.Contains("data"))){
          eeBadScFilter=1;
    	  double ecalBadCalibReducedExtraFilter=1;
          if(!(PrimaryVertexFilter==1 && globalSuperTightHalo2016Filter == 1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter && NVtx>0&& eeBadScFilter==1 && ecalBadCalibReducedExtraFilter==1)) continue;}
         
      }
    h_selectBaselineYields_->Fill("MetCleaning",wt);
    
    if(dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3 )
      {
        h_selectBaselineYields_->Fill("dPhi1 & dPhi2 >= 0.3",wt);
      }
    else continue;
    if(phoMatchingJetIndx>=0 && (Jets_v1[phoMatchingJetIndx].Pt())/(bestPhoton.Pt()) < 1.0) continue;
    if(phoMatchingJetIndx<0) continue;
    if(higMET)
      { 
    	if (bestPhoton.Pt()<100)continue;
    	if(MET<250) continue;                                                                                              
      }
    if(Debug)
      cout<<"===loading tree entry === "<<jentry<<"\t"<<"  ===after preselection  ==="<<endl;
    h_selectBaselineYields_v1->Fill("Pre-Selection",wt);
    h_Njets[1]->Fill(hadJets.size(),wt);
    h_Nbjets[1]->Fill(BTags,wt);
    h_MET_[1]->Fill(MET,wt);
    h_PhotonPt[1]->Fill(bestPhoton.Phi(),wt);
    h_Mt_PhoMET[1]->Fill(mTPhoMET,wt);
    h_dPhi_PhoMET[1]->Fill(dPhi_PhoMET,wt);
    h_St[1]->Fill(ST,wt);
    h_HT[1]->Fill(HT,wt);

    // nEvents_Selec[9]+=wt;    
    // float madMinPhotonDeltaR=0.0;
    // //    h_madminPhotonDeltaR_->Fill(madMinPhotonDeltaR,wt);
    // //cout<<hasGenPromptPhoton<<endl;
    // h_hasGenPromptPhoton->Fill(hasGenPromptPhoton,wt);
    // double mindr_Pho_genlep=getGenLep(bestPhoton);
    // //h_phoPt_promptPho->Fill(
    // // photon-pT distribution checks
    // if(s_sample.Contains("TTJets-HT") && madHT<600) continue;    
    
    // if(s_sample.Contains("TTJets-inc") && madHT>600) continue;
    // if(!genphocheck)
    //   {
    //     genphomatch_before++;
    //     double mindr_Pho_genlep=getGenLep(bestPhoton);
    //     if( s_sample.Contains("TTG") )
    //       {
    //         if(!hasGenPromptPhoton)
    //           {		
    // 		h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);
    //             if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
    //           }
    //         else if(hasGenPromptPhoton)
    //           {
    // 		h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
    //             if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )) {//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt); 
    // 		  continue;}
    // 		else
    // 		  {
    // 		    if(madMinPhotonDeltaR >= 0.5)
    // 		      h_selectBaselineYields_v1->Fill("mindR(q/g, #gamma)",wt);
    // 		    if(mindr_Pho_genlep >=0.5)
    // 		      h_selectBaselineYields_v1->Fill("mindR(l, #gamma)",wt);
    // 		  }
    //           }
    //       }//Gen prompt                                                                                                                                                                                  
    //     if(s_sample.Contains("WG"))
    //       {
    //         if(!hasGenPromptPhoton)
    //           {
    // 		h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);

    //             if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
    //           }
    //         else if(hasGenPromptPhoton)
    //           {
    // 		h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
    //             if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )){ //h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt); 
    // 		  continue;}
    // 		else
    //               {
    //                 if(madMinPhotonDeltaR >= 0.5)
    //                   h_selectBaselineYields_v1->Fill("mindR(q/g, #gamma)",wt);
    //                 if(mindr_Pho_genlep >=0.5)
    // 		      h_selectBaselineYields_v1->Fill("mindR(l, #gamma)",wt);
    //               }

    //           }
    //       }//Gen prompt    
                                                                                                                                                                               
    //     if(s_sample.Contains("WJets"))
    //       {
    //         if(!hasGenPromptPhoton)
    //           {
    // 		h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);

    // 		if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
    //           }
    //         else if(hasGenPromptPhoton)
    //           {
    // 		h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
    //             if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)){//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);   
    // 		  continue;}
    // 		else
    //               {
    //                 if(madMinPhotonDeltaR >= 0.5)
    //                   h_selectBaselineYields_v1->Fill("pass_mindR(q/g, #gamma)",wt);
    //                 if(mindr_Pho_genlep >=0.5)
    // 		      h_selectBaselineYields_v1->Fill("pass_mindR(l, #gamma)",wt);
    //               }

    //           }
    //       }
    //     if(s_sample.Contains("TTJets_HT") || s_sample.Contains("TTJets-HT")||s_sample.Contains("TTJets-inc")|| s_sample.Contains("TTJets_inc") || s_sample.Contains("TTJets2_v17")||s_sample.Contains("TTJets"))
    //       {
    //         if(!hasGenPromptPhoton)
    //           {
    // 		h_selectBaselineYields_v1->Fill("No gen prompt #gamma",wt);
    //             if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;
    //           }
    //         else if(hasGenPromptPhoton)
    //           {
    // 		h_selectBaselineYields_v1->Fill("Gen prompt #gamma",wt);
    //             if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)){// h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt); 
    // 		  continue;}
    // 		else
    //               {
    //                 if(madMinPhotonDeltaR >= 0.5)
    //                   h_selectBaselineYields_v1->Fill("pass_mindR(q/g, #gamma)",wt);
    //                 if(mindr_Pho_genlep >=0.5)
    //                   h_selectBaselineYields_v1->Fill("pass_mindR(l, #gamma)",wt);
    //               }


    //           }
    //       }
    //     if(hasGenPromptPhoton && (s_sample.Contains("GJets")))
    //       {
    //         if(!(madMinPhotonDeltaR>0.4)) continue;

    //       }
    //     if(hasGenPromptPhoton && (s_sample.Contains("QCD")))
    //       {
    //         if((madMinPhotonDeltaR>0.4 && hasGenPromptPhoton)) continue;
    //       }
    // 	if(hasGenPromptPhoton && ((s_sample.Contains("ZG"))|| (s_sample.Contains("ZNuNuG"))))
    //       {
    //         if(!(madMinPhotonDeltaR>0.5)) continue;
    //       }
    //     if(hasGenPromptPhoton && ((s_sample.Contains("ZJets"))|| (s_sample.Contains("ZNuNuJets"))))
    //       {
    //         if(!(madMinPhotonDeltaR<=0.5)) continue;
    //       }
    //     genphomatch_after++;

    //   }
    // h_selectBaselineYields_->Fill("after bkg cmp",wt);
    // nEvents_Selec[7]+=wt;
    // out_nEventsTags[7]="overlap bkg";
    // // /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // // //////////////////////////////   Evaluating the response from BDT  /////////////////////////////////////////
    
    // met = MET;
    // njets = nHadJets;
    // btags = BTags;
    // mTPhoMET = mTPhoMET;
    // dPhi_PhoMET = dPhi_phoMET;
    // dPhi_MetJet =dPhi_METjet1;
    // st = ST;
    // // Double_t mvaValue = reader1->EvaluateMVA( "BDT_100trees_2maxdepth method");
    // // h_mvaResponse->Fill(mvaValue,wt);
    // if(Debug)
    //   cout<<"===load tree entry before SR ==="<<"\t"<<jentry<<endl;
    
    // h_Njets[1]->Fill(nHadJets,wt);
    // h_Nbjets[1]->Fill(BTags,wt);
    // h_MET_[1]->Fill(MET,wt);
    // h_PhotonPt[1]->Fill(bestPhoton.Pt(),wt);
    // h_Mt_PhoMET[1]->Fill(mTPhoMET,wt);
    // h_dPhi_PhoMET[1]->Fill(dPhi_PhoMET,wt);
    // h_St[1]->Fill(ST,wt);
    // h_HT[1]->Fill(Ht,wt);
    // h_njets_vs_ST[1]->Fill(nHadJets,ST,wt);
    // h_njets_vs_HT[1]->Fill(nHadJets,Ht,wt);
    // h_ST_vs_ptPho[1]->Fill(ST,bestPhoton.Pt(),wt);
    // // h_mindr_Pho_genlep[1]->Fill(mindr_Pho_genlep);
    // // if (NElectrons == 0 && NMuons == 0 ) h_selectBaselineYields_->Fill("veto electron & Muon",wt);
    // // else continue;
    // // //{
    // // if(isoElectronTracks==0 && isoMuonTracks ==0 && isoPionTracks==0)
    // // 	{
    // // 	  int Sbin_prev=getBinNoV6_WithOnlyBLSelec(BTags,nHadJets); 
    // // 	  h_Sbins_LL[30]->Fill(Sbin_prev,wt);
    // // 	}
    // // else continue;
    // //}
    // nsurVived+=wt;
     } //loop over entries

  if(Debug)
    cout<<"filling the branches in tree"<<endl;
   //  cout<<nocut<<"\t"<<nSurvived<<"\t"<<bkg_comp<<endl;
  cout<<"Alpana-check"<<"\t"<<"events not falling in any LL CR/SR"<<"\t"<<counter<<endl;

   cout<<"Survived preselection ===: "<<"\t"<<nsurVived<<endl;
  cout << "============   Electron     ===============" <<endl;

  cout<<"applied weights "<<" "<<wt<<endl;
  cout<<"CR electron :==  "<<"\t"<<nCR_elec<<"\t"<<endl; 
  cout<<"SR electron :==  "<<"\t"<<nSR_elec<<"\t"<<endl;
  cout<<"SR e- : Fail acceptance "<<"\t"<<FailAccept_Elec<<"\t"<<endl;
  cout<<"SR e- : Fail Id " <<"\t"<<FailId_Elec<<"\t"<<endl;
  cout<< "SR e- : Fail Iso" <<"\t"<<FailIso_Elec<<"\t"<<endl;

  cout << "============   Muon     ===============" <<endl;

  // cout<<"applied weights "<<" "<<wt<<endl;
  // cout<<"CR muon :==  "<<"\t"<<nCR_mu<<"\t"<<endl;
  // cout<<"SR muon :==  "<<"\t"<<nSR_mu<<"\t"<<endl;
  // cout<<"SR mu : Fail acceptance "<<"\t"<<FailAccept_Mu<<"\t"<<endl;
  // cout<<"SR mu : Fail Id " <<"\t"<<FailId_Mu<<"\t"<<endl;
  // cout<< "SR mu : Fail Iso" <<"\t"<<FailIso_Mu<<"\t"<<wt*FailIso_Mu<<endl;

  // for (int i =0;i<49;i++)
  //   {
  //     cout<<out_nEventsTags[i]<<" :"<<"\t"<<"=====: "<<"\t"<<nEvents_Selec[i]<<endl;
      
  //     cout<<" "<<"\t"<<endl;
      
  //   }

  // for (int i =0;i<49;i++)
  //   {
  //     cout<<nEvents_Selec[i]<<endl;
  //   }  
  // cout<<"Alpana-check"<<"\t"<<"events not falling in any LL CR/SR"<<"\t"<<counter<<endl;
  // outfile->cd();
  // skim_tree->Write();
  // outfile->Close();
  // cout<<"outFile: "<<outfileName<<" written!!"<<endl;

}


int AnalyzeLightBSM::getBinNoV6_WithOnlyBLSelec(int nbjets, int nHadJets)
{
  
  int sBin=-100,m_i=0;
  if(nbjets==0 ){
    if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
    else if(nHadJets==5 || nHadJets==6){ sBin=7;}
    else if(nHadJets>=7)               { sBin=13;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)     { sBin=19;}
    else if(nHadJets==5 || nHadJets==6){ sBin=25;}
    else if(nHadJets>=7)               { sBin=31;}
  }
  if(sBin==0){
    for(int i=0;i<METLowEdge_v3.size()-1;i++){
      if(METLowEdge_v3[i]<99.99) continue;
      int sBin1=sBin;
      m_i++;
      if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
        break; }
      else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 7;
        break; }
    }
  }
  else if(sBin==7 || sBin==13 || sBin==19 || sBin==25 || sBin==31){
    int sBin1=sBin;
    for(int i=0;i<METLowEdge_v3_1.size()-1;i++){
      if(METLowEdge_v3_1[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge_v3_1[i] && MET < METLowEdge_v3_1[i+1]){ sBin = sBin+m_i;
        break;}
      else if(MET >= METLowEdge_v3_1[METLowEdge_v3_1.size()-1])  { sBin = sBin+6;
        break; }
    }
  }

  // else if(sBin==37){
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }
  //     // else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }                                                                                                                    

  //   }
  // }

  // else if(sBin==44){
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 52   ;break; }
  //     // else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 51   ;break; }                                                                                                                    

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
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3[i]<99.99) continue;
  //     int sBin1=sBin;
  //     m_i++;
  //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
  // 	break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 7;
  // 	break; }
  //   }
  // }
  // else if(sBin==8 || sBin==15 || sBin==22 || sBin==29 || sBin==36){
  //   int sBin1=sBin;
  //   for(int i=0;i<METLowEdge_v3_1.size()-1;i++){
  //     if(METLowEdge_v3_1[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_v3_1[i] && MET < METLowEdge_v3_1[i+1]){ sBin = sBin+m_i;
  // 	break;}
  //     else if(MET >= METLowEdge_v3_1[METLowEdge_v3_1.size()-1])  { sBin = sBin+6;
  // 	break; }
  //   }
  // }

  // else if(sBin==37){
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 44   ;break; }

  // }
  // }

  // else if(sBin==44){
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3[i]<99.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 52   ;break; }
  //     // else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 51   ;break; }

  //   }
  //}
  return sBin;
}
int AnalyzeLightBSM::getBinNoV7_le(int nbjets, int nHadJets){
  int sBin=-100,m_i=0;
  if(nbjets==0){
    if(nHadJets==2) { if(MET<150)sBin=1; else if (MET>150) sBin=2;}
    else if(nHadJets==3)     { if(MET<150)sBin=3; else if (MET>150) sBin=4;}
    else if(nHadJets==4)     { if(MET<150)sBin=5; else if (MET>150) sBin=6;}
    else if((nHadJets==5 || nHadJets==6)){ if(MET<150)sBin=7; else if (MET>150) sBin=8;}
    else if(nHadJets>=7)   { if(MET<150)sBin=9; else if (MET>150) sBin=10;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)      {if(MET<150)sBin=11; else if (MET>150) sBin=12;}
    else if((nHadJets==5 || nHadJets==6)){ if(MET<150)sBin=13; else if (MET>150) sBin=14;}
    else if(nHadJets>=7)   { if(MET<150)sBin=15; else if (MET>150) sBin=16;}
  }
  return sBin;
}

int AnalyzeLightBSM::getBinNoV1_le(int nbjets, int nHadJets){
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
int AnalyzeLightBSM::getBinNoV7(int nHadJets, int btags){
  int sBin=-100,m_i=0;
  if(btags==0){
    if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
    else if(nHadJets==5 || nHadJets==6){ sBin=6;}
    else if(nHadJets>=7)               { sBin=11;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)     { sBin=16;}
    else if(nHadJets==5 || nHadJets==6){ sBin=21;}
    else if(nHadJets>=7)               { sBin=26;}
  }
  if(sBin==0){
    for(int i=0;i<METLowEdge1.size()-1;i++){
      if(METLowEdge1[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge1[i] && MET < METLowEdge1[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge1[METLowEdge1.size()-1])  { sBin = 6         ;break; }
    }
  }
  else{
    for(int i=0;i<METLowEdge2.size()-1;i++){
      if(METLowEdge2[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge2[i] && MET < METLowEdge2[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge2[METLowEdge2.size()-1])  { sBin = sBin+5   ;break; }
    }
  }
  return sBin;
}
TLorentzVector AnalyzeLightBSM::getBestPhoton(int pho_ID){
  // TLorentzVector v1;
  // vector<TLorentzVector> goodPho;
  // for(int iPhoton=0;iPhoton<Photons->size();iPhoton++){
  //   if( ((*Photons_fullID)[iPhoton]) && ((*Photons_hasPixelSeed)[iPhoton]<0.001) ) goodPho.push_back( (*Photons)[iPhoton] );
  // }

  // if(goodPho.size()==0) return v1;
  // sortTLorVec(&goodPho);
  // return goodPho[0];
  vector<TLorentzVector> goodPho;
  vector<int> goodPhoIndx;
  for(int iPho=0;iPho<Photons_;iPho++){
    if(((*Photons_hasPixelSeed)[iPho]<0.001) && ( (*Photons_fullID)[iPho]))// && ((*Photons_hasPixelSeed)[iPho]<0.001) &&( pho_ID==0 || (pho_ID==1 && (((*Photons_cutBasedIDFall17V2)[iPho]==1 || (*Photons_cutBasedIDFall17V2)[iPho]==2))) || (pho_ID==2 && (*Photons_cutBasedIDFall17V2)[iPho]==2) || (pho_ID==3 && (*Photons_mvaValuesIDFall17V2)[iPho]>-0.02) || (pho_ID==4 && (*Photons_mvaValuesIDFall17V2)[iPho]>0.42))) ) {
      {
      goodPho.push_back(Photons_v1[iPho] );
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
  if(highPtIndx==-100){TLorentzVector v0;return v0;}
  else return goodPho[highPtIndx];

  
}

// TLorentzVector AnalyzeLightBSM::getBestPhoton_tightID(){
//   vector<TLorentzVector> goodPho;
//   vector<int> goodPhoIndx;
//   for(int iPho=0;iPho<Photons->size();iPho++){
//     if((*Photons_hadTowOverEM)[iPho]<0.02148 && (*Photons_sigmaIetaIeta)[iPho]<0.00996 &&  (*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001)) {
//       goodPho.push_back( (*Photons)[iPho] );
//       goodPhoIndx.push_back(iPho);
//     }
//   }

//   int highPtIndx=-100;
//   for(int i=0;i<goodPho.size();i++){
//     if(i==0) highPtIndx=0;
//     else if( (goodPho[highPtIndx].Pt()) < (goodPho[i].Pt()) ){highPtIndx=i;}
//   }

//   if(highPtIndx>=0){
//     bestPhotonIndxAmongPhotons = goodPhoIndx[highPtIndx];
//   }
//   else bestPhotonIndxAmongPhotons = -100;
//   if(highPtIndx==-100){TLorentzVector v0;return v0;}
//   else return goodPho[highPtIndx];
// }



double AnalyzeLightBSM::getGendRLepPho(int pho_ID){//MC only
  TLorentzVector genPho1,genLep1;
  int leadGenPhoIdx=-100;
  // vector<TLorentzVector> goodPho;
  // for(int iPhoton=0;iPhoton<Photons->size();iPhoton++){
  //   if( ((*Photons_fullID)[iPhoton]) && ((*Photons_hasPixelSeed)[iPhoton]<0.001) ) goodPho.push_back((*Photons)[iPhoton]);
  // }
  // if(goodPho.size()!=0) 
  genPho1 =getBestPhoton(pho_ID);
  
  for(int i=0;i<GenParticles_;i++){
     if(GenParticles_v1[i].Pt()!=0){
       if( (abs((*GenParticles_PdgId)[i])==11 || abs((*GenParticles_PdgId)[i])==13 || abs((*GenParticles_PdgId)[i])==15 ) && (abs((*GenParticles_ParentId)[i])<=25) && (abs((*GenParticles_ParentId)[i])!=15) ){
	 if(genLep1.Pt() < (GenParticles_v1[i]).Pt()) genLep1 = (GenParticles_v1[i]);
       }
     }
  }//for
  if(genPho1.Pt() > 0. && genLep1.Pt() > 0.) return genLep1.DeltaR(genPho1);
  else return 1000.0;
}

 double AnalyzeLightBSM::getdR_GenPho_RecoPho(TLorentzVector bestPhoton)
 {
   
   TLorentzVector genPho;
   int leadGenPhoIdx=-100;
   int minDR = 9999;
   vector<TLorentzVector> v_genPho;
   for (int igen=0; igen<GenParticles_; igen++)
     {
       if(GenParticles_v1[igen].Pt()!=0){
	 if((abs((*GenParticles_PdgId)[igen])==22) && ((abs((*GenParticles_ParentId)[igen])<=25) || ((*GenParticles_ParentId)[igen]==2212) ) && (*GenParticles_Status)[igen]==1){
	   genPho = (GenParticles_v1[igen]);
	   v_genPho.push_back(genPho);
	 }
       }
     }
   return MinDr(bestPhoton,v_genPho);   
 }
double AnalyzeLightBSM::getGendRElecPho(int pho_ID){//MC only                                                                                             
  TLorentzVector genPho1,genLep1;
  int leadGenPhoIdx=-100;
  genPho1 =getBestPhoton(pho_ID);
  for(int i=0;i<GenParticles_;i++){
    if(GenParticles_v1[i].Pt()!=0){
      if( abs((*GenParticles_PdgId)[i])==11 && (abs((*GenParticles_ParentId)[i])<=25 )&& (abs((*GenParticles_ParentId)[i])!=15)){
	if(genLep1.Pt() < (GenParticles_v1[i]).Pt()) genLep1 = (GenParticles_v1[i]);
      }
    }
  }//for                                                                                                                                       
  if(genPho1.Pt() > 0. && genLep1.Pt() > 0.) return genLep1.DeltaR(genPho1);
  else return 1000.0;

}
double AnalyzeLightBSM::getGenLep1(TLorentzVector bestPhoton, bool flag){                                     
  vector<TLorentzVector> v_genLep2;
  TLorentzVector genMu1, genEle1;
  if(flag)
    {
      for(int i=0 ; i < GenElectrons_; i++)
	{
	  if(GenElectrons_v1[i].Pt()!=0)
	    {
	      genEle1 = (GenElectrons_v1[i]);
	      v_genLep2.push_back(genEle1);
	    }

	}
    }
  else
    {
      for(int i=0 ; i < GenMuons_; i++)
	{
	  if(GenMuons_v1[i].Pt()!=0)
	    {
	      genMu1 = (GenMuons_v1[i]);
	      v_genLep2.push_back(genMu1);
	    }
	}
    }
  return MinDr(bestPhoton,v_genLep2);
}


double AnalyzeLightBSM::getGenLep(TLorentzVector bestPhoton){
  vector<TLorentzVector> v_genLep2;
  TLorentzVector genMu1, genEle1;
  // if(flag)
  //   {
      for(int i=0 ; i < GenElectrons_; i++)
        {
          if(GenElectrons_v1[i].Pt()!=0)
            {
              genEle1 = (GenElectrons_v1[i]);
              v_genLep2.push_back(genEle1);
            }

        }
  //   }
  // else
  //   {
      for(int i=0 ; i < GenMuons_; i++)
        {
          if(GenMuons_v1[i].Pt()!=0)
            {
              genMu1 = (GenMuons_v1[i]);
              v_genLep2.push_back(genMu1);
            }
        }
      //  }
  return MinDr(bestPhoton,v_genLep2);
}

double AnalyzeLightBSM::getGenRecodRLep(TLorentzVector recoLep1){//MC only                                                                                                
  TLorentzVector genPho1,genLep1;
  int leadGenPhoIdx=-100;
  for(int i=0;i<GenParticles_;i++){
    if(GenParticles_v1[i].Pt()!=0){
      if( (abs((*GenParticles_PdgId)[i])==11 || abs((*GenParticles_PdgId)[i])==13 || abs((*GenParticles_PdgId)[i])==15 ) && (abs((*GenParticles_ParentId)[i])<=25) && (abs((*GenParticles_ParentId)[i])!=15) ){
	if(genLep1.Pt() < (GenParticles_v1[i]).Pt()) genLep1 = (GenParticles_v1[i]);
      }
    }
  }//for                                                                                                                                          
  if(recoLep1.Pt() > 0. && genLep1.Pt() > 0.) return genLep1.DeltaR(recoLep1);
  else return 1000.0;
}

std::vector<int> AnalyzeLightBSM::dR_recoPho_GenParticle(TLorentzVector recoPho){//MC only                                                                     
  TLorentzVector genPho1,genLep1;
  int leadGenPhoIdx=-100, pdgID =-99999, parentID=-99999, count=0;
  int flag = 0;
  std::vector<int> array;//={};
  for(int i=0;i<GenParticles_;i++){
    if(GenParticles_v1[i].Pt()!=0){
      //if(abs((*GenParticles_PdgId)[i])==22 || abs((*GenParticles_ParentId)[i])==22) continue;
      if(recoPho.DeltaR(GenParticles_v1[i])<0.2 )
	{
	  //cout<<"You got to be kidding me"<<endl;	  
	  //h_parID_matchRecoPho->Fill((*GenParticles_PdgId)[i],wt);
	  //h_parentID_vsPartId_matchPho->Fill((*GenParticles_PdgId)[i],(*GenParticles_ParentId)[i],wt);
	  if(abs((*GenParticles_PdgId)[i])==22)
	    {
	      //h_phopT_BeamRemenant->Fill(recoPho.Pt(),wt);
	      continue;
	    }
	  if(abs((*GenParticles_ParentId)[i])==22) continue;	  
	  pdgID =(*GenParticles_PdgId)[i];
	  count++;
	  flag = 1;	  
	  //cout<<"You got to be kidding me"<<"\t"<<flag<<"\t"<<pdgID<<"\t"<<(*GenParticles_ParentId)[i]<<"\t"<<(*GenParticles_Status)[i]<<endl;
	}
    }
  }
  array.push_back(pdgID);
  array.push_back(flag);
  //cout<<"after kidding me"<<"\t"<<array[0]<<"\t"<<array[1]<<endl;
  return array;//pdgID, flag;
 
}
// double AnalyzeLightBSM::getMVArespone(float MET, float NJets, float BTags, float HT)
// {
//   TMVA::Tools::Instance();

//   TMVA::Reader *reader1 = new TMVA::Reader("Color:Silent");
//   float met=0.0,njets=0.0,btags=0.0,ht=0.0;
//   reader1->AddVariable( "MET", &met );
//   reader1->AddVariable( "NJets", &njets );
//   reader1->AddVariable( "BTags", &btags );
//   reader1->AddVariable( "HT", &ht );
//   reader1->BookMVA( "BDT_200trees_2maxdepth method", "./TMVAClassification_BDT_200trees_2maxdepth.weights.xml" );
//   met = MET;
//   njets= NJets;
//   btags = BTags;
//   ht = HT;
//   Double_t mvaValue = reader1->EvaluateMVA( "BDT_200trees_2maxdepth method");
//   return mvaValue;
//   delete reader1;
// }
