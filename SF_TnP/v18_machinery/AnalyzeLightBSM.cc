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

   if (argc <5 ) {
     cerr << "Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "which year dataset" <<" "<<"which Process"<< " "<<"Which pho_ID"<<endl;
     return -1;
   }
   const char *inputFileList = argv[1];
   const char *outFileName   = argv[2];
   const char *data          = argv[3];
   const char *sample=argv[4];
   const char *phoID = argv[5];
   //TString pho_ID = phoID;





   AnalyzeLightBSM ana(inputFileList, outFileName, data,sample, phoID);
   cout << "dataset " << data << " " << endl;
   cout<<"Which pho_ID: "<<"\t"<<phoID<<endl;
   ana.EventLoop(data,inputFileList,sample,outFileName,phoID);
   Tools::Instance();
   return 0;
 }

 void AnalyzeLightBSM::EventLoop(const char *data,const char *inputFileList, const char *sample , const char *outFileName,  const char* phoID) {
   cout<<"inside event loop"<<endl;
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   cout << "nentries " << nentries << endl;
   cout << "Analyzing dataset " << data << " " << endl;
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
   cout<<"which photon ID "<<"\t"<<"case"<<"\t"<<pho_ID<<endl;
   // TString pho_D = phoID;
   bool FR_elec_flag=false;
   float counter=0.0;
   TString s_data=data;
  char* s_cross = new char[100];
   sprintf(s_cross,"%s.%s",data,sample);
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
   bool applyTrgEff=true;
   bool applyHEMveto=true;
   bool applyL1TrigFire_cal=true;
   bool applyL1TrigFire_prob=true;
   bool applyPUwt=false;
   bool apply_pixelveto=false;
   bool check_flag_wBL =false;
   if(s_sample.Contains("UL")){
   if(s_data.Contains("2016preVFP")){ lumiInfb=19.5;deepCSVvalue = 0.6001; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01;}// APV
   if(s_data.Contains("2016postVFP")) { lumiInfb=16.5; deepCSVvalue = 0.5847; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01;} //2016

   if(s_data.Contains("2017")) {lumiInfb=41.48; deepCSVvalue = 0.4506;}
   if(s_data.Contains("2018")){ lumiInfb=59.83;deepCSVvalue = 0.4168;}
   if(s_data.Contains("signal"))lumiInfb= 137.19;

   if(!s_sample.Contains("data"))
     {
       cout<<" analyzing MC bkgs" <<endl;
       if(s_data.Contains("2016preVFP") && highdphi) {lumiInfb=19.5; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01; deepCSVvalue = 0.6001;}
       if(s_data.Contains("2016postVFP") && highdphi) {lumiInfb=16.5; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01; deepCSVvalue = 0.5847;}
      
       if(s_data.Contains("2017") && highdphi) {lumiInfb=41.48; p0=1.82e+02; p1=6.336e+01; p2=9.171e-01; deepCSVvalue = 0.4506;}
       if(s_data.Contains("2018") && highdphi) {lumiInfb=59.83; p0=1.787e+02; p1=6.657e+01; p2=9.47e-01; deepCSVvalue = 0.4168;}

       if(s_data.Contains("2016preVFP") && !highdphi) {lumiInfb=19.5; p0=1.586e+02; p1=6.83e+01; p2=9.28e-01; deepCSVvalue = 0.6001;}
       if(s_data.Contains("2016postVFP") && highdphi) {lumiInfb=16.5;p0=1.586e+02; p1=6.83e+01; p2=9.28e-01; deepCSVvalue = 0.5847;}

       if(s_data.Contains("2017") && !highdphi) {lumiInfb=41.48; p0=1.82e+02; p1=6.336e+01; p2=9.171e-01; deepCSVvalue = 0.4506;}
       if(s_data.Contains("2018") && !highdphi) {lumiInfb=59.83; p0=1.787e+02; p1=6.657e+01; p2=9.47e-01; deepCSVvalue = 0.4168;}

       if(s_data.Contains("FastSim")) lumiInfb=137.19;
       cout<<"Trigger efficiency flag "<<applyTrgEff<<" p0 "<<p0<<"  p1  "<<p1<<"  p2 "<< p2<<endl;

     }
   }
   else
     {

       if(!s_sample.Contains("data"))
	 {
	   if(s_data.Contains("2016")) {lumiInfb=35.922;deepCSVvalue = 0.6321;}
	   //        if(s_data.Contains("2017")) {lumiInfb=41.529;deepCSVvalue = 0.4941;}                                                                                
	   if(s_data.Contains("2017")) {lumiInfb=41.9;deepCSVvalue = 0.4941;}
	   //              if(s_data.Contains("2018")) {lumiInfb=59.74;deepCSVvalue = 0.4184;}                                                                           
	   if(s_data.Contains("2018")) {lumiInfb=58.703;deepCSVvalue = 0.4184;}
	   //if(s_data.Contains("FastSim") && s_data.Contains("2016")) lumiInfb=137.19;
	 }
       if(s_sample.Contains("data"))
	 {
	   if(s_data.Contains("2016")) {deepCSVvalue = 0.6321;}
	   if(s_data.Contains("2017")) {deepCSVvalue = 0.4941;}
	   if(s_data.Contains("2018")) {deepCSVvalue = 0.4184;}
	 }

     }
   // if(!s_data.Contains("data"))
   //   {

   //     if(s_data.Contains("2016") && highdphi) {lumiInfb=35.922; p0=1.66539e+02; p1=8.13254e+01; p2=9.71152e-01; };//deepCSVvalue = 0.6321;}
   //     if(s_data.Contains("2017") && highdphi) {lumiInfb=41.529; p0=1.86744e+02; p1=6.74978e+01; p2=9.65333e-01;};// deepCSVvalue = 0.4941;}
   //     if(s_data.Contains("2018") && highdphi) {lumiInfb=59.74; p0=1.89868e+02; p1=6.60434e+01; p2=9.79618e-01;};// deepCSVvalue = 0.4184;}

   //     if(s_data.Contains("2016") && !highdphi) {lumiInfb=35.922; p0=1.67229e+02; p1=8.52729e+01; p2=8.29784e-01; deepCSVvalue = 0.6321;}
   //     if(s_data.Contains("2017") && !highdphi) {lumiInfb=41.529; p0=1.67641e+02; p1=1.21487e+02; p2=9.23864e-01; deepCSVvalue = 0.4941;}
   //     if(s_data.Contains("2018") && !highdphi) {lumiInfb=59.74; p0=1.45529e+02; p1=1.08431e+02; p2=9.27220e-01; deepCSVvalue = 0.4184;}

   //     if(s_data.Contains("FastSim") && s_data.Contains("2016")) lumiInfb=137.19;                                                                                    
   //    }
   if(s_sample.Contains("data") && s_sample.Contains("UL"))
     {
       cout<<" analyzing data" <<endl;

       if(s_data.Contains("2016preVFP")) {deepCSVvalue = 0.6001;}
       if(s_data.Contains("2016postVFP")) {deepCSVvalue = 0.5847;}

       if(s_data.Contains("2017")) {deepCSVvalue = 0.4506;}
       if(s_data.Contains("2018")) {deepCSVvalue = 0.4168;}
     }

   cout<<"Calculating the weight for process "<<s_cross<<" deepCSVvalue  "<<deepCSVvalue<<"\t"<<"sdata  "<<s_data<<"\t"<<"s_sample "<<s_sample<<endl;
  cout<<"Trigger efficiency flag "<<applyTrgEff<<" p0 "<<p0<<"  p1  "<<p1<<"  p2 "<< p2<<endl;
   std::string s_process = s_cross;
   TString s_Process = s_process;
   double cross_section = getCrossSection(s_process);
   cout<<cross_section<<"\t"<<"analyzed process"<<"\t"<<s_cross<<"\t"<<s_process<<endl;
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
   // const char* filetag[42]={"TTGJets_2018","TTGJets_2017","TTGJets_2016postVFP","Run2_TTGJets","WGJets_2018","WGJets_2017","WGJets_2016postVFP","Run2_WGJets","TTJets_2018","TTJets_2017","TTJets_2016postVFP","Run2_TTJets","WJets_2018","WJets_2017","WJets_2016postVFP","Run2_WJets","WGJets+WJets_2018","WGJets+WJets_2017","WGJets+WJets_2016postVFP","Run2_WGJets+WJets","TTGJets+TTJets_2018","TTGJets+TTJets_2017","TTGJets+TTJets_2016postVFP","Run2_TTGJets+TTJets","W+TTBar_2018","W+TTBar_2017","W+TTBar_2016postVFP","W+TTBar_FullRun2","TTGJets_2016preVFP","WGJets_2016preVFP","TTJets_2016preVFP","WJets_2016preVFP","WGJets+WJets_2016preVFP","TTGJets+TTJets_2016preVFP","W+TTBar_2016preVFP","TTGJets_2016","WGJets_2016","TTJets_2016","WJets_2016","WGJets+WJets_2016","TTGJets+TTJets_2016","W+TTBar_2016"};

   //const char* filetag[8]={"TTGJets_2018","TTGJets_2017","TTGJets_2016","Run2_TTGJets","WGJets_2018","WGJets_2017","WGJets_2016","Run2_WGJets"};
  //  char* hname = new char [200];
  //  sprintf(hname, "Electron_FR_TFbins_v3_phopt_qmulti_phoID_loose_09Jan24.root");
  //  // char* hname1 = new char [200];
  //  // sprintf(hname1, "Lepton_LL_TFbins_v2_nJetsBjets_PhoPt_phoID_loose_09Jan24.root");//"TF_Electron_LLEstimation_binsV0_phoID_%s_08Aug23.root",phoID);
  //   char* hname2 = new char [200];
  //  sprintf(hname2, "Electron_FR_TFbins_v3_phopt_qmulti_phoID_loose_09Jan24.root");// "TF_allin1_LLEstimation_binsV0_phoID_%s_08Aug23.root",phoID);
  //  TFile* f_TFbins= new TFile(hname2);
  //  // TFile* f_TFbins_v1 =new TFile(hname1);
  //  // TFile* f_TFbins_v2= new TFile(hname2);
  //   char* histname = new char[2000];
  //  TH1F* h_TF;
  //  // TH1F* h_TF1;
  //  // TH1F* h_TF2;
  // // auto sample1="";
  // // if(s_sample.Contains("Autumn18.WGJets_MonoPhoton_PtG-130")|| s_sample.Contains("Autumn18.WGJets_MonoPhoton_PtG-40to130"))
  // //   sample1 = "WGJets";
  // // else
  // //   sample1 = sample;
  
  //  sprintf(histname,"h_TFbins_LL_W+TTBar_%s",data);//FullRun2");//h_TFbins_LL_W+TTBar_FullRun2");//2018",data);                                                                  
  // h_TF = (TH1F*)f_TFbins->Get(histname);
  // // h_TF1 = (TH1F*)f_TFbins_v1->Get(histname);
  // // h_TF2 = (TH1F*)f_TFbins_v2->Get(histname);
  // cout<<"Reading TF for lost electron estimation:  "<<"\t"<<histname<<endl;
  // for(int i=0; i<h_TF->GetNbinsX();i++)
  //   {cout<<"TFBIns_v0"<<"\t"<<i<<"\t"<<h_TF->GetBinContent(i)<<endl;}

  // for(int i=0; i<h_TF1->GetNbinsX();i++)
  //   {cout<<"TFBIns_v1"<<"\t"<<i<<"\t"<<h_TF1->GetBinContent(i)<<endl;}
  // for(int i=0; i<h_TF2->GetNbinsX();i++)
  //   {cout<<"TFBIns_v2"<<"\t"<<i<<"\t"<<h_TF2->GetBinContent(i)<<endl;}
   int coutt=0;
 double  count_nelec=0.0;
   //nentries=1000;
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
 	cout<<"===load tree entry ==="<<jentry<<endl;
       if (ientry < 0) break;
       nb = fChain->GetEntry(jentry);   nbytes += nb;
       if(Debug)
         cout<<"===load tree entry ==="<<jentry<<endl;
       wt =1.0;
       h_selectBaselineYields_v1->Fill("no cuts",wt);
       if(s_Process.Contains("2018.WGJets_MonoPhoton_PtG-40to130UL")|| s_Process.Contains("2018.WGJets_MonoPhoton_PtG-130UL")|| s_Process.Contains("2016preVFP.WGJets_MonoPhoton_PtG-40to130UL") ||s_Process.Contains("2016preVFP.WGJets_MonoPhoton_PtG-130UL") || s_Process.Contains("2017.WGJets_MonoPhoton_PtG-40to130UL")|| s_Process.Contains("2017.WGJets_MonoPhoton_PtG-130UL") ||s_Process.Contains("2016postVFP.WGJets_MonoPhoton_PtG-130UL")||s_Process.Contains("2016postVFP.WGJets_MonoPhoton_PtG-40to130UL"))//|| s_Process.Contains("2018.TTGJets_incUL")) //(s_data.Contains("2016") || s_data.Contains("2017") || s_data.Contains("2018") || s_sample.Contains("WGJets"))
       	{	  
 	  if(jentry==0)
       	    cout<<cross_section<<"\t"<<"analyzed process"<<"\t"<<s_process<<endl;
       	  wt = cross_section*lumiInfb*1000.0;//)/nentries; // Weight*lumiInfb*1000.0; //(cross_section*lumiInfb*1000.0)/nentries;
 	  if(Debug &&jentry<10) cout<<"Alpana urgent check "<<wt<<"\t"<<cross_section<<"\t"<<lumiInfb<<endl;
       	}
       else// if (s_data.Contains("2016") || s_data.Contains("2017") || s_data.Contains("2018"))// || s_sample.Contains("WGJets"))
       	{
       	  wt = Weight*lumiInfb*1000.0;
       	}
     if (jentry<10) cout<<"Printing weight befor apply trig eff "<<wt<<"\t"<<cross_section<<"\t"<<lumiInfb<<endl;
      if(s_sample.Contains("data")) wt =1;
      if(jentry<10)
	{
	  cout<<"Event "<<jentry<<" before event trig  "<<wt<<"\t"<<Weight<<" TriggerPass->size() "<<TriggerPass->size()<<endl;
	}
      bool tighte_trgpass=true;
      if(s_sample.Contains("data") && !s_sample.Contains("UL")&& applyTrgEff)
	{	if(TriggerPass->size()!=148) continue;
	  if((*TriggerPass)[36]==1 && s_data.Contains("2016"))
	    tighte_trgpass=true;
	  else if((*TriggerPass)[40]==1 && !s_data.Contains("2016"))//((*TriggerPass)[36]==1 || (*TriggerPass)[40]==1) && !s_data.Contains("2016"))
	    {
	      if((*TriggerPrescales)[40]!=1)     cout<<jentry<<" :  TriggerPrescales = "<<(*TriggerPrescales)[40]<<endl;
	      tighte_trgpass=true;
	    }
	  else tighte_trgpass=false;
	  if(jentry<10) cout<<"Alp  "<<(*TriggerVersion)[46] <<"\t"<<(*TriggerVersion)[51]<<"\t"<<endl;//(*TriggerVersion)[183]<<"\t"<<(*TriggerVersion)[164]<<"\t"<<(*TriggerVersion)[166]<<"\t"<<(*TriggerVersion)[170] <<"\t"<<(*TriggerVersion)[172] <<"\t"<<(*TriggerVersion)[174] <<"\t"<<endl; 
	  //if(tighte_trgpass==false)  continue;
	}
      


      if(s_data.Contains("2016") || s_data.Contains("2017")){ //L1 prefire issue
	wt= wt*NonPrefiringProb;
      }
     
      if(!s_sample.Contains("data") && applyPUwt)
	wt = wt*puWeight;

      
      // if(!s_sample.Contains("data") && !s_sample.Contains("signal") && applyTrgEff )
      // 	{
      // 	  wt = wt * (((TMath::Erf((MET - p0)/p1)+1)/2.0)*p2);
      // 	}
      if(Debug)
	cout<<"alpa"<<"\t"<<jentry<<"\t"<<endl;
      if(jentry<10)
	{    
	  cout<<"Event "<<jentry<<" event weight"<<"\t"<<wt<<"\t"<<cross_section<<"\t"<<Weight<<"lumiInfb  "<<lumiInfb<<endl;
	  cout<<"Event "<<jentry<<"\t"<<wt<<"\t"<<"MET  "<<MET<<"\t"<<(((TMath::Erf((MET - p0)/p1)+1)/2.0)*p2)<<endl;
	}
           //      h_selectBaselineYields_->Fill("No cuts, evt in 1/fb",wt);
      nocut++;
      if(Debug)
    	cout<<"Before filling lorentz vectors"<<endl;
      //filling the lorentz vector variables 

      // if(!s_sample.Contains("UL")){      
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
        HLTElectronObjects_v1.clear();
        TAPElectronTracks_v1.clear();
      	Electrons_v1 = (*Electrons);
      	Photons_v1 =  (*Photons);
	Muons_v1 = (*Muons);
	Electrons_ =Electrons_v1.size();
	Photons_ =Photons_v1.size();
	Muons_= Muons_v1.size();
	//Taus_v1 = (*Taus);
	Jets_v1 = (* Jets);
	Jets_ = Jets_v1.size();
	if(s_sample.Contains("data"))  {
	  HLTElectronObjects_v1 = (*HLTElectronObjects);
	  HLTElectronObjects_ = HLTElectronObjects_v1.size();
	}
	TAPElectronTracks_v1 = (*TAPElectronTracks);
	TAPElectronTracks_ =TAPElectronTracks_v1.size();
	if(!s_sample.Contains("data"))  {
	  GenParticles_v1 = (*GenParticles);
	  GenParticles_ =GenParticles_v1.size();
	  GenMuons_v1 = (*GenMuons);
	  GenMuons_ = GenMuons_v1.size();
	  GenTaus_v1= (*GenTaus);
	  GenTaus_ = GenTaus_v1.size();
	  GenJets_v1 = (*GenJets);
	  GenJets_ = GenJets_v1.size();
	  GenElectrons_v1 = (*GenElectrons);
	  GenElectrons_ = GenElectrons_v1.size();
	}
	//}
      // if(s_sample.Contains("UL")){
      // 	//set objetc pointer
      // 	Electrons_v1.clear();
      // 	GenElectrons_v1.clear();
      // 	GenParticles_v1.clear();
      // 	GenMuons_v1.clear();
      // 	GenTaus_v1.clear();
      // 	GenJets_v1.clear();//_v1.clear();
      // 	Electrons_v1.clear();
      // 	Photons_v1.clear();
      // 	Muons_v1.clear();
      // 	Taus_v1.clear();
      // 	Jets_v1.clear();
      // 	HLTElectronObjects_v1.clear();
      // 	TAPElectronTracks_v1.clear();
      // 	if(Debug)
      // 	  cout<<jentry<<"\t"<<"After setting pointer  lorentz vectors"<<endl;
      // 	Electrons_v1 = getLorentzVector(Electrons_,Electrons_fCoordinates_fPt,Electrons_fCoordinates_fEta,Electrons_fCoordinates_fPhi,Electrons_fCoordinates_fE);
      // 	Muons_v1 = getLorentzVector(Muons_,Muons_fCoordinates_fPt,Muons_fCoordinates_fEta,Muons_fCoordinates_fPhi,Muons_fCoordinates_fE);
      // 	Photons_v1 = getLorentzVector(Photons_,Photons_fCoordinates_fPt,Photons_fCoordinates_fEta,Photons_fCoordinates_fPhi,Photons_fCoordinates_fE);
      // 	Jets_v1 = getLorentzVector(Jets_,Jets_fCoordinates_fPt,Jets_fCoordinates_fEta,Jets_fCoordinates_fPhi,Jets_fCoordinates_fE);
      // 	if(s_sample.Contains("data"))
      // 	   HLTElectronObjects_v1 = getLorentzVector(HLTElectronObjects_,HLTElectronObjects_fCoordinates_fPt,HLTElectronObjects_fCoordinates_fEta,HLTElectronObjects_fCoordinates_fPhi,HLTElectronObjects_fCoordinates_fE);
      // 	TAPElectronTracks_v1 = getLorentzVector(TAPElectronTracks_,TAPElectronTracks_fCoordinates_fPt,TAPElectronTracks_fCoordinates_fEta,TAPElectronTracks_fCoordinates_fPhi,TAPElectronTracks_fCoordinates_fE);
	
      // 	if(!s_sample.Contains("data")) {
      // 	  GenElectrons_v1 = getLorentzVector(GenElectrons_,GenElectrons_fCoordinates_fPt,GenElectrons_fCoordinates_fEta,GenElectrons_fCoordinates_fPhi,GenElectrons_fCoordinates_fE);
      // 	  GenMuons_v1 = getLorentzVector(GenMuons_,GenMuons_fCoordinates_fPt,GenMuons_fCoordinates_fEta,GenMuons_fCoordinates_fPhi,GenMuons_fCoordinates_fE);
      // 	  GenTaus_v1 = getLorentzVector(GenTaus_,GenTaus_fCoordinates_fPt,GenTaus_fCoordinates_fEta,GenTaus_fCoordinates_fPhi,GenTaus_fCoordinates_fE);
      // 	  GenJets_v1 = getLorentzVector(GenJets_,GenJets_fCoordinates_fPt,GenJets_fCoordinates_fEta,GenJets_fCoordinates_fPhi,GenJets_fCoordinates_fE);
      // 	  GenParticles_v1 = getLorentzVector(GenParticles_,GenParticles_fCoordinates_fPt,GenParticles_fCoordinates_fEta,GenParticles_fCoordinates_fPhi,GenParticles_fCoordinates_fE);
      // 	}

      // if(Debug)
      // 	cout<<jentry<<"\t"<<"After filling  reco Jets"<<endl;

      // }
      
      int count=0;
      int branch_size = GenParticles_v1.size();
      int ele_branch=Electrons_v1.size();
      if(Debug)
        cout<<"Electrons "<<"\t"<<Electrons_<<"\t"<<Muons_<<"\t"<<Jets_<<endl;

      int pho_branch = Photons_v1.size();
      if(Debug)
    	        cout<<"Photons "<<"\t"<<Photons_<<endl;
      //      h_selectBaselineYields_v1->Fill("no cuts",wt);
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
    
     //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     ///////////////////////////////////   Making gen level collections of all leptons  ////////////////////////////////////////////////
     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if(Debug)
      cout<<"===load tree entry check1 - after filling gen level at entry ==="<<"\t"<<jentry<<endl;           
           
      int nGenMu=0,nGenEle=0,nGenTau=0,nGenHad=0,nGenLep=0,nGenEle_tau=0,nGenMu_tau=0,nGentau_lep=0,nGentau_had=0,nGentau_had1=0,nGenTau1=0,nGenEle1=0,nGenEle_tau1=0,nGenMu1=0,nGenMu_tau1=0,o=0,nelec_reco=0,nmu_reco=0;
      TLorentzVector genPho1,genEle1,genMu1,genTau1,recElec, recMu;//,genHad1,genLep1,gentau_lep1,gentau_had1,recElec, recMu;
      vector<TLorentzVector> v_genEle1, v_genPho1, v_genMu1,v_genTau1,v_genHad1,v_gentau_lep1,v_gentau_had1,v_gentau_had2,v_genTau2,v_genEle2,v_genMu2,v_genLep2,v_genLep1,v_recEle,v_recMu;
      v_genEle2.clear();
      v_genEle1.clear();
      v_genMu1.clear();
      v_genMu2.clear();
      v_genTau1.clear();
      v_genTau2.clear();
      bool HEMaffected=false, extracondition=false;
      if(s_data.Contains("2018") && applyHEMveto){
	for(int i=0; i<Electrons->size();i++)
	  if((*Electrons)[i].Pt() >30 && (*Electrons)[i].Eta() > -3.0 && (*Electrons)[i].Eta() < -1.4 && (*Electrons)[i].Phi() > -1.57 && (*Electrons)[i].Phi() < -0.87)   {HEMaffected = true; break;}
	for(int i=0; i<Jets->size();i++)
	  if((*Jets)[i].Pt() > 30 && (*Jets)[i].Eta() > -3.2 && (*Jets)[i].Eta() < -1.2 && (*Jets)[i].Phi() > -1.77 && (*Jets)[i].Phi() < -0.67 && DeltaPhi((*Jets)[i].Pt(),METPhi)<0.5) {HEMaffected = true; break;}
	if(HEMaffected == true) continue;
      }
      h_selectBaselineYields_v1->Fill("HET veto ",wt);
      if(!s_sample.Contains("FastSim")){
	if(s_sample.Contains("data"))
	  {
	    if(s_data.Contains("2016"))
	      ecalBadCalibReducedExtraFilter=1;
	    if(!(NVtx > 0 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadPFMuonFilter && eeBadScFilter==1 && ecalBadCalibReducedExtraFilter==1)) continue;
	  }
      
	if(!s_sample.Contains("data"))
	  {
	    eeBadScFilter = 1;
	    if(!(NVtx > 0 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadPFMuonFilter && eeBadScFilter==1 && ecalBadCalibReducedExtraFilter==1)) continue;
	    }

      }

      if(PFCaloMETRatio >=  5) continue;
      if(MET/CaloMET > 2.0) continue;
      
      h_selectBaselineYields_v1->Fill("MET filters ",wt);
      int hasEle=0,hasEle2=0, hasPho=0, hasPho_px=0,hastagEle=0,hasprobeEle=0;
      TLorentzVector bestPhoton=getBestPhoton(pho_ID);
      TLorentzVector bestPhoton_woFullID = getPhoton_withoutfullID();
      bool bestPhoHasPxlSeed=true, noPhoHasPxlSeed=true;
      if(bestPhotonIndxAmongPhotons>=0){
        if((*Photons_hasPixelSeed)[bestPhotonIndxAmongPhotons]<0.001) bestPhoHasPxlSeed=false;

        if(!bestPhoHasPxlSeed && bestPhoton.Pt()>100)
          {
            hasPho=1;hasPho_px=0;
          }

        else if(bestPhoHasPxlSeed && bestPhoton.Pt()>100) {hasPho_px=true;hasPho=0;}// npho_px++;}                                                                   
        else
          {
            hasPho_px=0;
            hasPho=0;
          }
      }

      double mindr_Pho_genlep=getGenLep(bestPhoton);

      if(!s_sample.Contains("data")){
	if(Debug)
	  cout<<"===loading tree entry === "<<jentry<<"\t"<<"  ===after preselection  ==="<<endl;
	h_madminPhotonDeltaR_beforeStitching->Fill(madMinPhotonDeltaR,wt);
	//remove the overlapped events from different MC samples as discussed in these slides - https://indico.cern.ch/event/1240842/contributions/5238493/attachments/2582794/4457078/ImpactOf_DiffPhotonIDs_ElectronFaking_Photon_28012023.pdf                                                                                                                                  

	if((s_sample.Contains("TTJets_HT")||s_sample.Contains("TTJets-HT")) && madHT<600) continue;
	if((s_sample.Contains("TTJets_inc")|| s_sample.Contains("TTJets_SingleLept") || s_sample.Contains("TTJets_DiLept")) && madHT>600) continue;
	if(!genphocheck)      {        genphomatch_before++;
	  //double mindr_Pho_genlep=getGenLep(bestPhoton);                                                                                                                                                                                                                                                                                                                        
	  if( s_sample.Contains("TTG") )
	    {  if(!hasGenPromptPhoton)              {
                h_selectBaselineYields->Fill("No gen prompt #gamma",wt);
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;          }
	      else if(hasGenPromptPhoton)              {
		h_selectBaselineYields->Fill("Gen prompt #gamma",wt);
		if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )) {//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);                                                                                                                                                                                                                                            
		  if(madMinPhotonDeltaR<0.5)                  h_selectBaselineYields->Fill("madMinPhotonDeltaR <0.5",wt);
                  if(mindr_Pho_genlep<0.5)                  h_selectBaselineYields->Fill("mindr_Pho_genlep<0.5",wt);
                  continue;}
                else
                  {                 if(madMinPhotonDeltaR >= 0.5)                     h_selectBaselineYields->Fill("mindR(q/g, #gamma)",wt);
                    if(mindr_Pho_genlep >=0.5)                        h_selectBaselineYields->Fill("mindR(l, #gamma)",wt);}
	      }
	    }
	  if(s_sample.Contains("WG"))
	    {            if(!hasGenPromptPhoton){                                 h_selectBaselineYields->Fill("No gen prompt #gamma",wt);
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;              }
	      else if(hasGenPromptPhoton)              {        h_selectBaselineYields->Fill("Gen prompt #gamma",wt);
		if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )){//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);                                                                                                                                                                                                                                             
                  if(madMinPhotonDeltaR<0.5)                    h_selectBaselineYields->Fill("madMinPhotonDeltaR <0.5",wt);
                  if(mindr_Pho_genlep<0.5)                    h_selectBaselineYields->Fill("mindr_Pho_genlep<0.5",wt);
                  continue;}
                else
                  {           if(madMinPhotonDeltaR >= 0.5)                     h_selectBaselineYields->Fill("mindR(q/g, #gamma)",wt);
                    if(mindr_Pho_genlep >=0.5)                        h_selectBaselineYields->Fill("mindR(l, #gamma)",wt);
                  }}}
	  if(s_sample.Contains("WJets"))
	    { if(!hasGenPromptPhoton){h_selectBaselineYields->Fill("No gen prompt #gamma",wt);
		if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;}
	      else if(hasGenPromptPhoton){  h_selectBaselineYields->Fill("Gen prompt #gamma",wt);
                if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)){//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);                                                                                                                                                                                                                                             
                  continue;}
                else{
		  if(madMinPhotonDeltaR >= 0.5)                      h_selectBaselineYields->Fill("pass_mindR(q/g, #gamma)",wt);
		  if(mindr_Pho_genlep >=0.5)                        h_selectBaselineYields->Fill("pass_mindR(l, #gamma)",wt);
                }}}

	  if(s_sample.Contains("TTJets_HT") || s_sample.Contains("TTJets-HT")||s_sample.Contains("TTJets-inc")|| s_sample.Contains("TTJets_inc") || s_sample.Contains("TTJets2_v17")||s_sample.Contains("TTJets") || s_sample.Contains("TTJets_SingleLept") || s_sample.Contains("TTJets_DiLept"))
	    {            if(!hasGenPromptPhoton){
                h_selectBaselineYields->Fill("No gen prompt #gamma",wt);
                if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;}
	      else if(hasGenPromptPhoton){
                h_selectBaselineYields->Fill("Gen prompt #gamma",wt);
                if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)){// h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);                                                                                                                                                                                                                                            
                  continue;}
                else{
		  if(madMinPhotonDeltaR >= 0.5)                      h_selectBaselineYields->Fill("pass_mindR(q/g, #gamma)",wt);
		  if(mindr_Pho_genlep >=0.5)                      h_selectBaselineYields->Fill("pass_mindR(l, #gamma)",wt);}                          }}
	  if(hasGenPromptPhoton && (s_sample.Contains("GJets")))
	    {            if(!(madMinPhotonDeltaR>0.4)) continue;}

	  if(hasGenPromptPhoton && (s_sample.Contains("QCD")))
	    {            if((madMinPhotonDeltaR>0.4 && hasGenPromptPhoton)) continue;}
	  if(hasGenPromptPhoton && ((s_sample.Contains("ZG"))|| (s_sample.Contains("ZNuNuG")) || s_sample.Contains("ZNuNuGJets")))          {
	    if(!(madMinPhotonDeltaR>0.5)) continue;}

	  if(hasGenPromptPhoton && ((s_sample.Contains("ZJets"))|| (s_sample.Contains("ZNuNuJets"))))
	    {            if(!(madMinPhotonDeltaR<=0.5)) continue;}
	  if(hasGenPromptPhoton && (s_sample.Contains("DYJets")))
	    {
	      if(!(madMinPhotonDeltaR<=0.5 || mindr_Pho_genlep<=0.5) ) continue;

	    }
	  genphomatch_after++;
	}
	h_madminPhotonDeltaR_preSelection->Fill(madMinPhotonDeltaR,wt);
      }
      h_selectBaselineYields_v1->Fill("overlap removal",wt);

      if(s_sample.Contains("data") && !s_sample.Contains("UL")&& applyTrgEff){
	if(tighte_trgpass==false)  continue;
      }
      h_selectBaselineYields_v1->Fill("tight trigger",wt);
      if(NMuons==0) h_selectBaselineYields_v1->Fill("0 #mu",wt);
      else continue;


      // //apply iso track veto                                                                                                                                          
      if(isoMuonTracks ==0 && isoPionTracks==0)
	{
	  h_selectBaselineYields_v1->Fill("Iso #mu & pion track",wt);
	}
      else continue;

      int leadGenPhoIdx=-100, mu_index=-100;
      int pass_accep_elec=0,fail_accep_elec=0,fail_isoEle=0,pass_isoElec=0,fail_IdElec=0,pass_IdElec=0;
      int pass_accep_mu=0,fail_accep_mu=0,fail_isoMu=0,pass_isoMu=0,fail_IdMu=0,pass_IdMu=0;
      // if(!s_sample.Contains("data")){
      // 	for(int i=0;i<GenParticles_;i++){
      // 	  if(GenParticles_v1[i].Pt()!=0){
      // 	    if((abs((*GenParticles_PdgId)[i])==22) && ((abs((*GenParticles_ParentId)[i])<=100) || ((*GenParticles_ParentId)[i]==2212) ) && (*GenParticles_Status)[i]==1 )
      // 	      {
      // 		leadGenPhoIdx = i;
      // 		genPho1 = (GenParticles_v1[i]);
      // 		v_genPho1.push_back(genPho1);
      // 	      }
      // 	  }
      // 	}
      // 	for(int i=0 ; i < GenElectrons_; i++)
      // 	  {
      // 	    if(GenElectrons_v1[i].Pt()!=0)// && (*GenElectrons)[i].Pt()!=inf)                                                                                                                                          
      // 	      {
      // 		nGenEle1++;
      // 		genEle1 = (GenElectrons_v1[i]);
      // 		v_genEle2.push_back(genEle1);
      // 		v_genLep2.push_back(genEle1);
      // 		//h_GenpT[3]->Fill(genEle1.Pt()); h_GenEta[3]->Fill(genEle1.Eta());
      // 	      }
      // 	  }
      // 	bool hadtau = false;
      // 	for(int i=0 ; i < GenMuons_; i++)
      //     {
      //       if(GenMuons_v1[i].Pt()!=0)// && (*GenMuons)[i].Pt()!=inf)                                                                                                
      // 	      {
      // 		nGenMu1++;
      //           genMu1 = (GenMuons_v1[i]);
      //           v_genMu2.push_back(genMu1);
      //           v_genLep2.push_back(genMu1);
      //           //h_GenpT[4]->Fill(genMu1.Pt()); h_GenEta[4]->Fill(genMu1.Eta());                                                                                     
      //         }
      //     }

      // 	for(int i=0 ; i < GenTaus_; i++)
      //     {
      //       if(GenTaus_v1[i].Pt()!=0) // && (*GenTaus)[i].Pt()!=inf)                                                                                                               
      // 	      {
      // 	      nGenTau1++;
      //           genTau1 = (GenTaus_v1[i]);
      //           v_genTau2.push_back(genTau1);
      //           if((*GenTaus_had)[i])
      //             nGentau_had1++;

      //         }
      //     }

      // 	if(v_genEle2.size()==0 ) continue;
      // }
      // sortTLorVec(&v_genEle2);
      // sortTLorVec(&v_genLep2);
      //      h_selectBaselineYields_v1->Fill("gen e=0",wt);
      if(Debug)
    	cout<<"=== load entry "<<jentry<<"\t"<<"GenTaus_  "<<GenTaus_<<"\t"<<"had tau decay"<<"\t"<<nGentau_had1<<endl;

      bool hasEMobj_pho=false, hasEMobj_pho_px=false,hasEMobj_elec=false; //flags - true if a good photon or good electron is present in the event 
      TLorentzVector emobj_pho, emobj_pho_px,emobj_elec;
      // int  hasPho=0;
      // int hasPho_px=0;
      // int hasEle=0;
      double mt_ele=0.0;
      int total_lost_el = 0,cr_el=0,sr_el,e_index=-1,nlep=0, NgenElec=0,leadGenPhoIdx1=0;
      bool Elec_passEtacut=false,Elec_passpTcut=false,Elec_passAccep = false,Elec_failAccep= false, Elec_passId= false, Elec_failId= false,Elec_passIso = false,Elec_failIso = false, genElec_passEtacut=false,genElec_passpTcut=false,genElec_passAccep = false, genElec_passId= false, genElec_passIso = false,genMu_passEtacut=false,genMu_passpTcut=false,genMu_passAccep = false, genMu_passId= false, genMu_passIso = false;
      bool Mu_passEtacut=false,Mu_passpTcut=false, Mu_passAccep= false,Mu_failAccep= false, Mu_passId= false, Mu_failId= false,Mu_passIso = false,Mu_failIso = false;
      int fail_isoElec=0;
      nlep=0; nelec_reco=0;
      sortTLorVec(&Electrons_v1);
      sortTLorVec(&Muons_v1);
      int e_index2=-1, two_e=0, one_e=0, one_e2=0,two_e2=0,two_e3=0;
      //    TLorentzVector //lep2;
      TLorentzVector tmp,tmp2;
      //double mt_ele=0;
      vector<TLorentzVector>reco_e;
      if(NElectrons==2 || NElectrons==1)
	{
	  for(int i=0 ; i<Electrons_; i++)
	    {
	      if((*Electrons_passIso)[i]==1 )
		{
		  nlep++;

		  if(nlep==1)
		    {
		      e_index=i;
		      tmp=Electrons_v1[e_index]; //1st electron                                                                                                     
		      one_e++;
		    }
		  if(nlep==2)
		    {
		      e_index2=i;
		      tmp2=Electrons_v1[e_index2]; // 2nd electron                                                                                                  
		      two_e++;
		    }
		}
	    }
	}

      if(nlep>0) reco_e.push_back(tmp);
      else if(nlep > 1) reco_e.push_back(tmp2);
      else if(nlep>2 || nlep==0) continue;
      h_selectBaselineYields_v1->Fill("NElectrons = 1 or 2 ",wt);
      h_nRecoElec->Fill(NElectrons,wt);
      h_ngenE->Fill(v_genEle2.size(),wt);
      //     if(tmp2.DeltaR(tmp)<0.2)continue;                                                                                                                       
      //if(nlep==0) continue;
      //      h_selectBaselineYields_v1->Fill("NElectrons = 1 or 2 ",wt);
      int minDR3indx=-100, minDR4indx=-100 , minDR5indx=-100, minDR6indx=-100;
      double minDR3=9999, minDR4=9999 , minDR5=9999, minDR6=9999;
      int par=0;
      vector<TLorentzVector> v_genz,v_hlt;
      TLorentzVector v_genz_,v_hlt_;
      if(!s_sample.Contains("data") )
      {

	for(int i=0;i<GenParticles_v1.size();i++){
	  if((abs((*GenParticles_ParentId)[i])==23 || abs((*GenParticles_ParentId)[i])==24) && (*GenParticles_Status)[i]==1  && abs((*GenParticles_PdgId)[i])==11 )
	    {
	      double dR=tmp.DeltaR(GenParticles_v1[i]);
	      if(dR<minDR3){minDR3=dR;minDR3indx=i;}
	      if(nlep==2){
		double dR2=tmp2.DeltaR(GenParticles_v1[i]);
		if(dR2<minDR4){minDR4=dR2;minDR4indx=i;}
	      }
	      if(nlep==1){
	      	double dR2=bestPhoton.DeltaR(GenParticles_v1[i]);
	      	if(dR2<minDR4){minDR4=dR2;minDR4indx=i;}
	      }
	    }
	}
	if((abs((*GenParticles_ParentId)[minDR4indx])==23 || abs((*GenParticles_ParentId)[minDR4indx])==24) && (*GenParticles_Status)[minDR4indx]==1 && abs((*GenParticles_PdgId)[minDR4indx])==11 )
	  minDR5=tmp.DeltaR(GenParticles_v1[minDR4indx]);
	if(nlep==2 && (abs((*GenParticles_ParentId)[minDR3indx])==23 || abs((*GenParticles_ParentId)[minDR3indx])==24) && (*GenParticles_Status)[minDR3indx]==1 && abs((*GenParticles_PdgId)[minDR3indx])==11 )
	  {
	    if(nlep==2) minDR6=tmp2.DeltaR(GenParticles_v1[minDR3indx]);
	    if(nlep==1) minDR6=bestPhoton.DeltaR(GenParticles_v1[minDR3indx]);
	  }

	//          if(hasPho==2) minDR6=bestPhoton.DeltaR(GenParticles_v1[minDR3indx]);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          

	if(minDR3indx==minDR4indx)
	  {
	    double minDR3_=minDR3, minDR4_=minDR4;
	    if(minDR3_>minDR4_)
	      {
		minDR3=9999, minDR3indx=-100;
		for(int i=0;i<GenParticles_v1.size();i++)
		  {//tag matching                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
		    if((abs((*GenParticles_ParentId)[i])==23 || abs((*GenParticles_ParentId)[i])==24) && (*GenParticles_Status)[i]==1 && abs((*GenParticles_PdgId)[i])==11  && i!=minDR4indx){
		      double dR=tmp.DeltaR(GenParticles_v1[i]);
		      if(dR<minDR3){minDR3=dR;minDR3indx=i;}
		    }
		  }
		if((abs((*GenParticles_ParentId)[minDR4indx])==23 || abs((*GenParticles_ParentId)[minDR4indx])==24) && (*GenParticles_Status)[minDR4indx]==1 && abs((*GenParticles_PdgId)[minDR4indx])==11 )
		  minDR5=tmp.DeltaR(GenParticles_v1[minDR4indx]);
		if(nlep==2 && (abs((*GenParticles_ParentId)[minDR3indx])==23 || abs((*GenParticles_ParentId)[minDR3indx])==24) && (*GenParticles_Status)[minDR3indx]==1 && abs((*GenParticles_PdgId)[minDR3indx])==11 )
		  minDR6=tmp2.DeltaR(GenParticles_v1[minDR3indx]);
		if(nlep==1 && (abs((*GenParticles_ParentId)[minDR3indx])==23 || abs((*GenParticles_ParentId)[minDR3indx])==24) && (*GenParticles_Status)[minDR3indx]==1 && abs((*GenParticles_PdgId)[minDR3indx])==11 )
		  minDR6=bestPhoton.DeltaR(GenParticles_v1[minDR3indx]);

	      }
	    if(minDR4_>minDR3_)
	      {
		minDR4=9999, minDR4indx=-100;
		for(int i=0;i<GenParticles_v1.size();i++)
		  {//tag matching                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
		    if((abs((*GenParticles_ParentId)[i])==23 || abs((*GenParticles_ParentId)[i])==24) && (*GenParticles_Status)[i]==1 && abs((*GenParticles_PdgId)[i])==11 && i!=minDR3indx){
		      if(nlep==2){
			double dR=tmp2.DeltaR(GenParticles_v1[i]);
			if(dR<minDR4){minDR4=dR;minDR4indx=i;}
		      }
		      if(nlep==1){
			double dR2=bestPhoton.DeltaR(GenParticles_v1[i]);
			if(dR2<minDR4){minDR4=dR2;minDR4indx=i;}
		      }
		    }
		  }
		if((abs((*GenParticles_ParentId)[minDR4indx])==23 || abs((*GenParticles_ParentId)[minDR4indx])==24) && (*GenParticles_Status)[minDR4indx]==1 && abs((*GenParticles_PdgId)[minDR4indx])==11 )
		  minDR5=tmp.DeltaR(GenParticles_v1[minDR4indx]);
		if((abs((*GenParticles_ParentId)[minDR3indx])==23 || abs((*GenParticles_ParentId)[minDR3indx])==24) && (*GenParticles_Status)[minDR3indx]==1 && abs((*GenParticles_PdgId)[minDR3indx])==11 )
		  {
		    if(nlep==2) minDR6=tmp2.DeltaR(GenParticles_v1[minDR3indx]);
		    if(nlep==1) minDR6=bestPhoton.DeltaR(GenParticles_v1[minDR3indx]);
		  }
	      }

	  }

	if(minDR3indx!=minDR4indx){
	  if(minDR3>=0) v_genz.push_back(GenParticles_v1[minDR3indx]);
	  if(minDR4>=0) v_genz.push_back(GenParticles_v1[minDR4indx]);
	  //  cout<<jentry<<" : "<<hasPho<<" , "<<minDR3<<" , "<<minDR4<<endl;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
	}
	if(tmp2.DeltaR(tmp)<0.2)continue;
	h_selectBaselineYields_v1->Fill("dR(tmp,tmp2)>0.2 ",wt);
	if(minDR3indx==minDR4indx && minDR3indx!=-100)
	  cout<<"3        "<<jentry<<" : "<<minDR3<<" , "<<minDR4<<" , "<<minDR5<<" , "<<minDR6<<endl;
	    // for(int i=0;i<GenElectrons_v1.size();i++){
	    //   v_genz.push_back(GenElectrons_v1[i]);
	    // }
	if(v_genz.size()==0 || nlep==0) continue;
	h_selectBaselineYields_v1->Fill("gen e!=0 ",wt);
      }

      float ratio =0.0, ratio1=0.0, mindr_genElecPho=-999, mindr_=-9999;
      int count_genEle=0,count_recEle=0;
      vector<TLorentzVector> goodPho_n;
      vector<int> goodPhoIndx_n;
      for(int iPho=0;iPho<Photons_;iPho++){
    	if((*Photons_fullID)[iPho] && ((*Photons_hasPixelSeed)[iPho]<0.001)) {
    	  goodPho_n.push_back( Photons_v1[iPho] );
    	  goodPhoIndx_n.push_back(iPho);
    	}
      }
      int nPhotons = goodPho_n.size();
      
      if(s_sample.Contains("data"))//  && !s_data.Contains("2018"))//_v18_2016") && s_data.Contains("data_v18_2017") )                                                                                     
        {
          if(Debug)  cout<<" printing before HLT loop "<<"\t"<<HLTElectronObjects_<<endl;
          for(int i=0;i<HLTElectronObjects_;i++)
            {//tag matching                                                                                                                                                                                
              double dR=tmp.DeltaR(HLTElectronObjects_v1[i]);
              if(dR<minDR3){minDR3=dR;minDR3indx=i;}
              // if(nlep>1){
	      double dR2=tmp2.DeltaR(HLTElectronObjects_v1[i]);
		if(dR2<minDR4){minDR4=dR2;minDR4indx=i;}
              //  }
              // else
              //   {
              //     double dR2=bestPhoton.DeltaR(HLTElectronObjects_v1[i]);
              //     if(dR2<minDR4){minDR4=dR2;minDR4indx=i;}

              //   }
            }
          if(Debug)  cout<<" printing before HLT loop "<<"\t"<<HLTElectronObjects_<<"\t"<<HLTElectronObjects_v1.size()<<"\t"<<minDR3indx<<"\t"<<minDR4indx<<"\t"<<minDR3<<"\t"<<minDR4<<"\t"<<nlep<<"\t"<<NElectrons<<"\t"<<Photons_<<"\t"<<(*TriggerPass)[51]<<endl;
          if(HLTElectronObjects_!=0){
	    minDR5=tmp.DeltaR(HLTElectronObjects_v1[minDR4indx]);
	    minDR6=tmp2.DeltaR(HLTElectronObjects_v1[minDR3indx]);
	    //saving HLT electrons matching with reco electron                                                                                                                                       
	    if(minDR3>=0) v_hlt.push_back(HLTElectronObjects_v1[minDR3indx]); // confirm these conditions again                                                                                        
	    if(minDR4>=0) v_hlt.push_back(HLTElectronObjects_v1[minDR4indx]);
	    }
	  h_nHLTE->Fill(HLTElectronObjects_,wt);
          //h_selectBaselineYields_v1->Fill("NHLT- Electrons != 0 ",wt);
          if(v_hlt.size()==0 || nlep==0) continue;                                                                                                                              h_selectBaselineYields_v1->Fill("gen e!=0 ",wt);                          
	  //          h_selectBaselineYields_v1->Fill("NHLT- Electrons != 0 ",wt);
	}
      
      if(!s_data.Contains("2018")&& s_sample.Contains("data")) {  if (v_hlt.size()==0) continue;}
     if(Debug)  cout<<" printing before nlep ==1 loop "<<"\t"<<nlep<<endl;
      //      if(Debug
      // 0 index for tag object and 1 index for probe                                                                                                      
      if(nlep==1 && hasPho==1)
	{
	  double mindr=100;
	  double mindr2=100;
	  if(!s_sample.Contains("data")) {mindr=tmp.DeltaR(v_genz[0]); mindr2=0.1; mindr2=bestPhoton.DeltaR(v_genz[1]);}//mindr=tmp.DeltaR(v_genz[0]); mindr2=0.1; mindr2=bestPhoton.DeltaR(v_genz[1]);//mindr=MinDr(tmp,v_genEle2); mindr2 =MinDr(bestPhoton,v_genEle2);}// mindr=tmp.DeltaR(v_genEle2[1]); mindr2=0.1; mindr2=bestPhoton.DeltaR(v_genEle2[0]);}
	  if(s_sample.Contains("data"))
	    {
	      //    mindr=MinDr(tmp,*HLTElectronObjects);                                                                                                         
	      if(Debug)  cout<<" printing in nlep ==1 loop "<<"\t"<<nlep<<endl;

	      if(s_data.Contains("2018"))mindr=0.1;
	      else mindr=tmp.DeltaR(v_hlt[0]);//MinDr(tmp,HLTElectronObjects_v1);//v_hlt);//tmp.DeltaR(v_hlt[0]);//MinDr(tmp,v_hlt);//tmp.DeltaR(v_hlt[0]);	      
	      // if(v_hlt.size()>1)
	      // 	mindr2=bestPhoton.DeltaR(v_hlt[1]);//0.1; // this should be checked wrt v_hlt[1] if exist Think about it
	      // else
		mindr2=0.1;
	    }
	  if(s_data.Contains("2018"))
	    { mindr=0.1, mindr2=0.1;} // should be removed once you have updated data samples
	  h_dR_HLTvsRecoElec1->Fill(mindr,wt);
	  h_dR_HLTvsRecoElec2->Fill(mindr2,wt);
	  h_dR_HLTvsRecoElec1_categ3->Fill(mindr2,wt);
	  h_dR_HLTvsRecoElec2_categ3->Fill(mindr,wt);
	  h_selectBaselineYields_v1->Fill("lep==1 & haspho",wt);
	  if(Debug)  cout<<" printing in nlep ==1 loop "<<"\t"<<nlep<<"\t"<<mindr<<"\t"<<mindr2<<endl;
	  //mindr2=0.1; //second HLT electron does not match with photon 

	  if(tmp.Pt()>40 && mindr<0.2 && (*Electrons_tightID)[e_index] && mindr2<0.2 && hasPho==1) // why by hand making mindR2 to be <0.2 - because in data no information about second electron which might be faking as photon.
	    {
	      h_selectBaselineYields_v1->Fill("lep==1 & haspho & pass dR cond",wt);
	      hasEle=1;
	      v_lep2=tmp; //tag                                                                                                                                      
	      if(hasPho==1)v_lep1=bestPhoton; //probe                                                                                                                
	      if(!s_sample.Contains("data")) v_genz_=v_genz[0];
	      if(s_sample.Contains("data") && !s_data.Contains("2018")) v_hlt_=v_hlt[0];
	      one_e2++;
	    }
	  // if(s_data.Contains("data"))h_minDr_bestphoremEle->Fill(mindr,wt);
	  // if(!s_data.Contains("data"))
	  //   {
	  //     h_minDr_bestphoremEle->Fill(mindr,wt);
	  //   }
	}

      if(nlep==2)
	{
	  if(Debug)  cout<<" printing in nlep ==2 loop "<<"\t"<<nlep<<"\t"<<endl;

	  int a=0;
	  double mindr=100;
	  double mindr2=100;
	  double mindr3=100;
	  double mindr4=100;
	  h_selectBaselineYields_v1->Fill("lep==2 & !haspho ",wt);
	  if(!s_sample.Contains("data"))
	    {
	      // mindr=tmp2.DeltaR(v_genEle2[1]);
	      // mindr2=tmp.DeltaR(v_genEle2[0]);
	      // mindr=MinDr(tmp2,v_genz);//
	      // mindr2=MinDr(tmp,v_genz);
	      mindr=tmp2.DeltaR(v_genz[1]);
	      mindr2=tmp.DeltaR(v_genz[0]);

	      mindr3=0.1;
	      mindr4=0.1;
	      //h_minDr_bestphoremEle->Fill(mindr2,wt);
	    }
	  if(s_sample.Contains("data"))
	    {
	      if(!s_data.Contains("2018")){
	      // mindr2=tmp.DeltaR(v_hlt[0]);
	      // mindr=tmp2.DeltaR(v_hlt[1]);
		//mindr2=MinDr(tmp,HLTElectronObjects_v1);mindr = MinDr(tmp2,HLTElectronObjects_v1);//v_hlt); mindr=MinDr(tmp2,v_hlt);
		mindr2=tmp.DeltaR(v_hlt[0]);
		mindr=tmp2.DeltaR(v_hlt[1]);

	      }
	      else
		{
		  mindr2=0.1;
		  mindr=0.1;
		}
	      mindr3=0.1;
	      mindr4=0.1;
	      //h_minDr_bestphoremEle->Fill(mindr2,wt);
	    }
	  // if (s_data.Contains("2018")){
	  //   mindr2=0.1;
	  //   mindr=0.1;
	  //   mindr3=0.1;
	  //   mindr4=0.1;

	  // }
	  h_dR_HLTvsRecoElec1_categ1->Fill(mindr2,wt);
	  h_dR_HLTvsRecoElec2_categ1->Fill(mindr,wt);
	  h_dR_HLTvsRecoElec1_categ3->Fill(mindr2,wt);
	  h_dR_HLTvsRecoElec2_categ3->Fill(mindr,wt);

	  if(tmp.Pt()>100 && tmp2.Pt()>100 && a==0 && mindr<0.2 && mindr2<0.2 && (*Electrons_tightID)[e_index] && (*Electrons_tightID)[e_index2])
	    {
	      //      v_lep1=tmp; v_lep2=tmp2; hasEle=2;wt=2*wt;a++;                                                                                                    
	      h_selectBaselineYields_v1->Fill("lep==2 & pass two tight e ",wt);
	      if(jentry%2==0){v_lep1=tmp; v_lep2=tmp2; hasEle=2;wt=2*wt; if(s_sample.Contains("data") && !s_data.Contains("2018") ) v_hlt_=v_hlt[1]; if(!s_sample.Contains("data"))  v_genz_=v_genz[1];a++;}
	      else{v_lep2=tmp; v_lep1=tmp2; hasEle=2;wt=2*wt; if(s_sample.Contains("data") && !s_data.Contains("2018")) v_hlt_=v_hlt[0]; if(!s_sample.Contains("data"))  v_genz_=v_genz[0]; a++; }
	      two_e2++;
	      //            two_e3++;                                                                                                                                   
	      //if(jentry==18163)                                                                                                                                       
	      //cout<<jentry<<" : "<<v_lep2.DeltaR(v_genz[0])<<" , "<<v_lep1.DeltaR(v_genz[1])<<endl;                                                                   
	      // h_dR_HLTvsRecoElec1_categ1->Fill(mindr2,wt);
	      // h_dR_HLTvsRecoElec2_categ1->Fill(mindr,wt);
 
	    }
	  
	  if(tmp.Pt()>100 && tmp2.Pt()>40 && a==0 && mindr<0.2 && (*Electrons_tightID)[e_index2] && mindr4<0.2) // there should be additiona mindR condition for second electron - check the plot and update
	    {
	      hasEle=2; v_lep2=tmp2; v_lep1=tmp;
	      a++;
	      h_selectBaselineYields_v1->Fill("lep==2 & pass only tight e ",wt);

	      if(!s_sample.Contains("data")) v_genz_=v_genz[1];
	      if(s_sample.Contains("data") && !s_data.Contains("2018")) v_hlt_=v_hlt[1];
	      two_e2++;
	      h_dR_HLTvsRecoElec1_categ2->Fill(mindr2,wt);
	      h_dR_HLTvsRecoElec2_categ2->Fill(mindr,wt);
	      // h_dR_HLTvsRecoElec1_categ3->Fill(mindr2,wt);
              // h_dR_HLTvsRecoElec2_categ3->Fill(mindr,wt);

	    }

	  if(tmp2.Pt()>100 && tmp.Pt()>40 && a==0 && mindr2<0.2 && (*Electrons_tightID)[e_index] && mindr3<0.2)
	    {
	      hasEle=2; v_lep2=tmp; v_lep1=tmp2;
	      a++;
	      h_selectBaselineYields_v1->Fill("lep==2 & pass only tight e ",wt);
	      if(!s_sample.Contains("data")) v_genz_=v_genz[0];
	      if(s_sample.Contains("data")&& !s_data.Contains("2018")) v_hlt_=v_hlt[0];
	      two_e3++;
	      two_e2++;
	      // h_dR_HLTvsRecoElec1_categ3->Fill(mindr2,wt);
	      // h_dR_HLTvsRecoElec2_categ3->Fill(mindr,wt);

	    }


	}
      if(Debug)  cout<<" outside nlep loops "<<"\t"<<nlep<<"\t"<<endl;

      bool bestEMObjIsEle=false, bestEMObjIsEle_px=false,bestEMobj=false; int probe_ele=0,tag_ele=0; //int hastagEle=0,hasprobeEle=0,probe_ele=0,tag_ele=0;
      TLorentzVector bestEMObj,tagEMObj;     

      if(hasEle==2 || hasEle==1){ hastagEle=1;}
      if (hasEle==2 && hasPho==0 && hastagEle==1) {bestEMObjIsEle=true; bestEMObjIsEle_px=false;hasprobeEle=1; tagEMObj = v_lep2; bestEMObj = v_lep1; bestEMobj=true; probe_ele++;}
      if(hasEle==1 &&  hasPho==1 && hastagEle==1) {bestEMObjIsEle=false;  bestEMObjIsEle_px=false; hasprobeEle=1; tagEMObj = v_lep2; bestEMObj = bestPhoton; bestEMobj=true; probe_ele++;}
      if(hasprobeEle==0) continue;
      else
	tag_ele++;
      h_selectBaselineYields_v1->Fill("After no probe object found ",wt);

      if(Debug)  cout<<" printing after selecting em objects "<<"\t"<<nlep<<"\t"<<endl;

      if(!bestEMobj) continue;
      h_selectBaselineYields_v1->Fill("1 tag & 1 probe ",wt);
      if(bestEMObj.DeltaR(tagEMObj) < 0.2) continue;
      h_selectBaselineYields_v1->Fill("dR(tag, probe)>0.2 ",wt);
      /////======Pixel Power supply issue ===============//////////////                                                                                                //update these 
      if(apply_pixelveto){
	if(s_data.Contains("2017") && ((bestEMObj.Eta() > 1.0 && bestEMObj.Eta() < 1.5 && bestEMObj.Phi() > 2.7) || (tagEMObj.Eta() > 1.0 && tagEMObj.Eta() < 1.5 && tagEMObj.Phi() > 2.7))) continue;
	if(s_data.Contains("2018") && ((bestEMObj.Eta() > 1.0 && bestEMObj.Eta() < 1.5 && bestEMObj.Phi() > 0.4 && bestEMObj.Phi() < 0.8) || (tagEMObj.Eta() > 1.0 && tagEMObj.Eta() < 1.5 && tagEMObj.Phi() > 0.4 && tagEMObj.Phi() < 0.8))) continue;
      }

      h_selectBaselineYields_v1->Fill("after pixel veto first",wt);

      bool fakePhoton=0,showEnt=0;
      //    oldfake=true;                                                                                                                                               
      double minDr=99999;
      bool notrack=true;
      int ntracks=0, minDrindx=-100, track_idx=-100;
      for(int i=0;i<TAPElectronTracks_;i++){
	if(TAPElectronTracks_v1[i].Pt()>5)
	  {
	    double dR2=tagEMObj.DeltaR(TAPElectronTracks_v1[i]);
	    if(dR2<minDr){minDr=dR2;minDrindx=i;}
	    if(dR2<0.2) ntracks++;
	  }
      }
      
      if(minDr<0.2)
	{
	  notrack=false;
	  track_idx=minDrindx;
	}
      if(notrack || ntracks==0) continue;
      h_selectBaselineYields_v1->Fill("!no track && ntracks",wt);

      if(apply_pixelveto){
	if(s_data.Contains("2017") && (TAPElectronTracks_v1[track_idx].Eta() > 0.0 && TAPElectronTracks_v1[track_idx].Eta() < 1.5 && TAPElectronTracks_v1[track_idx].Phi() > 2.7)) continue;
	if(s_data.Contains("2018") && (TAPElectronTracks_v1[track_idx].Eta() > 0.0 && TAPElectronTracks_v1[track_idx].Eta() < 1.5 && TAPElectronTracks_v1[track_idx].Phi() > 0.4 && TAPElectronTracks_v1[track_idx].Phi() < 0.8)) continue;
      }

      h_selectBaselineYields_v1->Fill("after pixel veto for track matched to tag e",wt);

      h_tagpT_vsTrackPT->Fill(TAPElectronTracks_v1[track_idx].Pt(),tagEMObj.Pt(),wt);
      h_nTracks->Fill(ntracks,wt);
      h_tagEle_pT->Fill(tagEMObj.Pt(),wt);
      h_tagEle_Eta->Fill(tagEMObj.Eta(),wt);
      h_tagEle_Phi->Fill(tagEMObj.Phi(),wt);
      h_trackEle_pT->Fill(TAPElectronTracks_v1[track_idx].Pt(),wt);
      h_trackEle_Eta->Fill(TAPElectronTracks_v1[track_idx].Eta(),wt);
      h_trackEle_Phi->Fill(TAPElectronTracks_v1[track_idx].Phi(),wt);
      
      if(!bestEMObjIsEle && bestEMobj)
	{
	  if(!s_sample.Contains("data"))
	    {

	  
	    //fail_realphoton++;
	    fakePhoton=true;
	  
	  
	  if(!fakePhoton) continue;
	    }
	  if(ntracks>1) continue;
	  h_selectBaselineYields_v1->Fill("ntracsk<1",wt);

	}


           
      // *******************  Selecting Jet objects ********************************//
      int minDRindx=-100,minDR2indx=-100,phoMatchingJetIndx=-100,hadJetID=-999,bJet1Idx=-100,nHadJets=0, nbjets=0,photonMatchingJetIndx=-100;
      double qmulti=0,leadjet_qmulti=-1,leadjet_Pt=-1,leadbjet_tag=-1;
      double minDR=99999,minDR2=9999,ST=0,Ht=0,remJetPt=0;
       vector<TLorentzVector> hadJets,bjets,remJets;
      hadJets.clear();
      bjets.clear();
      remJets.clear();
      vector<int> hadJets_hadronFlavor;
      vector<bool> hadJets_HTMask;
      vector<double> hadJets_bJetTagDeepCSVBvsAll;
      vector<TLorentzVector> nonbjets;
      vector<int> jetMatchindx;
      TLorentzVector Jet_matched;
      TLorentzVector leadJet;
      bool recoJetMatch_recoPho=false, genJetMatch_recoPho=false;
      if(Debug)
        cout<<"===load tree entry check2 at entry ==="<<"\t"<<jentry<<endl;
      for(int i=0;i<Jets_;i++){
	if( (Jets_v1[i].Pt() > 30.0) && (abs(Jets_v1[i].Eta()) <= 2.4) )
	  {
            {
	      double dR=bestEMObj.DeltaR(Jets_v1[i]);
	      double dR2=tagEMObj.DeltaR((*Jets)[i]);
	      if(dR2<minDR2){minDR2=dR2;minDR2indx=i;}
	      if(dR<minDR){minDR=dR;minDRindx=i;}
	    }
	  }
      }
      // if( minDR > 0.3 )  minDRindx   = -100;
      // if( minDR2 > 0.3 )  minDR2indx   = -100;

      int q=0;
      for(int i=0;i<Jets_;i++){
	if( (Jets_v1[i].Pt() > 30.0) && (abs(Jets_v1[i].Eta()) <= 2.4) ){
	  //q++;
	  //if((*Jets_ID)[i]) h2_mindr_jetmatchphoratio->Fill(minDR,Jets_v1[minDRindx].Pt()/bestEMObj.Pt(),wt);
	  if( (Jets_v1[i].Pt() > 30.0) && (abs(Jets_v1[i].Eta()) <= 2.4) ){
	    //	    	    	  if( !(minDR < 0.3 && i==minDRindx) ){
	    //	    if((i!=minDRindx) && ){//&& (i!=minDR2indx) ){
	    if( !(minDR < 0.3 && i==minDRindx) && !(minDR2 < 0.3 && i==minDR2indx) ){
	      hadJetID=(*Jets_ID)[i];
	      if(hadJetID)
		{
		  q++;
		  hadJets.push_back(Jets_v1[i]);
		  hadJets_hadronFlavor.push_back((*Jets_hadronFlavor)[i]);
		  hadJets_HTMask.push_back((*Jets_HTMask)[i]);
		  hadJets_bJetTagDeepCSVBvsAll.push_back((*Jets_bJetTagDeepCSVBvsAll)[i]);
		  if(q==1) leadjet_qmulti=(*Jets_chargedMultiplicity)[q];
		  if(q==1) leadjet_Pt=Jets_v1[q].Pt();
		  if((*Jets_bJetTagDeepCSVBvsAll)[i] > deepCSVvalue){
		    bjets.push_back(Jets_v1[i]); bJet1Idx = i;}
		}
	    }
	  else if(minDR<0.3 && i==minDRindx)
	    Jet_matched = Jets_v1[minDRindx];
	  }  
	}
      }
      //    cout<<jentry<<" ====== break ========"<<endl;
      if( minDR<0.3)
	{
	  photonMatchingJetIndx=minDRindx;
	  qmulti=(*Jets_chargedMultiplicity)[photonMatchingJetIndx];
	  leadbjet_tag=(*Jets_bJetTagDeepCSVBvsAll)[photonMatchingJetIndx];
	}
     
      if( photonMatchingJetIndx>=0 ){
	if( (Jets_v1[photonMatchingJetIndx].Pt()) > 1.1*(bestEMObj.Pt()) ){
	  if( ((Jets_v1[photonMatchingJetIndx] - bestEMObj).Pt())>30){
	    hadJets.push_back( Jets_v1[photonMatchingJetIndx] - bestEMObj );
	    remJetPt=(Jets_v1[photonMatchingJetIndx] - bestEMObj).Pt();
	    remJets.push_back( Jets_v1[photonMatchingJetIndx] - bestEMObj );
	  }
	}
      }


      for(int i=0;i<hadJets.size();i++)
	{
	  if( (abs(hadJets[i].Eta()) < 2.4) ){ST=ST+(hadJets[i].Pt());}
	  if( (abs(hadJets[i].Eta()) < 2.4) ){nHadJets++;}
	}
    
      for(int i=0;i<bjets.size();i++){
	if( (abs(bjets[i].Eta()) < 2.4) ){nbjets++;}
      }
      if( minDR<0.3 ){
	ST=ST+bestEMObj.Pt();
      }
     
      
     
      sortTLorVec(&hadJets);
      //      if(minDRindx<0) continue;
    if(hadJets.size()==0) continue;
    if(Debug)
      cout<<"===load tree entry ===  "<<"\t"<<jentry<<"\t"<<"No of B-Jets ===  "<<bjets.size()<<endl;
    BTags = nbjets;
     double genmindr=99999, recojetmindr=99999;
    if(Debug)
      cout<<"===load tree entry === "<<"\t"<<jentry<<" hadJets.size()  "<< hadJets.size()<<endl;    
    if(hadJets.size()<BTags)
      cout<<"Entry: "<<jentry<<"\t"<<"njets "<<nHadJets<<"\t"<<" hadJets.size() "<<hadJets.size()<<endl;
    // // ********* This is to account into all visible energy: Adding photon matching with Jet ******************//
    TLorentzVector Met; 
    Met.SetPtEtaPhiE(MET,0,METPhi,0);
    double mTPhoMET = sqrt(2*(bestEMObj.Pt())*MET*(1-cos(DeltaPhi(METPhi,bestEMObj.Phi()))));
    double dPhi_PhoMET= abs(Met.DeltaPhi(bestEMObj));//DeltaPhi(METPhi,bestPhoton.Phi()));
     
    if(Debug)
      cout<<"===load tree entry ==="<<"\t"<<jentry<<" mTPhoMET "<<mTPhoMET<<" METPhi  "<<bestPhoton.Eta()<<" Photons_ "<<Photons_<<" hadJets.size() " <<"\t"<<hadJets.size()<<" NElectrons "<<NElectrons<< "  Jets  "<<Jets_<<"\t"<< "Jets size  "<<Jets_v1.size()<<endl;
    if(jentry<10)
      cout<<"madMinPhotonDeltaR  "<<madMinPhotonDeltaR<<endl;
    h_madminPhotonDeltaR_noSelection->Fill(madMinPhotonDeltaR,wt);
			    		    
    // TLorentzVector Met;
    // Met.SetPtEtaPhiE(MET,0,METPhi,0);
    Double_t deta_jet_pho= 0.0,deta_jet_met=0.0,deta_met_pho=0.0;
    if(Debug)
      cout<<"METPhi "<<METPhi<<endl;
    double mT= 0.0, dPhi_METjet1=5, dPhi_METjet2=5, dPhi_phojet1=5, dPhi_phojet2=5, dPhi_phoMET=5;
    dPhi_METjet1 = abs(Met.DeltaPhi(hadJets[0]));
    dPhi_METjet2 = abs(Met.DeltaPhi(hadJets[1]));
     
    if(Debug)
      cout<<"===load tree entry ==="<<"\t"<<jentry<<" hadJets.size() "<<hadJets.size()<<endl;
    // h_selectBaselineYields_->Fill("all",wt);
    
    if(Debug)
      cout<<"===load tree entry ==="<<"\t"<<jentry<<" bestPhoHasPxlSeed "<<bestPhoHasPxlSeed<<endl;
    
    //cout<<"Event : " <<jentry <<"\t"<<(*Photons_mvaValuesID)[bestPhotonIndxAmongPhotons] << " (*Photons_mvaValuesID)[iPho]"<< "\t"<<bestPhotonIndxAmongPhotons<< " Photons _ "<<Photons_ <<"\t Photons_mvaValuesID->size() "<<Photons_mvaValuesID->size()<<"\t"<<Photons_v1.size()<<"\t"<<Photons_fullID->size()<<"\t"<<Photons_cutBasedID->size()<<endl;
    if (nHadJets>=2)h_selectBaselineYields_v1->Fill("njets>=2",wt);
    else continue;
    if(MET<200) 
      h_selectBaselineYields_v1->Fill("MET<200",wt);
    else continue;
    if(ST>200)
      h_selectBaselineYields_v1->Fill("ST>200",wt);
    else continue;
    // //recommended MET filtser for UL - taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Analysis_Recommendations_should
    // //check nvtx filter - less efficiency for 2016 - https://indico.cern.ch/event/1057110/#27-met-filters-performance-stu
    // if(PFCaloMETRatio >=  5) continue;
    // if(s_sample.Contains("data"))
    //   {
    // 	if(s_data.Contains("2016"))
    // 	  ecalBadCalibReducedExtraFilter=1;
    // 	if(!(NVtx > 0 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadPFMuonFilter && eeBadScFilter==1 && ecalBadCalibReducedExtraFilter==1)) continue;

    //   }
    // if(!s_sample.Contains("data"))
    //   {
    // 	eeBadScFilter = 1;
    // 	if(!(NVtx > 0 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 && HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter==1 && BadPFMuonFilter && eeBadScFilter==1 && ecalBadCalibReducedExtraFilter==1)) continue;
    //   }

    // if (s_data.Contains("2017") || s_data.Contains("2018"))
    //   if(!(PrimaryVertexFilter==1 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 &&HBHEIsoNoiseFilter==1 && EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter==1 && BadPFMuonDzFilter==1 && eeBadScFilter==1 && ecalBadCalibFilter==1 && NVtx>0)) continue;
    // if (s_data.Contains("2016")){
    //   if(!(PrimaryVertexFilter==1 && globalSuperTightHalo2016Filter==1 && HBHENoiseFilter==1 &&HBHEIsoNoiseFilter==1 &&EcalDeadCellTriggerPrimitiveFilter == 1 && BadPFMuonFilter==1 && BadPFMuonDzFilter==1 && eeBadScFilter==1 )) continue;
    // }

    h_selectBaselineYields_v1->Fill("EventCleaning",wt);
    if(highdphi){
    if(dPhi_METjet1 > 0.3 && dPhi_METjet2 > 0.3 )
      {
        h_selectBaselineYields_v1->Fill("dPhi1 & dPhi2 >= 0.3",wt);
      }
    else continue;
    }
    // if(photonMatchingJetIndx>=0 && (Jets_v1[photonMatchingJetIndx].Pt())/(bestEMObj.Pt()) < 1.0) continue;
    // h_selectBaselineYields_v1->Fill("jet pT/em obj pT<1",wt);
    // if(photonMatchingJetIndx<0) continue;
    // h_selectBaselineYields_v1->Fill("No jet matched with EM obj",wt);
    // if(higMET)
    //   { 
    // 	if (bestPhoton.Pt()<=100)continue;
    // 	if(MET<=200) continue;                                                                                              
    //   }
     
 //    //    FillHistogram_Kinematics(2,nHadJets,BTags,bestPhoton.Pt(),mTPhoMET,dPhi_PhoMET,ST,1);
//     double mindr_Pho_genlep=getGenLep(bestPhoton);
    
//     if(!s_sample.Contains("data")){
//       if(Debug)
//       cout<<"===loading tree entry === "<<jentry<<"\t"<<"  ===after preselection  ==="<<endl;
//     h_madminPhotonDeltaR_beforeStitching->Fill(madMinPhotonDeltaR,wt);
//     //remove the overlapped events from different MC samples as discussed in these slides - https://indico.cern.ch/event/1240842/contributions/5238493/attachments/2582794/4457078/ImpactOf_DiffPhotonIDs_ElectronFaking_Photon_28012023.pdf

//     if((s_sample.Contains("TTJets_HT")||s_sample.Contains("TTJets-HT")) && madHT<600) continue;        
//     if((s_sample.Contains("TTJets_inc")|| s_sample.Contains("TTJets_SingleLept") || s_sample.Contains("TTJets_DiLept")) && madHT>600) continue;
//     if(!genphocheck)      {        genphomatch_before++;
//       //double mindr_Pho_genlep=getGenLep(bestPhoton);
//         if( s_sample.Contains("TTG") )
//           {  if(!hasGenPromptPhoton)              {		
//     		h_selectBaselineYields->Fill("No gen prompt #gamma",wt);
//                 if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;	    }
//             else if(hasGenPromptPhoton)              {
// 	      h_selectBaselineYields->Fill("Gen prompt #gamma",wt);
// 	      if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )) {//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt); 
// 		if(madMinPhotonDeltaR<0.5)		    h_selectBaselineYields->Fill("madMinPhotonDeltaR <0.5",wt);
// 		  if(mindr_Pho_genlep<0.5)		    h_selectBaselineYields->Fill("mindr_Pho_genlep<0.5",wt);		
//     		  continue;}
//     		else
//     		  {    		    if(madMinPhotonDeltaR >= 0.5)    		      h_selectBaselineYields->Fill("mindR(q/g, #gamma)",wt);
//     		    if(mindr_Pho_genlep >=0.5)    		      h_selectBaselineYields->Fill("mindR(l, #gamma)",wt);}
// 	    }
// 	  }                                             
//         if(s_sample.Contains("WG"))
//           {            if(!hasGenPromptPhoton){                  		h_selectBaselineYields->Fill("No gen prompt #gamma",wt);
//                 if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;              }
//             else if(hasGenPromptPhoton)              {	      h_selectBaselineYields->Fill("Gen prompt #gamma",wt);
// 	      if(!(madMinPhotonDeltaR >= 0.5 && mindr_Pho_genlep >=0.5 )){//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt); 
// 		  if(madMinPhotonDeltaR<0.5)                    h_selectBaselineYields->Fill("madMinPhotonDeltaR <0.5",wt);
//                   if(mindr_Pho_genlep<0.5)                    h_selectBaselineYields->Fill("mindr_Pho_genlep<0.5",wt);
//     		  continue;}
//     		else
//                   {           if(madMinPhotonDeltaR >= 0.5)                     h_selectBaselineYields->Fill("mindR(q/g, #gamma)",wt);
//                     if(mindr_Pho_genlep >=0.5)    		      h_selectBaselineYields->Fill("mindR(l, #gamma)",wt);
//                   }}}
//         if(s_sample.Contains("WJets"))
//           { if(!hasGenPromptPhoton){h_selectBaselineYields->Fill("No gen prompt #gamma",wt);
// 	      if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;}              
//             else if(hasGenPromptPhoton){  h_selectBaselineYields->Fill("Gen prompt #gamma",wt);
//                 if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)){//h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt);   
//     		  continue;}
//     		else{                  
//                     if(madMinPhotonDeltaR >= 0.5)                      h_selectBaselineYields->Fill("pass_mindR(q/g, #gamma)",wt);
//                     if(mindr_Pho_genlep >=0.5)    		      h_selectBaselineYields->Fill("pass_mindR(l, #gamma)",wt);
// 		}}}
                        
//         if(s_sample.Contains("TTJets_HT") || s_sample.Contains("TTJets-HT")||s_sample.Contains("TTJets-inc")|| s_sample.Contains("TTJets_inc") || s_sample.Contains("TTJets2_v17")||s_sample.Contains("TTJets") || s_sample.Contains("TTJets_SingleLept") || s_sample.Contains("TTJets_DiLept"))
//           {            if(!hasGenPromptPhoton){           
//     		h_selectBaselineYields->Fill("No gen prompt #gamma",wt);
//                 if(jentry==0)cout<<"**********processing "<<s_sample<<" with non-prompt Gen photon"<<endl;}              
//             else if(hasGenPromptPhoton){             
//     		h_selectBaselineYields->Fill("Gen prompt #gamma",wt);
//                 if(!(madMinPhotonDeltaR < 0.5 || mindr_Pho_genlep < 0.5)){// h_phoPt_promptPho_rejected->Fill(bestPhoton.Pt(),wt); 
//     		  continue;}
//     		else{                  
//                     if(madMinPhotonDeltaR >= 0.5)                      h_selectBaselineYields->Fill("pass_mindR(q/g, #gamma)",wt);
//                     if(mindr_Pho_genlep >=0.5)                      h_selectBaselineYields->Fill("pass_mindR(l, #gamma)",wt);}                          }}
//         if(hasGenPromptPhoton && (s_sample.Contains("GJets")))
//           {            if(!(madMinPhotonDeltaR>0.4)) continue;}
          
//         if(hasGenPromptPhoton && (s_sample.Contains("QCD")))
//           {            if((madMinPhotonDeltaR>0.4 && hasGenPromptPhoton)) continue;}         
//     	if(hasGenPromptPhoton && ((s_sample.Contains("ZG"))|| (s_sample.Contains("ZNuNuG")) || s_sample.Contains("ZNuNuGJets")))          {
// 	  if(!(madMinPhotonDeltaR>0.5)) continue;}
          
//         if(hasGenPromptPhoton && ((s_sample.Contains("ZJets"))|| (s_sample.Contains("ZNuNuJets"))))
//           {            if(!(madMinPhotonDeltaR<=0.5)) continue;}          
// 	if(hasGenPromptPhoton && (s_sample.Contains("DYJets")))
// 	  {
// 	    if(!(madMinPhotonDeltaR<=0.5 || mindr_Pho_genlep<=0.5) ) continue;

// 	  }
//         genphomatch_after++;
//     }
//     h_madminPhotonDeltaR_preSelection->Fill(madMinPhotonDeltaR,wt); 
// 	 }
//     h_selectBaselineYields_v1->Fill("overlap removal",wt);
// '
    FillHistogram_Kinematics(1,nHadJets,BTags,bestEMObj.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestEMObj.Eta(),bestEMObj.Phi(),bestEMObj.E(),METPhi,qmulti, leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
			    
     	 //    h_madminPhotonDeltaR_preSelection->Fill(madMinPhotonDeltaR,wt);
    if(Debug)
      cout<<" just before photon identification - prompt/non-prompt ========   ===="<<jentry<<endl;
    ////////////////////////////////////////////////////////////////////////////////////
    //// ======= ======= photon identification - prompt/non-prompt ======== ////////////
    ////// based on the following studies   ///////////////////////////////////////////
    // https://indico.cern.ch/event/1240842/contributions/5238493/attachments/2582794/4457078/ImpactOf_DiffPhotonIDs_ElectronFaking_Photon_28012023.pdf
    // removing it for now
    TLorentzVector zvec;
    zvec=(bestEMObj + tagEMObj);
    double invariantmass = (bestEMObj + tagEMObj).M();
    //                                                                                                                                                                
    h_InvariantMass_noselec->Fill(invariantmass,wt);
        if(!(invariantmass>=80 && invariantmass<=100) )continue;
	h_selectBaselineYields_v1->Fill(" inv mass cut ",wt);
    if(Debug)
      cout<<"just before 1 Electron CR for FAKE RATE ====== "<<jentry<<endl;
    ////////////////////////////////////////////////////////////////////////////////////
    //  ====  ========== ======== 1 Electron CR + 0 photon  =====================  ///
    ///////////////////////////////////////////////////////////////////////////////////
    
    bool elec_CR= false; 
    bool pho_SR= false; 
    
    if(bestEMobj && bestEMObjIsEle && !hasPho ){
      if(isoMuonTracks !=0 || isoPionTracks!=0) continue; // veto muon/pions from 1 electron CR
      if(NMuons!=0 || NElectrons==0) continue;     
      if(Debug)
	cout<<"in Electron CR ====== "<<jentry<<endl;
      mt_ele=sqrt(2*bestEMObj.Pt()*MET*(1-cos(DeltaPhi(METPhi,bestEMObj.Phi()))));
      //double mTElecMET=sqrt(2*(Electrons_v1[e_index].Pt())*MET*(1-cos(DeltaPhi(METPhi,Electrons_v1[e_index].Phi()))));
      //if(mt_ele>100) { continue;}//h_selectBaselineYields_CR->Fill("e-CR: mT<100",wt);continue;} // remove signal contamination 
      elec_CR = true;      
      int searchBin = getBinNoV6_WithOnlyBLSelec(nHadJets,BTags);
      if(nHadJets>=7 && Debug)
	cout<<"Entry: "<<jentry<<"\t"<<searchBin<<"\t"<<nHadJets<<"\t"<<BTags<<"\t"<<MET<<"\t"<<bestPhoton.Pt()<<"\t"<<wt<<endl;
      //if (bestEMObj.Pt()<=40) continue;
      h_selectBaselineYields_v1->Fill("Elec CR ",wt);
      //cout<<" Debugging CR "<<"\t"<<hasPho<<"\t"<<Electrons_v1[e_index].Pt()<<"\t"<<bestEMObj.Pt()<<"\t"<<MET<<"\t"<<nHadJets<<"\t"<<ST<<"\t"<<bestPhoton.Pt()<<"\t"<<hasPho_px<<"\t"<<isoPionTracks<<"\t"<<isoMuonTracks<<"\t"<<NMuons<<"\t"<<NElectrons<<"\t"<<isoElectronTracks<<"\t"<<(*Electrons_passIso)[e_index]<<"\t"<<mt_ele<<"\t"<<mTElecMET<<endl;
      //double wt_LL = wt*h_TF->GetBinContent(TFbins+1);
      //FillTFBins_Valid(6, nHadJets, BTags,wt_LL, wt_LL1, wt_LL2, bestPhoton.Pt(),ST);//
      // h_TFbins_ElecLL_validation[6]->Fill(TFbins+1,wt_LL);
      // h_TFbins_ElecLL_validation_v1[6]->Fill(TFbins+1,wt_LL);
      // h_Sbins_LL_Validation[6]->Fill(searchBin,wt_LL);
      //      cout<<"TFbins " <<TFbins<<"\t"<<h_TF->GetBinContent(TFbins+1)<<"\t"<<wt_LL<<endl;
      if(elec_CR) {    
	h_invariantMass_noCut[2]->Fill(invariantmass,wt);
	// if(!(invariantmass>=80 && invariantmass<=100) )continue; 
	// h_selectBaselineYields_v1->Fill("Elec CR - inv mass cut ",wt);

	FillHistogram_Kinematics(2,nHadJets,BTags,bestEMObj.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestEMObj.Eta(),bestEMObj.Phi(),bestEMObj.E(),METPhi,qmulti, leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
	FillHistogram_Kinematics_varBin(2,nHadJets, BTags, bestEMObj.Pt(),ST,qmulti,wt);
	h_tagEle_pT_Elec_CR->Fill(tagEMObj.Pt(),wt);
	h_tagEle_Eta_Elec_CR->Fill(tagEMObj.Eta(),wt);
	h_tagEle_Phi_Elec_CR->Fill(tagEMObj.Phi(),wt);
	h_trackEle_pT_Elec_CR->Fill(TAPElectronTracks_v1[track_idx].Pt(),wt);
	h_trackEle_Eta_Elec_CR->Fill(TAPElectronTracks_v1[track_idx].Eta(),wt);
	h_trackEle_Phi_Elec_CR->Fill(TAPElectronTracks_v1[track_idx].Phi(),wt);
	h_tagEle_EtaVsPhi_Elec_CR->Fill(tagEMObj.Eta(),tagEMObj.Phi(),wt);
	h_tagEle_PtVsPhi_Elec_CR->Fill(tagEMObj.Pt(),tagEMObj.Phi(),wt);
	h_tagEle_PtVsPhi_Elec_CR->Fill(tagEMObj.Pt(),tagEMObj.Eta(),wt);
	
	h_tagEle_EtaVsPhi_Elec_CR->Fill(TAPElectronTracks_v1[track_idx].Eta(),TAPElectronTracks_v1[track_idx].Phi(),wt);
	h_tagEle_PtVsPhi_Elec_CR->Fill(TAPElectronTracks_v1[track_idx].Pt(),TAPElectronTracks_v1[track_idx].Phi(),wt);
	h_tagEle_PtVsPhi_Elec_CR->Fill(TAPElectronTracks_v1[track_idx].Pt(),TAPElectronTracks_v1[track_idx].Eta(),wt);
	
	// h2_minDr_tagEle_pt[2]->Fill();
	// h_minDr_bestphovsTagEle[2]->Fill();
	if(BTags==0)
	  FR_nbtagBins[2]->Fill(1,wt);
	else if(BTags>=1)
	  FR_nbtagBins[2]->Fill(2,wt);
	h_invariantMass[2]->Fill(invariantmass,wt);
	//Should be moved before the preselections once verified
	// if(HEMaffected == true) continue;
      
	// FillHistogram_Kinematics(4,nHadJets,BTags,bestEMObj.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestEMObj.Eta(),bestEMObj.Phi(),bestEMObj.E(),METPhi,qmulti,leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
	// FillHistogram_Kinematics_varBin(4,nHadJets, BTags, bestEMObj.Pt(),ST,qmulti,wt);

	//   //applying l1 prefire issue correction
	// // if(!extracondition) {
        //   // FillHistogram_Kinematics(6,nHadJets,BTags,bestEMObj.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestEMObj.Eta(),bestEMObj.Phi(),bestEMObj.E(),METPhi,qmulti,leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
        //   // FillHistogram_Kinematics_varBin(6,nHadJets, BTags, bestEMObj.Pt(),ST,qmulti,wt);
	  
	//   //}
	// double wt1= wt;
        // if(applyL1TrigFire_prob && (s_data.Contains("2016") ||  s_data.Contains("2017") ))
        //   {
        //     wt1 =wt*NonPrefiringProb;
        //   }

	// FillHistogram_Kinematics(8,nHadJets,BTags,bestEMObj.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestEMObj.Eta(),bestEMObj.Phi(),bestEMObj.E(),METPhi,qmulti,leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt1);
	// FillHistogram_Kinematics_varBin(8,nHadJets, BTags, bestEMObj.Pt(),ST,qmulti,wt1);
	// //applying validation
	// int TFbins = getBinNo_v0FR(bestEMObj.Pt(), qmulti, 0.1);
	// double wt_LL = wt;//*h_TF->GetBinContent(TFbins+1);
	// FillTFBins_Valid(8,nHadJets,BTags,bestEMObj.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestEMObj.Eta(),bestEMObj.Phi(),bestEMObj.E(),METPhi,qmulti,leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt_LL,wt_LL,wt_LL);
	
	// if (!extracondition)
	//   {	  FillHistogram_Kinematics(6,nHadJets,BTags,bestEMObj.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestEMObj.Eta(),bestEMObj.Phi(),bestEMObj.E(),METPhi,qmulti,leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
	//     FillHistogram_Kinematics_varBin(6,nHadJets, BTags, bestEMObj.Pt(),ST,qmulti,wt);
	//   }
      
      }
    }


         ////////////////////////////////////////////////////////////////////////////////////
    //  ====  ========== ======== 1 photon SR =====================  ///
    /////////////////////////////////////////////////////////////////////////////////// 
    bool fakeElectron=false, proc=false, proc1=false;
    if(hasPho && bestEMobj && !bestEMObjIsEle ){ //&& !s_sample.Contains("data")){
      // if(isoMuonTracks !=0 || isoPionTracks!=0 ) continue; // veto muon/pions/electrons 
      // if(NMuons!=0 ) continue;
      // h_selectBaselineYields_v1->Fill("Pho SR ",wt);
      if(!s_sample.Contains("data")){
      // 	//if(nGenMu1==0 && nGenEle=1=0 && v_genTau2.size()==0) continue; //reject w->qq events or all hadronic decay of TTbar                                                    
      //  // if(nGenMu1>0) continue;
      //  // if(nGenEle1!=2) continue;
      
      // 	//	if(((s_sample.Contains("WGJets") || s_sample.Contains("WJets"))&& nGentau_had1>0) ||( (s_sample.Contains("TTGJets") || s_sample.Contains("TTJets")) && nGentau_had1>1)) continue;
	//      if(nGenEle1==0) continue;      
      proc = v_genz.size()==1;
      proc1 = v_genz.size()>1;
      fakeElectron =(proc && (bestPhoton.DeltaR(v_genz[0])<0.2 && bestPhoton.Pt()/v_genz[0].Pt()<1.2 && bestPhoton.Pt()/v_genz[0].Pt()>0.8))||(proc1 && ((bestPhoton.DeltaR(v_genz[0])<0.2 && bestPhoton.Pt()/v_genz[0].Pt()<1.2 && bestPhoton.Pt()/v_genz[0].Pt()>0.8)|| (v_genz.size()>1 && (bestPhoton.DeltaR(v_genz[1])<0.2 && bestPhoton.Pt()/v_genz[1].Pt()<1.2 && bestPhoton.Pt()/v_genz[1].Pt()>0.8))));
      fakeElectron= true;
      if(fakeElectron) pho_SR= true;
      else continue;
      h_selectBaselineYields_v1->Fill("Pho SR fake photon condition ",wt);

       }
      pho_SR= true;
      //if(bestEMObj.Pt()<=40) continue;
      //cout<<" Debugging SR "<<"\t"<<hasPho<<"\t"<<"\t"<<bestEMObj.Pt()<<"\t"<<MET<<"\t"<<nHadJets<<"\t"<<ST<<"\t"<<bestPhoton.Pt()<<"\t"<<hasPho_px<<"\t"<<isoPionTracks<<"\t"<<isoMuonTracks<<"\t"<<NMuons<<"\t"<<NElectrons<<"\t"<<isoElectronTracks<<"\t"<<mt_ele<<"\t"<<endl;

      //if(ntracks>1) continue; //check this condition
      h_selectBaselineYields_v1->Fill("Pho SR nTrack conditions ",wt);

      if(pho_SR) {
	h_invariantMass_noCut[3]->Fill(invariantmass,wt);
        // if(!(invariantmass>=80 && invariantmass<=100) )continue;

	// h_selectBaselineYields_v1->Fill("Pho SR inv mass cut",wt);

	FillHistogram_Kinematics(3,nHadJets,BTags,bestEMObj.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestEMObj.Eta(),bestEMObj.Phi(),bestEMObj.E(),METPhi,qmulti, leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
	FillHistogram_Kinematics_varBin(3,nHadJets, BTags, bestEMObj.Pt(),ST,qmulti,wt);

	h_tagEle_pT_Pho_SR->Fill(tagEMObj.Pt(),wt);
        h_tagEle_Eta_Pho_SR->Fill(tagEMObj.Eta(),wt);
        h_tagEle_Phi_Pho_SR->Fill(tagEMObj.Phi(),wt);
        h_trackEle_pT_Pho_SR->Fill(TAPElectronTracks_v1[track_idx].Pt(),wt);
	h_trackEle_Eta_Pho_SR->Fill(TAPElectronTracks_v1[track_idx].Eta(),wt);
        h_trackEle_Phi_Pho_SR->Fill(TAPElectronTracks_v1[track_idx].Phi(),wt);
        h_tagEle_EtaVsPhi_Pho_SR->Fill(tagEMObj.Eta(),tagEMObj.Phi(),wt);
	h_tagEle_PtVsPhi_Pho_SR->Fill(tagEMObj.Pt(),tagEMObj.Phi(),wt);
        h_tagEle_PtVsPhi_Pho_SR->Fill(tagEMObj.Pt(),tagEMObj.Eta(),wt);

        h_tagEle_EtaVsPhi_Pho_SR->Fill(TAPElectronTracks_v1[track_idx].Eta(),TAPElectronTracks_v1[track_idx].Phi(),wt);
        h_tagEle_PtVsPhi_Pho_SR->Fill(TAPElectronTracks_v1[track_idx].Pt(),TAPElectronTracks_v1[track_idx].Phi(),wt);
        h_tagEle_PtVsPhi_Pho_SR->Fill(TAPElectronTracks_v1[track_idx].Pt(),TAPElectronTracks_v1[track_idx].Eta(),wt);
 
        if(BTags==0)
          FR_nbtagBins[3]->Fill(1,wt);
        else if(BTags>=1)
          FR_nbtagBins[3]->Fill(2,wt);
        h_invariantMass[3]->Fill(invariantmass,wt);

        // // //Should be moved before the preselections once verified                                                                                                     
	// if(HEMaffected == true) continue;
      
	// FillHistogram_Kinematics(5,nHadJets,BTags,bestEMObj.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestEMObj.Eta(),bestEMObj.Phi(),bestEMObj.E(),METPhi,qmulti,leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
	// FillHistogram_Kinematics_varBin(5,nHadJets, BTags, bestEMObj.Pt(),ST,qmulti,wt);

	
	// // if(!L1prefire_issue) {
	// //   FillHistogram_Kinematics(7,nHadJets,BTags,bestEMObj.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestEMObj.Eta(),bestEMObj.Phi(),bestEMObj.E(),METPhi,qmulti,leadjet_qmulti , leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
	// //   FillHistogram_Kinematics_varBin(7,nHadJets, BTags, bestEMObj.Pt(),ST,qmulti,wt);
	// // }
	// double wt1= wt;
	// if(applyL1TrigFire_prob && (s_data.Contains("2016") ||  s_data.Contains("2017") ))
	//   {
	//     wt1 =wt*NonPrefiringProb;
	//   }
	// FillHistogram_Kinematics(9,nHadJets,BTags,bestEMObj.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestEMObj.Eta(),bestEMObj.Phi(),bestEMObj.E(),METPhi,qmulti,leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt1);
	// FillHistogram_Kinematics_varBin(9,nHadJets, BTags, bestEMObj.Pt(),ST,qmulti,wt1);		

	// if (!extracondition)
        //   {       FillHistogram_Kinematics(7,nHadJets,BTags,bestEMObj.Pt(),mTPhoMET,dPhi_PhoMET,ST,bestEMObj.Eta(),bestEMObj.Phi(),bestEMObj.E(),METPhi,qmulti,leadjet_qmulti, leadjet_Pt,leadbjet_tag,minDR,Jet_matched, hadJets, hadJets[0], NVtx,mindr_Pho_genlep,wt);
        //     FillHistogram_Kinematics_varBin(7,nHadJets, BTags, bestEMObj.Pt(),ST,qmulti,wt);
        //   }


      }
    }
			    
     
     
    nsurVived+=wt;
     }			     //loop over entries

     
  if(Debug)
    cout<<"filling the branches in tree"<<endl;
   //  cout<<nocut<<"\t"<<nSurvived<<"\t"<<bkg_comp<<endl;
  cout<<"Alpana-check"<<"\t"<<"events not falling in any LL CR/SR"<<"\t"<<counter<<endl;
  cout<<"Survived preselection ===: "<<"\t"<<nsurVived<<"\t"<< count_nelec<<endl;
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
int AnalyzeLightBSM::getBinNo_v0FR(double pho_pt,double qmulti, double minDRindx){
  int sBin=0,m_i=1,sBin1=0,n_i=0; 
  for(int i=0;i<BestPhotonPtBinLowEdge.size()-1;i++){
    if(BestPhotonPtBinLowEdge[i]<99.99) continue; 
    //if(qmulti>=10) continue;
    if(i!=0)
      m_i++;  
    if(pho_pt >= BestPhotonPtBinLowEdge[i] && pho_pt < BestPhotonPtBinLowEdge[i+1])
      {
	sBin = sBin+((m_i-1)*4);

	break;

      }
    else if(pho_pt >= BestPhotonPtBinLowEdge[BestPhotonPtBinLowEdge.size()-1])
      {
	sBin = 40;
	break;
      }
  }

  if(sBin%4==0)// || sBin==4 || sBin==8 || sBin==12 || sBin==16 || sBin==20 )
    {
      //cout<<"bestEMObj : "<<pho_pt<<" , Qmult : "<<qmulti<<" , sBin = "<<sBin1<<endl;

      for(int i=0;i<QMultLowedge.size()-1;i++){
	n_i++;
	if(qmulti>=QMultLowedge[i] && qmulti<QMultLowedge[i+1]) {sBin1=sBin+n_i; break;}
	else if(qmulti>=QMultLowedge[QMultLowedge.size()-1]){sBin1=sBin1+(QMultLowedge.size()-1); break;}
	//else sBin1=-999;
      }
    }
  
  //   if((qmulti>=7 && qmulti<10) && pho_pt>=100 && pho_pt <120)
  // if(pho_pt>=200 && pho_pt <230)
  //   cout<<"bestEMObj : "<<pho_pt<<" , Qmult : "<<qmulti<<" , sBin = "<<sBin1<<endl;
  return sBin1;  
}
int AnalyzeLightBSM::getBinNo_v1FR(double pho_pt, int nHadJets){
  int sBin=0,m_i=1,sBin1=0,n_i=0; 
  for(int i=0;i<BestPhotonPtBinLowEdge.size()-1;i++){
    if(BestPhotonPtBinLowEdge[i]<99.99) continue; 
    if (i!=0)
      m_i++;  
    if(pho_pt >= BestPhotonPtBinLowEdge[i] && pho_pt < BestPhotonPtBinLowEdge[i+1])
      {
	sBin = sBin+((m_i-1)*3);

	break;

      }
    else if(pho_pt >= BestPhotonPtBinLowEdge[BestPhotonPtBinLowEdge.size()-1])
      {
	sBin = 30;
	break;
      }
  }

  if(sBin%3==0)// || sBin==4 || sBin==8 || sBin==12 || sBin==16 || sBin==20 )
    {
      for(int i=0;i<nJetsLowedge.size()-1;i++){
	n_i++;
	if(nHadJets>=nJetsLowedge[i] && nHadJets<nJetsLowedge[i+1]) {sBin1=sBin+n_i; break;}
	else if(nHadJets>=nJetsLowedge[nJetsLowedge.size()-1]){sBin1=sBin+(nJetsLowedge.size()-1); break;}
	//else sBin1=-999;
      }
    }
  
  //   if((nHadJets>=7 && nHadJets<10) && pho_pt>=100 && pho_pt <120)
  // if(pho_pt>=200 && pho_pt <230)
  //   cout<<"bestEMObj : "<<pho_pt<<" , Qmult : "<<nHadJets<<" , sBin = "<<sBin1<<endl;
  return sBin1;  
}

int AnalyzeLightBSM::getBinNo_v2FR(double qmulti, int nHadJets, int nbjets){
  int sBin=-100,m_i=0;
  if(nbjets==0){
    if(nHadJets>=2 && nHadJets<=4)      { sBin=1;}
    else if((nHadJets==5 || nHadJets==6)){ sBin=5;}
    else if(nHadJets>=7)   { sBin=9;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)      { sBin=13;}
    else if((nHadJets==5 || nHadJets==6)){ sBin=17;}
    else if(nHadJets>=7)   { sBin=21;}
  }
  
  for(int i=0;i<QMultLowedge.size()-1;i++){
    m_i++;
    if(qmulti>=QMultLowedge[i] && qmulti<QMultLowedge[i+1]) {sBin=sBin+m_i; break;}
    else if(qmulti>=QMultLowedge[QMultLowedge.size()-1]){sBin=sBin+(QMultLowedge.size()-1); break;}
    //else sBin1=-999;                                                                                                                                            
  }
  return sBin;
}

int AnalyzeLightBSM::getBinNoV6_WithOnlyBLSelec(int nHadJets,int nbjets)
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
    for(int i=0;i<METLowEdge_v3.size()-1;i++){
      if(METLowEdge_v3[i]<199.99) continue;
      int sBin1=sBin;
      m_i++;
      if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
	break; }
      else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 7         ;
        break; }
    }
  }
  else if(sBin==7 || sBin==33 || sBin==39){
    int sBin1=sBin;
    for(int i=0;i<METLowEdge_v3_1.size()-1;i++){
      if(METLowEdge_v3_1[i]<199.99) continue;
      m_i++;
      if(MET >= METLowEdge_v3_1[i] && MET < METLowEdge_v3_1[i+1]){ sBin = sBin+m_i;break;}
      else if(MET >= METLowEdge_v3_1[METLowEdge_v3_1.size()-1])  { sBin = sBin+6; break; }
    }
  }

  else 
    {
      for(int i=0;i<METLowEdge_v3_2.size()-1;i++){
	if(METLowEdge_v3_2[i]<199.99) continue;
	m_i++;
	if(MET >= METLowEdge_v3_2[i] && MET < METLowEdge_v3_2[i+1]){ sBin = sBin+m_i;break; }
	else if(MET >= METLowEdge_v3_2[METLowEdge_v3_2.size()-1])  { sBin = sBin+5; break; }
      }
    }
  // if(sBin==0){
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3[i]<199.99) continue;
  //     int sBin1=sBin;
  //     m_i++;
  //     if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
  //       break; }
  //     else if(MET >= METLowEdge_v3[METLowEdge_v3.size()])  { sBin = 6;
  //       break; }
  //   }
  //   else if (sBin==7){
  //     for(int i=0;i<METLowEdge_v3.size()-1;i++){
  // 	if(METLowEdge_v3[i]<199.99) continue;
  // 	int sBin1=sBin;
  // 	m_i++;
  // 	if(MET >= METLowEdge_v3[i] && MET < METLowEdge_v3[i+1]){ sBin = sBin+m_i;
  // 	  break; }
  // 	else if(MET >= METLowEdge_v3[METLowEdge_v3.size()-1])  { sBin = 6;
  // 	  break; }

  //   }
    
  // }
  // else if(sBin==7 || sBin==13 || sBin==19 || sBin==25 || sBin==31){
  //   int sBin1=sBin;
  //   for(int i=0;i<METLowEdge_v3.size()-1;i++){
  //     if(METLowEdge_v3_1[i]<199.99) continue;
  //     m_i++;
  //     if(MET >= METLowEdge_v3_1[i] && MET < METLowEdge_v3_1[i+1]){ sBin = sBin+m_i;
  //       break;}
  //     else if(MET >= METLowEdge_v3_1[METLowEdge_v3_1.size()-1])  { sBin = sBin+6;
  //       break; }
  //   }
  // }

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
int AnalyzeLightBSM::getBinNoV7_le(int nHadJets, int nbjets){
  int sBin=-100,m_i=0;
  if(nbjets==0){
    if(nHadJets==2) { if(MET<300)sBin=1; else if (MET>300) sBin=2;}
    else if(nHadJets==3)     { if(MET<300)sBin=3; else if (MET>300) sBin=4;}
    else if(nHadJets==4)     { if(MET<300)sBin=5; else if (MET>300) sBin=6;}
    else if((nHadJets==5 || nHadJets==6)){ if(MET<300)sBin=7; else if (MET>300) sBin=8;}
    else if(nHadJets>=7)   { if(MET<300)sBin=9; else if (MET>300) sBin=10;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)      {if(MET<300)sBin=11; else if (MET>300) sBin=12;}
    else if((nHadJets==5 || nHadJets==6)){ if(MET<300)sBin=13; else if (MET>300) sBin=14;}
    else if(nHadJets>=7)   { if(MET<300)sBin=15; else if (MET>300) sBin=16;}
  }
  return sBin;
}

int AnalyzeLightBSM::getBinNoV16_le(int nHadJets, int nbjets, double photon_pT){
  int sBin=-100,m_i=0;
  if(nbjets==0){
    if(nHadJets==2) { if(photon_pT<100)sBin=1; else if (photon_pT>100) sBin=2;}
    else if(nHadJets==3)     { if(photon_pT<100)sBin=3; else if (photon_pT>100) sBin=4;}
    else if(nHadJets==4)     { if(photon_pT<100)sBin=5; else if (photon_pT>100) sBin=6;}
    else if((nHadJets==5 || nHadJets==6)){ if(photon_pT<100)sBin=7; else if (photon_pT>100) sBin=8;}
    else if(nHadJets>=7)   { if(photon_pT<100)sBin=9; else if (photon_pT>100) sBin=10;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)      {if(photon_pT<100)sBin=11; else if (photon_pT>100) sBin=12;}
    else if((nHadJets==5 || nHadJets==6)){ if(photon_pT<100)sBin=13; else if (photon_pT>100) sBin=14;}
    else if(nHadJets>=7)   { if(photon_pT<100)sBin=15; else if (photon_pT>100) sBin=16;}
  }
  return sBin;
}


int AnalyzeLightBSM::getBinNoV1_le( int nHadJets, int nbjets){
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
    else if(nHadJets==5 || nHadJets==6){ sBin=4;}
    else if(nHadJets>=7)               { sBin=8;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)     { sBin=12;}
    else if(nHadJets>=5){ sBin=16;}
    //    else if(nHadJets>=7)               { sBin=26;}
  }
  //  if(sBin==-99)//
  //  cout<<sBin<<"\t"<<nHadJets<<"\t"<<btags<<endl;
  if(sBin==0){
    for(int i=0;i<METLowEdge_lowMET.size()-1;i++){
      if(METLowEdge_lowMET[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge_lowMET[i] && MET < METLowEdge_lowMET[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge_lowMET[METLowEdge_lowMET.size()-1])  { sBin = 4         ;break; }
    }
  }
  else{
    for(int i=0;i<METLowEdge_lowMET.size()-1;i++){
      if(METLowEdge_lowMET[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge_lowMET[i] && MET < METLowEdge_lowMET[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge_lowMET[METLowEdge_lowMET.size()-1])  { sBin = sBin+4   ;break; }
    }
  }
  return sBin;
}
int AnalyzeLightBSM::getBinNoV7_highMET(int nHadJets, int btags){
  int sBin=-100,m_i=0;
  if(btags==0){
    if(nHadJets>=2 && nHadJets<=4)     { sBin=0;}
    else if(nHadJets==5 || nHadJets==6){ sBin=4;}
    else if(nHadJets>=7)               { sBin=8;}
  }
  else{
    if(nHadJets>=2 && nHadJets<=4)     { sBin=12;}
    else if(nHadJets>=5){ sBin=16;}
    //    else if(nHadJets>=7)               { sBin=26;}                                                                                                                     
  }
  if(sBin==0){
    for(int i=0;i<METLowEdge_highMET.size()-1;i++){
      if(METLowEdge_highMET[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge_highMET[i] && MET < METLowEdge_highMET[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge_highMET[METLowEdge_highMET.size()-1])  { sBin = 4         ;break; }
    }
  }
  else{
    for(int i=0;i<METLowEdge_highMET.size()-1;i++){
      if(METLowEdge_highMET[i]<99.99) continue;
      m_i++;
      if(MET >= METLowEdge_highMET[i] && MET < METLowEdge_highMET[i+1]){ sBin = sBin+m_i;break; }
      else if(MET >= METLowEdge_highMET[METLowEdge_highMET.size()-1])  { sBin = sBin+4   ;break; }
    }
  }
  return sBin;
}

TLorentzVector AnalyzeLightBSM::getBestPhoton(int pho_ID){
  vector<TLorentzVector> goodPho;
  vector<int> goodPhoIndx;
  //  cout<<"\t"<<(*Photons_mvaValuesID).size()<<endl;
  for(int iPho=0;iPho<Photons_;iPho++){
    //if(((*Photons_hasPixelSeed)[iPho]<0.001) && ( (*Photons_fullID)[iPho]))
    //    cout<<(*Photons_hasPixelSeed)[iPho] <<"\t (*Photons_hasPixelSeed)[iPho]  "<<" (*Photons_fullID)[iPho]"<<"\t"<<(*Photons_fullID)[iPho]<<"\t"<<Photons_v1[iPho].Eta()<<endl;
    if((*Photons_fullID)[iPho])//abs(Photons_v1[iPho].Eta())<2.4 && ((*Photons_fullID)[iPho] && pho_ID==0 && (*Photons_hasPixelSeed)[iPho]<0.001))// && ( ((*Photons_fullID)[iPho] && pho_ID==0) || (pho_ID==1 &&(((*Photons_cutBasedID)[iPho]==1 || (*Photons_cutBasedID)[iPho]==2))) || (pho_ID==2 && (*Photons_cutBasedID)[iPho]==2) || (pho_ID==3 && (*Photons_mvaValuesID)[iPho]>-0.02) || (pho_ID==4 && (*Photons_mvaValuesID)[iPho]>0.42)))
      {
      goodPho.push_back(Photons_v1[iPho] );
      goodPhoIndx.push_back(iPho);
      //      cout<<(*Photons_mvaValuesID)[iPho] << " (*Photons_mvaValuesID)[iPho]"<<endl;
    }
  }
  //  cout<<(*Photons_hasPixelSeed)[iPho] <<"\t(*Photons_hasPixelSeed)[iPho]  "<<endl;
  //cout<<" goodPho.size()" <<"\t"<<  goodPho.size()<<endl;
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
vector <TLorentzVector> AnalyzeLightBSM::getLorentzVector(int size, Float_t Pt_size[],Float_t Eta_size[],Float_t Phi_size[],Float_t E_size[])
{
  //  cout<<"filling the vector of lorentz vector"<< "\t"<< "vector <TLorentzVector> AnalyzeLightBSM::getLorentzVector"<<endl;
  vector<TLorentzVector> Temp;
  for(int i=0; i<size;i++)
    {
      TLorentzVector P1;
      P1.SetPtEtaPhiE(Pt_size[i],Eta_size[i],Phi_size[i],E_size[i]);
      //cout<< "Pt "<<Pt_size[i]<<" Eta "<<Eta_size[i]<< " Phi "<<Phi_size[i]<< "E "<<E_size[i]<<endl;
      Temp.push_back(P1);
    }
  //cout<< "Out lorentz vectors size "<<Temp.size() <<endl;
  return Temp;
}

void AnalyzeLightBSM::FillHistogram_Kinematics(int i, int Njets, int btags, double pho_Pt, double mt_phoMET, double dPhi, double ST,float Eta , float Phi, double E, float METPhi,double qmulti,double leadJets_qmulti, double leadJet_Pt, double leadbjet_tag, double mindR, TLorentzVector Jet_matched, vector<TLorentzVector>hadJets, TLorentzVector leadJet,int nvrtx,double mindr_Pho_genElec, double wt){
  //cout<<"Alps "<<i<<"\t"<<Njets<<"\t"<<btags<<"\t"<<pho_Pt<<"\t"<<mt_phoMET<<"\t"<<dPhi<<"\t"<<ST<<endl;
  
  h_Njets[i]->Fill(Njets,wt);
  h_Nbjets[i]->Fill(btags,wt);
  h_MET_[i]->Fill(MET,wt);
  h_PhotonPt[i]->Fill(pho_Pt,wt);
  h_Mt_PhoMET[i]->Fill(mt_phoMET,wt);
  h_dPhi_PhoMET[i]->Fill(dPhi,wt);
  h_St[i]->Fill(ST,wt);
  h_HT[i]->Fill(HT,wt);
  h_Photon_Eta[i]->Fill(Eta,wt);
  h_Photon_Phi[i]->Fill(Phi,wt);
  
  h_Photon_E[i]->Fill(E,wt);
  h_MET_Phi[i]->Fill(METPhi,wt);
  h_qmulti_1[i]->Fill(qmulti,wt);
  h_qmultiVsEmobjPT[i]->Fill(pho_Pt,qmulti,wt);
  h_qmultiVsnJets[i]->Fill(Njets,qmulti,wt);
  h_ST_vs_EMObjPt[i]->Fill( pho_Pt,ST,wt);
  h_Emobj_PtvsEta[i]->Fill(pho_Pt,Eta,wt);
  h_Emobj_PtvsPhi[i]->Fill(pho_Pt,Phi,wt);
  h_Emobj_EtavsPhi[i]->Fill(Eta,Phi,wt);
  //  cout<<"just after eta vs phi  "<<endl;
  h_nvrtx[i]->Fill(nvrtx,wt);
  h_minDR_Jets_EMObject[i]->Fill(mindR,wt);
  h_Phi_matchedJet[i]->Fill(Jet_matched.Phi(),wt);
  h_Pt_matchedJet[i]->Fill(Jet_matched.Pt(),wt);
  h_Eta_matchedJet[i]->Fill(Jet_matched.Eta(),wt);
  h_EtavsPhi_matchedJet[i]->Fill(Jet_matched.Eta(),Jet_matched.Phi(),wt);
  h_PtvsEta_matchedJet[i]->Fill(Jet_matched.Pt(),Jet_matched.Eta(),wt);
  h_PtvsPhi_matchedJet[i]->Fill(Jet_matched.Pt(),Jet_matched.Phi(),wt);
  h_minDR_Jets_vs_Em_Pt[i]->Fill(Jet_matched.Pt(),mindR,wt);
  h_btaggervalue_vs_qmulti[i]->Fill(qmulti,leadbjet_tag,wt);
  h_btaggervalue_vs_minDR_Jets[i]->Fill(mindR,leadbjet_tag,wt);
  h_minDR_Jets_vsqMulti[i]->Fill(qmulti,mindR,wt);
  h_HT5HT[i]->Fill(HT5/HT,wt);
  h_Emobje_pt_vs_Jet_Pt[i]->Fill(pho_Pt,Jet_matched.Pt(),wt);
  for(int j=0; j<hadJets.size();j++){
    if (j>=4) continue;
    h_Phi_leadJet[j][i]->Fill(hadJets[j].Phi(),wt);
    h_Pt_leadJet[j][i]->Fill(hadJets[j].Pt(),wt);
    h_Eta_leadJet[j][i]->Fill(hadJets[j].Eta(),wt);
    h_EtavsPhi_leadJet[j][i]->Fill(hadJets[j].Eta(),hadJets[j].Phi(),wt);
    h_PtvsPhi_leadJet[j][i]->Fill(hadJets[j].Pt(),hadJets[j].Phi(),wt);
    h_PtvsEta_leadJet[j][i]->Fill(hadJets[j].Pt(),hadJets[j].Eta(),wt);
    TLorentzVector Met;
    Met.SetPtEtaPhiE(MET,0,METPhi,0);
    double dPhi_METjet = abs(Met.DeltaPhi(hadJets[j])); 
    h_HT5HT_vsdPhi_METJet[j][i]->Fill(dPhi_METjet,HT5/HT,wt);
    h_dPhi_METJet[j][i]->Fill(dPhi_METjet,wt);
    h_dPhi_METJet_vsMET[j][i]->Fill(MET,dPhi_METjet,wt);
  }

    
  h_nBjets_vs_qmulti[i]->Fill(btags,qmulti,wt);
  h_qmultiVs_MET[i]->Fill(MET,qmulti,wt);
  h_qmultiVs_ST[i]->Fill(ST,qmulti,wt);
  h_leadJets_qmulti[i]->Fill(leadJets_qmulti,wt);
  h_leadJet_Pt[i]->Fill(leadJet_Pt,wt);
  h_leadbjet_tag[i]->Fill(leadbjet_tag,wt);
  h_leadbjet_tag_vs_leadQmulti[i]->Fill(leadbjet_tag,leadJets_qmulti,wt);
  h_leadbjet_tag_vsQmulti[i]->Fill(leadbjet_tag,qmulti,wt);
  h_leadjet_ptvsqmulti[i]->Fill(leadJet_Pt,qmulti,wt);

  int TFbins_v2 = getBinNoV1_le(Njets,btags);
  int TFbins_v3 = getBinNoV16_le(Njets,btags,pho_Pt);
  int searchBin = getBinNoV6_WithOnlyBLSelec(Njets,btags);
  int TFbins_v4 = getBinNoV7_le(Njets,btags);
  int TFbins_v5 = getBinNo_v0FR(pho_Pt, qmulti, 0.1);
  int TFbins_v6 = getBinNo_v1FR(pho_Pt,Njets);
  int TFbins_v7 = getBinNo_v2FR(qmulti,Njets,btags);

  h_TFbins_LL_v4[i]->Fill(TFbins_v4,wt);
  h_TFbins_LL_v2[i]->Fill(TFbins_v2,wt);
  
  h_Sbins_LL[i]->Fill(searchBin,wt);
  h_TFbins_LL_v3[i]->Fill(TFbins_v3,wt);
  h_TFbins_LL_v5[i]->Fill(TFbins_v5,wt);
  h_TFbins_LL_v6[i]->Fill(TFbins_v6,wt);
  h_TFbins_LL_v7[i]->Fill(TFbins_v7,wt);
  

  
  //cout<<"Alps "<<i<<"\t"<<"TFbins_v1  "<<TFbins_v1<<" TFbins_v2 "<< TFbins_v2<<"\t MET "<<MET<<"\t NJets  "<<Njets<<"  btags "<<btags<<"\t"<<h_TFbins_LL_v1[i]->GetBinContent(TFbins_v1)<<"\t"<<wt<<"\t"<<searchBin<<endl;
  //cout<<h_Njets[i]->GetMean()<<endl;

}
void AnalyzeLightBSM::FillTFBins_Valid(int i, int Njets, int btags, double pho_Pt, double mt_phoMET, double dPhi, double ST,float Eta , float Phi, double E, float METPhi,double qmulti,double leadJets_qmulti, double leadJet_Pt, double leadbjet_tag, double mindR, TLorentzVector Jet_matched, vector<TLorentzVector>hadJets, TLorentzVector leadJet,int nvrtx,double mindr_Pho_genElec, double wt,double wt1, double wt2){
  //cout<<"enetring the validation loopjust after EM phi  "<<endl;

  h_Mt_PhoMET_validation[i]->Fill(mt_phoMET,wt);
  h_dPhi_PhoMET_validation[i]->Fill(dPhi,wt);
  h_St_validation[i]->Fill(ST,wt);
  h_Photon_Eta_validation[i]->Fill(Eta,wt);
  h_Photon_Phi_validation[i]->Fill(Phi,wt);
  //  cout<<"just after EM phi  "<<endl;

  h_Photon_E_validation[i]->Fill(E,wt);
  h_MET_Phi_validation[i]->Fill(METPhi,wt);
  h_qmulti_1_validation[i]->Fill(qmulti,wt);
  h_nvrtx_validation[i]->Fill(nvrtx,wt);
  //cout<<"just after n vertex "<<endl;

  h_minDR_Jets_EMObject_validation[i]->Fill(mindR,wt);
  h_Phi_matchedJet_validation[i]->Fill(Jet_matched.Phi(),wt);
  h_Pt_matchedJet_validation[i]->Fill(Jet_matched.Pt(),wt);
  h_Eta_matchedJet_validation[i]->Fill(Jet_matched.Eta(),wt);
  //  cout<<"just after matched jet  "<<endl;

  h_HT5HT_validation[i]->Fill(HT5/HT,wt);
  for(int j=0; j<hadJets.size();j++){
    if (j>=4) continue;
    h_Phi_leadJet_validation[j][i]->Fill(hadJets[j].Phi(),wt);
    h_Pt_leadJet_validation[j][i]->Fill(hadJets[j].Pt(),wt);
    h_Eta_leadJet_validation[j][i]->Fill(hadJets[j].Eta(),wt);
    TLorentzVector Met;
    Met.SetPtEtaPhiE(MET,0,METPhi,0);
    double dPhi_METjet = abs(Met.DeltaPhi(hadJets[j]));
    h_dPhi_METJet_validation[j][i]->Fill(dPhi_METjet,wt);
  }
  //cout<<"just after jet loop  "<<endl;

  h_leadJets_qmulti_validation[i]->Fill(leadJets_qmulti,wt);
  h_leadJet_Pt_validation[i]->Fill(leadJet_Pt,wt);
  h_leadbjet_tag_validation[i]->Fill(leadbjet_tag,wt);
  //cout<<"just after lead jet  "<<endl;

  //old ones
  int TFbins_v2 = getBinNo_v0FR(pho_Pt, qmulti, 0.1);
  //int TFbins_v2 = getBinNoV1_le(Njets,btags);
  int TFbins_v3 = getBinNoV16_le(Njets,btags,pho_Pt);
  int TFbins_v4 = getBinNoV7_le(Njets,btags);
  int searchBin = getBinNoV6_WithOnlyBLSelec(Njets,btags);
  h_TFbins_ElecLL_validation[i]->Fill(TFbins_v2,wt);
  h_TFbins_ElecLL_validation_v1[i]->Fill(TFbins_v2,wt);
  h_Sbins_LL_Validation[i]->Fill(searchBin,wt);
  h_Njets_validation[i]->Fill(Njets,wt);
  h_Nbjets_validation[i]->Fill(btags,wt);
  h_MET_validation[i]->Fill(MET,wt);
  h_PhotonPt_validation[i]->Fill(pho_Pt,wt);
  h_St_validation[i]->Fill(ST,wt);

  h_TFbins_ElecLL_validation_TFbins_v2[i]->Fill(TFbins_v3,wt1);
  h_TFbins_ElecLL_validation_TFbins_v2_v1[i]->Fill(TFbins_v3,wt1);
  h_Sbins_LL_Validation_TFbins_V2[i]->Fill(searchBin,wt1);
  h_Njets_validation_TFbins_v2[i]->Fill(Njets,wt1);
  h_Nbjets_validation_TFbins_v2[i]->Fill(btags,wt1);
  h_MET_validation_TFbins_v2[i]->Fill(MET,wt1);
  h_PhotonPt_validation_TFbins_v2[i]->Fill(pho_Pt,wt1);
  h_St_validation_TFbins_v2[i]->Fill(ST,wt1);

  h_TFbins_ElecLL_validation_TFbins_v3[i]->Fill(TFbins_v4,wt2);
  h_TFbins_ElecLL_validation_TFbins_v3_v1[i]->Fill(TFbins_v4,wt2);
  h_Sbins_LL_Validation_TFbins_V3[i]->Fill(searchBin,wt2);
  h_Njets_validation_TFbins_v3[i]->Fill(Njets,wt2);
  h_Nbjets_validation_TFbins_v3[i]->Fill(btags,wt2);
  h_MET_validation_TFbins_v3[i]->Fill(MET,wt2);
  h_PhotonPt_validation_TFbins_v3[i]->Fill(pho_Pt,wt2);
  h_St_validation_TFbins_v3[i]->Fill(ST,wt2);

  
}
void AnalyzeLightBSM::FillHistogram_Kinematics_varBin(int i, int Njets, int btags, double pho_Pt, double ST,double qmulti, double wt){
  //cout<<"Alps "<<i<<"\t"<<Njets<<"\t"<<btags<<"\t"<<pho_Pt<<"\t"<<mt_phoMET<<"\t"<<dPhi<<"\t"<<ST<<endl;                                        
  h_Njets_Varbin[i]->Fill(Njets,wt);
  h_Nbjets_Varbin[i]->Fill(btags,wt);
  h_MET_Varbin[i]->Fill(MET,wt);
  h_PhotonPt_Varbin[i]->Fill(pho_Pt,wt);
  h_St_Varbin[i]->Fill(ST,wt);
  h_qmulti_Varbin[i]->Fill(qmulti,wt);
    h_qmultiVs_phopT_Varbin[i]->Fill(pho_Pt,qmulti,wt);
  //h_HT_Varbin[i]->Fill(HT,wt);
  //cout<<h_Njets[i]->GetMean()<<endl;                                                                                                            

    
}
int AnalyzeLightBSM::Photons_OriginType(){
  int flag_photon = -1;
  if(!(*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons) && !(*Photons_electronFakes).at(bestPhotonIndxAmongPhotons) && hasGenPromptPhoton)
    flag_photon=0;
  else if( (*Photons_nonPrompt).at(bestPhotonIndxAmongPhotons)==1)
    flag_photon=1;
  else if ((*Photons_electronFakes).at(bestPhotonIndxAmongPhotons))
    flag_photon=2;
  else
    flag_photon = -1;
  return flag_photon;
}

TLorentzVector AnalyzeLightBSM::getPhoton_withoutfullID(){
  vector<TLorentzVector> nogoodPho;
  vector<int> nogoodPhoIndx;
  for(int iPho=0;iPho<Photons_;iPho++){
    if( !(*Photons_fullID)[iPho] ) 
      {
	if(abs(Photons_v1[iPho].Eta())<2.4){
	nogoodPho.push_back( Photons_v1[iPho] );
	nogoodPhoIndx.push_back(iPho);
	}
      }
  }

  int highPtIndx=-100;
  for(int i=0;i<nogoodPho.size();i++){
    if(i==0) highPtIndx=0;
    else if( (nogoodPho[highPtIndx].Pt()) < (nogoodPho[i].Pt()) ){highPtIndx=i;}
  }
  int eIndxAmongPhotons=-100;
  if(highPtIndx>=0){
    eIndxAmongPhotons = nogoodPhoIndx[highPtIndx];
  }
  else eIndxAmongPhotons = -100;

  if(highPtIndx==-100){TLorentzVector v0;return v0;}
  else return nogoodPho[highPtIndx];
}

