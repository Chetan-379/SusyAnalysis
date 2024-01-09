/////////////////////////////////////////////////////////////////////////
//   This class has been automatically generated 
//   (at Sat Nov  4 05:03:17 2023 by ROOT version 6.18/04)
//   from TTree PreSelection/PreSelection
/////////////////////////////////////////////////////////////////////////

#ifndef NtupleVarsTProxy_h
#define NtupleVarsTProxy_h

#define R__BRANCHPROXY_GENERATOR_VERSION 2

// ROOT headers needed by the proxy
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TPad.h>
#include <TH1.h>
#include <TSelector.h>
#include <TBranchProxy.h>
#include <TBranchProxyDirector.h>
#include <TBranchProxyTemplate.h>
#include <TFriendProxy.h>
using namespace ROOT::Internal;
using ROOT::Detail::TBranchProxy;

// Header needed by this particular proxy
#include "Math/GenVector/LorentzVector.h"
#include <vector>
#include "Math/GenVector/PtEtaPhiE4D.h"
#include <TLorentzVector.h>
#include "TMath.h"

using namespace std;
class NtupleVarsTProxy : public TSelector {
public :
   TTree          *fChain;         //!pointer to the analyzed TTree or TChain
   TBranchProxyDirector fDirector; //!Manages the proxys

   // Optional User methods
   TClass         *fClass;    // Pointer to this class's description

   // Wrapper class for each unwounded class
   struct TStlPx_ROOT__Math__PtEtaPhiE4D_float_
   {
      TStlPx_ROOT__Math__PtEtaPhiE4D_float_(TBranchProxyDirector* director,const char *top,const char *mid=0) :
         ffPrefix        (top,mid),
         obj             (director, top, mid),
         fPt             (director, ffPrefix, "fPt"),
         fEta            (director, ffPrefix, "fEta"),
         fPhi            (director, ffPrefix, "fPhi"),
         fE              (director, ffPrefix, "fE")
      {};
      TStlPx_ROOT__Math__PtEtaPhiE4D_float_(TBranchProxyDirector* director, TBranchProxy *parent, const char *membername, const char *top=0, const char *mid=0) :
         ffPrefix        (top,mid),
         obj             (director, parent, membername, top, mid),
         fPt             (director, ffPrefix, "fPt"),
         fEta            (director, ffPrefix, "fEta"),
         fPhi            (director, ffPrefix, "fPhi"),
         fE              (director, ffPrefix, "fE")
      {};
      ROOT::Internal::TBranchProxyHelper ffPrefix;
      InjecTBranchProxyInterface();
      const ROOT::Math::PtEtaPhiE4D<float>& At(UInt_t i) {
         static ROOT::Math::PtEtaPhiE4D<float> default_val;
         if (!obj.Read()) return default_val;
         ROOT::Math::PtEtaPhiE4D<float> *temp = (ROOT::Math::PtEtaPhiE4D<float> *)( obj.GetProxy()->GetStlStart(i) );
         if (temp) return *temp; else return default_val;
      }
      const ROOT::Math::PtEtaPhiE4D<float>& operator[](Int_t i) { return At(i); }
      const ROOT::Math::PtEtaPhiE4D<float>& operator[](UInt_t i) { return At(i); }
      Int_t GetEntries() { return obj.GetPtr()->size(); }
      const vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > >* operator->() { return obj.GetPtr(); }
      operator vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > >*() { return obj.GetPtr(); }
      TObjProxy<vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > > > obj;

      TStlFloatProxy   fPt;
      TStlFloatProxy   fEta;
      TStlFloatProxy   fPhi;
      TStlFloatProxy   fE;
   };
   struct TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__
   {
      TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__(TBranchProxyDirector* director,const char *top,const char *mid=0) :
         ffPrefix                               (top,mid),
         obj                                    (director, top, mid),
         fCoordinates                           (director, obj.GetProxy(), "fCoordinates", ffPrefix, "fCoordinates")
      {};
      TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__(TBranchProxyDirector* director, TBranchProxy *parent, const char *membername, const char *top=0, const char *mid=0) :
         ffPrefix                               (top,mid),
         obj                                    (director, parent, membername, top, mid),
         fCoordinates                           (director, obj.GetProxy(), "fCoordinates", ffPrefix, "fCoordinates")
      {};
      ROOT::Internal::TBranchProxyHelper      ffPrefix;
      InjecTBranchProxyInterface();
      const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> >& At(UInt_t i) {
         static ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > default_val;
         if (!obj.Read()) return default_val;
         ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > *temp = & obj.GetPtr()->at(i);
         if (temp) return *temp; else return default_val;
      }
      const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> >& operator[](Int_t i) { return At(i); }
      const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> >& operator[](UInt_t i) { return At(i); }
      Int_t GetEntries() { return obj.GetPtr()->size(); }
      const vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > >* operator->() { return obj.GetPtr(); }
      operator vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > >*() { return obj.GetPtr(); }
      TObjProxy<vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > > > obj;

      TStlPx_ROOT__Math__PtEtaPhiE4D_float_   fCoordinates;
   };

   // Proxy for each of the branches, leaves and friends of the tree
   TUIntProxy                                                         RunNum;
   TUIntProxy                                                         LumiBlockNum;
   TULong64Proxy                                                      EvtNum;
   TIntProxy                                                          BadChargedCandidateFilter;
   TBoolProxy                                                         BadPFMuonDzFilter;
   TIntProxy                                                          BadPFMuonFilter;
   TIntProxy                                                          BTagsDeepCSV;
   TIntProxy                                                          BTagsDeepCSVJECdown;
   TIntProxy                                                          BTagsDeepCSVJECup;
   TIntProxy                                                          BTagsDeepCSVJERdown;
   TIntProxy                                                          BTagsDeepCSVJERup;
   TFloatProxy                                                        CaloMET;
   TFloatProxy                                                        CaloMETPhi;
   TFloatProxy                                                        CrossSection;
   TIntProxy                                                          CSCTightHaloFilter;
   TFloatProxy                                                        DeltaPhi1;
   TFloatProxy                                                        DeltaPhi1_AK8;
   TFloatProxy                                                        DeltaPhi1JECdown;
   TFloatProxy                                                        DeltaPhi1JECup;
   TFloatProxy                                                        DeltaPhi1JERdown;
   TFloatProxy                                                        DeltaPhi1JERup;
   TFloatProxy                                                        DeltaPhi2;
   TFloatProxy                                                        DeltaPhi2_AK8;
   TFloatProxy                                                        DeltaPhi2JECdown;
   TFloatProxy                                                        DeltaPhi2JECup;
   TFloatProxy                                                        DeltaPhi2JERdown;
   TFloatProxy                                                        DeltaPhi2JERup;
   TFloatProxy                                                        DeltaPhi3;
   TFloatProxy                                                        DeltaPhi3JECdown;
   TFloatProxy                                                        DeltaPhi3JECup;
   TFloatProxy                                                        DeltaPhi3JERdown;
   TFloatProxy                                                        DeltaPhi3JERup;
   TFloatProxy                                                        DeltaPhi4;
   TFloatProxy                                                        DeltaPhi4JECdown;
   TFloatProxy                                                        DeltaPhi4JECup;
   TFloatProxy                                                        DeltaPhi4JERdown;
   TFloatProxy                                                        DeltaPhi4JERup;
   TFloatProxy                                                        DeltaPhiMin_AK8;
   TIntProxy                                                          ecalBadCalibFilter;
   TIntProxy                                                          EcalDeadCellBoundaryEnergyFilter;
   TIntProxy                                                          EcalDeadCellTriggerPrimitiveFilter;
   TIntProxy                                                          eeBadScFilter;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   Electrons;
   TStlPx_ROOT__Math__PtEtaPhiE4D_float_                              fCoordinates;
   TStlSimpleProxy<vector<int> >                                      Electrons_charge;
   TStlSimpleProxy<vector<float> >                                    Electrons_iso;
   TStlSimpleProxy<vector<bool> >                                     Electrons_mediumID;
   TStlSimpleProxy<vector<float> >                                    Electrons_MTW;
   TStlSimpleProxy<vector<bool> >                                     Electrons_passIso;
   TStlSimpleProxy<vector<bool> >                                     Electrons_tightID;
   TFloatProxy                                                        fixedGridRhoFastjetAll;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   GenElectrons;
   TFloatProxy                                                        GenHT;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   GenJets;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   GenJetsAK15;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   GenJetsAK8;
   TStlSimpleProxy<vector<int> >                                      GenJetsAK8_multiplicity;
   TStlSimpleProxy<vector<float> >                                    GenJetsAK8_softDropMass;
   TFloatProxy                                                        GenMET;
   TFloatProxy                                                        GenMETPhi;
   TFloatProxy                                                        GenMHT;
   TFloatProxy                                                        GenMHTPhi;
   TFloatProxy                                                        GenMT2_AK8;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   GenMuons;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   GenParticles;
   TStlSimpleProxy<vector<int> >                                      GenParticles_Charge;
   TStlSimpleProxy<vector<int> >                                      GenParticles_ParentId;
   TStlSimpleProxy<vector<int> >                                      GenParticles_ParentIdx;
   TStlSimpleProxy<vector<int> >                                      GenParticles_PdgId;
   TStlSimpleProxy<vector<int> >                                      GenParticles_Status;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   GenTaus;
   TStlSimpleProxy<vector<bool> >                                     GenTaus_had;
   TIntProxy                                                          globalSuperTightHalo2016Filter;
   TIntProxy                                                          globalTightHalo2016Filter;
   TBoolProxy                                                         hasGenPromptPhoton;
   TIntProxy                                                          HBHEIsoNoiseFilter;
   TIntProxy                                                          HBHENoiseFilter;
   TIntProxy                                                          hfNoisyHitsFilter;
   TFloatProxy                                                        HT;
   TFloatProxy                                                        HT5;
   TFloatProxy                                                        HT5JECdown;
   TFloatProxy                                                        HT5JECup;
   TFloatProxy                                                        HT5JERdown;
   TFloatProxy                                                        HT5JERup;
   TFloatProxy                                                        HTJECdown;
   TFloatProxy                                                        HTJECup;
   TFloatProxy                                                        HTJERdown;
   TFloatProxy                                                        HTJERup;
   TIntProxy                                                          isoElectronTracks;
   TIntProxy                                                          isoMuonTracks;
   TIntProxy                                                          isoPionTracks;
   TBoolProxy                                                         JetID;
   TBoolProxy                                                         JetIDAK15;
   TBoolProxy                                                         JetIDAK8;
   TBoolProxy                                                         JetIDAK8JECdown;
   TBoolProxy                                                         JetIDAK8JECup;
   TBoolProxy                                                         JetIDAK8JERdown;
   TBoolProxy                                                         JetIDAK8JERup;
   TBoolProxy                                                         JetIDJECdown;
   TBoolProxy                                                         JetIDJECup;
   TBoolProxy                                                         JetIDJERdown;
   TBoolProxy                                                         JetIDJERup;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   Jets;
   TStlSimpleProxy<vector<float> >                                    Jets_axismajor;
   TStlSimpleProxy<vector<float> >                                    Jets_axisminor;
   TStlSimpleProxy<vector<float> >                                    Jets_bDiscriminatorCSV;
   TStlSimpleProxy<vector<float> >                                    Jets_bJetTagDeepCSVBvsAll;
   TStlSimpleProxy<vector<float> >                                    Jets_bJetTagDeepCSVprobb;
   TStlSimpleProxy<vector<float> >                                    Jets_bJetTagDeepCSVprobbb;
   TStlSimpleProxy<vector<float> >                                    Jets_bJetTagDeepCSVprobc;
   TStlSimpleProxy<vector<float> >                                    Jets_bJetTagDeepCSVprobudsg;
   TStlSimpleProxy<vector<float> >                                    Jets_bJetTagDeepFlavourprobb;
   TStlSimpleProxy<vector<float> >                                    Jets_bJetTagDeepFlavourprobbb;
   TStlSimpleProxy<vector<float> >                                    Jets_bJetTagDeepFlavourprobc;
   TStlSimpleProxy<vector<float> >                                    Jets_bJetTagDeepFlavourprobg;
   TStlSimpleProxy<vector<float> >                                    Jets_bJetTagDeepFlavourproblepb;
   TStlSimpleProxy<vector<float> >                                    Jets_bJetTagDeepFlavourprobuds;
   TStlSimpleProxy<vector<float> >                                    Jets_chargedEmEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    Jets_chargedHadronEnergyFraction;
   TStlSimpleProxy<vector<int> >                                      Jets_chargedHadronMultiplicity;
   TStlSimpleProxy<vector<int> >                                      Jets_chargedMultiplicity;
   TStlSimpleProxy<vector<float> >                                    Jets_electronEnergyFraction;
   TStlSimpleProxy<vector<int> >                                      Jets_electronMultiplicity;
   TStlSimpleProxy<vector<int> >                                      Jets_hadronFlavor;
   TStlSimpleProxy<vector<float> >                                    Jets_hfEMEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    Jets_hfHadronEnergyFraction;
   TStlSimpleProxy<vector<bool> >                                     Jets_HTMask;
   TStlSimpleProxy<vector<bool> >                                     Jets_ID;
   TStlSimpleProxy<vector<float> >                                    Jets_jecFactor;
   TStlSimpleProxy<vector<float> >                                    Jets_jecUnc;
   TStlSimpleProxy<vector<float> >                                    Jets_jerFactor;
   TStlSimpleProxy<vector<float> >                                    Jets_jerFactorDown;
   TStlSimpleProxy<vector<float> >                                    Jets_jerFactorUp;
   TStlSimpleProxy<vector<bool> >                                     Jets_LeptonMask;
   TStlSimpleProxy<vector<bool> >                                     Jets_MHTMask;
   TStlSimpleProxy<vector<int> >                                      Jets_multiplicity;
   TStlSimpleProxy<vector<float> >                                    Jets_muonEnergyFraction;
   TStlSimpleProxy<vector<int> >                                      Jets_muonMultiplicity;
   TStlSimpleProxy<vector<float> >                                    Jets_neutralEmEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    Jets_neutralHadronEnergyFraction;
   TStlSimpleProxy<vector<int> >                                      Jets_neutralHadronMultiplicity;
   TStlSimpleProxy<vector<int> >                                      Jets_neutralMultiplicity;
   TStlSimpleProxy<vector<int> >                                      Jets_origIndex;
   TStlSimpleProxy<vector<int> >                                      Jets_partonFlavor;
   TStlSimpleProxy<vector<float> >                                    Jets_photonEnergyFraction;
   TStlSimpleProxy<vector<int> >                                      Jets_photonMultiplicity;
   TStlSimpleProxy<vector<float> >                                    Jets_pileupId;
   TStlSimpleProxy<vector<float> >                                    Jets_ptD;
   TStlSimpleProxy<vector<float> >                                    Jets_qgLikelihood;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   JetsAK15;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_axismajor;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_axisminor;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_chargedEmEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_chargedHadronEnergyFraction;
   TStlSimpleProxy<vector<int> >                                      JetsAK15_chargedHadronMultiplicity;
   TStlSimpleProxy<vector<int> >                                      JetsAK15_chargedMultiplicity;
   TStlSimpleProxy<vector<int> >                                      JetsAK15_constituentsIndex;
   TStlSimpleProxy<vector<int> >                                      JetsAK15_constituentsIndexCounts;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_DeepMassDecorrelTagbbvsLight;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_DeepMassDecorrelTagHbbvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_DeepMassDecorrelTagTvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_DeepMassDecorrelTagWvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_DeepMassDecorrelTagZbbvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_DeepMassDecorrelTagZHbbvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_DeepMassDecorrelTagZvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_DeepTagHbbvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_DeepTagTvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_DeepTagWvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_DeepTagZbbvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_DeepTagZvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_doubleBDiscriminator;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_ecfC2b1;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_ecfC2b2;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_ecfD2b1;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_ecfD2b2;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_ecfM2b1;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_ecfM2b2;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_ecfN2b1;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_ecfN2b2;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_electronEnergyFraction;
   TStlSimpleProxy<vector<int> >                                      JetsAK15_electronMultiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_girth;
   TStlSimpleProxy<vector<int> >                                      JetsAK15_hadronFlavor;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_hfEMEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_hfHadronEnergyFraction;
   TStlSimpleProxy<vector<bool> >                                     JetsAK15_ID;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_jecFactor;
   TStlSimpleProxy<vector<int> >                                      JetsAK15_multiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_muonEnergyFraction;
   TStlSimpleProxy<vector<int> >                                      JetsAK15_muonMultiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_neutralEmEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_neutralHadronEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_neutralHadronMultiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_neutralMultiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_NsubjettinessTau1;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_NsubjettinessTau2;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_NsubjettinessTau3;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_NsubjettinessTau4;
   TStlSimpleProxy<vector<int> >                                      JetsAK15_NumBhadrons;
   TStlSimpleProxy<vector<int> >                                      JetsAK15_NumChadrons;
   TStlSimpleProxy<vector<int> >                                      JetsAK15_partonFlavor;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_pfMassIndependentDeepDoubleBvLJetTagsProbHbb;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_photonEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_photonMultiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_ptD;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_softDropMass;
   TStlSimpleProxy<vector<float> >                                    JetsAK15_softDropMassBeta1;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   JetsAK15_subjets;
   TStlSimpleProxy<vector<int> >                                      JetsAK15_subjetsCounts;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   JetsAK8;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_axismajor;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_axisminor;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_chargedEmEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_chargedHadronEnergyFraction;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_chargedHadronMultiplicity;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_chargedMultiplicity;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_constituentsIndex;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_constituentsIndexCounts;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_DeepMassDecorrelTagbbvsLight;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_DeepMassDecorrelTagHbbvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_DeepMassDecorrelTagTvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_DeepMassDecorrelTagWvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_DeepMassDecorrelTagZbbvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_DeepMassDecorrelTagZHbbvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_DeepMassDecorrelTagZvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_DeepTagHbbvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_DeepTagTvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_DeepTagWvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_DeepTagZbbvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_DeepTagZvsQCD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_doubleBDiscriminator;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_ecfN2b1;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_ecfN2b2;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_ecfN3b1;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_ecfN3b2;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_electronEnergyFraction;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_electronMultiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_girth;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_hadronFlavor;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_hfEMEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_hfHadronEnergyFraction;
   TStlSimpleProxy<vector<bool> >                                     JetsAK8_ID;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_jecFactor;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_jecUnc;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_jerFactor;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_jerFactorDown;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_jerFactorUp;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_multiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_muonEnergyFraction;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_muonMultiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_neutralEmEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_neutralHadronEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_neutralHadronMultiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_neutralMultiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_NsubjettinessTau1;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_NsubjettinessTau2;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_NsubjettinessTau3;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_NsubjettinessTau4;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_NumBhadrons;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_NumChadrons;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_origIndex;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_partonFlavor;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_photonEnergyFraction;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_photonMultiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_ptD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_softDropMass;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   JetsAK8_subjets;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_subjetsCounts;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_subjets_axismajor;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_subjets_axisminor;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_subjets_jecFactor;
   TStlSimpleProxy<vector<int> >                                      JetsAK8_subjets_multiplicity;
   TStlSimpleProxy<vector<float> >                                    JetsAK8_subjets_ptD;
   TStlSimpleProxy<vector<float> >                                    JetsAK8JECdown_jerFactor;
   TStlSimpleProxy<vector<int> >                                      JetsAK8JECdown_origIndex;
   TStlSimpleProxy<vector<float> >                                    JetsAK8JECup_jerFactor;
   TStlSimpleProxy<vector<int> >                                      JetsAK8JECup_origIndex;
   TStlSimpleProxy<vector<int> >                                      JetsAK8JERdown_origIndex;
   TStlSimpleProxy<vector<int> >                                      JetsAK8JERup_origIndex;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   JetsConstituents;
   TStlSimpleProxy<vector<float> >                                    JetsConstituents_dxy;
   TStlSimpleProxy<vector<float> >                                    JetsConstituents_dxysig;
   TStlSimpleProxy<vector<float> >                                    JetsConstituents_dz;
   TStlSimpleProxy<vector<float> >                                    JetsConstituents_dzsig;
   TStlSimpleProxy<vector<int> >                                      JetsConstituents_PdgId;
   TStlSimpleProxy<vector<float> >                                    JetsConstituents_PuppiWeight;
   TStlSimpleProxy<vector<float> >                                    JetsJECdown_jerFactor;
   TStlSimpleProxy<vector<int> >                                      JetsJECdown_origIndex;
   TStlSimpleProxy<vector<float> >                                    JetsJECup_jerFactor;
   TStlSimpleProxy<vector<int> >                                      JetsJECup_origIndex;
   TStlSimpleProxy<vector<int> >                                      JetsJERdown_origIndex;
   TStlSimpleProxy<vector<int> >                                      JetsJERup_origIndex;
   TFloatProxy                                                        madHT;
   TIntProxy                                                          madMinDeltaRStatus;
   TFloatProxy                                                        madMinPhotonDeltaR;
   TFloatProxy                                                        MET;
   TStlSimpleProxy<vector<float> >                                    METDown;
   TFloatProxy                                                        METPhi;
   TStlSimpleProxy<vector<float> >                                    METPhiDown;
   TStlSimpleProxy<vector<float> >                                    METPhiUp;
   TFloatProxy                                                        METSignificance;
   TStlSimpleProxy<vector<float> >                                    METUp;
   TFloatProxy                                                        MHT;
   TFloatProxy                                                        MHTJECdown;
   TFloatProxy                                                        MHTJECup;
   TFloatProxy                                                        MHTJERdown;
   TFloatProxy                                                        MHTJERup;
   TFloatProxy                                                        MHTPhi;
   TFloatProxy                                                        MHTPhiJECdown;
   TFloatProxy                                                        MHTPhiJECup;
   TFloatProxy                                                        MHTPhiJERdown;
   TFloatProxy                                                        MHTPhiJERup;
   TFloatProxy                                                        MJJ_AK8;
   TFloatProxy                                                        Mmc_AK8;
   TFloatProxy                                                        MT_AK8;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   Muons;
   TStlSimpleProxy<vector<int> >                                      Muons_charge;
   TStlSimpleProxy<vector<float> >                                    Muons_iso;
   TStlSimpleProxy<vector<bool> >                                     Muons_mediumID;
   TStlSimpleProxy<vector<float> >                                    Muons_MTW;
   TStlSimpleProxy<vector<bool> >                                     Muons_passIso;
   TStlSimpleProxy<vector<bool> >                                     Muons_tightID;
   TIntProxy                                                          nAllVertices;
   TIntProxy                                                          NElectrons;
   TIntProxy                                                          NJets;
   TIntProxy                                                          NJetsISR;
   TIntProxy                                                          NJetsISRJECdown;
   TIntProxy                                                          NJetsISRJECup;
   TIntProxy                                                          NJetsISRJERdown;
   TIntProxy                                                          NJetsISRJERup;
   TIntProxy                                                          NJetsJECdown;
   TIntProxy                                                          NJetsJECup;
   TIntProxy                                                          NJetsJERdown;
   TIntProxy                                                          NJetsJERup;
   TIntProxy                                                          NMuons;
   TFloatProxy                                                        NonPrefiringProb;
   TFloatProxy                                                        NonPrefiringProbDown;
   TFloatProxy                                                        NonPrefiringProbECAL;
   TFloatProxy                                                        NonPrefiringProbECALDown;
   TFloatProxy                                                        NonPrefiringProbECALUp;
   TFloatProxy                                                        NonPrefiringProbMuon;
   TFloatProxy                                                        NonPrefiringProbMuonDown;
   TFloatProxy                                                        NonPrefiringProbMuonUp;
   TFloatProxy                                                        NonPrefiringProbUp;
   TFloatProxy                                                        NumEvents;
   TIntProxy                                                          NumInteractions;
   TIntProxy                                                          NVtx;
   TStlSimpleProxy<vector<float> >                                    PDFweights;
   TFloatProxy                                                        PFCaloMETRatio;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   Photons;
   TStlSimpleProxy<vector<int> >                                      Photons_cutBasedID;
   TStlSimpleProxy<vector<bool> >                                     Photons_electronFakes;
   TStlSimpleProxy<vector<bool> >                                     Photons_fullID;
   TStlSimpleProxy<vector<float> >                                    Photons_genMatched;
   TStlSimpleProxy<vector<float> >                                    Photons_hadTowOverEM;
   TStlSimpleProxy<vector<bool> >                                     Photons_hasPixelSeed;
   TStlSimpleProxy<vector<float> >                                    Photons_isEB;
   TStlSimpleProxy<vector<float> >                                    Photons_mvaValuesID;
   TStlSimpleProxy<vector<bool> >                                     Photons_nonPrompt;
   TStlSimpleProxy<vector<float> >                                    Photons_passElectronVeto;
   TStlSimpleProxy<vector<float> >                                    Photons_pfChargedIso;
   TStlSimpleProxy<vector<float> >                                    Photons_pfChargedIsoRhoCorr;
   TStlSimpleProxy<vector<float> >                                    Photons_pfGammaIso;
   TStlSimpleProxy<vector<float> >                                    Photons_pfGammaIsoRhoCorr;
   TStlSimpleProxy<vector<float> >                                    Photons_pfNeutralIso;
   TStlSimpleProxy<vector<float> >                                    Photons_pfNeutralIsoRhoCorr;
   TStlSimpleProxy<vector<float> >                                    Photons_sigmaIetaIeta;
   TIntProxy                                                          PrimaryVertexFilter;
   TStlSimpleProxy<vector<float> >                                    PSweights;
   TFloatProxy                                                        PuppiMET;
   TStlSimpleProxy<vector<float> >                                    PuppiMETDown;
   TFloatProxy                                                        PuppiMETPhi;
   TStlSimpleProxy<vector<float> >                                    PuppiMETPhiDown;
   TStlSimpleProxy<vector<float> >                                    PuppiMETPhiUp;
   TStlSimpleProxy<vector<float> >                                    PuppiMETUp;
   TFloatProxy                                                        puSysDown;
   TFloatProxy                                                        puSysUp;
   TFloatProxy                                                        puWeight;
   TStlSimpleProxy<vector<float> >                                    ScaleWeights;
   TStlSimpleProxy<vector<float> >                                    SignalParameters;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   TAPElectronTracks;
   TStlSimpleProxy<vector<float> >                                    TAPElectronTracks_dxypv;
   TStlSimpleProxy<vector<bool> >                                     TAPElectronTracks_leptonMatch;
   TStlSimpleProxy<vector<float> >                                    TAPElectronTracks_mT;
   TStlSimpleProxy<vector<float> >                                    TAPElectronTracks_pfRelIso03chg;
   TStlSimpleProxy<vector<float> >                                    TAPElectronTracks_trkiso;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   TAPMuonTracks;
   TStlSimpleProxy<vector<float> >                                    TAPMuonTracks_dxypv;
   TStlSimpleProxy<vector<bool> >                                     TAPMuonTracks_leptonMatch;
   TStlSimpleProxy<vector<float> >                                    TAPMuonTracks_mT;
   TStlSimpleProxy<vector<float> >                                    TAPMuonTracks_pfRelIso03chg;
   TStlSimpleProxy<vector<float> >                                    TAPMuonTracks_trkiso;
   TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__   TAPPionTracks;
   TStlSimpleProxy<vector<float> >                                    TAPPionTracks_dxypv;
   TStlSimpleProxy<vector<bool> >                                     TAPPionTracks_leptonMatch;
   TStlSimpleProxy<vector<float> >                                    TAPPionTracks_mT;
   TStlSimpleProxy<vector<float> >                                    TAPPionTracks_pfRelIso03chg;
   TStlSimpleProxy<vector<float> >                                    TAPPionTracks_trkiso;
   TStlSimpleProxy<vector<int> >                                      TriggerPass;
   TStlSimpleProxy<vector<int> >                                      TriggerPrescales;
   TStlSimpleProxy<vector<int> >                                      TriggerVersion;
   TFloatProxy                                                        TrueNumInteractions;
   TFloatProxy                                                        Weight;


   //   NtupleVarsTProxy(TTree *tree=0) : 
   NtupleVarsTProxy(TTree *tree=0) : 
      fChain(0),
      fDirector(tree,-1),
      // fClass                (TClass::GetClass("NtupleVarsTProxy")),
      RunNum                                                            (&fDirector,"RunNum"),
      LumiBlockNum                                                      (&fDirector,"LumiBlockNum"),
      EvtNum                                                            (&fDirector,"EvtNum"),
      BadChargedCandidateFilter                                         (&fDirector,"BadChargedCandidateFilter"),
      BadPFMuonDzFilter                                                 (&fDirector,"BadPFMuonDzFilter"),
      BadPFMuonFilter                                                   (&fDirector,"BadPFMuonFilter"),
      BTagsDeepCSV                                                      (&fDirector,"BTagsDeepCSV"),
      BTagsDeepCSVJECdown                                               (&fDirector,"BTagsDeepCSVJECdown"),
      BTagsDeepCSVJECup                                                 (&fDirector,"BTagsDeepCSVJECup"),
      BTagsDeepCSVJERdown                                               (&fDirector,"BTagsDeepCSVJERdown"),
      BTagsDeepCSVJERup                                                 (&fDirector,"BTagsDeepCSVJERup"),
      CaloMET                                                           (&fDirector,"CaloMET"),
      CaloMETPhi                                                        (&fDirector,"CaloMETPhi"),
      CrossSection                                                      (&fDirector,"CrossSection"),
      CSCTightHaloFilter                                                (&fDirector,"CSCTightHaloFilter"),
      DeltaPhi1                                                         (&fDirector,"DeltaPhi1"),
      DeltaPhi1_AK8                                                     (&fDirector,"DeltaPhi1_AK8"),
      DeltaPhi1JECdown                                                  (&fDirector,"DeltaPhi1JECdown"),
      DeltaPhi1JECup                                                    (&fDirector,"DeltaPhi1JECup"),
      DeltaPhi1JERdown                                                  (&fDirector,"DeltaPhi1JERdown"),
      DeltaPhi1JERup                                                    (&fDirector,"DeltaPhi1JERup"),
      DeltaPhi2                                                         (&fDirector,"DeltaPhi2"),
      DeltaPhi2_AK8                                                     (&fDirector,"DeltaPhi2_AK8"),
      DeltaPhi2JECdown                                                  (&fDirector,"DeltaPhi2JECdown"),
      DeltaPhi2JECup                                                    (&fDirector,"DeltaPhi2JECup"),
      DeltaPhi2JERdown                                                  (&fDirector,"DeltaPhi2JERdown"),
      DeltaPhi2JERup                                                    (&fDirector,"DeltaPhi2JERup"),
      DeltaPhi3                                                         (&fDirector,"DeltaPhi3"),
      DeltaPhi3JECdown                                                  (&fDirector,"DeltaPhi3JECdown"),
      DeltaPhi3JECup                                                    (&fDirector,"DeltaPhi3JECup"),
      DeltaPhi3JERdown                                                  (&fDirector,"DeltaPhi3JERdown"),
      DeltaPhi3JERup                                                    (&fDirector,"DeltaPhi3JERup"),
      DeltaPhi4                                                         (&fDirector,"DeltaPhi4"),
      DeltaPhi4JECdown                                                  (&fDirector,"DeltaPhi4JECdown"),
      DeltaPhi4JECup                                                    (&fDirector,"DeltaPhi4JECup"),
      DeltaPhi4JERdown                                                  (&fDirector,"DeltaPhi4JERdown"),
      DeltaPhi4JERup                                                    (&fDirector,"DeltaPhi4JERup"),
      DeltaPhiMin_AK8                                                   (&fDirector,"DeltaPhiMin_AK8"),
      ecalBadCalibFilter                                                (&fDirector,"ecalBadCalibFilter"),
      EcalDeadCellBoundaryEnergyFilter                                  (&fDirector,"EcalDeadCellBoundaryEnergyFilter"),
      EcalDeadCellTriggerPrimitiveFilter                                (&fDirector,"EcalDeadCellTriggerPrimitiveFilter"),
      eeBadScFilter                                                     (&fDirector,"eeBadScFilter"),
      Electrons                                                         (&fDirector,"Electrons"),
      fCoordinates                                                      (&fDirector,"Electrons.fCoordinates"),
      Electrons_charge                                                  (&fDirector,"Electrons_charge"),
      Electrons_iso                                                     (&fDirector,"Electrons_iso"),
      Electrons_mediumID                                                (&fDirector,"Electrons_mediumID"),
      Electrons_MTW                                                     (&fDirector,"Electrons_MTW"),
      Electrons_passIso                                                 (&fDirector,"Electrons_passIso"),
      Electrons_tightID                                                 (&fDirector,"Electrons_tightID"),
      fixedGridRhoFastjetAll                                            (&fDirector,"fixedGridRhoFastjetAll"),
      GenElectrons                                                      (&fDirector,"GenElectrons"),
      GenHT                                                             (&fDirector,"GenHT"),
      GenJets                                                           (&fDirector,"GenJets"),
      GenJetsAK15                                                       (&fDirector,"GenJetsAK15"),
      GenJetsAK8                                                        (&fDirector,"GenJetsAK8"),
      GenJetsAK8_multiplicity                                           (&fDirector,"GenJetsAK8_multiplicity"),
      GenJetsAK8_softDropMass                                           (&fDirector,"GenJetsAK8_softDropMass"),
      GenMET                                                            (&fDirector,"GenMET"),
      GenMETPhi                                                         (&fDirector,"GenMETPhi"),
      GenMHT                                                            (&fDirector,"GenMHT"),
      GenMHTPhi                                                         (&fDirector,"GenMHTPhi"),
      GenMT2_AK8                                                        (&fDirector,"GenMT2_AK8"),
      GenMuons                                                          (&fDirector,"GenMuons"),
      GenParticles                                                      (&fDirector,"GenParticles"),
      GenParticles_Charge                                               (&fDirector,"GenParticles_Charge"),
      GenParticles_ParentId                                             (&fDirector,"GenParticles_ParentId"),
      GenParticles_ParentIdx                                            (&fDirector,"GenParticles_ParentIdx"),
      GenParticles_PdgId                                                (&fDirector,"GenParticles_PdgId"),
      GenParticles_Status                                               (&fDirector,"GenParticles_Status"),
      GenTaus                                                           (&fDirector,"GenTaus"),
      GenTaus_had                                                       (&fDirector,"GenTaus_had"),
      globalSuperTightHalo2016Filter                                    (&fDirector,"globalSuperTightHalo2016Filter"),
      globalTightHalo2016Filter                                         (&fDirector,"globalTightHalo2016Filter"),
      hasGenPromptPhoton                                                (&fDirector,"hasGenPromptPhoton"),
      HBHEIsoNoiseFilter                                                (&fDirector,"HBHEIsoNoiseFilter"),
      HBHENoiseFilter                                                   (&fDirector,"HBHENoiseFilter"),
      hfNoisyHitsFilter                                                 (&fDirector,"hfNoisyHitsFilter"),
      HT                                                                (&fDirector,"HT"),
      HT5                                                               (&fDirector,"HT5"),
      HT5JECdown                                                        (&fDirector,"HT5JECdown"),
      HT5JECup                                                          (&fDirector,"HT5JECup"),
      HT5JERdown                                                        (&fDirector,"HT5JERdown"),
      HT5JERup                                                          (&fDirector,"HT5JERup"),
      HTJECdown                                                         (&fDirector,"HTJECdown"),
      HTJECup                                                           (&fDirector,"HTJECup"),
      HTJERdown                                                         (&fDirector,"HTJERdown"),
      HTJERup                                                           (&fDirector,"HTJERup"),
      isoElectronTracks                                                 (&fDirector,"isoElectronTracks"),
      isoMuonTracks                                                     (&fDirector,"isoMuonTracks"),
      isoPionTracks                                                     (&fDirector,"isoPionTracks"),
      JetID                                                             (&fDirector,"JetID"),
      JetIDAK15                                                         (&fDirector,"JetIDAK15"),
      JetIDAK8                                                          (&fDirector,"JetIDAK8"),
      JetIDAK8JECdown                                                   (&fDirector,"JetIDAK8JECdown"),
      JetIDAK8JECup                                                     (&fDirector,"JetIDAK8JECup"),
      JetIDAK8JERdown                                                   (&fDirector,"JetIDAK8JERdown"),
      JetIDAK8JERup                                                     (&fDirector,"JetIDAK8JERup"),
      JetIDJECdown                                                      (&fDirector,"JetIDJECdown"),
      JetIDJECup                                                        (&fDirector,"JetIDJECup"),
      JetIDJERdown                                                      (&fDirector,"JetIDJERdown"),
      JetIDJERup                                                        (&fDirector,"JetIDJERup"),
      Jets                                                              (&fDirector,"Jets"),
      Jets_axismajor                                                    (&fDirector,"Jets_axismajor"),
      Jets_axisminor                                                    (&fDirector,"Jets_axisminor"),
      Jets_bDiscriminatorCSV                                            (&fDirector,"Jets_bDiscriminatorCSV"),
      Jets_bJetTagDeepCSVBvsAll                                         (&fDirector,"Jets_bJetTagDeepCSVBvsAll"),
      Jets_bJetTagDeepCSVprobb                                          (&fDirector,"Jets_bJetTagDeepCSVprobb"),
      Jets_bJetTagDeepCSVprobbb                                         (&fDirector,"Jets_bJetTagDeepCSVprobbb"),
      Jets_bJetTagDeepCSVprobc                                          (&fDirector,"Jets_bJetTagDeepCSVprobc"),
      Jets_bJetTagDeepCSVprobudsg                                       (&fDirector,"Jets_bJetTagDeepCSVprobudsg"),
      Jets_bJetTagDeepFlavourprobb                                      (&fDirector,"Jets_bJetTagDeepFlavourprobb"),
      Jets_bJetTagDeepFlavourprobbb                                     (&fDirector,"Jets_bJetTagDeepFlavourprobbb"),
      Jets_bJetTagDeepFlavourprobc                                      (&fDirector,"Jets_bJetTagDeepFlavourprobc"),
      Jets_bJetTagDeepFlavourprobg                                      (&fDirector,"Jets_bJetTagDeepFlavourprobg"),
      Jets_bJetTagDeepFlavourproblepb                                   (&fDirector,"Jets_bJetTagDeepFlavourproblepb"),
      Jets_bJetTagDeepFlavourprobuds                                    (&fDirector,"Jets_bJetTagDeepFlavourprobuds"),
      Jets_chargedEmEnergyFraction                                      (&fDirector,"Jets_chargedEmEnergyFraction"),
      Jets_chargedHadronEnergyFraction                                  (&fDirector,"Jets_chargedHadronEnergyFraction"),
      Jets_chargedHadronMultiplicity                                    (&fDirector,"Jets_chargedHadronMultiplicity"),
      Jets_chargedMultiplicity                                          (&fDirector,"Jets_chargedMultiplicity"),
      Jets_electronEnergyFraction                                       (&fDirector,"Jets_electronEnergyFraction"),
      Jets_electronMultiplicity                                         (&fDirector,"Jets_electronMultiplicity"),
      Jets_hadronFlavor                                                 (&fDirector,"Jets_hadronFlavor"),
      Jets_hfEMEnergyFraction                                           (&fDirector,"Jets_hfEMEnergyFraction"),
      Jets_hfHadronEnergyFraction                                       (&fDirector,"Jets_hfHadronEnergyFraction"),
      Jets_HTMask                                                       (&fDirector,"Jets_HTMask"),
      Jets_ID                                                           (&fDirector,"Jets_ID"),
      Jets_jecFactor                                                    (&fDirector,"Jets_jecFactor"),
      Jets_jecUnc                                                       (&fDirector,"Jets_jecUnc"),
      Jets_jerFactor                                                    (&fDirector,"Jets_jerFactor"),
      Jets_jerFactorDown                                                (&fDirector,"Jets_jerFactorDown"),
      Jets_jerFactorUp                                                  (&fDirector,"Jets_jerFactorUp"),
      Jets_LeptonMask                                                   (&fDirector,"Jets_LeptonMask"),
      Jets_MHTMask                                                      (&fDirector,"Jets_MHTMask"),
      Jets_multiplicity                                                 (&fDirector,"Jets_multiplicity"),
      Jets_muonEnergyFraction                                           (&fDirector,"Jets_muonEnergyFraction"),
      Jets_muonMultiplicity                                             (&fDirector,"Jets_muonMultiplicity"),
      Jets_neutralEmEnergyFraction                                      (&fDirector,"Jets_neutralEmEnergyFraction"),
      Jets_neutralHadronEnergyFraction                                  (&fDirector,"Jets_neutralHadronEnergyFraction"),
      Jets_neutralHadronMultiplicity                                    (&fDirector,"Jets_neutralHadronMultiplicity"),
      Jets_neutralMultiplicity                                          (&fDirector,"Jets_neutralMultiplicity"),
      Jets_origIndex                                                    (&fDirector,"Jets_origIndex"),
      Jets_partonFlavor                                                 (&fDirector,"Jets_partonFlavor"),
      Jets_photonEnergyFraction                                         (&fDirector,"Jets_photonEnergyFraction"),
      Jets_photonMultiplicity                                           (&fDirector,"Jets_photonMultiplicity"),
      Jets_pileupId                                                     (&fDirector,"Jets_pileupId"),
      Jets_ptD                                                          (&fDirector,"Jets_ptD"),
      Jets_qgLikelihood                                                 (&fDirector,"Jets_qgLikelihood"),
      JetsAK15                                                          (&fDirector,"JetsAK15"),
      JetsAK15_axismajor                                                (&fDirector,"JetsAK15_axismajor"),
      JetsAK15_axisminor                                                (&fDirector,"JetsAK15_axisminor"),
      JetsAK15_chargedEmEnergyFraction                                  (&fDirector,"JetsAK15_chargedEmEnergyFraction"),
      JetsAK15_chargedHadronEnergyFraction                              (&fDirector,"JetsAK15_chargedHadronEnergyFraction"),
      JetsAK15_chargedHadronMultiplicity                                (&fDirector,"JetsAK15_chargedHadronMultiplicity"),
      JetsAK15_chargedMultiplicity                                      (&fDirector,"JetsAK15_chargedMultiplicity"),
      JetsAK15_constituentsIndex                                        (&fDirector,"JetsAK15_constituentsIndex"),
      JetsAK15_constituentsIndexCounts                                  (&fDirector,"JetsAK15_constituentsIndexCounts"),
      JetsAK15_DeepMassDecorrelTagbbvsLight                             (&fDirector,"JetsAK15_DeepMassDecorrelTagbbvsLight"),
      JetsAK15_DeepMassDecorrelTagHbbvsQCD                              (&fDirector,"JetsAK15_DeepMassDecorrelTagHbbvsQCD"),
      JetsAK15_DeepMassDecorrelTagTvsQCD                                (&fDirector,"JetsAK15_DeepMassDecorrelTagTvsQCD"),
      JetsAK15_DeepMassDecorrelTagWvsQCD                                (&fDirector,"JetsAK15_DeepMassDecorrelTagWvsQCD"),
      JetsAK15_DeepMassDecorrelTagZbbvsQCD                              (&fDirector,"JetsAK15_DeepMassDecorrelTagZbbvsQCD"),
      JetsAK15_DeepMassDecorrelTagZHbbvsQCD                             (&fDirector,"JetsAK15_DeepMassDecorrelTagZHbbvsQCD"),
      JetsAK15_DeepMassDecorrelTagZvsQCD                                (&fDirector,"JetsAK15_DeepMassDecorrelTagZvsQCD"),
      JetsAK15_DeepTagHbbvsQCD                                          (&fDirector,"JetsAK15_DeepTagHbbvsQCD"),
      JetsAK15_DeepTagTvsQCD                                            (&fDirector,"JetsAK15_DeepTagTvsQCD"),
      JetsAK15_DeepTagWvsQCD                                            (&fDirector,"JetsAK15_DeepTagWvsQCD"),
      JetsAK15_DeepTagZbbvsQCD                                          (&fDirector,"JetsAK15_DeepTagZbbvsQCD"),
      JetsAK15_DeepTagZvsQCD                                            (&fDirector,"JetsAK15_DeepTagZvsQCD"),
      JetsAK15_doubleBDiscriminator                                     (&fDirector,"JetsAK15_doubleBDiscriminator"),
      JetsAK15_ecfC2b1                                                  (&fDirector,"JetsAK15_ecfC2b1"),
      JetsAK15_ecfC2b2                                                  (&fDirector,"JetsAK15_ecfC2b2"),
      JetsAK15_ecfD2b1                                                  (&fDirector,"JetsAK15_ecfD2b1"),
      JetsAK15_ecfD2b2                                                  (&fDirector,"JetsAK15_ecfD2b2"),
      JetsAK15_ecfM2b1                                                  (&fDirector,"JetsAK15_ecfM2b1"),
      JetsAK15_ecfM2b2                                                  (&fDirector,"JetsAK15_ecfM2b2"),
      JetsAK15_ecfN2b1                                                  (&fDirector,"JetsAK15_ecfN2b1"),
      JetsAK15_ecfN2b2                                                  (&fDirector,"JetsAK15_ecfN2b2"),
      JetsAK15_electronEnergyFraction                                   (&fDirector,"JetsAK15_electronEnergyFraction"),
      JetsAK15_electronMultiplicity                                     (&fDirector,"JetsAK15_electronMultiplicity"),
      JetsAK15_girth                                                    (&fDirector,"JetsAK15_girth"),
      JetsAK15_hadronFlavor                                             (&fDirector,"JetsAK15_hadronFlavor"),
      JetsAK15_hfEMEnergyFraction                                       (&fDirector,"JetsAK15_hfEMEnergyFraction"),
      JetsAK15_hfHadronEnergyFraction                                   (&fDirector,"JetsAK15_hfHadronEnergyFraction"),
      JetsAK15_ID                                                       (&fDirector,"JetsAK15_ID"),
      JetsAK15_jecFactor                                                (&fDirector,"JetsAK15_jecFactor"),
      JetsAK15_multiplicity                                             (&fDirector,"JetsAK15_multiplicity"),
      JetsAK15_muonEnergyFraction                                       (&fDirector,"JetsAK15_muonEnergyFraction"),
      JetsAK15_muonMultiplicity                                         (&fDirector,"JetsAK15_muonMultiplicity"),
      JetsAK15_neutralEmEnergyFraction                                  (&fDirector,"JetsAK15_neutralEmEnergyFraction"),
      JetsAK15_neutralHadronEnergyFraction                              (&fDirector,"JetsAK15_neutralHadronEnergyFraction"),
      JetsAK15_neutralHadronMultiplicity                                (&fDirector,"JetsAK15_neutralHadronMultiplicity"),
      JetsAK15_neutralMultiplicity                                      (&fDirector,"JetsAK15_neutralMultiplicity"),
      JetsAK15_NsubjettinessTau1                                        (&fDirector,"JetsAK15_NsubjettinessTau1"),
      JetsAK15_NsubjettinessTau2                                        (&fDirector,"JetsAK15_NsubjettinessTau2"),
      JetsAK15_NsubjettinessTau3                                        (&fDirector,"JetsAK15_NsubjettinessTau3"),
      JetsAK15_NsubjettinessTau4                                        (&fDirector,"JetsAK15_NsubjettinessTau4"),
      JetsAK15_NumBhadrons                                              (&fDirector,"JetsAK15_NumBhadrons"),
      JetsAK15_NumChadrons                                              (&fDirector,"JetsAK15_NumChadrons"),
      JetsAK15_partonFlavor                                             (&fDirector,"JetsAK15_partonFlavor"),
      JetsAK15_pfMassIndependentDeepDoubleBvLJetTagsProbHbb             (&fDirector,"JetsAK15_pfMassIndependentDeepDoubleBvLJetTagsProbHbb"),
      JetsAK15_photonEnergyFraction                                     (&fDirector,"JetsAK15_photonEnergyFraction"),
      JetsAK15_photonMultiplicity                                       (&fDirector,"JetsAK15_photonMultiplicity"),
      JetsAK15_ptD                                                      (&fDirector,"JetsAK15_ptD"),
      JetsAK15_softDropMass                                             (&fDirector,"JetsAK15_softDropMass"),
      JetsAK15_softDropMassBeta1                                        (&fDirector,"JetsAK15_softDropMassBeta1"),
      JetsAK15_subjets                                                  (&fDirector,"JetsAK15_subjets"),
      JetsAK15_subjetsCounts                                            (&fDirector,"JetsAK15_subjetsCounts"),
      JetsAK8                                                           (&fDirector,"JetsAK8"),
      JetsAK8_axismajor                                                 (&fDirector,"JetsAK8_axismajor"),
      JetsAK8_axisminor                                                 (&fDirector,"JetsAK8_axisminor"),
      JetsAK8_chargedEmEnergyFraction                                   (&fDirector,"JetsAK8_chargedEmEnergyFraction"),
      JetsAK8_chargedHadronEnergyFraction                               (&fDirector,"JetsAK8_chargedHadronEnergyFraction"),
      JetsAK8_chargedHadronMultiplicity                                 (&fDirector,"JetsAK8_chargedHadronMultiplicity"),
      JetsAK8_chargedMultiplicity                                       (&fDirector,"JetsAK8_chargedMultiplicity"),
      JetsAK8_constituentsIndex                                         (&fDirector,"JetsAK8_constituentsIndex"),
      JetsAK8_constituentsIndexCounts                                   (&fDirector,"JetsAK8_constituentsIndexCounts"),
      JetsAK8_DeepMassDecorrelTagbbvsLight                              (&fDirector,"JetsAK8_DeepMassDecorrelTagbbvsLight"),
      JetsAK8_DeepMassDecorrelTagHbbvsQCD                               (&fDirector,"JetsAK8_DeepMassDecorrelTagHbbvsQCD"),
      JetsAK8_DeepMassDecorrelTagTvsQCD                                 (&fDirector,"JetsAK8_DeepMassDecorrelTagTvsQCD"),
      JetsAK8_DeepMassDecorrelTagWvsQCD                                 (&fDirector,"JetsAK8_DeepMassDecorrelTagWvsQCD"),
      JetsAK8_DeepMassDecorrelTagZbbvsQCD                               (&fDirector,"JetsAK8_DeepMassDecorrelTagZbbvsQCD"),
      JetsAK8_DeepMassDecorrelTagZHbbvsQCD                              (&fDirector,"JetsAK8_DeepMassDecorrelTagZHbbvsQCD"),
      JetsAK8_DeepMassDecorrelTagZvsQCD                                 (&fDirector,"JetsAK8_DeepMassDecorrelTagZvsQCD"),
      JetsAK8_DeepTagHbbvsQCD                                           (&fDirector,"JetsAK8_DeepTagHbbvsQCD"),
      JetsAK8_DeepTagTvsQCD                                             (&fDirector,"JetsAK8_DeepTagTvsQCD"),
      JetsAK8_DeepTagWvsQCD                                             (&fDirector,"JetsAK8_DeepTagWvsQCD"),
      JetsAK8_DeepTagZbbvsQCD                                           (&fDirector,"JetsAK8_DeepTagZbbvsQCD"),
      JetsAK8_DeepTagZvsQCD                                             (&fDirector,"JetsAK8_DeepTagZvsQCD"),
      JetsAK8_doubleBDiscriminator                                      (&fDirector,"JetsAK8_doubleBDiscriminator"),
      JetsAK8_ecfN2b1                                                   (&fDirector,"JetsAK8_ecfN2b1"),
      JetsAK8_ecfN2b2                                                   (&fDirector,"JetsAK8_ecfN2b2"),
      JetsAK8_ecfN3b1                                                   (&fDirector,"JetsAK8_ecfN3b1"),
      JetsAK8_ecfN3b2                                                   (&fDirector,"JetsAK8_ecfN3b2"),
      JetsAK8_electronEnergyFraction                                    (&fDirector,"JetsAK8_electronEnergyFraction"),
      JetsAK8_electronMultiplicity                                      (&fDirector,"JetsAK8_electronMultiplicity"),
      JetsAK8_girth                                                     (&fDirector,"JetsAK8_girth"),
      JetsAK8_hadronFlavor                                              (&fDirector,"JetsAK8_hadronFlavor"),
      JetsAK8_hfEMEnergyFraction                                        (&fDirector,"JetsAK8_hfEMEnergyFraction"),
      JetsAK8_hfHadronEnergyFraction                                    (&fDirector,"JetsAK8_hfHadronEnergyFraction"),
      JetsAK8_ID                                                        (&fDirector,"JetsAK8_ID"),
      JetsAK8_jecFactor                                                 (&fDirector,"JetsAK8_jecFactor"),
      JetsAK8_jecUnc                                                    (&fDirector,"JetsAK8_jecUnc"),
      JetsAK8_jerFactor                                                 (&fDirector,"JetsAK8_jerFactor"),
      JetsAK8_jerFactorDown                                             (&fDirector,"JetsAK8_jerFactorDown"),
      JetsAK8_jerFactorUp                                               (&fDirector,"JetsAK8_jerFactorUp"),
      JetsAK8_multiplicity                                              (&fDirector,"JetsAK8_multiplicity"),
      JetsAK8_muonEnergyFraction                                        (&fDirector,"JetsAK8_muonEnergyFraction"),
      JetsAK8_muonMultiplicity                                          (&fDirector,"JetsAK8_muonMultiplicity"),
      JetsAK8_neutralEmEnergyFraction                                   (&fDirector,"JetsAK8_neutralEmEnergyFraction"),
      JetsAK8_neutralHadronEnergyFraction                               (&fDirector,"JetsAK8_neutralHadronEnergyFraction"),
      JetsAK8_neutralHadronMultiplicity                                 (&fDirector,"JetsAK8_neutralHadronMultiplicity"),
      JetsAK8_neutralMultiplicity                                       (&fDirector,"JetsAK8_neutralMultiplicity"),
      JetsAK8_NsubjettinessTau1                                         (&fDirector,"JetsAK8_NsubjettinessTau1"),
      JetsAK8_NsubjettinessTau2                                         (&fDirector,"JetsAK8_NsubjettinessTau2"),
      JetsAK8_NsubjettinessTau3                                         (&fDirector,"JetsAK8_NsubjettinessTau3"),
      JetsAK8_NsubjettinessTau4                                         (&fDirector,"JetsAK8_NsubjettinessTau4"),
      JetsAK8_NumBhadrons                                               (&fDirector,"JetsAK8_NumBhadrons"),
      JetsAK8_NumChadrons                                               (&fDirector,"JetsAK8_NumChadrons"),
      JetsAK8_origIndex                                                 (&fDirector,"JetsAK8_origIndex"),
      JetsAK8_partonFlavor                                              (&fDirector,"JetsAK8_partonFlavor"),
      JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb              (&fDirector,"JetsAK8_pfMassIndependentDeepDoubleBvLJetTagsProbHbb"),
      JetsAK8_photonEnergyFraction                                      (&fDirector,"JetsAK8_photonEnergyFraction"),
      JetsAK8_photonMultiplicity                                        (&fDirector,"JetsAK8_photonMultiplicity"),
      JetsAK8_ptD                                                       (&fDirector,"JetsAK8_ptD"),
      JetsAK8_softDropMass                                              (&fDirector,"JetsAK8_softDropMass"),
      JetsAK8_subjets                                                   (&fDirector,"JetsAK8_subjets"),
      JetsAK8_subjetsCounts                                             (&fDirector,"JetsAK8_subjetsCounts"),
      JetsAK8_subjets_axismajor                                         (&fDirector,"JetsAK8_subjets_axismajor"),
      JetsAK8_subjets_axisminor                                         (&fDirector,"JetsAK8_subjets_axisminor"),
      JetsAK8_subjets_jecFactor                                         (&fDirector,"JetsAK8_subjets_jecFactor"),
      JetsAK8_subjets_multiplicity                                      (&fDirector,"JetsAK8_subjets_multiplicity"),
      JetsAK8_subjets_ptD                                               (&fDirector,"JetsAK8_subjets_ptD"),
      JetsAK8JECdown_jerFactor                                          (&fDirector,"JetsAK8JECdown_jerFactor"),
      JetsAK8JECdown_origIndex                                          (&fDirector,"JetsAK8JECdown_origIndex"),
      JetsAK8JECup_jerFactor                                            (&fDirector,"JetsAK8JECup_jerFactor"),
      JetsAK8JECup_origIndex                                            (&fDirector,"JetsAK8JECup_origIndex"),
      JetsAK8JERdown_origIndex                                          (&fDirector,"JetsAK8JERdown_origIndex"),
      JetsAK8JERup_origIndex                                            (&fDirector,"JetsAK8JERup_origIndex"),
      JetsConstituents                                                  (&fDirector,"JetsConstituents"),
      JetsConstituents_dxy                                              (&fDirector,"JetsConstituents_dxy"),
      JetsConstituents_dxysig                                           (&fDirector,"JetsConstituents_dxysig"),
      JetsConstituents_dz                                               (&fDirector,"JetsConstituents_dz"),
      JetsConstituents_dzsig                                            (&fDirector,"JetsConstituents_dzsig"),
      JetsConstituents_PdgId                                            (&fDirector,"JetsConstituents_PdgId"),
      JetsConstituents_PuppiWeight                                      (&fDirector,"JetsConstituents_PuppiWeight"),
      JetsJECdown_jerFactor                                             (&fDirector,"JetsJECdown_jerFactor"),
      JetsJECdown_origIndex                                             (&fDirector,"JetsJECdown_origIndex"),
      JetsJECup_jerFactor                                               (&fDirector,"JetsJECup_jerFactor"),
      JetsJECup_origIndex                                               (&fDirector,"JetsJECup_origIndex"),
      JetsJERdown_origIndex                                             (&fDirector,"JetsJERdown_origIndex"),
      JetsJERup_origIndex                                               (&fDirector,"JetsJERup_origIndex"),
      madHT                                                             (&fDirector,"madHT"),
      madMinDeltaRStatus                                                (&fDirector,"madMinDeltaRStatus"),
      madMinPhotonDeltaR                                                (&fDirector,"madMinPhotonDeltaR"),
      MET                                                               (&fDirector,"MET"),
      METDown                                                           (&fDirector,"METDown"),
      METPhi                                                            (&fDirector,"METPhi"),
      METPhiDown                                                        (&fDirector,"METPhiDown"),
      METPhiUp                                                          (&fDirector,"METPhiUp"),
      METSignificance                                                   (&fDirector,"METSignificance"),
      METUp                                                             (&fDirector,"METUp"),
      MHT                                                               (&fDirector,"MHT"),
      MHTJECdown                                                        (&fDirector,"MHTJECdown"),
      MHTJECup                                                          (&fDirector,"MHTJECup"),
      MHTJERdown                                                        (&fDirector,"MHTJERdown"),
      MHTJERup                                                          (&fDirector,"MHTJERup"),
      MHTPhi                                                            (&fDirector,"MHTPhi"),
      MHTPhiJECdown                                                     (&fDirector,"MHTPhiJECdown"),
      MHTPhiJECup                                                       (&fDirector,"MHTPhiJECup"),
      MHTPhiJERdown                                                     (&fDirector,"MHTPhiJERdown"),
      MHTPhiJERup                                                       (&fDirector,"MHTPhiJERup"),
      MJJ_AK8                                                           (&fDirector,"MJJ_AK8"),
      Mmc_AK8                                                           (&fDirector,"Mmc_AK8"),
      MT_AK8                                                            (&fDirector,"MT_AK8"),
      Muons                                                             (&fDirector,"Muons"),
      Muons_charge                                                      (&fDirector,"Muons_charge"),
      Muons_iso                                                         (&fDirector,"Muons_iso"),
      Muons_mediumID                                                    (&fDirector,"Muons_mediumID"),
      Muons_MTW                                                         (&fDirector,"Muons_MTW"),
      Muons_passIso                                                     (&fDirector,"Muons_passIso"),
      Muons_tightID                                                     (&fDirector,"Muons_tightID"),
      nAllVertices                                                      (&fDirector,"nAllVertices"),
      NElectrons                                                        (&fDirector,"NElectrons"),
      NJets                                                             (&fDirector,"NJets"),
      NJetsISR                                                          (&fDirector,"NJetsISR"),
      NJetsISRJECdown                                                   (&fDirector,"NJetsISRJECdown"),
      NJetsISRJECup                                                     (&fDirector,"NJetsISRJECup"),
      NJetsISRJERdown                                                   (&fDirector,"NJetsISRJERdown"),
      NJetsISRJERup                                                     (&fDirector,"NJetsISRJERup"),
      NJetsJECdown                                                      (&fDirector,"NJetsJECdown"),
      NJetsJECup                                                        (&fDirector,"NJetsJECup"),
      NJetsJERdown                                                      (&fDirector,"NJetsJERdown"),
      NJetsJERup                                                        (&fDirector,"NJetsJERup"),
      NMuons                                                            (&fDirector,"NMuons"),
      NonPrefiringProb                                                  (&fDirector,"NonPrefiringProb"),
      NonPrefiringProbDown                                              (&fDirector,"NonPrefiringProbDown"),
      NonPrefiringProbECAL                                              (&fDirector,"NonPrefiringProbECAL"),
      NonPrefiringProbECALDown                                          (&fDirector,"NonPrefiringProbECALDown"),
      NonPrefiringProbECALUp                                            (&fDirector,"NonPrefiringProbECALUp"),
      NonPrefiringProbMuon                                              (&fDirector,"NonPrefiringProbMuon"),
      NonPrefiringProbMuonDown                                          (&fDirector,"NonPrefiringProbMuonDown"),
      NonPrefiringProbMuonUp                                            (&fDirector,"NonPrefiringProbMuonUp"),
      NonPrefiringProbUp                                                (&fDirector,"NonPrefiringProbUp"),
      NumEvents                                                         (&fDirector,"NumEvents"),
      NumInteractions                                                   (&fDirector,"NumInteractions"),
      NVtx                                                              (&fDirector,"NVtx"),
      PDFweights                                                        (&fDirector,"PDFweights"),
      PFCaloMETRatio                                                    (&fDirector,"PFCaloMETRatio"),
      Photons                                                           (&fDirector,"Photons"),
      Photons_cutBasedID                                                (&fDirector,"Photons_cutBasedID"),
      Photons_electronFakes                                             (&fDirector,"Photons_electronFakes"),
      Photons_fullID                                                    (&fDirector,"Photons_fullID"),
      Photons_genMatched                                                (&fDirector,"Photons_genMatched"),
      Photons_hadTowOverEM                                              (&fDirector,"Photons_hadTowOverEM"),
      Photons_hasPixelSeed                                              (&fDirector,"Photons_hasPixelSeed"),
      Photons_isEB                                                      (&fDirector,"Photons_isEB"),
      Photons_mvaValuesID                                               (&fDirector,"Photons_mvaValuesID"),
      Photons_nonPrompt                                                 (&fDirector,"Photons_nonPrompt"),
      Photons_passElectronVeto                                          (&fDirector,"Photons_passElectronVeto"),
      Photons_pfChargedIso                                              (&fDirector,"Photons_pfChargedIso"),
      Photons_pfChargedIsoRhoCorr                                       (&fDirector,"Photons_pfChargedIsoRhoCorr"),
      Photons_pfGammaIso                                                (&fDirector,"Photons_pfGammaIso"),
      Photons_pfGammaIsoRhoCorr                                         (&fDirector,"Photons_pfGammaIsoRhoCorr"),
      Photons_pfNeutralIso                                              (&fDirector,"Photons_pfNeutralIso"),
      Photons_pfNeutralIsoRhoCorr                                       (&fDirector,"Photons_pfNeutralIsoRhoCorr"),
      Photons_sigmaIetaIeta                                             (&fDirector,"Photons_sigmaIetaIeta"),
      PrimaryVertexFilter                                               (&fDirector,"PrimaryVertexFilter"),
      PSweights                                                         (&fDirector,"PSweights"),
      PuppiMET                                                          (&fDirector,"PuppiMET"),
      PuppiMETDown                                                      (&fDirector,"PuppiMETDown"),
      PuppiMETPhi                                                       (&fDirector,"PuppiMETPhi"),
      PuppiMETPhiDown                                                   (&fDirector,"PuppiMETPhiDown"),
      PuppiMETPhiUp                                                     (&fDirector,"PuppiMETPhiUp"),
      PuppiMETUp                                                        (&fDirector,"PuppiMETUp"),
      puSysDown                                                         (&fDirector,"puSysDown"),
      puSysUp                                                           (&fDirector,"puSysUp"),
      puWeight                                                          (&fDirector,"puWeight"),
      ScaleWeights                                                      (&fDirector,"ScaleWeights"),
      SignalParameters                                                  (&fDirector,"SignalParameters"),
      TAPElectronTracks                                                 (&fDirector,"TAPElectronTracks"),
      TAPElectronTracks_dxypv                                           (&fDirector,"TAPElectronTracks_dxypv"),
      TAPElectronTracks_leptonMatch                                     (&fDirector,"TAPElectronTracks_leptonMatch"),
      TAPElectronTracks_mT                                              (&fDirector,"TAPElectronTracks_mT"),
      TAPElectronTracks_pfRelIso03chg                                   (&fDirector,"TAPElectronTracks_pfRelIso03chg"),
      TAPElectronTracks_trkiso                                          (&fDirector,"TAPElectronTracks_trkiso"),
      TAPMuonTracks                                                     (&fDirector,"TAPMuonTracks"),
      TAPMuonTracks_dxypv                                               (&fDirector,"TAPMuonTracks_dxypv"),
      TAPMuonTracks_leptonMatch                                         (&fDirector,"TAPMuonTracks_leptonMatch"),
      TAPMuonTracks_mT                                                  (&fDirector,"TAPMuonTracks_mT"),
      TAPMuonTracks_pfRelIso03chg                                       (&fDirector,"TAPMuonTracks_pfRelIso03chg"),
      TAPMuonTracks_trkiso                                              (&fDirector,"TAPMuonTracks_trkiso"),
      TAPPionTracks                                                     (&fDirector,"TAPPionTracks"),
      TAPPionTracks_dxypv                                               (&fDirector,"TAPPionTracks_dxypv"),
      TAPPionTracks_leptonMatch                                         (&fDirector,"TAPPionTracks_leptonMatch"),
      TAPPionTracks_mT                                                  (&fDirector,"TAPPionTracks_mT"),
      TAPPionTracks_pfRelIso03chg                                       (&fDirector,"TAPPionTracks_pfRelIso03chg"),
      TAPPionTracks_trkiso                                              (&fDirector,"TAPPionTracks_trkiso"),
      TriggerPass                                                       (&fDirector,"TriggerPass"),
      TriggerPrescales                                                  (&fDirector,"TriggerPrescales"),
      TriggerVersion                                                    (&fDirector,"TriggerVersion"),
      TrueNumInteractions                                               (&fDirector,"TrueNumInteractions"),
      Weight                                                            (&fDirector,"Weight")
      { }

   ~NtupleVarsTProxy();
   Int_t   Version() const {return 1;}
   void    Init(::TTree *tree);

   // Functions used in analysis
   double DeltaPhi(double, double);
   double DeltaR(double eta1, double phi1, double eta2, double phi2);
   void   sortTLorVec(vector<TLorentzVector> *);
   double TransMass(double phi1, double phi2, double pt1, double pt2);
   double MinDr(TLorentzVector v1,vector<TLorentzVector> v2);
   double MinDr2(vector<TLorentzVector> v1,TLorentzVector v2);
   double getCrossSection(std::string process_name);
   //ClassDef(NtupleVarsTProxy,0)

   //ClassDef(NtupleVarsTProxy,0);
};  //ends class NtupleVarsTProxy ...  

#endif // NtupleVarsTProxy_h

//================ 
/*
#ifdef __MAKECINT__
#pragma link C++ class NtupleVarsTProxy::TStlPx_ROOT__Math__PtEtaPhiE4D_float_-;
#pragma link C++ class NtupleVarsTProxy::TStlPx_ROOT__Math__LorentzVector_ROOT__Math__PtEtaPhiE4D_float__-;
#pragma link C++ class vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > >;
#pragma link C++ class NtupleVarsTProxy;
#endif
*/

//== == == == == == == == == == == ==  == == == == == == == == == == == == == 
// Functions related to Ntuple reading

#ifdef NtupleVarsTProxy_cxx

NtupleVarsTProxy::~NtupleVarsTProxy() {
   // destructor. Clean up helpers.
}

void NtupleVarsTProxy::Init(TTree *tree)
{
  if (tree == 0) return;
  std::cout << "NtupleVarsTProxy::Init(TTree *tree) tree->GetEntries()" << tree->GetEntries() << std::endl;
  // Set branch addresses
  fChain = tree;
  fDirector.SetTree(fChain);
  std::cout << "NtupleVarsTProxy::Init(TTree *tree) fChain->GetEntries() " << fChain->GetEntries() << std::endl;
}

#endif // #ifdef NtupleVarsTProxy_cxx
