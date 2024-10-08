#!/bin/bash
./analyzeTProxytBSM runFiles/runList_TTGJets_Tune_Summer20UL18.txt Summer20UL18_TTGJets_Tune.root 2018 TTG LL MVA
./analyzeTProxytBSM runFiles/runList_ZNuNuGJets_MonoPhoton_Summer20UL18.txt Summer20UL18_ZNuNuGJets_MonoPhoton.root 2018 ZG LL MVA
./analyzeTProxytBSM runFiles/runList_ZJetsToNuNu_HT_Summer20UL18.txt Summer20UL18_ZJetsToNuNu.root 2018 ZNuNuJets LL MVA
./analyzeTProxytBSM runFiles/runList_WJetsToLNu_HT_Summer20UL18.txt Summer20UL18_WJetsToLNu_HT.root 2018 WJets LL MVA
./analyzeTProxytBSM runFiles/runList_QCD_HT_Summer20UL18.txt Summer20UL18_QCD_HT.root 2018 QCD LL MVA
./analyzeTProxytBSM runFiles/runList_GJets_DR-0p4_HT_Summer20UL18.txt Summer20UL18_GJets_DR-0p4_HT.root 2018 GJets LL MVA

./analyzeTProxytBSM WGJets_40to130.txt Summer20UL18_WGJets_MonoPhoton_40to130.root 2018 WGJets_MonoPhoton_PtG-40to130UL LL MVA
./analyzeTProxytBSM WGJets_130.txt Summer20UL18_WGJets_MonoPhoton_130.root 2018 WGJets_MonoPhoton_PtG-130UL LL MVA
hadd -f Summer20UL18_WGJets_MonoPhoton.root Summer20UL18_WGJets_MonoPhoton_40to130.root Summer20UL18_WGJets_MonoPhoton_130.root

./analyzeTProxytBSM runFiles/runList_TTJets_HT_Summer20UL18.txt Summer20UL18_TTJets_HT.root 2018 TTJets LL MVA
./analyzeTProxytBSM runFiles/runList_TTJets_Leptons_Summer20UL18.txt Summer20UL18_TTJets_Leptons.root 2018 TTJets_Leptons LL MVA
hadd -f Summer20UL18_TTJets.root Summer20UL18_TTJets_HT.root Summer20UL18_TTJets_Leptons.root
