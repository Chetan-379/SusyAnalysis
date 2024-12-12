#!/bin/bash

path="/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/root_files"
hadd -f $path/Summer20UL16APV_TTJets.root $path/Summer20UL16APV_TTJets_HT.root $path/Summer20UL16APV_TTJets_Leptons.root
hadd -f $path/Summer20UL16APV_WGJets_MonoPhoton.root $path/Summer20UL16APV_WGJets_MonoPhoton_40to130.root $path/Summer20UL16APV_WGJets_MonoPhoton_130Inf.root
hadd -f $path/Summer20UL16APV_TTGJets_process_combined.root $path/Summer20UL16APV_TTJets.root $path/Summer20UL16APV_TTGJets_Tune.root
hadd -f $path/Summer20UL16APV_WGJets_process_combined.root $path/Summer20UL16APV_WJetsToLNu_HT.root $path/Summer20UL16APV_WGJets_MonoPhoton.root
hadd -f $path/Summer20UL16APV_All_bkg_Combined.root $path/Summer20UL16APV_TTJets.root $path/Summer20UL16APV_TTGJets_Tune.root $path/Summer20UL16APV_WJetsToLNu_HT.root $path/Summer20UL16APV_WGJets_MonoPhoton.root

hadd -f $path/Summer20UL16_TTJets.root $path/Summer20UL16_TTJets_HT.root $path/Summer20UL16_TTJets_Leptons.root
hadd -f $path/Summer20UL16_WGJets_MonoPhoton.root $path/Summer20UL16_WGJets_MonoPhoton_40to130.root $path/Summer20UL16_WGJets_MonoPhoton_130Inf.root
hadd -f $path/Summer20UL16_TTGJets_process_combined.root $path/Summer20UL16_TTJets.root $path/Summer20UL16_TTGJets_Tune.root
hadd -f $path/Summer20UL16_WGJets_process_combined.root $path/Summer20UL16_WJetsToLNu_HT.root $path/Summer20UL16_WGJets_MonoPhoton.root
hadd -f $path/Summer20UL16_All_bkg_Combined.root $path/Summer20UL16_TTJets.root $path/Summer20UL16_TTGJets_Tune.root $path/Summer20UL16_WJetsToLNu_HT.root $path/Summer20UL16_WGJets_MonoPhoton.root

hadd -f $path/Summer20UL17_TTJets.root $path/Summer20UL17_TTJets_HT.root $path/Summer20UL17_TTJets_Leptons.root
hadd -f $path/Summer20UL17_WGJets_MonoPhoton.root $path/Summer20UL17_WGJets_MonoPhoton_40to130.root $path/Summer20UL17_WGJets_MonoPhoton_130Inf.root
hadd -f $path/Summer20UL17_TTGJets_process_combined.root $path/Summer20UL17_TTJets.root $path/Summer20UL17_TTGJets_Tune.root
hadd -f $path/Summer20UL17_WGJets_process_combined.root $path/Summer20UL17_WJetsToLNu_HT.root $path/Summer20UL17_WGJets_MonoPhoton.root
hadd -f $path/Summer20UL17_All_bkg_Combined.root $path/Summer20UL17_TTJets.root $path/Summer20UL17_TTGJets_Tune.root $path/Summer20UL17_WJetsToLNu_HT.root $path/Summer20UL17_WGJets_MonoPhoton.root

hadd -f $path/Summer20UL18_TTJets.root $path/Summer20UL18_TTJets_HT.root $path/Summer20UL18_TTJets_Leptons.root
hadd -f $path/Summer20UL18_WGJets_MonoPhoton.root $path/Summer20UL18_WGJets_MonoPhoton_40to130.root $path/Summer20UL18_WGJets_MonoPhoton_130Inf.root
hadd -f $path/Summer20UL18_TTGJets_process_combined.root $path/Summer20UL18_TTJets.root $path/Summer20UL18_TTGJets_Tune.root
hadd -f $path/Summer20UL18_WGJets_process_combined.root $path/Summer20UL18_WJetsToLNu_HT.root $path/Summer20UL18_WGJets_MonoPhoton.root
hadd -f $path/Summer20UL18_All_bkg_Combined.root $path/Summer20UL18_TTJets.root $path/Summer20UL18_TTGJets_Tune.root $path/Summer20UL18_WJetsToLNu_HT.root $path/Summer20UL18_WGJets_MonoPhoton.root

hadd -f $path/Summer20ULFullRun2_TTJets.root $path/Summer20UL18_TTJets.root $path/Summer20UL17_TTJets.root $path/Summer20UL16_TTJets.root $path/Summer20UL16APV_TTJets.root
hadd -f $path/Summer20ULFullRun2_TTGJets_Tune.root $path/Summer20UL18_TTGJets_Tune.root $path/Summer20UL17_TTGJets_Tune.root $path/Summer20UL16_TTGJets_Tune.root $path/Summer20UL16APV_TTGJets_Tune.root
hadd -f $path/Summer20ULFullRun2_WJetsToLNu_HT.root $path/Summer20UL18_WJetsToLNu_HT.root $path/Summer20UL17_WJetsToLNu_HT.root $path/Summer20UL16_WJetsToLNu_HT.root $path/Summer20UL16APV_WJetsToLNu_HT.root
hadd -f $path/Summer20ULFullRun2_WGJets_MonoPhoton.root $path/Summer20UL18_WGJets_MonoPhoton.root $path/Summer20UL17_WGJets_MonoPhoton.root $path/Summer20UL16_WGJets_MonoPhoton.root $path/Summer20UL16APV_WGJets_MonoPhoton.root
hadd -f $path/Summer20ULFullRun2_All_bkg_Combined.root $path/Summer20UL18_All_bkg_Combined.root $path/Summer20UL17_All_bkg_Combined.root $path/Summer20UL16_All_bkg_Combined.root $path/Summer20UL16APV_All_bkg_Combined.root
