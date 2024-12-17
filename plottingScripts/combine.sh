# #!/bin/bash
path="/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/root_files_test"
#hadd -f $path/Summer20UL16APV_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_TTJets_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_TTJets_inc_PhoIdloose_phopt40_MET200.root

# hadd -f $path/Summer20UL_total2016_singleTop_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_singleTop_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_singleTop_PhoIdloose_phopt40_MET200.root
# hadd -f $path/FullRun2_singleTop_PhoIdloose_phopt40_MET200.root  $path/Summer20UL16APV_singleTop_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_singleTop_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_singleTop_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_singleTop_PhoIdloose_phopt40_MET200.root

hadd -f $path/Summer20UL16APV_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_WGJets_MonoPhoton_PtG-40to130_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_WGJets_MonoPhoton_PtG-130_PhoIdloose_phopt40_MET200.root


#hadd -f $path/Summer20UL16_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_TTJets_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_TTJets_inc_PhoIdloose_phopt40_MET200.root #checked

#hadd -f $path/Summer20UL17_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_TTJets_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_TTJets_inc_PhoIdloose_phopt40_MET200.root #checked
#hadd -f $path/Summer20UL18_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_TTJets_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_TTJets_inc_PhoIdloose_phopt40_MET200.root #checked

hadd -f $path/Summer20UL18_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_WGJets_MonoPhoton_PtG-40to130_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_WGJets_MonoPhoton_PtG-130_PhoIdloose_phopt40_MET200.root #checked 

hadd -f $path/Summer20UL17_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_WGJets_MonoPhoton_PtG-40to130_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_WGJets_MonoPhoton_PtG-130_PhoIdloose_phopt40_MET200.root #checked 

hadd -f $path/Summer20UL16_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_WGJets_MonoPhoton_PtG-40to130_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_WGJets_MonoPhoton_PtG-130_PhoIdloose_phopt40_MET200.root #checked 

hadd -f $path/FullRun2_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_TTJets_PhoIdloose_phopt40_MET200.root  #checked 

hadd -f $path/FullRun2_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_WGJets_PhoIdloose_phopt40_MET200.root #checked 

hadd -f $path/FullRun2_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root #checked 


hadd -f $path/Summer20UL18_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root  $path/Summer20UL18_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root #checked 

hadd -f $path/Summer20UL16_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root  $path/Summer20UL16_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root #checked 

hadd -f $path/Summer20UL16APV_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root  $path/Summer20UL16APV_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root

hadd -f $path/Summer20UL17_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root  $path/Summer20UL17_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root #checked 
hadd -f $path/Summer20UL18_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root  $path/Summer20UL18_TTGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_TTJets_PhoIdloose_phopt40_MET200.root #checked 

hadd -f $path/Summer20UL16_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root  $path/Summer20UL16_TTGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_TTJets_PhoIdloose_phopt40_MET200.root #checked 

hadd -f $path/Summer20UL16APV_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root  $path/Summer20UL16APV_TTGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_TTJets_PhoIdloose_phopt40_MET200.root

hadd -f $path/Summer20UL17_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root  $path/Summer20UL17_TTGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_TTJets_PhoIdloose_phopt40_MET200.root #checked 


hadd -f $path/FullRun2_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root #checked 

hadd -f	$path/FullRun2_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root #checked 

hadd -f $path/Summer20UL18_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_singleTop_PhoIdloose_phopt40_MET200.root #checked 


hadd -f	$path/Summer20UL17_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_singleTop_PhoIdloose_phopt40_MET200.root #checked 

hadd -f	$path/Summer20UL16_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_singleTop_PhoIdloose_phopt40_MET200.root #checked 

hadd -f $path/Summer20UL16APV_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_singleTop_PhoIdloose_phopt40_MET200.root #checked                                                                                                            

hadd -f $path/FullRun2_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root #checked $path/FullRun2_singleTop_PhoIdloose_phopt40_MET200.root #checked 


hadd -f $path/Summer20UL_total2016_TTGJets_PhoIdloose_phopt40_MET200.root  $path/Summer20UL16APV_TTGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_TTGJets_PhoIdloose_phopt40_MET200.root #checked

hadd -f $path/Summer20UL_total2016_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_TTJets_PhoIdloose_phopt40_MET200.root #checked
hadd -f $path/Summer20UL_total2016_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_WGJets_PhoIdloose_phopt40_MET200.root #checked

hadd -f $path/Summer20UL_total2016_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root  $path/Summer20UL16APV_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root #checked

hadd -f $path/Summer20UL_total2016_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL_total2016_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL_total2016_WGJets_PhoIdloose_phopt40_MET200.root #checked
hadd -f $path/Summer20UL_total2016_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL_total2016_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL_total2016_TTGJets_PhoIdloose_phopt40_MET200.root #checked
hadd -f $path/Summer20UL_total2016_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root  $path/Summer20UL_total2016_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL_total2016_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root  #$path/Summer20UL_total2016_singleTop_PhoIdloose_phopt40_MET200.root #checked


hadd -f $path/FullRun2_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL17_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root
#hadd -f $path/FullRun2_WGJets_PhoIdloose_phopt40_MET200.root   $path/Summer20UL17_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_WGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_WGJets_PhoIdloose_phopt40_MET200.root
hadd -f $path/FullRun2_TTGJets_PhoIdloose_phopt40_MET200.root       $path/Summer20UL17_TTGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_TTGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_TTGJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_TTGJets_PhoIdloose_phopt40_MET200.root
hadd -f $path/FullRun2_TTJets_PhoIdloose_phopt40_MET200.root        $path/Summer20UL17_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL18_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16_TTJets_PhoIdloose_phopt40_MET200.root $path/Summer20UL16APV_TTJets_PhoIdloose_phopt40_MET200.root


# hadd -f out_Data_UL2016_Allruns_MET_phoID_loose_pt40_MET200_Muon.root out_Data_UL2016_Run2016*_MET_phoID_loose_pt40_MET200_Muon.*
# hadd -f out_Data_UL2017_Allruns_MET_phoID_loose_pt40_MET200_Muon.root out_Data_UL2017_Run2017*_MET_phoID_loose_pt40_MET200_Muon.*
# hadd -f out_Data_UL2018_Allruns_MET_phoID_loose_pt40_MET200_Muon.root out_Data_UL2018_Run2018*_MET_phoID_loose_pt40_MET200_Muon.*

# hadd -f out_Data_UL2016_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2016_Run2016*_MET_phoID_loose_pt40_MET200.*
# hadd -f out_Data_UL2017_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2017_Run2017*_MET_phoID_loose_pt40_MET200.*
# hadd -f out_Data_UL2018_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2018_Run2018*_MET_phoID_loose_pt40_MET200.*
# hadd -f out_Data_UL2016APV_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2016APV_Run2016*_MET_phoID_loose_pt40_MET200.*

# hadd -f out_Data_UL2018_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root out_Data_UL2018_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2018_Allruns_MET_phoID_loose_pt40_MET200_Muon.root
# hadd -f out_Data_UL2017_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root out_Data_UL2017_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2017_Allruns_MET_phoID_loose_pt40_MET200_Muon.root
# hadd -f out_Data_UL2016_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root out_Data_UL2016_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2016_Allruns_MET_phoID_loose_pt40_MET200_Muon.root

# hadd -f out_Data_$path/FullRun2_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root out_Data_UL2018_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root out_Data_UL2017_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root out_Data_UL2016_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root
# hadd -f	out_Data_$path/FullRun2_Allruns_MET_phoID_loose_pt40_MET200_Muon.root out_Data_UL2018_Allruns_MET_phoID_loose_pt40_MET200_Muon.root out_Data_UL2017_Allruns_MET_phoID_loose_pt40_MET200_Muon.root out_Data_UL2016_Allruns_MET_phoID_loose_pt40_MET200_Muon.root

# hadd -f	out_Data_$path/FullRun2_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2018_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2017_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2016_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2016APV_Allruns_MET_phoID_loose_pt40_MET200.root

# hadd -f	out_Data_UL20_total2016_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2016APV_Allruns_MET_phoID_loose_pt40_MET200.root out_Data_UL2016_Allruns_MET_phoID_loose_pt40_MET200.root

# hadd -f out_Data_UL20_total2016_Allruns_MET_phoID_loose_pt40_MET200_Lepton.root out_Data_UL20_total2016_Allruns_MET_phoID_loose_pt40_MET200_Muon.root out_Data_UL20_total2016_Allruns_MET_phoID_loose_pt40_MET200.root


