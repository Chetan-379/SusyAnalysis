#!/bin/bash                                                                                                                                      
path=/store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp
# # hadd -f out_Data_UL2018_Run2018D_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2018_Run2018D_v1_MET_pt100_MET200'`
# # hadd -f out_Data_UL2018_Run2018A_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2018_Run2018A_v2_MET_pt100_MET200'`
# # hadd -f out_Data_UL2018_Run2018B_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2018_Run2018B_v2_MET_pt100_MET200'`
# # hadd -f out_Data_UL2018_Run2018C_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2018_Run2018C_v1_MET_pt100_MET200'`

# # hadd -f out_Data_UL2017_Run2017B_MET_phoID_loose_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2017_Run2017B_v1_MET_pt100_MET200'`
# # hadd -f out_Data_UL2017_Run2017C_MET_phoID_loose_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2017_Run2017C_v1_MET_pt100_MET200'`
# # hadd -f out_Data_UL2017_Run2017D_MET_phoID_loose_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2017_Run2017D_v1_MET_pt100_MET200'`
# # hadd -f out_Data_UL2017_Run2017E_MET_phoID_loose_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2017_Run2017E_v1_MET_pt100_MET200'`
# # hadd -f out_Data_UL2017_Run2017F_MET_phoID_loose_pt40_MET200.root  `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2017_Run2017F_v1_MET_pt100_MET200'`
# # hadd -f out_Data_UL2016APV_Run2016B_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016APV_Run2016B-UL2016_HIPM-ver2-v2_MET_pt100_MET200'`
# # hadd -f out_Data_UL2016APV_Run2016C_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016APV_Run2016C-UL2016_HIPM-v2_MET_pt100_MET200'`
# # hadd -f out_Data_UL2016APV_Run2016D_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016APV_Run2016D-UL2016_HIPM-v2_MET_pt100_MET200'`
# # hadd -f out_Data_UL2016APV_Run2016E_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016APV_Run2016E-UL2016_HIPM-v2_MET_pt100_MET200'`
# # hadd -f out_Data_UL2016APV_Run2016F_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016APV_Run2016F-UL2016_HIPM-v2_MET_pt100_MET200'`

# # hadd -f out_Data_UL2016_Run2016F_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016_Run2016F-UL2016-v2_MET_pt100_MET200'`
# # hadd -f out_Data_UL2016_Run2016G_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016_Run2016G-UL2016-v2_MET_pt100_MET200'`
# # hadd -f out_Data_UL2016_Run2016H_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016_Run2016H-UL2016-v2_MET_pt100_MET200'`



hadd -f Summer20UL16APV_TTGJets_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16APV_TTGJets'`
hadd -f Summer20UL16APV_WGJets_MonoPhoton_PtG-40to130_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16APV_WGJets_MonoPhoton_PtG-40to130'`
hadd -f Summer20UL16APV_WGJets_MonoPhoton_PtG-130_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16APV_WGJets_MonoPhoton_PtG-130'`




# hadd -f Summer20UL17_TTGJets_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL17_TTGJets'`
hadd -f Summer20UL16_TTGJets_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16_TTGJets'`
# hadd -f Summer20UL18_TTGJets_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL18_TTGJets'`



# hadd -f Summer20UL17_WGJets_MonoPhoton_PtG-40to130_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL17_WGJets_MonoPhoton_PtG-40to130'`
hadd -f Summer20UL16_WGJets_MonoPhoton_PtG-40to130_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16_WGJets_MonoPhoton_PtG-40to130'`
# hadd -f Summer20UL18_WGJets_MonoPhoton_PtG-40to130_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL18_WGJets_MonoPhoton_PtG-40to130'`
# hadd -f Summer20UL17_WGJets_MonoPhoton_PtG-130_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL17_WGJets_MonoPhoton_PtG-130'`
hadd -f Summer20UL16_WGJets_MonoPhoton_PtG-130_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16_WGJets_MonoPhoton_PtG-130'`
# hadd -f Summer20UL18_WGJets_MonoPhoton_PtG-130_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL18_WGJets_MonoPhoton_PtG-130'`


# hadd -f Summer20UL17_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL17_WJetsToLNu_HT'`
hadd -f Summer20UL16_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16_WJetsToLNu_HT'`
# hadd -f Summer20UL18_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL18_WJetsToLNu_HT'`
hadd -f Summer20UL16APV_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16APV_WJetsToLNu_HT'`
# hadd -f Summer20UL17_TTJets_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL17_TTJets'`
hadd -f Summer20UL16_TTJets_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16_TTJets'`
# hadd -f Summer20UL18_TTJets_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL18_TTJets'`
hadd -f Summer20UL16APV_TTJets_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16APV_TTJets'`


# # hadd -f Summer20UL18_singleTop_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL18_ST'`
# # hadd -f Summer20UL17_singleTop_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL17_ST'`
# # hadd -f Summer20UL16_singleTop_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16_ST'`
# # hadd -f Summer20UL16APV_singleTop_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16APV_ST'`

# hadd -f Summer20UL18_ZLLGJets_MonoPhoton_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL18_ZLLGJets_MonoPhoton_PtG'`
# hadd -f Summer20UL17_ZLLGJets_MonoPhoton_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL17_ZLLGJets_MonoPhoton_PtG'`
# hadd -f Summer20UL16APV_ZLLGJets_MonoPhoton_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16APV_ZLLGJets_MonoPhoton_PtG'`
# hadd -f Summer20UL16_ZLLGJets_MonoPhoton_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16_ZLLGJets_MonoPhoton_PtG'`


# hadd -f Summer20UL18_DYJets_Mt50_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL18_DYJetsToLL_M-50_HT'`
# hadd -f Summer20UL17_DYJets_Mt50_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL17_DYJetsToLL_M-50_HT'`
# hadd -f Summer20UL16_DYJets_Mt50_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16_DYJetsToLL_M-50_HT'`
# hadd -f Summer20UL16APV_DYJets_Mt50_PhoIdloose_phopt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_Summer20UL16APV_DYJetsToLL_M-50_HT'`



# hadd -f out_Data_UL2018_Run2018A_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2018_Run2018A'`

# hadd -f out_Data_UL2018_Run2018B_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2018_Run2018B'`

# hadd -f out_Data_UL2018_Run2018C_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2018_Run2018C'`

# hadd -f out_Data_UL2018_Run2018D_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2018_Run2018D'`

# hadd -f out_Data_UL2017_Run2017B_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2017_Run2017B'`
# hadd -f out_Data_UL2017_Run2017C_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2017_Run2017C'`
# hadd -f out_Data_UL2017_Run2017D_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2017_Run2017D'`
# hadd -f out_Data_UL2017_Run2017E_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2017_Run2017E'`
# hadd -f out_Data_UL2017_Run2017F_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2017_Run2017F'`

# hadd -f out_Data_UL2016APV_Run2016B_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016APV_Run2016B'`
# hadd -f out_Data_UL2016APV_Run2016C_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016APV_Run2016C'`
# hadd -f out_Data_UL2016APV_Run2016D_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016APV_Run2016D'`
# hadd -f out_Data_UL2016APV_Run2016E_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016APV_Run2016E'`
# hadd -f out_Data_UL2016APV_Run2016F_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016APV_Run2016F'`

# hadd -f out_Data_UL2016_Run2016F_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016_Run2016F'`
# hadd -f out_Data_UL2016_Run2016G_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016_Run2016G'`
# hadd -f out_Data_UL2016_Run2016H_MET_phoID_loose_pt40_MET200.root `xrdfsls -u /store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp | grep 'phoID_loose_runList_UL2016_Run2016H'`





















