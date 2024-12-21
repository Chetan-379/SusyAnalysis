root -l -q 'splitRunList.C("inputfiles/runList_UL2018_Run2018A_v1_EGamma.txt",20,"2018A","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2018_Run2018B_v1_EGamma.txt",20,"2018B","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2018_Run2018C_v1_EGamma.txt",20,"2018C","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2018_Run2018D_v2_EGamma.txt",20,"2018D","dataUL","loose")'


root -l -q 'splitRunList.C("inputfiles/runList_UL2017_Run2017B_UL2017_v1_SingleElectron.txt",20,"2017B","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2017_Run2017C_UL2017_v1_SingleElectron.txt",20,"2017C","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2017_Run2017D_UL2017_v1_SingleElectron.txt",20,"2017D","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2017_Run2017E_UL2017_v1_SingleElectron.txt",20,"2017E","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2017_Run2017F_UL2017_v1_SingleElectron.txt",20,"2017F","dataUL","loose")'

root -l -q 'splitRunList.C("inputfiles/runList_UL2016APV_Run2016B_v1_SingleElectron.txt",20,"2016preVFPB","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2016APV_Run2016C_v1_SingleElectron.txt",20,"2016preVFPC","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2016APV_Run2016D_v1_SingleElectron.txt",20,"2016preVFPD","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2016APV_Run2016E_v1_SingleElectron.txt",20,"2016preVFPE","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2016APV_Run2016F_v1_SingleElectron.txt",20,"2016preVFPF","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2016_Run2016F_v1_SingleEelctron.txt",20,"2016postVFPF","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2016_Run2016G_v1_SingleEelctron.txt",20,"2016postVFPG","dataUL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_UL2016_Run2016H_v1_SingleEelctron.txt",20,"2016postVFPH","dataUL","loose")'

## ZLL G samples
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL18_ZLLGJets_MonoPhoton_PtG-15to130_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",25,"2018","ZLLGJets_MonoPhoton_PtG-15to130UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL17_ZLLGJets_MonoPhoton_PtG-15to130_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",25,"2017","ZLLGJets_MonoPhoton_PtG-15to130UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16_ZLLGJets_MonoPhoton_PtG-15to130_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",25,"2016postVFP","ZLLGJets_MonoPhoton_PtG-15to130UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL18_ZLLGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",25,"2018","ZLLGJets_MonoPhoton_PtG-130UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL17_ZLLGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",25,"2017","ZLLGJets_MonoPhoton_PtG-130UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16_ZLLGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",25,"2016postVFP","ZLLGJets_MonoPhoton_PtG-130UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16APV_ZLLGJets_MonoPhoton_PtG-130_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",25,"2016preVFP","ZLLGJets_MonoPhoton_PtG-130UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16APV_ZLLGJets_MonoPhoton_PtG-15to130_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",25,"2016preVFP","ZLLGJets_MonoPhoton_PtG-15to130UL","loose")'


# # # #### DY samples
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL17_DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2017","DYJetsToLL_M-50_HT-100to200UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL17_DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2017","DYJetsToLL_M-50_HT-200to400UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL17_DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2017","DYJetsToLL_M-50_HT-400to600UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL17_DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2017","DYJetsToLL_M-50_HT-600to800UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL17_DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2017","DYJetsToLL_M-50_HT-800to1200UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL17_DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2017","DYJetsToLL_M-50_HT-1200to2500UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL17_DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2017","DYJetsToLL_M-50_HT-2500toInfUL","loose")'

         
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL18_DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2018","DYJetsToLL_M-50_HT-100to200UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL18_DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2018","DYJetsToLL_M-50_HT-200to400UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL18_DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2018","DYJetsToLL_M-50_HT-400to600UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL18_DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2018","DYJetsToLL_M-50_HT-600to800UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL18_DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2018","DYJetsToLL_M-50_HT-800to1200UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL18_DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2018","DYJetsToLL_M-50_HT-1200to2500UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL18_DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2018","DYJetsToLL_M-50_HT-2500toInfUL","loose")'


         
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16_DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016postVFP","DYJetsToLL_M-50_HT-100to200UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16_DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016postVFP","DYJetsToLL_M-50_HT-200to400UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16_DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016postVFP","DYJetsToLL_M-50_HT-400to600UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16_DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016postVFP","DYJetsToLL_M-50_HT-600to800UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16_DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016postVFP","DYJetsToLL_M-50_HT-800to1200UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16_DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016postVFP","DYJetsToLL_M-50_HT-1200to2500UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16_DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016postVFP","DYJetsToLL_M-50_HT-2500toInfUL","loose")'

         
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16APV_DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016preVFP","DYJetsToLL_M-50_HT-100to200UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16APV_DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016preVFP","DYJetsToLL_M-50_HT-200to400UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16APV_DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016preVFP","DYJetsToLL_M-50_HT-400to600UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16APV_DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016preVFP","DYJetsToLL_M-50_HT-600to800UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16APV_DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016preVFP","DYJetsToLL_M-50_HT-800to1200UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16APV_DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016preVFP","DYJetsToLL_M-50_HT-1200to2500UL","loose")'
root -l -q 'splitRunList.C("inputfiles/runList_Summer20UL16APV_DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",25,"2016preVFP","DYJetsToLL_M-50_HT-2500toInfUL","loose")'



















# # root -l -q 'splitRunList.C("runList_UL2018_Run2018A_v2_MET.txt",20,"2018A","dataUL","loose")'
# # root -l -q 'splitRunList.C("runList_UL2018_Run2018B_v2_MET.txt",20,"2018B","dataUL","loose")'
# # root -l  -q 'splitRunList.C("runList_UL2018_Run2018C_v1_MET.txt",20,"2018C","dataUL","loose")'

# # root -l  -q 'splitRunList.C("runList_UL2018_Run2018D_v1_MET.txt",20,"2018D","dataUL","loose")'

# # root -l -q 'splitRunList.C("runList_UL2017_Run2017B_v1_MET.txt",20,"2017B","dataUL","loose")'
# # root -l -q 'splitRunList.C("runList_UL2017_Run2017C_v1_MET.txt",20,"2017C","dataUL","loose")'                                                             
# # root -l -q 'splitRunList.C("runList_UL2017_Run2017D_v1_MET.txt",20,"2017D","dataUL","loose")'                                                             
# # root -l -q 'splitRunList.C("runList_UL2017_Run2017E_v1_MET.txt",20,"2017E","dataUL","loose")'                                                             
# # root -l -q 'splitRunList.C("runList_UL2017_Run2017F_v1_MET.txt",20,"2017F","dataUL","loose")'                                                             

# # root -l -q 'splitRunList.C("runList_UL2016_Run2016F-UL2016-v2_MET.txt",1,"2016postVFPF","dataUL","loose")'
# # root -l -q 'splitRunList.C("runList_UL2016_Run2016G-UL2016-v2_MET.txt",1,"2016postVFPG","dataUL","loose")'
# # root -l -q 'splitRunList.C("runList_UL2016_Run2016H-UL2016-v2_MET.txt",1,"2016postVFPH","dataUL","loose")'


# # root -l -q 'splitRunList.C("runList_UL2016APV_Run2016B-UL2016_HIPM-ver2-v2_MET.txt",20,"2016preVFPB","dataUL","loose")'
# # root -l -q 'splitRunList.C("runList_UL2016APV_Run2016C-UL2016_HIPM-v2_MET.txt",20,"2016preVFPC","dataUL","loose")'
# # root -l -q 'splitRunList.C("runList_UL2016APV_Run2016D-UL2016_HIPM-v2_MET.txt",20,"2016preVFPD","dataUL","loose")'
# # root -l -q 'splitRunList.C("runList_UL2016APV_Run2016E-UL2016_HIPM-v2_MET.txt",20,"2016preVFPE","dataUL","loose")'
# # root -l -q 'splitRunList.C("runList_UL2016APV_Run2016F-UL2016_HIPM-v2_MET.txt",1,"2016preVFPF","dataUL","loose")'










# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL18_WGJets_MonoPhoton_PtG-40to130_v1.txt",1,"2018","WGJets_MonoPhoton_PtG-40to130UL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL18_WGJets_MonoPhoton_PtG-130_v1.txt",1,"2018","WGJets_MonoPhoton_PtG-130UL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL16_WGJets_MonoPhoton_PtG-40to130_v1.txt",1,"2016postVFP","WGJets_MonoPhoton_PtG-40to130UL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL16_WGJets_MonoPhoton_PtG-130_v1.txt",1,"2016postVFP","WGJets_MonoPhoton_PtG-130UL","loose")'

# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL17_WGJets_MonoPhoton_PtG-40to130_v1.txt",1,"2017","WGJets_MonoPhoton_PtG-40to130UL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL17_WGJets_MonoPhoton_PtG-130_v1.txt",1,"2017","WGJets_MonoPhoton_PtG-130UL","loose")'

# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL18_WJetsToLNu_HT_v1.txt",1,"2018","WJetsToLNu_HTUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL17_WJetsToLNu_HT_v1.txt",1,"2017","WJetsToLNu_HTUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL16_WJetsToLNu_HT_v1.txt",1,"2016postVFP","WJetsToLNu_HTUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL18_TTJets_inc_v1.txt",1,"2018","TTJets_incUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL17_TTJets_inc_v1.txt",1,"2017","TTJets_incUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL16_TTJets_inc_v1.txt",1,"2016postVFP","TTJets_incUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL18_TTJets_HT_v1.txt",1,"2018","TTJets_HTUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL17_TTJets_HT_v1.txt",1,"2017","TTJets_HTUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL16_TTJets_HT_v1.txt",1,"2016postVFP","TTJets_HTUL","loose")'

# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL18_TTGJets_inc_v1.txt",1,"2018","TTGJets_incUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL17_TTGJets_inc_v1.txt",1,"2017","TTGJets_incUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL16_TTGJets_inc_v1.txt",1,"2016postVFP","TTGJets_incUL","loose")'

# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL16APV_TTGJets_inc_v1.txt",1,"2016preVFP","TTGJets_incUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL16APV_TTJets_HT_v1.txt",1,"2016preVFP","TTJets_HTUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL16APV_TTJets_inc_v1.txt",1,"2016preVFP","TTJets_incUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL16APV_WJetsToLNu_HT_v1.txt",1,"2016preVFP","WJetsToLNu_HTUL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL16APV_WGJets_MonoPhoton_PtG-40to130_v1.txt",1,"2016preVFP","WGJets_MonoPhoton_PtG-40to130UL","loose")'
# # root -l -q 'splitRunList.C("runList_skimmed_Summer20UL16APV_WGJets_MonoPhoton_PtG-130_v1.txt",1,"2016preVFP","WGJets_MonoPhoton_PtG-130UL","loose")'






