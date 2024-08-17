#!/bin/bash                                                                                                                                    

path=/eos/user/s/seema/SUSYGMSB2023/SkimsUL_June2023 
for year in 16APV 16 17 18
do
    ls ${path}/skimmed_Summer20UL${year}_TTJets_HT-*.root > runList_TTJets_HT_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_TTJets_*Lep*.root > runList_TTJets_Leptons_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_WJetsToLNu_HT*.root > runList_WJetsToLNu_HT_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_WGJets_MonoPhoton*.root > runList_WGJets_MonoPhoton_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_ZJetsToNuNu_HT-*.root > runList_ZJetsToNuNu_HT_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_ZNuNuGJets_MonoPhoton*.root > runList_ZNuNuGJets_MonoPhoton_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_TTGJets_TuneC*.root > runList_TTGJets_Tune_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_QCD_HT*.root > runList_QCD_HT_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_GJets_DR-0p4_HT-*.root > runList_GJets_DR-0p4_HT_Summer20UL${year}.txt
	
done

 #ls ${path}/skimmed_${year}_TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8.root > runList_TTJets_${year}.txt

 #ls ${path}/SR_${year}.TTJets_HT*.root > runList_TTJets-HT_${year}_v18.txt
# ls ${path}/SR_${year}.TTJets_*Lept*.root > runList_TTJets-inc_${year}_v18.txt
# ls ${path}/SR_${year}.WGJets_MonoPhoton_PtG-*_TuneCUETP8M1_13TeV-madgraph_v18.root > runList_WGJets_${year}_v18.txt
# ls ${path}/SR_${year}.WJetsToLNu_HT-*.root > runList_WJets_${year}_v18.txt
# ls ${path}/SR_${year}.ZJetsToNuNu_HT*.root > runList_ZJets_${year}_v18.txt
# ls ${path}/SR_${year}.ZNuNuGJets_MonoPhoton_PtG-*.root > runList_ZGJets_${year}_v18.txt
# ls ${path}/SR_${year}.QCD_HT*.root > runList_QCD_${year}_v18.txt
# ls ${path}/SR_${year}.GJets_DR*.root > runList_GJets_DR_${year}_v18.txt
