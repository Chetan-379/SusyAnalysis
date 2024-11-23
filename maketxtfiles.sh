#!/bin/bash                                                                                                                                    

path=/eos/user/c/cagrawal/Samples/FR 
for year in 16APV 16 17 18
do
    ls ${path}/skimmed_Summer20UL${year}_TTJets_HT-*.root > runFiles/runList_TTJets_HT_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_TTJets_*Lep*.root > runFiles/runList_TTJets_Leptons_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_WJetsToLNu_HT*.root > runFiles/runList_WJetsToLNu_HT_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_TTGJets_TuneC*.root > runFiles/runList_TTGJets_Tune_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_WGJets_MonoPhoton_PtG-40to130*.root > runFiles/runList_WGJets_MonoPhoton_40to130_Summer20UL${year}.txt
    ls ${path}/skimmed_Summer20UL${year}_WGJets_MonoPhoton_PtG-130*.root > runFiles/runList_WGJets_MonoPhoton_130Inf_Summer20UL${year}.txt    
	
done

