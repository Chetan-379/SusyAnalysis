import os,sys,re,fileinput,string,shutil
from datetime import date
min_count=0
year=["Summer20UL16APV","Summer20UL16","Summer20UL17","Summer20UL18"]
#year=["Summer20UL18"]

bkg_samples=["TTJets_HT","TTJets_Leptons","TTGJets_Tune","WJetsToLNu_HT","WGJets_MonoPhoton_40to130","WGJets_MonoPhoton_130Inf"] #,"ZJetsToNuNu_HT","ZNuNuGJets_MonoPhoton","QCD_HT", "GJets_DR-0p4_HT"]
#bkg_samples=["ZJetsToNuNu_HT","ZNuNuGJets_MonoPhoton","QCD_HT","GJets_DR-0p4_HT"]
count=0;
for i in year:
    for j in bkg_samples:
        condorSubmit = "condor_submit Submit_condor/submitCondor_%s_%s_Sbin"%(i,j)
        print(condorSubmit)
        os.system(condorSubmit)
