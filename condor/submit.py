import os,sys,re,fileinput,string,shutil
from datetime import date
min_count=0
# year=["Fall17","Autumn18","Summer16v3"]
# number=["2017","2018","2016"]

year=["Summer20UL16APV","Summer20UL16","Summer20UL17","Summer20UL16"]
number=["2016APV","2016","2017","2018"]

#sample_tag=["GJets_DR","QCD_Jets","temp","temp","temp","temp","temp","temp"]                                                                                                
#bkg_samples=["GJets_DR","QCD","ZGJets","ZJets","TTGJets","TTJets-HT","WGJets","WJets","TTJets-inc"]
bkg_samples=["TTJets_HT","TTJets_Leptons","TTGJets_Tune","WJetsToLNu_HT","WGJets_MonoPhoton","ZJetsToNuNu_HT","ZNuNuGJets_MonoPhoton","QCD_HT"]

#bkg_samples=["ZGJets","ZJets","TTGJets","TTJets-HT","WGJets","WJets"]                                                                                                       
#bkg_samples=["ZGJets","GJets","ZJets","QCD"]                                                                                                                                
#bkg_samples=["TTJets-HT","TTJets-inc"]                                                                                                                                      
count=0;
for i in year:
    for j in bkg_samples:
        condorSubmit = "condor_submit condor_submit/submitCondor_%s_%s_Sbin"%(i,j)
        print condorSubmit
        os.system(condorSubmit)
