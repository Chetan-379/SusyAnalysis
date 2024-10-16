import os,sys,re,fileinput,string,shutil
from datetime import date
count=0
min_count=0
year=["Summer20UL16APV","Summer20UL16","Summer20UL17","Summer20UL18"]
number=["2016preVFP","2016postVFP","2017","2018"]
sample_tag=["TTJets_HT","TTJets_Leptons","TTG","WJets","WGJets_MonoPhoton_PtG-40to130UL","WGJets_MonoPhoton_PtG-130UL","ZNuNuJets","ZG","QCD","GJets"]
bkg_samples=["TTJets_HT","TTJets_Leptons","TTGJets_Tune","WJetsToLNu_HT","WGJets_MonoPhoton_40to130","WGJets_MonoPhoton_130Inf","ZJetsToNuNu_HT","ZNuNuGJets_MonoPhoton","QCD_HT","GJets_DR-0p4_HT"]
                                                                     
count1=0;
for i in year:
    count=0
    for j in bkg_samples:
        condorSubmit = "Submit_condor/submitCondor_%s_%s_Sbin"%(i,j)
        fname1 = "/afs/cern.ch/user/c/cagrawal/public/MyAnalysis/CMSSW_14_0_2/src/SusySoft2023Ana/runFiles/runList_%s_%s.txt"%(j,i)
        fname = "/afs/cern.ch/user/c/cagrawal/public/MyAnalysis/CMSSW_14_0_2/src/SusySoft2023Ana/condor/root_files_condor/%s_%s.root"%(i,j)
        jobName = "%s_%s"%(i,j)
        #out_tag="%s_%s" %(sample_tag[count1],year[count1])                                                             
        shutil.copyfile("proto_condor_submit",condorSubmit)
        for line in fileinput.FileInput(condorSubmit, inplace=1):
            line=line.replace("MIN", str(i))
            line=line.replace("MAX", str(j))
            line=line.replace("outfile", fname)
            line=line.replace("filelist",fname1)
            line=line.replace("year",number[count1])
            line=line.replace("sample",sample_tag[count])
            print(line.rstrip())
            
        print(condorSubmit)
        count+=1
    count1+=1
            
