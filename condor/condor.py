import os,sys,re,fileinput,string,shutil
from datetime import date
count=0
min_count=0
year=["Summer20UL16APV","Summer20UL16","Summer20UL17","Summer20UL18"]
number=["2016APV","2016","2017","2018"]
#sample_tag=["GJets_DR","QCD_Jets","temp","temp","temp","temp","temp","temp"]                                                                                                
#sample_tag=["GJets_DR","QCD","ZGJets","ZJets","TTGJets","TTJets-HT","WGJets","WJets","TTJets-inc"]
sample_tag=["TTJets_HT","TTJets_Leptons","TTGJets_Tune","WJetsToLNu_HT","WGJets_MonoPhoton","ZJetsToNuNu_HT","ZNuNuGJets_MonoPhoton","QCD_HT"]
bkg_samples=["TTJets_HT","TTJets_Leptons","TTGJets_Tune","WJetsToLNu_HT","WGJets_MonoPhoton","ZJetsToNuNu_HT","ZNuNuGJets_MonoPhoton","QCD_HT"]
#bkg_samples=["TTJets-HT","TTJets-inc"] #"GJets_DR","QCD"]                                                                                                                   
#sample_tag=["TTJets-HT","TTJets-inc"]                                                                                                                                       
count1=0;
for i in year:
    count=0
    for j in bkg_samples:
        #min_ = i                                                                                                                                                            
        #max_ = min_+ 100000                                                                                                                                                 
        # count=i                                                                                                                                                            
        # min_count= i+100000                                                                                                                                                
        condorSubmit = "condor_submit/submitCondor_%s_%s_Sbin"%(i,j)
        #fname1 = "/home/kalpana/t3store3/public/SUSY_GMSB_Aug_2020/CMSSW_11_1_0_pre3/src/Susy_Analysis_2020/runList/filelist/bkg_samples/v1/runList_%s_%s_v18.txt"%(j,i)
        fname1 = "/afs/cern.ch/user/c/cagrawal/public/MyAnalysis/CMSSW_14_0_2/src/SusySoft2023Ana/runList_%s_%s.txt"%(j,i)
        #fname = "/home/kalpana/t3store3/public/SUSY_GMSB_Aug_2020/CMSSW_11_1_0_pre3/src/Susy_Analysis_2020/Out_%s_%s_FinalOptimization_temp_v18.root"%(i,j)
        fname = "/afs/cern.ch/user/c/cagrawal/public/MyAnalysis/CMSSW_14_0_2_/src/SusySoft2023Ana/%s_%s.root"%(i,j)
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
            #line=line.replace("IN1",fname1)                                                                                                                                 
            #    line=line.replace("IN2",fname2)                                                                                                                             
            # line=line.replace("JOBNAME", jobName)
            print line.rstrip()

     # for line in fileinput.FileInput(condorSubmit, inplace=1):
     #        line=line.replace("MIN", str(i))
     #        line=line.replace("MAX", str(j))
     #        line=line.replace("outfile", fname)
     #        line=line.replace("filelist",fname1)
     #        line=line.replace("physicsList",number[count1])
     #        line=line.replace("sample_tag",sample_tag[count])
     #        #line=line.replace("IN1",fname1)                                                                                                                                 
     #        #    line=line.replace("IN2",fname2)                                                                                                                             
     #        line=line.replace("JOBNAME", jobName)
   
           

        print condorSubmit
        count+=1
    count1+=1
