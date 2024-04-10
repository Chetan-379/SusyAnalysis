# Fake rate background estimation

## Instructions to run the code -
```
1. git clone -b UL_FakeRate_bkgEstimation  https://github.com/alpana-hep/susy_2020.git .
2. make (run make after you change anything in any of the source/header file)
3. ./analyzeLightBSM <filelist> <outfile> <year <process>  <photon ID>
```
Note - if you are reading nevents & cross section from the 'map_crosssection_SMprocess_v1.txt' file then make sure to keep the <process> name similar to the saved in 'map_crosssection_SMprocess.txt' and should contain UL in it if you are analyzing UL ntuples.
<photon ID>: 'loose', 'medium', 'tight','mva_wp90','mva_p80'
(Note-  first three are cutbased ID recommended by Egamma group.)
<filelist>: see under the 'inputFiles' directory
<outfile>: as you want to name your file
<process>: MC smaples for which job is running or data for data files

To submit the condor jobs:

<executable> is 'worker2.sh' (change or add destination path in worker2.sh). If no path is added than it will store in the parent directory from where the jobs are submitted. 
spliRunlist.C - create condor files and submit the condor jobs (improtant to add the files which you want to transfer)

```
root -l -q 'splitRunList.C("runList_Summer20UL18_WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8.txt",25,"2018","WJetsToLNu_HT-100To200UL","loose")'

```
submitMany2.sh - submit multiple jobs at a time for MC
```
source submitMany2.sh
```
For data
```
source submitMany1.sh

```
clean*.sh - to clean the log files of the condor jobs

run 

```
source cleanupBatchfiles.sh

```

Example to run the script

```
./analyzeLightBSM inputfiles/runList_Summer20UL18_TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8.txtp out_Summer20UL18_TTJets_HT_1200to2500.root 2018  TTJets_HTUL  loose

```

To make the plots

```
cd plottingScripts
source working.sh

```

## Scale factor calculation
We will be using Single electron and Egamms dataset instead of MET and DY Jets (ZLL Gamma )samples
```
cd SF_TnP/
make
./analyzeLightBSM <filelist> <outfile> <year <process>  <photon ID>
```
Note - if you are reading nevents & cross section from the 'map_crosssection_SMprocess_v1.txt' file then make sure to keep the <process> name similar to the saved in 'map_crosssection_SMprocess.txt' and should contain UL in it if you are analyzing UL ntuples.
<photon ID>: 'loose', 'medium', 'tight','mva_wp90','mva_p80'
(Note-  first three are cutbased ID recommended by Egamma group.)
<filelist>: see under the 'inputFiles' directory
<outfile>: as you want to name your file
<process>: MC smaples for which job is running or data for data files

To submit the condor jobs:

<executable> is 'worker2.sh' (change or add destination path in worker2.sh). If no path is added than it will store in the parent directory from where the jobs are submitted.
spliRunlist.C - create condor files and submit the condor jobs (improtant to add the files which you want to transfer)

```
root -l -q 'splitRunList.C("runList_Summer20UL18_WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8.txt",25,"2018","WJetsToLNu_HT-100To200UL","loose")'

```
submitMany1.sh - submit multiple jobs at a time for MC & electron data
```
source submitMany1.sh
```
For WJets & TTJets
```
source submitMany2.sh

```
clean*.sh - to clean the log files of the condor jobs

run

```
source cleanupBatchfiles.sh

```

Example to run the script

```
./analyzeLightBSM inputfiles/runList_Summer20UL18_TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8.txtp out_Summer20UL18_TTJets_HT_1200to2500.root 2018  TTJets_HTUL  loose

```

To make the plots

```
cd plottingScripts
source working.sh

```
