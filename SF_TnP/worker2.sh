#!/bin/sh

executable=$1
inputFileTag=$2
outputFileTag=$3
#commitHash=$4
datasetName=$4
process=$5
#LL=$6
phoID=$6
currDir=$(pwd)
######################################
# SETUP CMSSW STUFF...
######################################

#uncomment only below para:
# source /cvmfs/cms.cern.ch/cmsset_default.sh
# export SCRAM_ARCH=el9_amd64_gcc12
# scram p CMSSW CMSSW_14_0_0_pre0
# cd CMSSW_14_0_0_pre0/src
cd /afs/cern.ch/user/c/cagrawal/public/MyAnalysis/CMSSW_14_0_2/src/Alpana_FR_machinery/SF_TnP

# export SCRAM_ARCH=slc7_amd64_gcc820
# scram p CMSSW CMSSW_11_1_0_pre3
# cd CMSSW_14_0_0_pre0/src

# eval `scramv1 runtime -sh`
pwd
# echo $CMSSW_RELEASE_BASE
# cd $currDir
# echo $currDir


######################################
# SETUP PRIVATE STUFF...
######################################
echo "RUNNING ANALYSIS"
pwd
echo $executable
echo $inputFileTag
./$executable $inputFileTag $outputFileTag $datasetName $process $phoID
echo "processed. ls"
ls
echo "COPYING OUTPUT"

#xrdcp -f skimmed_ntuple_$datasetName'_'$process'.root' root://cmseos.fnal.gov//store/user/kalpana/Susy_phoMet/SkimmedNtuples/
#xrdcp -f $outputFileTag root://cmseos.fnal.gov//store/user/kalpana/ul_rootop_Analys_May2/FR/SF_tnp 
mv -f $outputFileTag /eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/SF_TnP/unskimmed_root_files/ 

#xrdcp -f ${datasetName}'_'${outputFileTag} root://cmseos.fnal.gov//store/group/lpcsusyhad/SusyPhotonMET/Run2ProductionV17/background/skims/${outputFileTag}
#rm $outputFileTag
rm skimmed_ntuple_$datasetName'_'$process'.root'
#xrdcp -f $outputFileTag root://cmseos.fnal.gov//store/group/lpcsusyhad/SusyPhotonMET/Run2ProductionV17/background/$outputFileTag
#xrdcp -f $outputFileTag root://cmseos.fnal.gov//store/user/bkansal/GMSB_skims_TreesV17/for_bkg_estimation/lost_electron/new2/CR_$outputFileTag
#xrdcp -f $outputFileTag root://cmseos.fnal.gov//store/user/bkansal/GMSB_skims_TreesV18/SignalRegion/skims/SR_$outputFileTag
#xrdcp -f $outputFileTag root://cmseos.fnal.gov//store/user/bkansal/GMSB_skims_TreesV18/for_bkg_estimation/lost_electron/CR_$outputFileTag
#rm $outputFileTag
