#!/bin/bash
cd /afs/cern.ch/user/c/cagrawal/public/MyAnalysis/CMSSW_14_0_2/src/SusySoft2023Ana
# source /cvmfs/cms.cern.ch/cmsset_default.sh
# eval `scram runtime -sh`

source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=el9_amd64_gcc12

eval `scramv1 runtime -sh`


echo $PWD
cd /afs/cern.ch/user/c/cagrawal/public/MyAnalysis/CMSSW_14_0_2/src/SusySoft2023Ana

./analyzeTProxytBSM $1 $2 $3 $4 $5 $6
