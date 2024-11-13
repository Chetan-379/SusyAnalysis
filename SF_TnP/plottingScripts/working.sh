
# ## MC-Data comparisons 1 electron CR
## For SF                                                                                                          
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",0,0)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",0,1)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",0,2)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",0,3)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",0,4)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",0,5)'


root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",1,0)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",1,1)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",1,2)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",1,3)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",1,4)'
root -b -q 'CRcompare_DataVsMC.C("Results/data_mcComparisons/",1,5)'

## Calc FR in Data and MC as a function of different kinematics and processes
root -b -q 'DataMC_SRvsCR_kinem.C("Results/SRvsCR_KinemVar/data",1)'
root -b -q 'DataMC_SRvsCR_kinem.C("Results/SRvsCR_KinemVar/",0)'
## purity calculation
root -b -q 'v1_purityCalc.C("Results/PurityFact",0)'
# calculate SF
root -b -q 'CorrSF_calc.C("Results/SF_plots",0)'
