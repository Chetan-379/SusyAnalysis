## For fakerate
# root -b  'two2D_SRvsCR.C("Results/MET_200_EMpt100/twoDplots",0)'
# root -b  'two2D_SRvsCR.C("Results/MET_200_EMpt100/twoDplots/HEm_veto",1)'
# root -b  'two2D_SRvsCR.C("Results/MET_200_EMpt100/twoDplots/L1TrigProb",2)'


# root -b 'TF_varKinem_ratio.C("Results/MET_200_EMpt100/TFVarKineMatic",0)'
# root -b 'TF_varKinem_ratio.C("Results/MET_200_EMpt100/TFVarKineMatic/HEM_veto",1)'
# root -b 'TF_varKinem_ratio.C("Results/MET_200_EMpt100/TFVarKineMatic/L1TrigProb",2)'

# root -b -q 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar",0)'
# root -b -q 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar/HEM_veto",1)'
root -b -q 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar/",0)'

# root -b -q 'two2D_SRvsCR.C("Results/SRvsCR_KinemVar/2dPlots",0)'
# root -b -q 'two2D_SRvsCR.C("Results/SRvsCR_KinemVar/2dPlots/HEM_veto",1)'
root -b -q 'two2D_SRvsCR.C("Results/SRvsCR_KinemVar/2dPlots/",0)'

## Kinematics comparisons SR vs CR
#root -b 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar")'

## Transfer factors
#root -b 'SRvsCRStacked_LL_varRatio.C("Results/TransferFactors/TFbins_v1_nJets_BJets",1)'
root -b -q 'SRvsCRStacked_LL_varRatio.C("Results/TransferFactors/TFbins_v3_phopt_qmulti",3)'

## MC-Data comparisons 1 electron CR
                                                                                                          
root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",0,0)'
root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",0,1)'
root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",0,2)'
root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",0,3)'
root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",0,4)'
root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",0,5)'


# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,0)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,1)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,2)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,3)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,4)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,5)'


# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",2,0)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",2,1)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",2,2)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",2,3)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",2,4)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/",2,5)'

## validation on MC
root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/SRbins/Bins_inPhopT_Qmulti",1,1,1)'
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/SRbins/Bins_inPhopT_Qmulti/v1",1,1,1)'
# root -b 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/Bins_inPhopT_Qmulti/v1",1,1,1)'
root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/Bins_inPhopT_Qmulti/",1,1,1)'


## new SR bins
root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/SRbins/Bins_inPhopT_Qmulti",1,4,1)'
# extra plots
root -b 'extraPlots_CRvsSR.C("Results/ExtraPlots_validationChecks/SRvsCR",0)'
root -b 'extraPlots_valid.C("Results/ExtraPlots_validationChecks/Validation_MC",0)'
root -b 'extraPlots_2D.C("Results/ExtraPlots_validationChecks/2Dplots",0)'
