## For fakerate
# root -b  'two2D_SRvsCR.C("Results/MET_200_EMpt100/twoDplots",0)'
# root -b  'two2D_SRvsCR.C("Results/MET_200_EMpt100/twoDplots/HEm_veto",1)'
# root -b  'two2D_SRvsCR.C("Results/MET_200_EMpt100/twoDplots/L1TrigProb",2)'

# root -b 'TF_varKinem_ratio.C("Results/MET_200_EMpt100/TFVarKineMatic",0)'
# root -b 'TF_varKinem_ratio.C("Results/MET_200_EMpt100/TFVarKineMatic/HEM_veto",1)'
# root -b 'TF_varKinem_ratio.C("Results/MET_200_EMpt100/TFVarKineMatic/L1TrigProb",2)'

# root -b -q 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar",0)'
# root -b -q 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar/HEM_veto",1)'
# root -b -q 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar/L1TrigProb",2)'

# root -b -q 'two2D_SRvsCR.C("Results/SRvsCR_KinemVar/2dPlots",0)'
# root -b -q 'two2D_SRvsCR.C("Results/SRvsCR_KinemVar/2dPlots/HEM_veto",1)'
# root -b -q 'two2D_SRvsCR.C("Results/SRvsCR_KinemVar/2dPlots/L1TrigProb",2)'

# ## Kinematics comparisons SR vs CR
# root -b 'TF_varKinem_ratio.C("Results/SRvsCR_KinemVar")'

# ## Transfer factors
# root -b 'SRvsCRStacked_LL_varRatio.C("Results/TransferFactors/TFbins_v1_nJets_BJets",1)'
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

root -b -q 'DataMC_SRvsCR_kinem.C("Results/SRvsCR_KinemVar/data",1)'
root -b -q 'DataMC_SRvsCR_kinem.C("Results/SRvsCR_KinemVar/",0)'
root -b 'CorrSF_calc.C("Results/SF_plots",0)'

# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,0)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,1)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,2)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,3)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,4)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_veto/",1,5)'


# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_vetoAnd_L1TrigProb/",2,0)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_vetoAnd_L1TrigProb/",2,1)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_vetoAnd_L1TrigProb/",2,2)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_vetoAnd_L1TrigProb/",2,3)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_vetoAnd_L1TrigProb/",2,4)'
# root -b -q 'plotAlps_RatioPlots.C("Results/data_mcComparisons/AfterHEM_vetoAnd_L1TrigProb/",2,5)'

# ## validation on MC
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/SRbins/Bins_inPhopT_Qmulti",1,1,1)'
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/SRbins/Bins_inPhopT_Qmulti/v1",1,1,1)'
# root -b 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/Bins_inPhopT_Qmulti/v1",1,1,1)'
# root -b 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/Bins_inPhopT_Qmulti/",1,1,1)'
# # ## data MC comparisons
# # root -l 'plotAlps_RatioPlots.C("path to save plots",which lepton (1 for electron, 0 muon, 2 for e+mu),which year (0=2016 (pre+post VFP), 1=2017, 2=2018, 3=2016pre,4=2016post,5=fullrun2))'
# ## total lepton
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",2,0)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",2,1)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",2,2)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",2,3)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",2,4)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",2,5)'

# ## Electron                                                                                                                 
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",1,0)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",1,1)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",1,2)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",1,3)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",1,4)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",1,5)'

# ## Muon                                                                                                                 
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",0,0)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",0,1)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",0,2)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",0,3)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",0,4)'
# root -b -q 'plotAlps_RatioPlots.C("data_mcComparisons/",0,5)'

# ## showing TF as 1 SR and 1 CR (as it is in Bhumika's case)
                                                          
# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v1",1,1,1)'
# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v1",2,1,1)'
# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v1",3,1,1)'

# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v2",1,2,1)'
# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v2",2,2,1)'
# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v2",3,2,1)'

# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v3",1,3,1)'
# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v3",2,3,1)'
# root -b -q 'SRvsCRStacked_LL_varRatio.C("TransferFactors/SRvsCR/TFbins_v3",3,3,1)'


# ## calculating TF for lost leptons v1 - 8 bins

# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/",3,1,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/",2,1,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/",1,1,1)'
                                                                                                     
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/TFBins_v1_phoPt_nJets_Btags",3,2,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/TFBins_v1_phoPt_nJets_Btags",2,2,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/TFBins_v1_phoPt_nJets_Btags",1,2,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/TFBins_v2_MET_nJets_Btags",3,3,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/TFBins_v2_MET_nJets_Btags",2,3,1)'
# root -b -q 'DiffStacked_LL_varRatio.C("TransferFactors/TFBins_v2_MET_nJets_Btags",1,3,1)'

# ## TF vs different categories for lost lepton (failing acceptance, isolation and identification)
# root -b -q 'EveCateg_LL_ratio.C("TransferFactors/EventCateg_TF/TFbins_v1",1,1,1)'
# root -b -q 'EveCateg_LL_ratio.C("TransferFactors/EventCateg_TF/TFbins_v1",2,1,1)'
# root -b -q 'EveCateg_LL_ratio.C("TransferFactors/EventCateg_TF/TFbins_v2",1,2,1)'
# root -b -q 'EveCateg_LL_ratio.C("TransferFactors/EventCateg_TF/TFbins_v2",2,2,1)'
# root -b -q 'EveCateg_LL_ratio.C("TransferFactors/EventCateg_TF/TFbins_v3",1,3,1)'
# root -b -q 'EveCateg_LL_ratio.C("TransferFactors/EventCateg_TF/TFbins_v3",2,3,1)'


# ## TF in varbins for different kinematics
# root -b -q 'TF_varKinem_ratio.C("TF_inKinematics",1,3,1)'
# root -b -q 'TF_varKinem_ratio.C("TF_inKinematics",2,3,1)'
# root -b -q 'TF_varKinem_ratio.C("TF_inKinematics",3,3,1)'
# ## validation of TF on MC in search bins
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v1",1,1,1)' #electron for TFBins v1
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v1",2,1,1)' # muon
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v1",3,1,1)' # all leptons

                                                                                              
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v2",1,2,1)' #electron for TFBins v2                                                                          
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v2",2,2,1)'
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v2",3,2,1)'

                                                                                              
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v3",1,3,1)' #electron for TFBins v2                                                                          
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v3",2,3,1)'
# root -b -q 'Valid_SRBins_LL_wrRatio.C("Validation_MC/Sbins/TFBIns_v3",3,3,1)'


# ## validation of TF in MC - different kinematic variables
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v3",3,3,1)'
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v3",2,3,1)'
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v3",1,3,1)'

                                                                      
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v2",3,2,1)'
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v2",2,2,1)'
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v2",1,2,1)'

                                                                      
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v1",3,1,1)'
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v1",2,1,1)'
# root -b -q 'Valid_diffKinematics_wrRatio.C("Validation_MC/Kinematics_Valid/TFBins_v1",1,1,1)'



# ## prediction in data and MC expected comparisons in search bins
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,0,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,1,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,2,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,3,1)'                                                                                           
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,4,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",3,5,1)'                                                                              
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",2,0,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",2,1,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",2,2,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",2,3,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",2,4,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",2,5,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",1,0,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",1,1,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",1,2,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",1,3,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",1,4,1)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v1",1,5,1)'

                                                                                      
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",3,0,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",3,1,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",3,2,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",3,3,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",3,4,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",3,5,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",2,0,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",2,1,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",2,2,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",2,3,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",2,4,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",2,5,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",1,0,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",1,1,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",1,2,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",1,3,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",1,4,2)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v2",1,5,2)'


                                                                                      
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",3,0,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",3,1,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",3,2,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",3,3,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",3,4,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",3,5,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",2,0,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",2,1,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",2,2,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",2,3,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",2,4,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",2,5,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",1,0,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",1,1,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",1,2,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",1,3,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",1,4,3)'
# root -b -q 'plot_DataValida_SRbins.C("Validation_data/SRbins/TFbins_v3",1,5,3)'


# ## validation of TF in Data and MC - different kinematic variables
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",3,0,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",3,1,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",3,2,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",3,3,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",3,4,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",3,5,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",2,0,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",2,1,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",2,2,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",2,3,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",2,4,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",2,5,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",1,0,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",1,1,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",1,2,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",1,3,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",1,4,1)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v1",1,5,1)'

                                                                              
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",3,0,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",3,1,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",3,2,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",3,3,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",3,4,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",3,5,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",2,0,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",2,1,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",2,2,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",2,3,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",2,4,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",2,5,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",1,0,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",1,1,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",1,2,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",1,3,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",1,4,2)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v2",1,5,2)'

# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",3,0,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",3,1,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",3,2,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",3,3,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",3,4,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",3,5,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",2,0,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",2,1,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",2,2,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",2,3,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",2,4,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",2,5,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",1,0,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",1,1,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",1,2,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",1,3,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",1,4,3)'
# root -b -q 'plot_DataValida_kinem.C("Validation_data/Kinematics/TFbins_v3",1,5,3)'

