void get_yield()
{
  char* hist_name  = new char[200];
  vector<string> years;
  years = {"18","17","16","16APV","FullRun2"};

  vector<string> f;

  for (int iyear =0; iyear < years.size(); iyear++){
    cout << "yieds for the year 20" << years[iyear] << ": " << endl;
    char filename1[1000], filename2[1000], filename3[1000], filename4[1000], filename5[1000], filename6[1000];
    if (iyear != 4) {
      sprintf(filename1,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/Summer20UL%s_TTJets_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
      sprintf(filename2,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/Summer20UL%s_TTGJets_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
      sprintf(filename3,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/Summer20UL%s_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
      sprintf(filename4,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/Summer20UL%s_WGJets_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
      sprintf(filename5,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/Summer20UL%s_singleTop_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
      sprintf(filename6,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/Summer20UL%s_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
    }

    if (iyear == 4) {
      sprintf(filename1,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/%s_TTJets_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
      sprintf(filename2,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/%s_TTGJets_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
      sprintf(filename3,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/%s_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
      sprintf(filename4,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/%s_WGJets_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
      sprintf(filename5,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/%s_singleTop_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
      sprintf(filename6,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/%s_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
    }

    vector<string> f;
    f= {filename1, filename2, filename3, filename4, filename5, filename6};

    vector<string> histnames;
    //histnames = {"h_NhadJets_Pho_SR", "h_NhadJets_Elec_CR"};
    histnames = {"h_Sbins_LL_Pho_SR", "h_Sbins_LL_Validation_Elec_CR"};
    //histnames = {"h_Sbins_LL_newSbins_v7_Pho_SR", "h_Sbins_LL_newSbins_Validation_v7_Elec_CR"};
    
    
    for (int i_file =0; i_file <f.size(); i_file++){
      cout << "Yields for file " << f[i_file] << ": " << endl;
      TFile *root_file = new TFile(f[i_file].c_str()); 
      for(int i_hist=0; i_hist < histnames.size();i_hist++){
	sprintf(hist_name,"%s",histnames[i_hist].c_str());
	TH1F* resp = (TH1F*) root_file->Get(hist_name); 
	// if (i_hist == 0) cout << "Expected SR: " << resp->Integral() << endl;
	// if (i_hist == 1) cout << "CR: " << resp->Integral() << endl;
	if (i_hist == 0) cout << resp->Integral() << endl;
	if (i_hist == 1) cout << resp->Integral() << endl;
	
      }
      cout << endl;
    }
    cout << "\n\n";
  }
}
