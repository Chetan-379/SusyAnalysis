void verify()
{
  char* hist_name  = new char[200];
  vector<string> years;
  years = {"18","17","16","16APV","FullRun2"};

  //for (int iyear =0; iyear < years.size(); iyear++){
  //cout << "FOR Year 20" << years[iyear] << "====================================================== " << endl;
  char filename1[1000], filename2[1000], filename3[1000], filename4[1000], filename5[1000], filename6[1000];
    
    // if(iyear !=4){
    // sprintf(filename1,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/Summer20UL%s_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
    // sprintf(filename2,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/Alpana_FR_root_files/Summer20UL%s_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
    // }

    // else {
    //   sprintf(filename1,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/%s_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
    //   sprintf(filename2,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/Alpana_FR_root_files/%s_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root", years[iyear].c_str());
    // }
  
  sprintf(filename1,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/unskimmed_root_files/%s_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root", years[4].c_str());
  sprintf(filename2,"/eos/home-c/cagrawal/SusySoftPhoAna/FR_estimation/Alpana_FR_root_files/%s_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root", years[4].c_str());
    
      
    vector<string> histnames;
    
    //histnames = {"h_Sbins_LL_newSbins_Validation_v7_Elec_CR"};
    histnames = {"h_Sbins_LL_newSbins_v7_Pho_SR"};
    
    
 
    //for (int i_file =0; i_file <f.size(); i_file++){
    //cout << "Yields for file " << f[i_file] << ": " << endl;
    TFile *root_file_me = new TFile(filename1);
    TFile *root_file_Alpana = new TFile(filename2); 
    for(int i_hist=0; i_hist < histnames.size();i_hist++){
      sprintf(hist_name,"%s",histnames[i_hist].c_str());
      TH1F* resp1 = (TH1F*) root_file_me->Get(hist_name);
      TH1F* resp2 = (TH1F*) root_file_Alpana->Get(hist_name);
      //cout << "\n\naddress: " << resp2 << endl;

      if (resp1->GetNbinsX() != resp2->GetNbinsX()) break;

      cout << "Results from Alpana's root file======================================" << endl;
      for(int ibin=0; ibin<resp2->GetNbinsX(); ibin++) {
	//double difference = resp1 ->GetBinContent(ibin) - resp2 ->GetBinContent(ibin); 
	//if ((resp1 ->GetBinContent(ibin) != resp2 ->GetBinContent(ibin)) && abs(difference) > pow(10,-5)) {
	// cout << "my hist bin: " << resp1 -> GetBinContent(ibin) << endl;
	// cout << "Alpana hist bin: " << resp2-> GetBinContent(ibin) << endl;
	// cout << "difference: " << resp1->GetBinContent(ibin)-resp2->GetBinContent(ibin) << endl;
	// cout << endl;
	cout << resp2->GetBinContent(ibin) << endl;
	//}
      }

      cout << "Results from Chetan's root file======================================" << endl;
      for(int ibin=0; ibin<resp1->GetNbinsX(); ibin++) {
	//double difference = resp1 ->GetBinContent(ibin) - resp2 ->GetBinContent(ibin); 
	//if ((resp1 ->GetBinContent(ibin) != resp2 ->GetBinContent(ibin)) && abs(difference) > pow(10,-5)) {
	// cout << "my hist bin: " << resp1 -> GetBinContent(ibin) << endl;
	// cout << "Alpana hist bin: " << resp2-> GetBinContent(ibin) << endl;
	// cout << "difference: " << resp1->GetBinContent(ibin)-resp2->GetBinContent(ibin) << endl;
	// cout << endl;
	cout << resp1->GetBinContent(ibin) << endl;
	//}
      }

	
      // cout << "my hist integral: " << resp1->Integral() << endl;
      // cout << "alpana hist integral: " << resp2->Integral() << endl;	
	
    }
    //cout << endl;
    //}
}
