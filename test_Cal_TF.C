void test_Cal_TF()
{
  char* hist_name  = new char[200];                                                                                                   
  vector<string> f,tf1,tf2,tf3,tf4;
  vector<vector<string>> rootFileName;
  f = {"Summer20UL16APV_LL_All_Combined.root","Summer20UL16_LL_All_Combined.root","Summer20UL17_LL_All_Combined.root","Summer20UL18_LL_All_Combined.root"};
  
  vector<string> histnames1, histnames2, histnames3;
  histnames1 = {"lost_e_SR_binned", "lost_e_CR_binned"};
  histnames2 = {"lost_mu_SR_binned", "lost_mu_CR_binned"};
  histnames3 = {"LL_SR_binned", "LL_CR_binned"};

  tf1 = {"2016APV_lost_e_TF_NBJet.root","2016APV_lost_mu_TF_NBJet.root","2016APV_LL_TF_NBJet.root"};
  tf2 = {"2016_lost_e_TF_NBJet.root","2016_lost_mu_TF_NBJet.root","2016_LL_TF_NBJet.root"};
  tf3 = {"2017_lost_e_TF_NBJet.root","2017_lost_mu_TF_NBJet.root","2017_LL_TF_NBJet.root"};
  tf4 = {"2018_lost_e_TF_NBJet.root","2018_lost_mu_TF_NBJet.root","2018_LL_TF_NBJet.root"};

  rootFileName = {tf1, tf2, tf3, tf4};
  
  vector<vector<string>> cat;
  cat = {histnames1, histnames2, histnames3};
  
  for (int i_year =0; i_year < f.size(); i_year++){
    for (int icat =0; icat < cat.size(); icat++){
      vector<TH1F*> hist_list;
      TFile *root_file = new TFile(f[i_year].c_str()); 
      for(int i_hist=0; i_hist < cat[icat].size();i_hist++)
	{
	  sprintf(hist_name,"%s",cat[icat][i_hist].c_str());
	  TH1F* resp = (TH1F*) root_file->Get(hist_name); //reading hist from the TFile
	  hist_list.push_back(resp);	      
	}
      
      TFile *TF_NBJet[4][3];                    //root file to store TF histogram
      TF_NBJet[i_year][icat] = new TFile(rootFileName[i_year][icat].c_str(), "RECREATE"); 
      
      TH1F* h_total_SR = (TH1F*)hist_list[0]-> Clone();	 
      TH1F* h_total_CR = (TH1F*)hist_list[1]-> Clone();
      cout << "SR Integral of sample " << f[i_year] << ": " << hist_list[0]->Integral() << endl;
      
      TH1D* h_TF = (TH1D*) h_total_SR->Clone("combined_TF");
      h_TF->Divide(h_total_CR);
      TF_NBJet[i_year][icat] ->cd();
      h_TF->Write();
      cout << "TF stored in: " << rootFileName[i_year][icat] << "\n\n";

      if (i_year ==4){
	if (icat==0) {cout << "TFs for lost e:" << endl;
	  for (int ibin=0; ibin < h_TF->GetNbinsX(); ibin++){
	    cout << "TF in bin " << ibin << ": " << h_TF->GetBinContent(ibin) << endl;
	  }
	  cout << "\n\n";
	}

	if (icat==1) {cout << "TFs for lost mu:" << endl;
	  for (int ibin=0; ibin < h_TF->GetNbinsX(); ibin++){
	    cout << "TF in bin " << ibin << ": " << h_TF->GetBinContent(ibin) << endl;
	  }
	  cout << "\n\n";
	}

	if (icat==2) {cout << "TFs for LL combined categories:" << endl;
	  for (int ibin=0; ibin < h_TF->GetNbinsX(); ibin++){
	    cout << "TF in bin " << ibin << ": " << h_TF->GetBinContent(ibin) << endl;
	  }
	  cout << "\n\n";
	}

      }
    }    
  }
}
