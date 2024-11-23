void year_avg_Cal_TF()
{
  char* hist_name  = new char[200];                                                                                                   
  string f;
  vector<string> rootFileName;
  f = "root_files/AllYears_LL_All_Combined.root";
  vector<string> histnames1, histnames2, histnames3;
  histnames1 = {"lost_e_SR_binned", "lost_e_CR_binned"};
  histnames2 = {"lost_mu_SR_binned", "lost_mu_CR_binned"};
  histnames3 = {"LL_SR_binned", "LL_CR_binned"};

  rootFileName = {"year_avg_lost_e_TF_NBJet.root","year_avg_lost_mu_TF_NBJet.root","year_avg_LL_TF_NBJet.root"};
  
  vector<vector<string>> cat;
  cat = {histnames1, histnames2, histnames3};

  TFile *TF_NBJet[3];                    //root file to store TF histogram  
  for (int icat =0; icat < cat.size(); icat++){
    vector<TH1F*> hist_list;
    TFile *root_file = new TFile(f.c_str()); 
    for(int i_hist=0; i_hist < cat[icat].size();i_hist++)
      {
	sprintf(hist_name,"%s",cat[icat][i_hist].c_str());
	TH1F* resp = (TH1F*) root_file->Get(hist_name); //reading hist from the TFile
	hist_list.push_back(resp);	      
      }
    
    TF_NBJet[icat] = new TFile(rootFileName[icat].c_str(), "RECREATE");
    if (icat == 0){
      cout << "Total SR of lost e: " << hist_list[0]->Integral() << endl;
      cout << "Total CR of lost e: " << hist_list[1]->Integral() << "\n\n";     
    }

    if (icat == 1){
      cout << "Total SR of lost mu: " << hist_list[0]->Integral() << endl;
      cout << "Total CR of lost mu: " << hist_list[1]->Integral() << "\n\n";     
    }

    if (icat == 2){
      cout << "Total SR of LL: " << hist_list[0]->Integral() << endl;
      cout << "Total CR of LL: " << hist_list[1]->Integral() << "\n\n";     
    }
      
    TH1F* h_total_SR = (TH1F*)hist_list[0]-> Clone();	 
    TH1F* h_total_CR = (TH1F*)hist_list[1]-> Clone();
      
    TH1D* h_TF = (TH1D*) h_total_SR->Clone("combined_TF");
    h_TF->Divide(h_total_CR);
    TF_NBJet[icat] ->cd();
    h_TF->Write();
    cout << "TF stored in: " << rootFileName[icat] << "\n\n";    
  }      
}
