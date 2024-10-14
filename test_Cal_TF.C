void test_Cal_TF()
{
  char* hist_name  = new char[200];
                                                                                                   
  //define your files here
  vector<string> f,rootFileName;
  f = {"All_combined_Summer20UL18.root"};
  vector<string> histnames1, histnames2, histnames3;
  histnames1 = {"lost_e_SR_binned", "lost_e_CR_binned"};
  histnames2 = {"lost_mu_SR_binned", "lost_mu_CR_binned"};
  histnames3 = {"LL_SR_binned", "LL_CR_binned"};

  rootFileName = {"e_TF_NBJet.root","mu_TF_NBJet.root","LL_TF_NBJet.root"};

  vector<vector<string>> cat;
  cat = {histnames1, histnames2, histnames3};
  
  for (int icat =0; icat < cat.size(); icat++){
    vector<vector<TH1F*>> hist_file;
      vector<TH1F*> hist_list;
      TFile *root_file = new TFile(f[0].c_str()); 
      for(int i_hist=0; i_hist < cat[icat].size();i_hist++)
	{
	  sprintf(hist_name,"%s",cat[icat][i_hist].c_str());
	  TH1F* resp = (TH1F*) root_file->Get(hist_name); //reading hist from the TFile
	  hist_list.push_back(resp);	      
	}
    
    TFile *TF_NBJet[3];                    //root file to store TF histogram
    TF_NBJet[icat] = new TFile(rootFileName[icat].c_str(), "RECREATE"); 

    TH1F* h_total_SR = (TH1F*)hist_list[0]-> Clone();	 
    TH1F* h_total_CR = (TH1F*)hist_list[1]-> Clone();
    cout << "SR Integral of sample " << f[0] << ": " << hist_list[0]->Integral() << endl;
    
    cout << "Total Integral: " << h_total_SR->Integral() << endl;
    
    TH1D* h_TF = (TH1D*) h_total_SR->Clone("combined_TF");
    h_TF->Divide(h_total_CR);
    TF_NBJet[icat] ->cd();
    h_TF->Write();
    cout << "TF stored in: " << rootFileName[icat] << "\n\n";
  }	  
}
