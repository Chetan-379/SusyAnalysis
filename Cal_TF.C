void Cal_TF()
{
  char* hist_name  = new char[200];
                                                                                                   
  //define your files here
  vector<string> f;
  f = {"Summer20UL18_TTJets.root", "Summer20UL18_TTGJets_Tune.root", "Summer20UL18_WJetsToLNu_HT.root", "Summer20UL18_WGJets_MonoPhoton.root"};

  vector<string> histnames;
  //histnames = {"lost_e_SR_binned", "lost_e_CR_binned"};
  //histnames = {"lost_mu_SR_binned", "lost_mu_CR_binned"};
  histnames = {"LL_SR_binned", "LL_CR_binned"};
  vector<vector<TH1F*>> hist_file;
	   
  for(int i_file=0; i_file < f.size(); i_file++){              
    vector<TH1F*> hist_list;
    TFile *root_file = new TFile(f[i_file].c_str()); 
    for(int i_hist=0; i_hist < histnames.size();i_hist++) 	    
      {
	sprintf(hist_name,"%s",histnames[i_hist].c_str());
	cout<<"i_file "<<i_file<<"\t"<<i_hist<<"\t"<<root_file->GetName()<< "\t" << hist_name<<endl;
	TH1F* resp = (TH1F*) root_file->Get(hist_name); //reading hist from the TFile
	hist_list.push_back(resp);	      
      }
    hist_file.push_back(hist_list);
  }

  TFile *TF_NBJet;                    //root file to store TF histogram
  TF_NBJet = new TFile("LL_TF_NBJet.root", "RECREATE"); 
  //TF_NBJet = new TFile("mu_TF_NBJet.root", "RECREATE"); 
  TH1F* h_total_SR = (TH1F*)hist_file[0][0]-> Clone();	 
  TH1F* h_total_CR = (TH1F*)hist_file[0][1]-> Clone();
  cout << "SR Integral of sample " << f[0] << ": " << hist_file[0][0]->Integral() << endl;
  for (int ifile =1; ifile < f.size(); ifile++){
    cout << "SR Integral of sample " << f[ifile] << ": " << hist_file[ifile][0]->Integral() << endl;
    h_total_SR ->Add(hist_file[ifile][0]);	    
    h_total_CR ->Add(hist_file[ifile][1]);	    
  }
  
  cout << "Total Integral: " << h_total_SR->Integral() << endl;
  
  TH1D* h_TF = (TH1D*) h_total_SR->Clone("combined_TF");
  h_TF->Divide(h_total_CR);
  TF_NBJet ->cd();
  h_TF->Write();	
	  
}
