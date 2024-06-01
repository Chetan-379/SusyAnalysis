TString legend_text = "put whatever";
const int n_pl = 4;
bool logx = false;
TH1F* setLastBinAsOverFlow(TH1F*, int);
TH1F* setMyRange(TH1F*,double,double);
TH1F* DrawOverflow(TH1F*);
TH1F* DrawOverflow(TH1F* h,int xmin, int xrange){
  //function to paint the histogram h with an extra bin for overflows
  // This function paint the histogram h with an extra bin for overflows
   UInt_t nx    = h->GetNbinsX()+1;
   Double_t *xbins= new Double_t[nx+1];
   for (UInt_t i=0;i<nx;i++)
     xbins[i]=h->GetBinLowEdge(i+1);
   xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
   char *tempName= new char[strlen(h->GetName())+10];
   sprintf(tempName,"%swtOverFlow",h->GetName());
   h->GetXaxis()->SetLimits(xmin,xrange);
   // Book a temporary histogram having ab extra bin for overflows
   TH1F *htmp = new TH1F(tempName, h->GetTitle(), nx, xbins);
   htmp->GetXaxis()->SetRange(xmin,xrange);
   // Reset the axis labels
   htmp->SetXTitle(h->GetXaxis()->GetTitle());
   htmp->SetYTitle(h->GetYaxis()->GetTitle());
   // Fill the new hitogram including the extra bin for overflows
   for (UInt_t i=1; i<=nx; i++)
     htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
   // Fill the underflows
   htmp->Fill(h->GetBinLowEdge(1)-1, h->GetBinContent(0));
   // Restore the number of entries
   htmp->SetEntries(h->GetEntries());
   // FillStyle and color
   // htmp->SetFillStyle(h->GetFillStyle());
   // htmp->SetFillColor(h->GetFillColor());
   return htmp;
}
TH1F* setLastBinAsOverFlow(TH1F* h_hist, int xrange){
  //     h_hist = setMyRange(h_hist,0,xrange);
  //  h_hist->GetXaxis()->SetRangeUser(0,xrange);
  double lastBinCt =h_hist->GetBinContent(h_hist->GetNbinsX()),overflCt =h_hist->GetBinContent(h_hist->GetNbinsX());
  //  cout<<h_hist->GetNbinsX()<<"\t"<<lastBinCt<<"\t"<<overflCt<<endl;

  double lastBinErr=h_hist->GetBinError(h_hist->GetNbinsX()),  overflErr=h_hist->GetBinError(h_hist->GetNbinsX()+1); 
  if(lastBinCt!=0 && overflCt!=0)
    lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

  else if(lastBinCt==0 && overflCt!=0)
    lastBinErr = overflErr;
  else if(lastBinCt!=0 && overflCt==0)
    lastBinErr = lastBinErr;
  else lastBinErr=0;
  //h_temp->GetXaxis()->SetRangeUser(0,xrange);

  lastBinCt = lastBinCt+overflCt;
  //  cout<<lastBinCt<<endl;
  TH1F* h_temp = (TH1F*)h_hist->Clone();
  h_temp->SetBinContent(h_hist->GetNbinsX(),lastBinCt);
  h_temp->SetBinError(h_hist->GetNbinsX(),lastBinErr);
  //  h_temp->GetXaxis()->SetRangeUser(0,xrange);

  // h_hist = setMyRange(h_hist,0,xrange);
  //
  return h_temp;
}


TH1F* setMyRange(TH1F *h1,double xLow,double xHigh){
  //call it after setting last bin as overflow                                                                                                    
  double err=0;
  if(xHigh > 13000) return h1;
  if(xLow < -13000) return h1;

  // h1->Print("all");
  //  h1->GetXaxis()->SetRangeUser(xLow,xHigh);  
  int nMax=h1->FindBin(xHigh);
  h1->SetBinContent(nMax,h1->IntegralAndError(nMax,h1->GetNbinsX(),err));
  h1->SetBinError(nMax,err);

  //  cout<<nMax<<endl;
  //  h1->GetXaxis()->SetRangeUser(xLow,xHigh);
  for(int i=nMax+1;i<=h1->GetNbinsX()+1;i++){
    h1->SetBinContent(i,0);
    h1->SetBinError(i,0);
    //    cout<<":";
    //h1->GetXaxis()->SetRangeUser(xLow,xHigh); 
  }
  return h1;
}
//important function - change in this function if you want to change decoration of a plot
//void generate_1Dplot(vector<TH1F*> hist, char const *tag_name="", float energy=-1, int xmax=-1,int xmin=-1,char const *leg_head="",
//		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title="", const char* xtitile="", int rebin=-1){  
void generate_1Dplot(TH1F* hist[10][10], char const *tag_name="", float energy=-1, int xmax=-1,int xmin=-1,char const *leg_head="",
		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title="", const char* xtitile="", int rebin=-1){  
  

  TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name,900,750);//600,600,1200,1200);
  canvas_n1->Range(-60.25,-0.625,562.25,0.625);
  canvas_n1->SetFillColor(0);
  canvas_n1->SetBorderMode(0);
  canvas_n1->SetBorderSize(2);
  canvas_n1->SetRightMargin(0.045);
  canvas_n1->SetLeftMargin(0.12);
  canvas_n1->SetTopMargin(0.06);
  canvas_n1->SetBottomMargin(0.12);
  // auto *pad_1 = new TPad("pad_1","pad_1",0.,0.0,1.,0.32); pad_1->Draw();
  // pad_1->SetTopMargin(0.06);
  // pad_1->SetBottomMargin(0.3);
  // pad_1->SetRightMargin(0.025);
  // pad_1->SetLeftMargin(0.14);
  
  // auto *p1 = new TPad("p1","p1",0.,0.32,1.,1.);  p1->Draw();
  // p1->SetBottomMargin(0.04);
  // p1->SetRightMargin(0.025);
  // p1->SetLeftMargin(0.14);
  // p1->SetTopMargin(0.05);
  // p1->cd();
  gStyle->SetOptStat(0);
  
  vector<TString> legName;
  std::string leg_head_str = leg_head;
  double x = 0.15;
  double y = 0.90;
  TLegend *legend; //legend to be drawn on the plot - shift x,ys if you want to move this on the canvas
  legend = new TLegend(0.65,0.65,0.79,0.9);  
  legend->SetTextSize(0.030);
  legend->SetLineColor(kWhite);
  char* lhead = new char[100];

  // TLegend *legend1; //legend to be drawn on the plot - shift x,ys if you want to move this on the canvas
  // legend1 = new TLegend(0.65,0.55,0.79,0.8);  
  // legend1->SetTextSize(0.030);
  // legend1->SetLineColor(kWhite);
  //char* lhead = new char[100];
  
  sprintf(lhead,"#bf{%s} ",title); // legend header -- for example you are plotting this for WGJets then title string can be "WGJets"
  legend->SetHeader(lhead);
  legend->SetLineColor(kWhite);

  //legend1->SetHeader(lhead);
  //legend1->SetLineColor(kWhite);

  TLegendEntry* leg_entry;
  float x_label_size = 0.045;
  double ymin = 100000.0;
  double ymax = 0.0;


  for(int i = 0; i < 3; i++) {
    // if(DoRebin) {
    //  hist[i][j]->Rebin(2);

    // }
    //hist.at(i)= setLastBinAsOverFlow(hist.at(i),xrange);
     

    // normalize = true;

    // Int_t NEvnts;
    // float Evnt_frac;
    // Int_t no_of_events;
    
    // no_of_events = hist[i][j]->GetEntries();
    // NEvnts = hist[i][j]->Integral();
    // Evnt_frac = (no_of_events/(hist.at(0)->GetEntries()))*100;
    for (int j=0; j<3; j++){
      int overflCt= (hist[i][j]->GetBinContent(hist[i][j]->GetNbinsX()+1));
    
    //cout << no_of_events << "    (frac: " << setprecision(3) << Evnt_frac << "%)" << overflCt << endl;
    cout << "overflow: " << overflCt << endl;
    //cout << "histname: " << hist[i][j]->GetName() << endl;
    if(normalize) {
      	hist[i][j]->Scale(1.0/hist[i][j]->Integral());
      	hist[i][j]->GetYaxis()->SetTitle("Normalized");
    }
    else {
      hist[i][j]->GetYaxis()->SetTitle("Entries");
       }
    hist[i][j]->GetXaxis()->SetRangeUser(xmin,xmax); //to be given while calling this function
    //hist[i][j]->SetLineWidth(line_width[i]); //these are defined on top of this script    
    // hist[i][j]->SetLineStyle(line_style[i]);
    // hist[i][j]->SetLineColor(line_color[i]);
    hist[i][j]->SetTitle(" "); //you can change it if you want
    //setLastBinAsOverFlow(hist,0);
    //
    hist[i][j]->GetXaxis()->SetTitleSize(0.05);
    hist[i][j]->GetXaxis()->SetLabelSize(x_label_size);
    hist[i][j]->GetXaxis()->SetLabelSize(0.0450);
    //    hist[i][j]->GetXaxis()->SetRange(0,1500);
    hist[i][j]->GetYaxis()->SetTitleSize(0.05);
    hist[i][j]->GetYaxis()->SetLabelSize(0.05);
    hist[i][j]->GetYaxis()->SetTitleOffset(1.2);
    hist[i][j]->GetYaxis()->SetLabelSize(x_label_size);
    hist[i][j]->GetXaxis()->SetTitle(xtitile); //setting the title of X axis

    hist[i][j]->GetXaxis()->SetTitleSize(0.05);
    hist[i][j]->GetXaxis()->SetLabelSize(0.04);
    hist[i][j]->GetYaxis()->SetLabelSize(0.04);
    hist[i][j]->GetYaxis()->SetTitleSize(0.05);
    hist[i][j]->GetYaxis()->SetTitleOffset(1.1);
    hist[i][j]->GetXaxis()->SetTitleOffset(1.1);

    TLatex* textOnTop = new TLatex();
  textOnTop->SetTextSize(0.04);

  textOnTop->DrawLatexNDC(0.13,0.96,"CMS it{#bf{Simulation Preliminary}}");
  sprintf(lhead,"%s ",title);
  textOnTop->DrawLatexNDC(0.79,0.96,lhead);
  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.035);
    if(DoRebin) { //if rebin flag is on - this will reduce the bin size by half
     hist[i][j]->Rebin(rebin);
    }
    //setting up the legend style and all
       vector <string> Intg;
    char final_label[100];
    //int Entries = hist[i][j]->GetEntries().setprecision(1)
    int Entries = hist[i][j]->GetEntries();
      //sprintf(final_label,"Integral: %f",setprecision(2) << (hist[i][j]->Integral()));
    sprintf(final_label,"%s (%d)",legend_text.Data(), Entries); 

    
    //leg_entry[i] = legend->AddEntry(hist,(hist[i][j]->Integral().c_str()),"xsl");

    legName.push_back(hist[i][j]->GetName());
    leg_entry = legend->AddEntry(hist[i][j], final_label, "l");
    //leg_entry1[i] = legend1->AddEntry(hist,final_label,"l");
    //leg_entry[i] = legend->AddEntry(hist,(hist[i][j]->Integral().c_str()),"xsl");
    leg_entry->SetTextColor(hist[i][j]->GetLineColor());
    //leg_entry1[i]->SetTextColor(hist[i][j]->GetLineColor());
    if(hist[i][j]->GetMaximum() > ymax) ymax = hist[i][j]->GetMaximum();
    if(hist[i][j]->GetMinimum() < ymin) ymin = hist[i][j]->GetMinimum();
    //setLastBinAsOverFlow(hist,0);
     
  
  if(ymin == 0.0) ymin = 1e-3;
  if(ymin<0.0) ymin = 1e-4;
  //  if(ymax<=10) ymax=10;
  //for(int i = 0; i < (int)hist.size(); i++) {
    if(!normalize) hist[i][j]->GetYaxis()->SetRangeUser(0.001,ymax*10);
    else
      hist[i][j]->GetYaxis()->SetRangeUser(0.00001,ymax*60.0);
    //    p1->SetGrid();
    
     hist[i][j]->Draw("hist");
     cout << "hist Name: " << hist[0][1]->GetName() << endl;   	
    }
  }
  
  legend->Draw();
  //legend1->Draw();


  if(log_flag) 
    gPad->SetLogy();
  
  if(logx)
    gPad->SetLogx();
  
  gPad->Update();

 
  TLatex* textOnTop = new TLatex();
  textOnTop->SetTextSize(0.04);
  textOnTop->DrawLatexNDC(0.12,0.96,"CMS #it{#bf{Preliminary}}");
  
  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.04);
  float inlumi=energy;
  sprintf(en_lat,"#bf{%0.2f fb^{-1} (13 TeV)}",inlumi);
  textOnTop->DrawLatexNDC(0.7,0.96,en_lat);

 

  char* canvas_name = new char[1000];
  //c->Print(canvas_name);
  //saving the file
  if(save_canvas) {
    sprintf(canvas_name,"%s.png",tag_name);//.png",tag_name);//_wnormalize.png",tag_name);
    canvas_n1->SaveAs(canvas_name);   
    // sprintf(canvas_name,"%s.pdf",tag_name);
    // canvas_n1->SaveAs(canvas_name);
 
  }
  
}
const int nfiles=100,nBG=6;                                                                                                                                                              
//TFile *f[nfiles];


void single(string pathname)
{
  char* hname_GenPt = new char[200];
  char* hist_name  = new char[200];
  char* hist_name1 = new char[200];
  char* hist_name2 = new char[200];
  char* hist_name3 = new char[200];
  char* hist_name4 = new char[200];
  char* hist_name5 = new char[200];
  char* hist_name6 = new char[200];
  char* hist_name7 = new char[200];
  char* full_path = new char[2000];
  char* full_path1 = new char[2000];
  char* full_path2 = new char[2000];
  char* path2 = new char[2000];
  char* title= new char[2000];
  //string filetag;//=new char[20000];                                                                                                                                                                   
  char* full_path3 = new char[2000];
  char* full_path4 = new char[2000];
  char* full_path5 = new char[2000];
  char* full_path6 = new char[2000];
  char* full_path7 = new char[2000];
  char* full_path8 = new char[2000];
  char* full_path9 = new char[2000];
  char* full_path10 = new char[2000];
  char* full_path11= new char[2000];
  char *leg_head = new char[200];
  //define your files here
  vector<string> f;
  //f = {"Summer20UL18_TTJets_HT.root", "Summer20UL18_TTGJets_Tune.root", "Summer20UL18_TTJets_Leptons.root", "Summer20UL18_WJetsToLNu_HT.root", "Summer20UL18_WGJets_MonoPhoton.root", "Summer20UL18_ZJetsToNuNu_HT.root", "Summer20UL18_ZNuNuGJets_MonoPhoton.root", "Summer20UL18_QCD_HT.root"};
  f={"Summer20UL18_WJetsToLNu_HT.root"};
 
  //define your histograms to be read from here
  int n_files=f.size(); //you have n files in this example
 
  string histname[10][10];//, histname2[100], histname3[100];
  // char hname_NHadJets[100], hname_Jet_Pt[100], hname_Jet_Eta[100], hname_Jet_Phi[100], hname_Met[100], hname_PhoPt[100], hname_PhoEta[100], hname_PhoPhi[100];
  
  // Book your histograms & summary counters here
    vector<string> selection = {"ST", "Lep_veto", "Iso_Lep_Trk_veto"};
    vector<string> GenParticleId = {"Electron","Muon","Tau"};
  for (int i=0; i<selection.size();i++)
    {
      for (int j=0; j<GenParticleId.size(); j++){
	sprintf(hname_GenPt,"h_Gen%s_Pt_%s",GenParticleId[j].c_str(),selection[i].c_str());
      histname[i][j] = hname_GenPt;
      // sprintf(hname_Met, "h_MET_%s",selection[i].c_str());
      // histname2[i] = hname_Met; 
      // sprintf(hname_PhoPt, "h_Pho_Pt_%s",selection[i].c_str());
      // histname1[i] = hname_PhoPt; 
      // sprintf(hname_PhoEta, "h_Pho_Eta_%s",selection[i].c_str());
      // sprintf(hname_PhoPhi, "h_Pho_Phi_%s",selection[i].c_str());

      //cout << "HISTOGram: " << histname3[i] << endl;
      }
    }

  //vector<vector<string>> bigbaseline;
  vector<string> baseline1; //, baseline2, baseline3;

  //baseline1 = {"h_pho_pT0","h_pho_pT3","h_pho_pT4", "h_pho_pT5"};
  //  baseline1 = {histname1[0], histname1[1], histname1[2]}//, histname1[3], histname1[4], histname1[5], histname1[6]};
  //baseline2 = {"h_pho_eta0","h_pho_eta3"};s
  //baseline3 = {"h_pho_phi0","h_pho_phi3"};
  //baseline4 = {"h_Jet_pT2", "h_Jet_pT3"};
  //baseline5 = {"h_Jet_Eta2", "h_Jet_Eta3"};
  //baseline6 = {"h_NHadJets2", "h_NJets3"};
  // baseline2 = {histname2[0], histname2[1], histname2[2], histname2[3], histname2[4], histname2[5], histname2[6]}; 
  // baseline3 = {histname3[0], histname3[1], histname3[2], histname3[3], histname3[4], histname3[5], histname3[6]};
  //baseline4 = {"h_NHadJets0", "h_NHadJets3", "h_NHadJets4", "h_NHadJets5"};
  


  //string to be added to output file name - useful when you have different files and reading the same histograms from these
  //bigbaseline = {baseline1, baseline2, baseline3};

  // for (int bigi=0; bigi<bigbaseline.size(); bigi++)
  //   {
      vector<string> filetag;
      //filetag={"TTJets_2018","TTGJets_2018", "TTJets_Leptons_2018", "WJetsToLNu_2018", "WGJets_2018", "ZJetsToNuNu_2018", "ZNuNuGJets_2018", "QCD_2018"};
      filetag={"WJetsToLNu_2018"};
      //luminosity for each year - depends if you want to use it or - generate1Dplot uses this number and add it on the top
      vector<float>energyy;
      energyy={59.74,41.529};//,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,19.5,19.5,19.5,19.5,19.5,19.5,19.5,36.0,36.0,36.0,36.0,36.0,36.0,36.0};

      //rebin values
      vector<int >rebin = {2,2,2,2,2,2,2,2,2}; //keep it 1 if you don't want to change hist bins
      
      //x axis title for all plots///
      // vector<string>xtitle = {title[bigi], title[bigi]};
      //x axis range
      vector<int>xmax = {2000,2000,2000,2000,2000,2000,2000,2000};
      vector<int>xmin = {0,0,0,0,0,0,0,0};

      //looping over each files///  
      for(int i_file=0; i_file < f.size(); i_file++) //looping over each file
	{ TFile *root_file = new TFile(f[i_file].c_str()); 
	  TH1F* hist_list[10][10];
	  for(int i=0; i<GenParticleId.size();i++) 
	    {
	      for(int j=0; j<selection.size(); j++)
		{
	    //sprintf(hist_name,"%s",baseline[i_cut].c_str());
		sprintf(hist_name,"%s",histname[i][j].c_str());
	       cout<<"i_file "<<i_file<<"\t"<<i<<"\t"<<root_file->GetName()<< "\t" << hist_name<<endl;
	      TH1F* resp = (TH1F*)root_file->Get(hist_name); //reading hist from the TFile
	      hist_list[i][j] = resp;
	      //cout << "random: " << hist_list[i][j]-> Integral() << endl;
	      //}
	  float energy=energyy[0];
	  int xrange=0.0;

	   //x axis title for all plots///
	  //vector<string>diff_title;
	  string diff_title = "Pt";
	    //diff_title = { "Pho_pT" , "Pt_Miss", "NHadJets"};
	  vector<string>xtitle;
	  xtitle = {diff_title, diff_title, diff_title, diff_title, diff_title, diff_title, diff_title, diff_title};
	 
	  //path to save the files a jpg or pdf
	  vector<string> folder;
	  //folder = {"plots/TT/gen/", "plots/TTG/gen/", "plots/TTL/gen/",  "plots/WLNu/gen/", "plots/WG/gen/", "plots/ZNuNu/gen/", "plots/ZGNuNu/gen/", "plots/QCD/gen/"};
	  folder = {"plots/WLNu/"};
	  sprintf(full_path,"%s/%s%s_MET_comparisons_%s",pathname.c_str(),folder[i_file].c_str(),diff_title.c_str(),filetag[i_file].c_str());
	  //cout << "random: " << hist_list[i][j]-> Integral() << endl;
	  //calling generate_1Dplot which will take this vector of histograms and 
	  generate_1Dplot(hist_list,full_path,energy,xmax[i_file],xmin[i_file],leg_head,false,true,false,true,filetag[i_file].c_str(),xtitle[i_file].c_str(),rebin[i_file]);
		}
	    }
	  
	}
}






