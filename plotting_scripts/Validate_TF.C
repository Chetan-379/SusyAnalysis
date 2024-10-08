const int n_pl = 4;
bool logx = false;
//defining the legends for each plots
TString legend_text[10] = {"expected SR","predicted SR", "hist4", "pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};
int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
int line_color[9] = {kBlack, kRed, kGreen+2, kMagenta, kBlue, kAzure + 7 , kCyan + 1 , kGreen + 3 };
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
void generate_1Dplot(vector<TH1F*> hist, TH1* hist_ratio,char const *tag_name="", float energy=-1, int xmax=-1,int xmin=-1,char const *leg_head="",
		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title="", const char* xtitile="", int rebin=-1){  
  

  TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name,900,750);//600,600,1200,1200);
  canvas_n1->Range(-60.25,-0.625,562.25,0.625);
  canvas_n1->SetFillColor(0);
  canvas_n1->SetBorderMode(0);
  canvas_n1->SetBorderSize(2);
  // canvas_n1->SetRightMargin(0.045);
  // canvas_n1->SetLeftMargin(0.12);
  // canvas_n1->SetTopMargin(0.06);
  // canvas_n1->SetBottomMargin(0.12);
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

  //  auto *pad_1 = new TPad("pad_1","pad_1",0.,0.0,1.,0.32); pad_1->Draw();
  auto *pad_1 = new TPad("pad_1","pad_1",0.,0.0,1.,0.32); pad_1->Draw();
  pad_1->SetTopMargin(0.04);
  //pad_1->SetTopMargin(0.001);
  pad_1->SetBottomMargin(0.33);
  pad_1->SetRightMargin(0.035);
  pad_1->SetLeftMargin(0.13);
  //auto *p1 = new TPad("p1","p1",0.,0.32,1.,1.);  p1->Draw();
  auto *p1 = new TPad("p1","p1",0.,0.32,1.,1.);  p1->Draw();
  //p1->SetBottomMargin(0.026);
  p1->SetBottomMargin(0.040);  
  p1->SetRightMargin(0.035);
  p1->SetLeftMargin(0.13);
  p1->SetTopMargin(1);

  p1->cd();
  //p1->SetGridx();
  //pad_1->SetGridx();
  
  vector<TString> legName;
  std::string leg_head_str = leg_head;
  double x = 0.15;
  double y = 0.90;

  TLegend *legend; //legend to be drawn on the plot - shift x,ys if you want to move this on the canvas
  legend = new TLegend(0.18,0.82,0.72,0.88);  
  legend->SetTextSize(0.030);
  legend->SetLineColor(kWhite);
  legend->SetNColumns(2);
  char* lhead = new char[100];
  
  sprintf(lhead,"#bf{%s} ",title); // legend header -- for example you are plotting this for WGJets then title string can be "WGJets"
  legend->SetHeader(lhead);
  legend->SetLineColor(kWhite);


  TLegendEntry* leg_entry[11];
  TLegendEntry* leg_entry1[11];
  float x_label_size = 0.045;
  double ymin = 100000.0;
  double ymax = 0.0;
  double xrange = xmax;

  cout<<" hist.size() = "<<hist.size()<<endl;
 
  for(int i = 0; i < (int)hist.size(); i++) {
    if(normalize) {
      	hist.at(i)->Scale(1.0/hist.at(i)->Integral());
      	hist.at(i)->GetYaxis()->SetTitle("Normalized");
	
    }
    else {
      hist.at(i)->GetYaxis()->SetTitle("Entries");
       }
    hist.at(i)->GetXaxis()->SetRangeUser(xmin,xmax); //to be given while calling this function
    hist.at(i)->SetLineWidth(line_width[i]); //these are defined on top of this script    
    hist.at(i)->SetLineStyle(line_style[i]);
    hist.at(i)->SetLineColor(line_color[i]);
    hist.at(i)->SetTitle(" "); //you can change it if you want
    //setLastBinAsOverFlow(hist.at(i),0);
    //
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(x_label_size);
    hist.at(i)->GetXaxis()->SetLabelSize(0.0450);
    //    hist.at(i)->GetXaxis()->SetRange(0,1500);
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.0);
    hist.at(i)->GetYaxis()->SetLabelSize(x_label_size);
    hist.at(i)->GetXaxis()->SetTitle(xtitile); //setting the title of X axis

    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(0.04);
    hist.at(i)->GetYaxis()->SetLabelSize(0.04);
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.1);
    hist.at(i)->GetXaxis()->SetTitleOffset(1.1);

    // //labelling the bins
    // vector<string> binlabel{"(2,0)","(3,0)","(4,0)","(5-6,0)","(>=7,0)","(2-4,1)","(5-6,1)","(>=7,1)"};
    // for (int j=0; j < 8; j++){
    //   hist.at(i)->GetXaxis()-> SetBinLabel(j+2,binlabel[j].c_str());}

    TLatex* textOnTop = new TLatex();
    textOnTop->SetTextSize(0.04);
    
    //textOnTop->DrawLatexNDC(0.13,0.96,"CMS it{#bf{Simulation Preliminary}}");
    textOnTop->DrawLatexNDC(0.13,0.96,"it{#bf{Work in Progress}}");
    
    sprintf(lhead,"%s ",title);
    //textOnTop->DrawLatexNDC(0.79,0.96,lhead);
    textOnTop->DrawLatexNDC(0.79,0.4,lhead);
    char* en_lat = new char[500];
    textOnTop->SetTextSize(0.035);
    if(DoRebin) { //if rebin flag is on - this will reduce the bin size by half
      hist.at(i)->Rebin(rebin);
    }
    
    if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    
    
    //setting up the legend style and all
    vector <string> Intg;
    char final_label[100];
    float Integral = hist.at(i)->Integral();
    sprintf(final_label,"%s (%0.2f)",legend_text[i].Data(), Integral); 
    legName.push_back(hist.at(i)->GetName());
    leg_entry[i] = legend->AddEntry(hist.at(i), final_label, "l");
    //leg_entry1[i] = legend1->AddEntry(hist.at(i),final_label,"l");
    //leg_entry[i] = legend->AddEntry(hist.at(i),(hist.at(i)->Integral().c_str()),"xsl");
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor());
    //leg_entry1[i]->SetTextColor(hist.at(i)->GetLineColor());
    // if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    // if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    //setLastBinAsOverFlow(hist.at(i),0);
    
  }
  

  if(ymin == 0.0) ymin = 1e-3;
  if(ymin<0.0) ymin = 1e-4;
  //  if(ymax<=10) ymax=10;
  for(int i = 0; i < (int)hist.size(); i++) {
    //if(!normalize) hist.at(i)->GetYaxis()->SetRangeUser(ymin*10,ymax*(ymax/5));
    if(!normalize) hist.at(i)->GetYaxis()->SetRangeUser(0.001,ymax*100);
    else
      hist.at(i)->GetYaxis()->SetRangeUser(ymin,ymax*10);
    //    p1->SetGrid();
    
    if(!i) hist.at(i)->Draw();
    else   hist.at(i)->Draw("sames"); //overlaying the histograms
	
  }
  
  legend->Draw();

  if(log_flag) 
    gPad->SetLogy();
  
  if(logx)
    gPad->SetLogx();
  
  gPad->Update();

 
  TLatex* textOnTop = new TLatex();
  textOnTop->SetTextSize(0.04);
  //textOnTop->DrawLatexNDC(0.12,0.96,"CMS #it{#bf{Preliminary}}");
  textOnTop->DrawLatexNDC(0.12,0.91," #it{#bf{Work in Progress}}");
  
  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.04);
  float inlumi=energy;
  sprintf(en_lat,"#bf{%0.2f fb^{-1} (13 TeV)}",inlumi);
  textOnTop->DrawLatexNDC(0.7,0.91,en_lat);

  TLine *line1V7=new TLine( 8.0,0.001, 8.0,ymax*5);
  TLine *line2V7=new TLine(14.0,0.001, 14.0,ymax*5);
  TLine *line3V7=new TLine(19.0,0.001, 19.0,ymax*5);
  TLine *line4V7=new TLine(24.0,0.001, 24.0,ymax*5);
  TLine *line5V7=new TLine(29.0,0.001, 29.0,ymax*5);

  
  line1V7->Draw();      line2V7->Draw();  line3V7->Draw();
  line4V7->Draw();      line5V7->Draw();           

  TArrow *arrow1 = new TArrow( 1.0,ymax*5, 8.0,ymax*5,0.01,"<|>");
  TArrow *arrow2 = new TArrow( 8.0,ymax*5, 14.0,ymax*5,0.01,"<|>");
  TArrow *arrow3 = new TArrow(14.0,ymax*5, 19.0,ymax*5,0.01,"<|>");
  TArrow *arrow4 = new TArrow(19.0,ymax*5, 24.0,ymax*5,0.01,"<|>");
  TArrow *arrow5 = new TArrow(24.0,ymax*5, 29.0,ymax*5,0.01,"<|>");
  TArrow *arrow6 = new TArrow(29.0,ymax*5, 35.0,ymax*5,0.01,"<|>");
    
  arrow1->Draw(); arrow2->Draw(); arrow3->Draw();
  arrow4->Draw(); arrow5->Draw(); arrow6->Draw();
  
  TLatex Tl;
  Tl.SetTextSize(0.035);
  Tl.DrawLatex(3.5,ymax*7.0,"N^{ 0}_{ 2-4}");
  Tl.DrawLatex(9.5,ymax*7.0,"N^{ 0}_{ 5-6}");
  Tl.DrawLatex(15.5,ymax*7.0,"N^{ 0}_{ #geq7}");
  Tl.DrawLatex(19.5,ymax*7.0,"N^{ #geq1}_{ 2-4}");
  Tl.DrawLatex(25.5,ymax*7.0,"N^{ #geq1}_{ 5-6}");
  Tl.DrawLatex(30.5,ymax*7.0,"N^{ #geq1}_{ #geq7}");
 
  hist_ratio->SetLineWidth(2);
  hist_ratio->SetLineStyle(1);
  hist_ratio->SetMarkerSize(0.2);
  hist_ratio->SetLineColor(kBlack);
  hist_ratio->SetTitle(" ");  
  hist_ratio->GetXaxis()->SetTitleSize(0.13);
  hist_ratio->GetYaxis()->SetTitle("exp/pred");//TF = #frac{N_{SR}}{N_{CR}}");//(0#mu,1#gamma)}{(1#mu,1#gamma)}");
  hist_ratio->GetXaxis()->SetLabelSize(0.1);
  hist_ratio->GetYaxis()->SetRangeUser(-0.5,2.5);
  hist_ratio->GetXaxis()->SetRangeUser(xmin,xrange+4);
  // if(which_TFbins==1) //default 8 bins                                                                                                    
    //   hist_ratio->GetXaxis()->SetRangeUser(0,10);//xmin,xrange);                                                                                                 
    // else if(which_TFbins==2) // v2 TF bins including photon pT>100 and pT<100
    //   hist_ratio->GetXaxis()->SetRangeUser(0,18);
    // else if(which_TFbins==3) // v3 TF bins including MET<300 and MET>300                                                                                        
    //    hist_ratio->GetXaxis()->SetRangeUser(0,39);
    
    //    hist_ratio->GetXaxis()->SetLabelSize(0.0450);
  hist_ratio->GetYaxis()->SetTitleSize(0.13);
  hist_ratio->GetYaxis()->SetLabelSize(0.08);
  hist_ratio->GetYaxis()->SetTitleOffset(1.0);
  hist_ratio->SetMarkerSize(1.0);
  hist_ratio->SetMarkerStyle(20);
  hist_ratio->SetMarkerColor(kBlue);
  hist_ratio->GetXaxis()->SetTitle(xtitile);
  hist_ratio->GetYaxis()->SetNdivisions(505);
  
  hist_ratio->GetYaxis()->CenterTitle(true);
  hist_ratio->GetXaxis()->SetTitleSize(0.05);
  hist_ratio->GetXaxis()->SetLabelSize(0.12);
  hist_ratio->GetYaxis()->SetTitleSize(0.12);
  hist_ratio->GetYaxis()->SetNdivisions(505);
  
  hist_ratio->GetXaxis()->SetTitleOffset(1);
  hist_ratio->GetYaxis()->SetTitleOffset(0.41);
  hist_ratio->GetXaxis()->SetTitleSize(0.14);
  
  hist_ratio->GetYaxis()->SetLabelSize(0.12);

  //    hist_ratio->GetYaxis()->SetLabelSize(x_label_size);
  //hist_ratio->SetOptStat(0);
  
  pad_1->cd();
  gStyle->SetOptStat(0);

  
  
  TLine *l =new TLine(xmin,1.0,xrange+4,1.0);
  hist_ratio->Draw("");
  l->Draw("sames");
  TLine *l1 =new TLine(xmin,1.5,xrange+4,1.5);
  l1->SetLineStyle(7);
  l1->Draw("sames");
  TLine *l2 =new TLine(xmin,0.5,xrange+4,0.5);
  l2->SetLineStyle(7);
  
  l2->Draw("sames");
  
  char* canvas_name = new char[1000];
  //c->Print(canvas_name);
  
  if(save_canvas) {
    sprintf(canvas_name,"%s.png",tag_name);//.png",tag_name);//_wnormalize.png",tag_name);
    canvas_n1->SaveAs(canvas_name);   
    // sprintf(canvas_name,"%s.pdf",tag_name);
    // canvas_n1->SaveAs(canvas_name);
    
  }

  
}
const int nfiles=100,nBG=6;                                                                                                                                                              
//TFile *f[nfiles];


void Validate_TF(string pathname)
{
  char* hname = new char[200];
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
  vector<string> filetag;
  vector<string> folder;
  string identifier;
  
  int combine;
  combine =4;

  if (combine==0) {
    f = {"./root_files/Summer20UL18_TTJets.root", "./root_files/Summer20UL18_TTGJets_Tune.root", "./root_files/Summer20UL18_WJetsToLNu_HT.root", "./root_files/Summer20UL18_WGJets_MonoPhoton.root"};
    //f = {"./root_files/Summer20UL18_WGJets_MonoPhoton.root"};
    filetag={"TTJets_2018", "TTGJets_2018", "WJets_2018", "WGJets_2018"};
    //filetag={"WGJets_2018"};
    folder = {"plots/LL_plots/Lost_mu_had_tau/TT_/", "plots/LL_plots/Lost_mu_had_tau/TTG_/", "plots/LL_plots/Lost_mu_had_tau/WLNu_/", "plots/LL_plots/Lost_mu_had_tau/WG_/"};
    //folder = {"plots/LL_plots/Lost_mu_had_tau/WG_/"};
    identifier = "Single";
  }

  if (combine==2) {
    f = {"./root_files/TTGJets_process_combined_Summer20UL18.root", "./root_files/WGJets_process_combined_Summer20UL18.root"};
    filetag={"TTGJets_Process_2018", "WGJets_Process_2018"};
    folder = {"plots/LL_plots/combined/", "plots/LL_plots/combined/"};
    identifier = "Combined";
  }

  if (combine==4) {
    f= {"./root_files/All_combined_Summer20UL18.root"};
    filetag={"TTG+WG_2018"};
    folder = {"plots/LL_plots/All_cat_combined/combined/"};
    //folder = {"plots/LL_plots/Lost_mu_had_tau/combined/"};
    //folder = {"plots/LL_plots/Lost_e/combined/"};
    identifier = "Combined";
  }
  
  vector<string> histnames;
  //histnames = {"lost_e_SR_srch_binned", "lost_e_SR_srch_binned_pred", "lost_e_CR_srch_binned"};
  //histnames = {"lost_mu_SR_srch_binned", "lost_mu_SR_srch_binned_pred", "lost_mu_CR_srch_binned"};
  histnames = {"LL_SR_srch_binned", "LL_SR_srch_binned_pred", "LL_CR_srch_binned_pred"};
    
  //luminosity for each year - depends if you want to use it or - generate1Dplot uses this number and add it on the top
  vector<float>energyy;
  energyy={59.74,41.529};//,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,19.5,19.5,19.5,19.5,19.5,19.5,19.5,36.0,36.0,36.0,36.0,36.0,36.0,36.0};
  
  //rebin values
  vector<int >rebin = {2,2,2,2,2,2,2,2,2}; //keep it 1 if you don't want to change hist bins
  
  //x axis title for all plots///
  vector<int>xmax = {2000,2000,2000,2000,2000,2000,2000,2000};
  vector<int>xmin = {0,0,0,0,0,0,0,0};            
  vector<vector<TH1F*>> hist_file;
 
  //looping over each files///  
  for(int ibkg=0; ibkg < f.size(); ibkg++){ //looping over each file
    vector<TH1F*> hist_list;
    TFile *root_file = new TFile(f[ibkg].c_str()); 
    for(int i_hist=0; i_hist < 2;i_hist++) //looping over different histograms which should be overlayed on the same canvas and these histograms are saved in the same file	    
      {
	sprintf(hist_name,"%s",histnames[i_hist].c_str());
	cout<<"i_file "<< ibkg <<"\t"<<i_hist<<"\t"<<root_file->GetName()<< "\t" << hist_name<<endl;
	TH1F* resp = (TH1F*) root_file->Get(hist_name); //reading hist from the TFile
	hist_list.push_back(resp);	      
      }
    cout << "all samples root file combine SR:" << hist_list[0]->Integral() << endl;
	  
    float energy=energyy[0];
    int xrange=0.0;
    //vector<TH1D*> h_TF;
    
    //TH1D* h_CR = (TH1D*)hist_add[ibkg].at(1)->Clone();
    TH1D* h_TF = (TH1D*)hist_list.at(0)->Clone();
    h_TF->Divide(hist_list.at(1));
    //x axis title for all plots///
    vector<string>diff_title;
    diff_title = {"Bin_No."};
    vector<string>xtitle;
    xtitle = {diff_title[0], diff_title[0], diff_title[0], diff_title[0]};
    
    sprintf(full_path,"%s/%s%s_LL_Validate_binned_SR_CR_overlay_%s",pathname.c_str(),folder[ibkg].c_str(),identifier.c_str(),filetag[ibkg].c_str());
    
    
    // 	//calling generate_1Dplot which will take this vector of histograms and 
    //generate_1Dplot(h_added,h_TF,full_path,energy,xmax[0],xmin[0],leg_head,false,true,false,true,filetag.c_str(),xtitle[0].c_str(),rebin[0]);
    generate_1Dplot(hist_list, h_TF,full_path,energy,xmax[ibkg],xmin[ibkg],leg_head,false,true,false,true,filetag[ibkg].c_str(),xtitle[ibkg].c_str(),rebin[ibkg]);
  }
}









