const int n_pl = 4;
bool logx = false;
TString pls[n_pl] = {"FTFP_BERT_EMN","QGSP_FTFP_BERT_EMN","FTFP_BERT_EMM","QGSP_BERT"};
TString dRs[9] = {"dR < 0.56cm","dR < 1.0cm","dR < 2.0cm","dR < 3.0cm","dR < 5.0cm","dR < 8.0cm","dR < 12.0cm","dR < 18.0cm","dR < 20.0cm"};
//TString legend_text[11] = {"No cuts","skimmed","lep-veto","isotrk-veto","Pho-Pt>20","Njets>=2","Dphi-cut","MET>100","MET>250","ST>300","Pho-pt>100"};
TString legend_text[5] ={"#tau-had SR","lost #mu SR","lost e SR","(1l,1#gamma) CR","Failed Iso"};//,"Failed acceptance","1e CR"};// {"(0#mu,0e) SR","(1#mu,0e) CR","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};//{"No cuts","skimmed","lep-veto","isotrk-veto","Dphi-cut","MET>250","ST>300","Pho-pt>100"};
//TString legend_text[4] = {"(0#mu,0e) SR","(1#mu,0e) CR","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451"};
int line_width[12] = {2,2,2,2,2,2,2,2,2,2,2,2};
int line_style[12] = {1,1,1,1,1,1,1,1,1,1,1,1};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};
// /int line_color[n_pl+1] = {kBlack, kRed, kGreen+2, kBlue, kRed};
int line_color[9] = {kBlue,kRed,kGreen+2,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};//{9,kCyan+2,45,kMagenta,kGray+1,kRed,kBlue+2,kMagenta,kCyan};        
int line_color1[9]= {kBlue,kRed,kGray+1,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2};
int line_color2[9] = {kRed,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta};
//int line_color[9] = {kMagenta+2, kGray+2, kRed, kGreen+2, kMagenta, kRed - 3, kBlue + 2 , kCyan + 1 , kGreen + 3 };                                                 
vector<int> col={kRed,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
vector<int> Style={3008,1001,3008,1001};
//int line_color[11] = {kPink+1, kRed, kBlue,kGray+1 , kGreen+2, kMagenta, kYellow + 2 , kCyan+3,  kBlue + 2 ,kRed+2,kGreen + 3 };
void decorate(TH1F*,int,int );
 
void decorate(TH1F* hist,int i, int j){
  //  hist->SetLineColor(col[i]);
  //   hist->SetFillColor(col[i]);
   // vector<int> col;
   // vector<int> Style;
   // if(j==1){
   //   col={kGray+1,kPink-9,kViolet,kGreen+2,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
   //   Style= {1001,3023,3019,3244};

   // }
   // else if(j==2){
   //   col={kGray+1,kPink-9,kBlue,kGreen+2,kBlue,kGreen+2,kGray+1,kPink-9,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
   //   Style={1001,3023,3008,1001,3244};
   // }
   // else if(j==3){
   //   col={kGray,kPink-9,kViolet,kBlue,kGreen+2,kGreen+2,kBlue,kViolet,kGray,kViolet+2,kGreen-2,kYellow+1,kGray+2,kMagenta,kBlue+2,kMagenta,kCyan};
   //   Style={1001,3023,1001,3008,1001};
   // }
   //  // if(i!=4)
   //  //   {
   // hist->SetFillColor(col[i]);
   
   // hist->SetFillStyle(Style[i]);
    //   }
    // else 
    //   {
   // hist->SetFillColor(kGray+1);
   
   //hist->SetFillStyle(1001);
   //     }
    hist->SetLineWidth(3);

  //  if(i<nBG) {                                                                                                                                 
  //hist->SetFillColor(col[i]);
    // if(i!=0)
    //   {
    // 	//	hist->SetFillColor(col[i]);

    //  	hist->SetFillStyle(Style[i]);
    // //   }
    // // else
    // //   {
    // // 	hist->SetFillColor(kGray+2);

    // // 	hist->SetFillStyle(1001);
    // // 	}
    // hist->SetLineWidth(2);
    //}                                                                                                                                           
}

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

void generate_1Dplot(vector<TH1F*> hist, TH1* hist_ratio, char const *tag_name="", float energy=-1, int xrange=-1,int xmin=-1,char const *leg_head="",
		     bool normalize=false, bool log_flag=true, bool DoRebin=false, bool save_canvas=true, char const *title="", vector<string> legend_texts={"nil"}, int which_TFbins=-1, int which_Lept=-1){  
  

   TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name,950,850);
       canvas_n1->Range(-60.25,-0.625,562.25,0.625);
       canvas_n1->SetFillColor(0);
       canvas_n1->SetBorderMode(0);
       canvas_n1->SetBorderSize(2);
       auto *pad_1 = new TPad("pad_1","pad_1",0.,0.0,1.,0.32); pad_1->Draw();
       pad_1->SetTopMargin(0.04);
       pad_1->SetBottomMargin(0.33);
       pad_1->SetRightMargin(0.025);
       pad_1->SetLeftMargin(0.14);
         auto *p1 = new TPad("p1","p1",0.,0.32,1.,1.);  p1->Draw();
      p1->SetBottomMargin(0.026);
       p1->SetRightMargin(0.025);
       p1->SetLeftMargin(0.14);
       p1->SetTopMargin(0.1);
       p1->cd();


   // TCanvas *canvas_n1 = new TCanvas(tag_name, tag_name,900,750);//600,600,1200,1200);
   // canvas_n1->Range(-60.25,-0.625,562.25,0.625);
   // canvas_n1->SetFillColor(0);
   // canvas_n1->SetBorderMode(0);
   // canvas_n1->SetBorderSize(2);
   // auto *pad_1 = new TPad("pad_1","pad_1",0.,0.0,1.,0.32); pad_1->Draw();
   // pad_1->SetTopMargin(0.013);
   // pad_1->SetBottomMargin(0.3);
   // pad_1->SetRightMargin(0.025);
   // pad_1->SetLeftMargin(0.14);
   
   // auto *p1 = new TPad("p1","p1",0.,0.32,1.,1.);  p1->Draw();
   // p1->SetBottomMargin(0.01);
   // p1->SetRightMargin(0.025);
   // p1->SetLeftMargin(0.14);
   // p1->SetTopMargin(0.05);
   // p1->cd();
   // p1->SetGrid();
   
//   TCanvas *c = new TCanvas(tag_name, tag_name, 700, 600);
//   c->SetRightMargin(0.009);
//   c->SetLeftMargin(0.11);
//   c->SetTopMargin(0.05);
// c->SetLeftMargin(0.11);
//   c->SetTopMargin(0.05);
//   c->SetRightMargin(0.025);
//   c->SetLeftMargin(0.13);
//   c->SetTopMargin(0.06);
//    c->SetBottomMargin(0.12);
     THStack *hs_var=new THStack("var_Stack","");
    gStyle->SetOptStat(1111111);
  //   gStyle->SetOptStat(0);
  //double pvt_x_min = 0.6;
  double pvt_x_min = 0.75;
  double pvt_x_max = 0.99;
  double pvt_y_min = 0.9;
  //double pvt_dely = 0.18;
  double pvt_dely = 0.15;
  gStyle->SetOptStat(0);

  //gStyle->SetOptFit(0);
  vector<TString> legName;
  //TLegend *legend = new TLegend(0.65,0.95,0.99,0.75);
  std::string leg_head_str = leg_head;
  double x = 0.15;
  double y = 0.90;
  TLegend *legend;
  //legend = new TLegend(0.60,0.88,0.98,0.72);  
  legend = new TLegend(0.2,0.78,0.55,0.88);  
  legend->SetTextSize(0.045);
  legend->SetLineColor(kWhite);
  legend->SetNColumns(4);
  char* lhead = new char[100];

  sprintf(lhead,"#bf{%s} ",title);
  legend->SetHeader(lhead);
  legend->SetLineColor(kWhite);

  TLegendEntry* leg_entry[11];
  float x_label_size = 0.045;
  double ymin = 100000.0;
  double ymax = 0.0;
  cout<<" hist.size() = "<<hist.size()<<endl;
  // if(hist.size()==2)
  //   {  col[1]=col[1]-41;
  //     legend_text[1] = "(0e,0#mu) SR";
  //   }
  // else
  //   {
  //     col[1]=45;
  //     legend_text[1] = "Failed acceptance";
  //   }
  for(int i =0; i<(int)hist.size();i++) {
    // if(DoRebin) {
    //  hist.at(i)->Rebin(2);

    // }
    hist.at(i)= setLastBinAsOverFlow(hist.at(i),xrange);
     

    //    normalize = true;
    if(normalize) {
      	hist.at(i)->Scale(1.0/hist.at(i)->Integral());
      	hist.at(i)->GetYaxis()->SetTitle("Normalized");
    }
    else {
      hist.at(i)->GetYaxis()->SetTitle("Entries");
    }
    if(which_TFbins==1) //default 8 bins
      hist.at(i)->GetXaxis()->SetRangeUser(0,10);//xmin,xrange);
    else if(which_TFbins==2) // v2 TF bins including photon pT>100 and pT<100
      hist.at(i)->GetXaxis()->SetRangeUser(0,50);
    else if(which_TFbins==3) // v3 TF bins including MET<300 and MET>300
      hist.at(i)->GetXaxis()->SetRangeUser(0,50);
    hist.at(i)->SetLineWidth(line_width[i]);    
    hist.at(i)->SetLineStyle(line_style[i]);
    if(which_Lept==1)
      hist.at(i)->SetLineColor(line_color[i]);
    else if(which_Lept==2)
      hist.at(i)->SetLineColor(line_color1[i]);
    else  if(which_Lept==3)
      hist.at(i)->SetLineColor(line_color2[i]);

    hist.at(i)->SetTitle(" ");
    //setLastBinAsOverFlow(hist.at(i),0);
    //
    hist.at(i)->GetXaxis()->SetTitleSize(0.05);
    hist.at(i)->GetXaxis()->SetLabelSize(x_label_size);
    hist.at(i)->GetXaxis()->SetLabelSize(0.0450);
    //    hist.at(i)->GetXaxis()->SetTitle("#frac{N_{SR/CR}}{N_{SR}+N_{CR}}");
    hist.at(i)->GetYaxis()->SetTitleSize(0.05);
    hist.at(i)->GetYaxis()->SetLabelSize(0.05);
    hist.at(i)->GetYaxis()->SetTitleOffset(1.2);
    hist.at(i)->GetYaxis()->SetLabelSize(x_label_size);

    //new ones
               hist.at(i)->GetXaxis()->SetTitleSize(0.08);
    hist.at(i)->GetXaxis()->SetLabelSize(0.06);

    hist.at(i)->GetYaxis()->SetTitleSize(0.07);
    hist.at(i)->GetYaxis()->SetLabelSize(0.06);

    hist.at(i)->GetXaxis()->SetTitleOffset(1.8);
    hist.at(i)->GetXaxis()->SetLabelOffset(1.6);

    hist.at(i)->GetYaxis()->SetTitleOffset(0.9);

    decorate(hist.at(i),i, which_Lept);
    hist.at(i)->GetYaxis()->SetNdivisions(506);
    
    // if(DoRebin) {
    //  hist.at(i)->Rebin(2);
    //   //hist.at(i)->Rebin(1);
    // }

  //     double lastBinCt =hist.at(i)->GetBinContent(hist.at(i)->GetNbinsX()),overflCt =hist.at(i)->GetBinContent(hist.at(i)->GetNbinsX()+1);
  // double lastBinErr=hist.at(i)->GetBinError(hist.at(i)->GetNbinsX()),  overflErr=hist.at(i)->GetBinError(hist.at(i)->GetNbinsX()+1);
  // if(lastBinCt!=0 && overflCt!=0)
  //   lastBinErr = (lastBinCt+overflCt)* (sqrt( ((lastBinErr/lastBinCt)*(lastBinErr/lastBinCt)) + ((overflErr/overflCt)*(overflErr/overflCt)) ) );

  // else if(lastBinCt==0 && overflCt!=0)
  //   lastBinErr = overflErr;
  // else if(lastBinCt!=0 && overflCt==0)
  //   lastBinErr = lastBinErr;
  // else lastBinErr=0;

  // lastBinCt = lastBinCt+overflCt;
  // hist.at(i)->SetBinContent(hist.at(i)->GetNbinsX(),lastBinCt);
  // hist.at(i)->SetBinError(hist.at(i)->GetNbinsX(),lastBinErr);
  //  hist.at(i)->GetXaxis()->SetRange(1, hist.at(i)->GetNbinsX() + 1);
    /* hist.at(i)->GetXaxis()->SetRangeUser(x_min[energy],x_max[energy]); */
    //    hist.at(i)= DrawOverflow(hist.at(i));
    legName.push_back(hist.at(i)->GetName());
    leg_entry[i] = legend->AddEntry(hist.at(i),legend_texts[i].c_str(),"f");
    leg_entry[i]->SetTextColor(hist.at(i)->GetLineColor()+2);
    
    if(hist.at(i)->GetMaximum() > ymax) ymax = hist.at(i)->GetMaximum();
    if(hist.at(i)->GetMinimum() < ymin) ymin = hist.at(i)->GetMinimum();
    //setLastBinAsOverFlow(hist.at(i),0);
    
    

  }
  if(ymin == 0.0) ymin = 1e-3;
  if(ymin<0.0) ymin = 1e-4;
  //  if(ymax<=10) ymax=10;
  for(int i =0;i< (int)hist.size(); i++) {
    if(!normalize) hist.at(i)->GetYaxis()->SetRangeUser(0.001,100.0*ymax);
    else
      {  hist.at(i)->GetYaxis()->SetRangeUser(0.000001,ymax*60.0);
	//	hist.at(i)->GetXaxis()->SetRangeUser(0,xmax_[i]);
	
      }
    //    p1->SetGrid();
    hs_var->Add(hist.at(i));
    // hs_var->SetMinimum(0.00001);
    // hs_var->SetMaximum(ymax*60);

    	if(!i) hist.at(i)->Draw("HIST");
    else   hist.at(i)->Draw(" HISTsames");
	
  }
//  hs_var->SetMinimum(0.0);
//  hs_var->SetMaximum(ymax+0.5);


//   hs_var->Draw("BAR HIST");
//   hs_var->Draw("HIST");
//   hs_var->GetXaxis()->SetRangeUser(0,xrange);
//  hs_var->GetXaxis()->SetTitle(title);
// hs_var->GetXaxis()->SetTitleOffset(1.0);
//   gPad->Modified(); gPad->Update();
//   hs_var->GetXaxis()->SetTitle(title);
//   //  hs_var->GetYaxis()->SetTitleSize(
//   hs_var->GetYaxis()->SetTitle("#frac{N_{SR/CR}}{N_{SR}+N_{CR}}");//hs_var->GetYaxis()->SetTitle("Events");
//   hs_var->SetTitle(0);
//    hs_var->GetYaxis()->SetTitleOffset(1.2);
// hs_var->GetXaxis()->SetTitleSize(00.05);
//     hs_var->GetXaxis()->SetLabelSize(0.04);
//     hs_var->GetYaxis()->SetLabelSize(0.04);
//     hs_var->GetYaxis()->SetTitleSize(00.055);
//     hs_var->GetYaxis()->SetTitleOffset(1.0);

  // hs_var->SetMinimum(0.0);
  // hs_var->SetMaximum(1.5);
  
  
  // // hs_var->Draw("BAR HIST");
  // // hs_var->Draw("HIST");
  // if(which_TFbins==1) //default 8 bins                                                                                                                                     
  //   hs_var->GetXaxis()->SetRangeUser(0,10);//xmin,xrange);                                                                                                    
  // else if(which_TFbins==2) // v2 TF bins including photon pT>100 and pT<100                                                                                                
  //   hs_var->GetXaxis()->SetRangeUser(0,18);
  // else if(which_TFbins==3) // v3 TF bins including MET<300 and MET>300                                                                                                     
  //   hs_var->GetXaxis()->SetRangeUser(0,18);
  // hs_var->GetXaxis()->SetTitle(title);
  // hs_var->GetXaxis()->SetTitleOffset(1.2);
  // gPad->Modified(); gPad->Update();
  // hs_var->GetXaxis()->SetTitle(title);
  // //  hs_var->GetYaxis()->SetTitleSize(
  // hs_var->GetYaxis()->SetTitle("#frac{N_{SR/CR}}{N_{SR}+N_{CR}}");//hs_var->GetYaxis()->SetTitle("Events");
  // hs_var->SetTitle(0);
  // hs_var->GetYaxis()->SetTitleOffset(1.2);
  // hs_var->GetXaxis()->SetTitleSize(00.05);
  // hs_var->GetXaxis()->SetLabelSize(0.04);
  // hs_var->GetYaxis()->SetLabelSize(0.04);
  // hs_var->GetYaxis()->SetTitleSize(00.055);
  // hs_var->GetYaxis()->SetTitleOffset(1.0);


  legend->Draw();

    //    legend->Draw();


  if(log_flag) {
      gPad->SetLogy();
    }
  // if(logx)
  //   gPad->SetLogx();

  
  gPad->Update();

 
  TLatex* textOnTop = new TLatex();
  textOnTop->SetTextSize(0.054);
  textOnTop->DrawLatexNDC(0.146,0.925,"CMS #it{#bf{Preliminary}}");
  
  char* en_lat = new char[500];
  textOnTop->SetTextSize(0.054);
  float inlumi=energy;
  sprintf(en_lat,"#bf{%0.2f fb^{-1} (13 TeV)}",inlumi);
  textOnTop->DrawLatexNDC(0.72,0.925,en_lat);
  
  TArrow *arrow1 = new TArrow( 1.0,0.10, 2.0,0.1,0.01,"<|>");
  TArrow *arrow2 = new TArrow( 2.0,0.10,3.0,0.1,0.01,"<|>");
  TArrow *arrow3 = new TArrow(3.0,0.10,4.0,0.1,0.01,"<|>");
  TArrow *arrow4 = new TArrow(4.0,0.10, 5.0,0.1,0.01,"<|>");
  TArrow *arrow5 = new TArrow(5.0,0.10, 6.0,0.1,0.01,"<|>");
  TArrow *arrow6 = new TArrow(6.0,0.10, 7.0,0.1,0.01,"<|>");
  TArrow *arrow7 = new TArrow(7.0,0.1, 8.0,0.1,0.01,"<|>");
  TArrow *arrow8 = new TArrow(8.0,0.1, 9.0,0.1,0.01,"<|>");
  TArrow *arrow9 = new TArrow( 9.0,0.10, 10.0,0.1,0.01,"<|>");
  TArrow *arrow10 = new TArrow( 10.0,0.10,11.0,0.1,0.01,"<|>");
  TArrow *arrow11 = new TArrow(11.0,0.10,12.0,0.1,0.01,"<|>");
  TArrow *arrow12 = new TArrow(12.0,0.10, 13.0,0.1,0.01,"<|>");
  TArrow *arrow13 = new TArrow(13.0,0.10, 14.0,0.1,0.01,"<|>");
  TArrow *arrow14 = new TArrow(14.0,0.10, 15.0,0.1,0.01,"<|>");
  TArrow *arrow15 = new TArrow(15.0,0.1, 16.0,0.1,0.01,"<|>");
  TArrow *arrow16 = new TArrow(16.0,0.1, 17.0,0.1,0.01,"<|>");
  TArrow *arrow17 = new TArrow(17.0,0.1, 18.0,0.1,0.01,"<|>");

  if(which_TFbins==1){
  arrow1->Draw(); arrow2->Draw(); arrow3->Draw();
  arrow4->Draw(); arrow5->Draw(); arrow6->Draw();
  arrow7->Draw(); arrow8->Draw();
  }
  // else if (which_TFbins==2 || which_TFbins==3){
  //   arrow2->Draw(); arrow3->Draw();
  //   arrow4->Draw(); arrow5->Draw(); arrow6->Draw();
  //   arrow7->Draw(); arrow8->Draw(); arrow9->Draw();arrow10->Draw(); arrow11->Draw(); arrow12->Draw();arrow13->Draw(); arrow14->Draw(); arrow15->Draw(); arrow16->Draw(); 
  // }
  // TLatex Tl;
  // Tl.SetTextSize(0.055);
  // Tl.SetTextColor(kBlack);
  // TLatex Tl1;
  //  Tl1.SetTextSize(0.03);
  //   Tl1.SetTextColor(kBlack);
  // if(which_TFbins==1){
  // Tl.DrawLatex(1.0,0.2,"N^{ 0}_{ 2}");
  
  // Tl.DrawLatex(2.0,0.2,"N^{ 0}_{ 3}");
  // Tl.DrawLatex(3.0,0.2,"N^{ 0}_{ 4}");
  // Tl.DrawLatex(4.0,0.2,"N^{ 0}_{ 5,6}");
  // Tl.DrawLatex(5.0,0.2,"N^{ 0}_{#geq7}");
  // Tl.DrawLatex(6.0,0.2,"N^{ #geq1}_{ #geq2}");
  // Tl.DrawLatex(7.0,0.2,"N^{ #geq1}_{ 5,6}");
  // Tl.DrawLatex(8.0,0.2,"N^{ #geq1}_{ #geq7}");
  // }
  // else if (which_TFbins==2){
  //   Tl.SetTextSize(0.035);
  //   Tl.DrawLatex(1.0,0.2,"N^{ 0}_{ 2}");
  //   Tl.DrawLatex(2.0,0.2,"N^{ 0}_{ 2}");
  //   Tl.DrawLatex(3.0,0.2,"N^{ 0}_{ 3}");
  //   Tl.DrawLatex(4.0,0.2,"N^{ 0}_{ 3}");
  //   Tl.DrawLatex(5.0,0.2,"N^{ 0}_{ 4}");
  //   Tl.DrawLatex(6.0,0.2,"N^{ 0}_{ 4}");
  //   Tl.DrawLatex(7.0,0.2,"N^{ 0}_{ 5,6}");
  //   Tl.DrawLatex(8.0,0.2,"N^{ 0}_{ 5,6}");
  //   Tl.DrawLatex(9.0,0.2,"N^{ 0}_{>=7}");
  //   Tl.DrawLatex(10.0,0.2,"N^{ 0}_{>=7}");
  //   Tl.DrawLatex(11.0,0.2,"N^{ #geq1}_{ #geq2}");
  //   Tl.DrawLatex(12.0,0.2,"N^{ #geq1}_{ #geq2}");
  //   Tl.DrawLatex(13.0,0.2,"N^{ #geq1}_{ 5,6}");
  //   Tl.DrawLatex(14.0,0.2,"N^{ #geq1}_{ 5,6}");
  //   Tl.DrawLatex(15.0,0.2,"N^{ #geq1}_{ #geq7}");
  //   Tl.DrawLatex(16.0,0.2,"N^{ #geq1}_{ #geq7}");

  //   Tl1.DrawLatex(1.0,0.35,"P_{T}^{#gamma}<0.1");//N^{ 0}_{ 2}");
  //   Tl1.DrawLatex(2.0,0.35,"P_{T}^{#gamma}>0.1");//N^{ 0}_{ 2}");
  //   Tl1.DrawLatex(3.0,0.35,"P_{T}^{#gamma}<0.1");//N^{ 0}_{ 3}");
  //   Tl1.DrawLatex(4.0,0.35,"P_{T}^{#gamma}>0.1");//N^{ 0}_{ 3}");
  //   Tl1.DrawLatex(5.0,0.35,"P_{T}^{#gamma}<0.1");//N^{ 0}_{ 4}");
  //   Tl1.DrawLatex(6.0,0.35,"P_{T}^{#gamma}>0.1");//N^{ 0}_{ 4}");
  //   Tl1.DrawLatex(7.0,0.35,"P_{T}^{#gamma}<0.1");//N^{ 0}_{ 5,6}");
  //   Tl1.DrawLatex(8.0,0.35,"P_{T}^{#gamma}>0.1");//N^{ 0}_{ 5,6}");
  //   Tl1.DrawLatex(9.0,0.35,"P_{T}^{#gamma}<0.1");//N^{ 0}_{>=7}");
  //   Tl1.DrawLatex(10.0,0.35,"P_{T}^{#gamma}>0.1");//N^{ 0}_{>=7}");
  //   Tl1.DrawLatex(11.0,0.35,"P_{T}^{#gamma}<0.1");//N^{ #geq1}_{ #geq2}");
  //   Tl1.DrawLatex(12.0,0.35,"P_{T}^{#gamma}>0.1");//N^{ #geq1}_{ #geq2}");
  //   Tl1.DrawLatex(13.0,0.35,"P_{T}^{#gamma}<0.1");//N^{ #geq1}_{ 5,6}");
  //   Tl1.DrawLatex(14.0,0.35,"P_{T}^{#gamma}>0.1");//N^{ #geq1}_{ 5,6}");
  //   Tl1.DrawLatex(15.0,0.35,"P_{T}^{#gamma}<0.1");//N^{ #geq1}_{ #geq7}");
  //   Tl1.DrawLatex(16.0,0.35,"P_{T}^{#gamma}>0.1");//N^{ #geq1}_{ #geq7}");

  // }

  // else if (which_TFbins==3){
    
  //   Tl.SetTextSize(0.035);
  //   Tl.DrawLatex(1.0,0.2,"N^{ 0}_{ 2}");
  //   Tl.DrawLatex(2.0,0.2,"N^{ 0}_{ 2}");
  //   Tl.DrawLatex(3.0,0.2,"N^{ 0}_{ 3}");
  //   Tl.DrawLatex(4.0,0.2,"N^{ 0}_{ 3}");
  //   Tl.DrawLatex(5.0,0.2,"N^{ 0}_{ 4}");
  //   Tl.DrawLatex(6.0,0.2,"N^{ 0}_{ 4}");
  //   Tl.DrawLatex(7.0,0.2,"N^{ 0}_{ 5,6}");
  //   Tl.DrawLatex(8.0,0.2,"N^{ 0}_{ 5,6}");
  //   Tl.DrawLatex(9.0,0.2,"N^{ 0}_{>=7}");
  //   Tl.DrawLatex(10.0,0.2,"N^{ 0}_{>=7}");
  //   Tl.DrawLatex(11.0,0.2,"N^{ #geq1}_{ #geq2}");
  //   Tl.DrawLatex(12.0,0.2,"N^{ #geq1}_{ #geq2}");
  //   Tl.DrawLatex(13.0,0.2,"N^{ #geq1}_{ 5,6}");
  //   Tl.DrawLatex(14.0,0.2,"N^{ #geq1}_{ 5,6}");
  //   Tl.DrawLatex(15.0,0.2,"N^{ #geq1}_{ #geq7}");
  //   Tl.DrawLatex(16.0,0.2,"N^{ #geq1}_{ #geq7}");
  //   Tl1.DrawLatex(1.0,0.35,"MET<0.3");//N^{ 0}_{ 2}");                                                                                                                                                  
  //   Tl1.DrawLatex(2.0,0.35,"MET>0.3");//N^{ 0}_{ 2}");                                                                                                                                                  
  //   Tl1.DrawLatex(3.0,0.35,"MET<0.3");//N^{ 0}_{ 3}");                                                                                                                                                  
  //   Tl1.DrawLatex(4.0,0.35,"MET>0.3");//N^{ 0}_{ 3}");                                                                                                                                                  
  //   Tl1.DrawLatex(5.0,0.35,"MET<0.3");//N^{ 0}_{ 4}");                                                                                                                                                  
  //   Tl1.DrawLatex(6.0,0.35,"MET>0.3");//N^{ 0}_{ 4}");                                                                                                                                                  
  //   Tl1.DrawLatex(7.0,0.35,"MET<0.3");//N^{ 0}_{ 5,6}");                                                                                                                                                
  //   Tl1.DrawLatex(8.0,0.35,"MET>0.3");//N^{ 0}_{ 5,6}");                                                                                                                                                
  //   Tl1.DrawLatex(9.0,0.35,"MET<0.3");//N^{ 0}_{>=7}");                                                                                                                                                 
  //   Tl1.DrawLatex(10.0,0.35,"MET>0.3");//N^{ 0}_{>=7}");                                                                                                                                                
  //   Tl1.DrawLatex(11.0,0.35,"MET<0.3");//N^{ #geq1}_{ #geq2}");                                                                                                                                         
  //   Tl1.DrawLatex(12.0,0.35,"MET>0.3");//N^{ #geq1}_{ #geq2}");                                                                                                                                         
  //   Tl1.DrawLatex(13.0,0.35,"MET<0.3");//N^{ #geq1}_{ 5,6}");                                                                                                                                           
  //   Tl1.DrawLatex(14.0,0.35,"MET>0.3");//N^{ #geq1}_{ 5,6}");                                                                                                                                           
  //   Tl1.DrawLatex(15.0,0.35,"MET<0.3");//N^{ #geq1}_{ #geq7}");                                                                                                                                         
  //   Tl1.DrawLatex(16.0,0.35,"MET>0.3");//N^{ #geq1}_{ #geq7}");           
  

  // }


  
//     TArrow *arrow1 = new TArrow( 0.0,300, 7.0,300,0.01,"<|>");
//     TArrow *arrow2 = new TArrow( 7.0,300,13.0,300,0.01,"<|>");
//     TArrow *arrow3 = new TArrow(13.0,300,19.0,300,0.01,"<|>");
//     TArrow *arrow4 = new TArrow(19.0,300, 25.0,300,0.01,"<|>");
//     TArrow *arrow5 = new TArrow(25.0,300, 31.0,300,0.01,"<|>");
//     TArrow *arrow6 = new TArrow(31.0,300, 38.0,300,0.01,"<|>");

// //     TArrow *arrow1 = new TArrow( 1.0,0.10, 3.0,0.1,0.01,"<|>");
// //     TArrow *arrow2 = new TArrow( 3.0,0.10,5.0,0.1,0.01,"<|>");
// //     TArrow *arrow3 = new TArrow(5.0,0.10,7.0,0.1,0.01,"<|>");
// //     TArrow *arrow4 = new TArrow(7.0,0.10, 9.0,0.1,0.01,"<|>");
// //     TArrow *arrow5 = new TArrow(9.0,0.10, 11.0,0.1,0.01,"<|>");
// //     TArrow *arrow6 = new TArrow(11.0,0.10, 13.0,0.1,0.01,"<|>");
// // TArrow *arrow7 = new TArrow(13.0,0.1, 15.0,0.1,0.01,"<|>");
// //  TArrow *arrow8 = new TArrow(15.0,0.1, 17.0,0.1,0.01,"<|>");

//     arrow1->Draw(); arrow2->Draw(); arrow3->Draw();
//     arrow4->Draw(); arrow5->Draw(); arrow6->Draw();
//     // arrow7->Draw(); arrow8->Draw();
//     TLatex Tl;
//     Tl.SetTextSize(0.055);
    //    Tl.SetTextColor(kWhite);
    // Tl.DrawLatex(3.5,100,"N^{ 0}_{ 2-4}");
    // Tl.DrawLatex(10.5,100,"N^{ 0}_{ 5-6}");
    // Tl.DrawLatex(15.5,100,"N^{ 0}_{ #geq7}");
    // Tl.DrawLatex(21.5,100,"N^{ #geq1}_{ 2-4}");
    // Tl.DrawLatex(26.5,100,"N^{ #geq1}_{ 5-6}");
    // Tl.DrawLatex(33.5,100,"N^{ #geq1}_{ #geq7}");
    // Tl.DrawLatex(1.5,0.1,"N^{ 0}_{ 2}");
    // Tl.DrawLatex(3.5,0.1,"N^{ 0}_{ 3}");
    // Tl.DrawLatex(5.5,0.1,"N^{ 0}_{ 4}");
    // Tl.DrawLatex(7.5,0.1,"N^{ 0}_{ 5,6}");
    // Tl.DrawLatex(9.5,0.1,"N^{ 0}_{>=7}");
    // Tl.DrawLatex(11.5,0.1,"N^{ #geq1}_{ #geq2}");
    // Tl.DrawLatex(13.5,0.1,"N^{ #geq1}_{ 5,6}");
    // Tl.DrawLatex(15.5,0.1,"N^{ #geq1}_{ #geq7}");


  gPad->Modified();
                                                                                       
    hist_ratio->SetLineWidth(line_width[0]+1);
    hist_ratio->SetLineStyle(1);
    hist_ratio->SetLineColor(kBlue);
    hist_ratio->SetTitle(" ");
    hist_ratio->GetXaxis()->SetTitleSize(0.13);
    hist_ratio->GetYaxis()->SetTitle("TF = #frac{N_{SR}}{N_{CR}}");//(0#mu,1#gamma)}{(1#mu,1#gamma)}");
    hist_ratio->GetXaxis()->SetLabelSize(0.1);
    if(which_Lept==2)
      hist_ratio->GetYaxis()->SetRangeUser(0.0,3.5);
    else
      hist_ratio->GetYaxis()->SetRangeUser(0.0,0.03);
    
    if(which_TFbins==1) //default 8 bins                                                                                                    
      hist_ratio->GetXaxis()->SetRangeUser(0,10);//xmin,xrange);                                                                                                 
    else if(which_TFbins==2) // v2 TF bins including photon pT>100 and pT<100
      {  hist_ratio->GetXaxis()->SetRangeUser(0,50);
	hist_ratio->GetYaxis()->SetRangeUser(0.0,0.1);
      }
    else if(which_TFbins==3) // v3 TF bins including MET<300 and MET>300
      {
	hist_ratio->GetXaxis()->SetRangeUser(0,50);
	hist_ratio->GetYaxis()->SetRangeUser(0.0,0.1);
      }
    //    hist_ratio->GetXaxis()->SetLabelSize(0.0450);
    hist_ratio->GetYaxis()->SetTitleSize(0.13);
    hist_ratio->GetYaxis()->SetLabelSize(0.08);
    hist_ratio->GetYaxis()->SetTitleOffset(.4);
    hist_ratio->GetYaxis()->SetNdivisions(505);
    hist_ratio->GetYaxis()->CenterTitle(true);

    //new ones

    hist_ratio->GetXaxis()->SetTitleSize(0.05);
    hist_ratio->GetXaxis()->SetLabelSize(0.12);
    hist_ratio->GetYaxis()->SetTitleSize(0.125);
    hist_ratio->GetYaxis()->SetNdivisions(505);
    
    hist_ratio->GetXaxis()->SetTitleOffset(0.8);
    hist_ratio->GetYaxis()->SetTitleOffset(0.51);
    hist_ratio->GetXaxis()->SetTitleSize(0.14);
    
    hist_ratio->GetYaxis()->SetLabelSize(0.12);

    //    hist_ratio->GetYaxis()->SetLabelSize(x_label_size);
   pad_1->cd();
   //   pad_1->SetGrid();
   if(which_TFbins==1){
   TLine *l =new TLine(0,.01,10,.01);   
   hist_ratio->Draw("");
   l->Draw("sames");
   TLine *l1 =new TLine(0,0.02,10,0.02);
   l1->SetLineStyle(7);
   l1->Draw("sames");
   TLine *l2 =new TLine(0,0.015,10,0.015);
   l2->SetLineStyle(7);

      l2->Draw("sames");
   }

   else{

     TLine *l =new TLine(0,.01,50,.01);
   hist_ratio->Draw("");
   l->Draw("sames");
   TLine *l1 =new TLine(0,0.02,50,.02);
   l1->SetLineStyle(7);
   l1->Draw("sames");
   TLine *l2 =new TLine(0,0.015,50,0.015);
   l2->SetLineStyle(7);

   l2->Draw("sames");
   }
  char* canvas_name = new char[1000];
  //c->Print(canvas_name);
  
  if(save_canvas) {
    sprintf(canvas_name,"%s.png",tag_name);//.png",tag_name);//_wnormalize.png",tag_name);
     canvas_n1->SaveAs(canvas_name);   
     sprintf(canvas_name,"%s.pdf",tag_name);
    canvas_n1->SaveAs(canvas_name);
    
  }
  
}
const int nfiles=100,nBG=6;                                                                                                                                                              
TFile *f[nfiles];
TFile *f1[nfiles];


void SRvsCRStacked_LL_varRatio(string pathname, int which_TFBins)
{
  char* hname = new char[200];
    char* hname1 = new char[200];
char* hname2 = new char[200];
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
  sprintf(path2,"Results/dPhi_METHadJets");
  int n=0;
  char *dataset=new char[200];
  char *year =new char[200];
  //  float energyy[2]={};
  int n_files=42;
  char *string_png = new char[200];
  vector<string>baseline1;
  vector<string> legend_texts; 
 char *TFbins_str= new char[2000];
  if(which_TFBins==1)
    sprintf(TFbins_str,"TFbins_v1_nJetsBjets");
  else if (which_TFBins==2)
    sprintf(TFbins_str,"TFbins_v2_nJetsBjets_PhoPt");
  else if(which_TFBins==3)
    sprintf(TFbins_str,"TFbins_v3_phopt_qmulti");
  
    baseline1={"ProbL1Trig_Elec_CR","ProbL1Trig_Pho_SR"};//,"TauHad_SR","Mu_SR","Elec_SR","FailAcep_ElecSR","FailId_ElecSR","FailIso_ElecSR","Elec_SR","Elec_"};//
    legend_texts ={"(1e,0#gamma) CR","(0e,1#gamma) SR"};//,"#tau-had SR","lost #mu SR","lost e SR","(1l,1#gamma) CR","Failed Iso"};
    
    //sprintf(hname,"");
    sprintf(string_png,"Electron_FR_%s",TFbins_str);
    sprintf(hname,"%s_phoID_loose_09Jan24",string_png);
    cout<<string_png<<"\t"<<TFbins_str<<"\t"<<which_TFBins<<endl;
      f[0] = new TFile("Summer20UL18_TTGJets_PhoIdloose_phopt40_MET200.root");
      f[1] = new TFile("Summer20UL17_TTGJets_PhoIdloose_phopt40_MET200.root");
      f[2] = new TFile("Summer20UL16_TTGJets_PhoIdloose_phopt40_MET200.root");
      f[3]= new TFile("FullRun2_TTGJets_PhoIdloose_phopt40_MET200.root");

      f[4] = new TFile("Summer20UL18_WGJets_PhoIdloose_phopt40_MET200.root");
      f[5] = new TFile("Summer20UL17_WGJets_PhoIdloose_phopt40_MET200.root");
      f[6] = new TFile("Summer20UL16_WGJets_PhoIdloose_phopt40_MET200.root");
      f[7]= new TFile("FullRun2_WGJets_PhoIdloose_phopt40_MET200.root");
      f[8] = new TFile("Summer20UL18_TTJets_PhoIdloose_phopt40_MET200.root");
      f[9] = new TFile("Summer20UL17_TTJets_PhoIdloose_phopt40_MET200.root");
      f[10] = new TFile("Summer20UL16_TTJets_PhoIdloose_phopt40_MET200.root");
      f[11]= new TFile("FullRun2_TTJets_PhoIdloose_phopt40_MET200.root");
      f[12] = new TFile("Summer20UL18_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root");
      f[13] = new TFile("Summer20UL17_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root");
      f[14] = new TFile("Summer20UL16_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root");
      f[15]= new TFile("FullRun2_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root");
      f[16] = new TFile("Summer20UL18_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root");
      f[17] = new TFile("Summer20UL17_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root");
      f[18] = new TFile("Summer20UL16_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root");
      f[19] = new TFile("FullRun2_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root");
      f[20] = new TFile("Summer20UL18_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root");
      f[21] = new TFile("Summer20UL17_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root");
      f[22] = new TFile("Summer20UL16_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root");
      f[23] = new TFile("FullRun2_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root");
      f[24]= new TFile("Summer20UL18_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root");
      f[25]= new TFile("Summer20UL17_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root");
      f[26]= new TFile("Summer20UL16_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root");
      f[27]= new TFile("FullRun2_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root");
      f[28] = new TFile("Summer20UL16APV_TTGJets_PhoIdloose_phopt40_MET200.root");
      f[29] = new TFile("Summer20UL16APV_WGJets_PhoIdloose_phopt40_MET200.root");
      f[30] = new TFile("Summer20UL16APV_TTJets_PhoIdloose_phopt40_MET200.root");
      f[31] = new TFile("Summer20UL16APV_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root");
      f[32] = new TFile("Summer20UL16APV_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root");
      f[33] = new TFile("Summer20UL16APV_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root");
      f[34]= new TFile("Summer20UL16APV_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root");
      f[35] = new TFile("Summer20UL_total2016_TTGJets_PhoIdloose_phopt40_MET200.root");
      f[36] = new TFile("Summer20UL_total2016_WGJets_PhoIdloose_phopt40_MET200.root");
      f[37] = new TFile("Summer20UL_total2016_TTJets_PhoIdloose_phopt40_MET200.root");
      f[38] = new TFile("Summer20UL_total2016_WJetsToLNu_HT_PhoIdloose_phopt40_MET200.root");
      f[39] = new TFile("Summer20UL_total2016_combined_WGJets_WJets_PhoIdloose_phopt40_MET200.root");
      f[40] = new TFile("Summer20UL_total2016_combined_TTGJets_TTJets_PhoIdloose_phopt40_MET200.root");
      f[41]= new TFile("Summer20UL_total2016_WGJets_TTGJets_Allcombined_PhoIdloose_phopt40_MET200.root");
      
  vector<string> filetag=  {"TTGJets_2018","TTGJets_2017","TTGJets_2016postVFP","Run2_TTGJets","WGJets_2018","WGJets_2017","WGJets_2016postVFP","Run2_WGJets","TTJets_2018","TTJets_2017","TTJets_2016postVFP","Run2_TTJets","WJets_2018","WJets_2017","WJets_2016postVFP","Run2_WJets","WGJets+WJets_2018","WGJets+WJets_2017","WGJets+WJets_2016postVFP","Run2_WGJets+WJets","TTGJets+TTJets_2018","TTGJets+TTJets_2017","TTGJets+TTJets_2016postVFP","Run2_TTGJets+TTJets","W+TTBar_2018","W+TTBar_2017","W+TTBar_2016postVFP","W+TTBar_FullRun2","TTGJets_2016preVFP","WGJets_2016preVFP","TTJets_2016preVFP","WJets_2016preVFP","WGJets+WJets_2016preVFP","TTGJets+TTJets_2016preVFP","W+TTBar_2016preVFP","TTGJets_2016","WGJets_2016","TTJets_2016","WJets_2016","WGJets+WJets_2016","TTGJets+TTJets_2016","W+TTBar_2016"};
 vector<float> energyy={59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,19.5,19.5,19.5,19.5,19.5,19.5,19.5,36.0,36.0,36.0,36.0,36.0,36.0,36.0};


  // vector<string> filetag={"W+TTBar_2016preVFP","W+TTBar_2016postVFP","W+TTBar_2017","W+TTBar_2018","W+TTBar_FullRun2","TTGJets+TTJets_2016preVFP","TTGJets+TTJets_2016postVFP","TTGJets+TTJets_2017","TTGJets+TTJets_2018","TTGJets+TTJets_FullRun2","WGJets+WJets_2016preVFP","WGJets+WJets_2016postVFP","WGJets+WJets_2017","WGJets+WJets_2018","WGJets+WJets_FullRun2","TTGJets_2018","TTGJets_2017","TTGJets_2016postVFP","Run2_TTGJets","WGJets_2018","WGJets_2017","WGJets_2016postVFP","Run2_WGJets","TTJets_2018","TTJets_2017","TTJets_2016postVFP","Run2_TTJets","WJets_2018","WJets_2017","WJets_2016postVFP","Run2_WJets","WGJets+WJets_2018","WGJets+WJets_2017","WGJets+WJets_2016postVFP","Run2_WGJets+WJets","TTGJets+TTJets_2018","TTGJets+TTJets_2017","TTGJets+TTJets_2016postVFP","Run2_TTGJets+TTJets","W+TTBar_2018","W+TTBar_2017","W+TTBar_2016postVFP","W+TTBar_FullRun2","TTGJets_2016preVFP","WGJets_2016preVFP","TTJets_2016preVFP","WJets_2016preVFP","WGJets+WJets_2016preVFP","TTGJets+TTJets_2016preVFP","W+TTBar_2016preVFP","TTGJets_2016","WGJets_2016","TTJets_2016","WJets_2016","WGJets+WJets_2016","TTGJets+TTJets_2016","W+TTBar_2016"};
  // vector<float> energyy = {19.5,16.5,41.529,59.74,137.19,19.5,16.5,41.529,59.74,137.19,19.5,16.5,41.529,59.74,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,59.74,41.529,16.5,137.19,19.5,19.5,19.5,19.5,19.5,19.5,19.5,36.0,36.0,36.0,36.0,36.0,36.0,36.0};
//   if(which_Lept==3){
//     energyy = {19.5,16.5,41.529,59.74,137.19,19.5,16.5,41.529,59.74,137.19,19.5,16.5,41.529,59.74,137.19,36.0,36.0,36.0};
//     filetag={"W+TTBar_2016preVFP","W+TTBar_2016postVFP","W+TTBar_2017","W+TTBar_2018","W+TTBar_FullRun2","TTGJets+TTJets_2016preVFP","TTGJets+TTJets_2016postVF\
// P","TTGJets+TTJets_2017","TTGJets+TTJets_2018","TTGJets+TTJets_FullRun2","WGJets+WJets_2016preVFP","WGJets+WJets_2016postVFP","WGJets+WJets_2017","WGJets+WJets_2018","WGJets+WJets_FullRun2","WGJets+WJets_2016","TTGJets+TTJets_2016","W+TTBar_2016"};
//   }
  bool flag=false;
  //  const char* filetag[10]={"TTGJets","pMSSM_MCMC_70_90438","pMSSM_MCMC_106_19786","pMSSM_MCMC_473_54451","WJets","GJets","T5bbbbZG_10","T5bbbbZG_50","T5bbbbZG_200","T5bbbbZG_1500"};
   /* vector<TH1F*> hist_list_Njets; */
   /*    vector<TH1F*> hist_list_Bjets; */
   /*    vector<TH1F*> hist_list_MET; */
   /*    vector<TH1F*> hist_list_PhoPt; */
   /*    //vector<TH1F*> hist_list_Mt;                                                                                                                                    */
   /*    vector<TH1F*> hist_list_ST; */
   /*    vector<TH1F*> hist_list_HT; */

  //  TFile* fout = new TFile("TF_allin1_LLEstimation_electron.root","RECREATE");
  sprintf(hname1,"%s.root",hname); 
  TFile* fout = new TFile(hname1,"RECREATE");
  // sprintf(hname,"EventYields_TF_LL_muon_allProcess_binsV3_phoID_loose_09Jan24.txt");
  // std::ofstream file_;
  // file_.open(hname,ios::out);
  //n_files=1;  
  for(int i_file=0; i_file<n_files;i_file++)
    {      
      vector<TH1F*> hist_list_Njets;
      vector<TH1F*> hist_list_Bjets;
      vector<TH1F*> hist_list_MET;
      vector<TH1F*> hist_list_PhoPt;
      //vector<TH1F*> hist_list_Mt;
      vector<TH1F*> hist_list_ST;
      vector<TH1F*> hist_list_HT;
      vector<TH1F*> hist_BDTresp;
      // sprintf(hist_name,"h_TFbins_ElecLL_%s_%s",baseline1[4].c_str(),filetag[i_file]);
      // TH1F* hNjets_SR = new TH1F(hist_name,hist_name,30,0,30);
      // sprintf(hist_name,"h_TFbins_ElecCR_%s_%s",baseline1[5].c_str(),filetag[i_file]);
      // TH1F* hNjets_CR = new TH1F(hist_name,hist_name,30,0,30);

      // sprintf(hist_name,"h_TFbins_LL_%s_%s",baseline1[5].c_str(),filetag[i_file]);
      // TH1F* h_TFBins = new TH1F(hist_name,hist_name,30,0,30);
      // cout<<hist_name<<endl;
      // cout<<f[i_file]->GetName()<<endl;
      float scale_fac = 0.0;//hNjets_temp->Integral();
      //cout<<scale_fac<<"\t"<<"E_SR"<<endl;
      sprintf(hist_name,"h_TFbins_PhoSR_%s",filetag[i_file].c_str());
      TH1F* hNjets_SR = new TH1F(hist_name,hist_name,30,0,30);
      sprintf(hist_name,"h_TFbins_ElecCR_%s",filetag[i_file].c_str());
      TH1F* hNjets_CR = new TH1F(hist_name,hist_name,30,0,30);

      sprintf(hist_name,"h_TFbins_LL_%s",filetag[i_file].c_str());
      TH1F* h_TFBins = new TH1F(hist_name,hist_name,30,0,30);
      sprintf(hname1,"%s_%s.txt",hname,filetag[i_file].c_str());
      std::ofstream file_i;
      file_i.open(hname1,ios::out);

      for(int i_cut=0; i_cut<baseline1.size();i_cut++)
	{
	  if(which_TFBins==1)
	    sprintf(hist_name,"h_TFbins_ElecLL_v2_%s",baseline1[i_cut].c_str());
	  else if (which_TFBins==2)
	    sprintf(hist_name,"h_TFbins_ElecLL_v7_%s",baseline1[i_cut].c_str());
	  else if(which_TFBins==3)
	    sprintf(hist_name,"h_TFbins_ElecLL_v5_%s",baseline1[i_cut].c_str());
	  else if(which_TFBins==4)
	    sprintf(hist_name,"h_TFbins_ElecLL_v6_%s",baseline1[i_cut].c_str());
	  
 	  cout<<hist_name<<"\t"<<i_cut<<endl;
	  TH1F* h_TFbins = (TH1F*)f[i_file]->Get(hist_name);
	  //cout<<h_TFbins->Integral()<<"\t"<<"Elec"<<"\t"<<h_TFbins->GetNbinsX()<<endl;
	  TH1F* h_TFallin1 = (TH1F*)h_TFbins->Clone();
	  cout<<h_TFbins->Integral()<<"\t"<<" h_TFbins->Integral() "<<"\t"<<string_png<<endl;//<<
	  // if (which_Lept==3){
	  //   if(i_cut>=2)
	  //     hist_list_Njets.push_back(h_TFbins);
	  // }
	  // else
	  
	  hist_list_Njets.push_back(h_TFbins);
	  //if(i_cut>=3)
	  hist_list_Bjets.push_back(h_TFbins);
	}
      // }
      cout<<" hist_list_Njets.size() "<<hist_list_Njets.size()<<"\t hist_list_Bjets.size() "<<hist_list_Bjets.size()<<endl;
      //path to save the png file
      float energy=energyy[i_file];
      int xrange=0.0;
      //      sprintf(hist_name,"h_TFbins_ElecLL_%s",baseline1[4].c_str());
	TH1F* hNjets_ratio =(TH1F*)hist_list_Njets.at(1)->Clone();
	hNjets_ratio->GetXaxis()->SetTitle("Bin no.");
	hNjets_ratio->GetYaxis()->SetTitle(" ");
	hNjets_ratio->SetTitle("");
	TH1F* hNjets_ratio1 =(TH1F*)hist_list_Njets.at(0)->Clone();
	hNjets_ratio1->GetXaxis()->SetTitle("Bin no.");
	hNjets_ratio1->GetYaxis()->SetTitle(" ");
	hNjets_ratio1->SetTitle("");
	hNjets_ratio1->Add(hist_list_Njets.at(1));
	hNjets_ratio->Divide(hist_list_Njets.at(0));
	// hist_list_Njets.at(0)->Divide(hNjets_ratio1);
	// hist_list_Njets.at(1)->Divide(hNjets_ratio1);
	// hist_list_Njets.at(0)->GetXaxis()->SetRangeUser(0,10);
	// hist_list_Njets.at(1)->GetXaxis()->SetRangeUser(0,10);
	float  xbin_cr=0, xbin_sr=0;
	sprintf(full_path,"%s/%s_%s",pathname.c_str(),string_png,filetag[i_file].c_str());
	fout->cd();
	h_TFBins = (TH1F*)hNjets_ratio->Clone();
	sprintf(hist_name,"h_TFbins_LL_%s",filetag[i_file].c_str());
	h_TFBins->SetName(hist_name);
	h_TFBins->Write();
	for (int j =2; j<50;j++)
	  {
	    file_i<<j<<"\t"<<h_TFBins->GetBinContent(j)<<"\t"<<h_TFBins->GetBinError(j)<<"\n";
	    cout<<j<<"\t"<<h_TFBins->GetBinContent(j)<<"\t"<<h_TFBins->GetBinError(j)<<endl;
	  }
	generate_1Dplot(hist_list_Njets,hNjets_ratio,full_path,energy,40,0,leg_head,false,true,false,true,filetag[i_file].c_str(),legend_texts,which_TFBins, 3);//which_Lept);
    }
  fout->Close();

}





