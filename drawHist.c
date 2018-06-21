// Chad Harrington 1/25/2018 - Draw histogram

void setStyle();

TH1D* getpT5( TString fname ) {
  TFile* file = TFile::Open(fname);
  if (file == 0) return 0;
  TH1D* h_pt5 = 0;
  TIter nextkey(file->GetListOfKeys());
  TKey* key;
  while ( (key = (TKey*)nextkey()) ) {
    if ( (TString)key->GetTitle() == "pt(5)" ) { h_pt5 = (TH1D*) key->ReadObj(); break; }
  }
  return h_pt5;
}

void divide(double x, double xerr, double y, double yerr) {
  double q = x/y;
  double qerr = q * sqrt( xerr/x*xerr/x + yerr/y*yerr/y );

  cout << q << "\t" << qerr << endl;
}

//plot acceptance, genjet, recojet, jet0_genjet0

//void drawHist( TString fname = "plots_nlo.root", TString plotname = "Vpt" ) {
void drawHist( TString fname = "plots_DY0J.root", TString plotname = "Ljet" ) {
//void drawHist( TString fname = "plots_accpt.root", TString plotname = "convolute" ) {
//void drawHist( TString fname = "unfold.root", TString plotname = "jet_genjet" ) {

  TFile* file = TFile::Open(fname);
  TFile* file_1J = TFile::Open("plots_DY1J.root");
  TFile* file_2J = TFile::Open("plots_DY2J.root");
  vector< pair<TString, EColor> > flavors = { {"b", kRed},  {"c", (EColor)8}, {"udsg", kBlack}, {"incl", kBlue} };
  //vector< pair<TString, EColor> > flavors = { {"", kRed},  {"_low", (EColor)8}, {"udsg", kBlack}, {"_high", kBlue} };

  setStyle();

  gStyle->SetPadRightMargin(0.11);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.11);
  gStyle->SetPaintTextFormat("4.2f");

  if (plotname != "jet_genjet") {
    gStyle->SetPadRightMargin(0.03);
    gStyle->SetPadLeftMargin(0.14);
  }

  TCanvas* c = new TCanvas("c", "c", plotname == "jet_genjet" ? 900 : 600, 600);

  //gPad->SetLogy();
  //gPad->SetGridx();
  //gPad->SetGridy();

  float b_scale = 0.3, t_scale = 1 - b_scale;
  TPad* top = new TPad("top", "top", 0, b_scale, 1, 1);
  TPad* bottom = new TPad("bottom", "bottom", 0, 0, 1, b_scale);

  if ( plotname == "acceptance" ) {
    top->SetTopMargin(0.05);
    top->SetBottomMargin(0.05);
    top->Draw();
    bottom->SetBottomMargin(0.35);
    bottom->SetGridy();
    bottom->Draw();
    top->cd();
  }

//  TLegend* leg = new TLegend(.6,.7,.75,.88);
  TLegend* leg = new TLegend(.6,.2,.75,.35);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  int nbins = 100;
  TH1D* h = new TH1D("h", "h", nbins, 0, 5*nbins);

  vector<double> genbins =  {30, 40, 60, 80, 110, 200};
  vector<double> recobins = {30, 35, 40, 50, 60, 80, 110, 140, 200}; //Jetpt
  if (plotname == "Vpt") recobins = {0, 20, 30, 40, 50, 70, 90, 120, 200}; //Vpt

//  TH2D* h = new TH2D("h","h",genbins.size()-1,&genbins[0],recobins.size()-1,&recobins[0]);
  h->Draw();

  //TString xtitle = "Leading Jet p_{T} [GeV]";
  TString xtitle = "p_{T}^{gen} (GeV)";
  if      (plotname == "genjet")       xtitle = "Genjet p_{T} (GeV)";
  else if (plotname == "Ljet")          xtitle = "Leading Jet p_{T} (GeV)";
  else if (plotname == "Vpt")          xtitle = "V p_{T} (GeV)";
  else if (plotname == "convolute")    xtitle = "p_{T} (GeV)";
  else if (plotname == "jet_genjet") {

    xtitle = "genjet p_{T} (GeV)";
    h->GetYaxis()->SetTitle("recojet p_{T} (GeV)");

    h->GetXaxis()->SetLabelSize(0);
    h->GetYaxis()->SetLabelSize(0);
    h->GetXaxis()->SetTickLength(0);
    h->GetYaxis()->SetTickLength(0);
    h->GetZaxis()->SetLabelSize(0.03);

    //h_response
    //c_jet_genjet_rebin
    //b_jet_genjet_rebin

    TH2D* h2  = (TH2D*) file->Get( "h_response" );

//    h2->Scale( 1 / h2->Integral() );
    h2->Draw("sametextcolz");

    TLine l;  l.SetLineStyle(2);
    TLatex t;  t.SetTextSize(0.04);  t.SetTextFont(42);
    double binmax = recobins[recobins.size()-1];
    double binmin = recobins[0];
    for (int i=0, n=genbins.size(); i<n; i++) {
      t.DrawLatex( genbins[i]-2, binmin-10, to_string((int)genbins[i]).data() );
      if (i==0 || i==n-1) continue;
      l.DrawLine( genbins[i], binmin, genbins[i], binmax );
    }
    for (int i=0, n=recobins.size(); i<n; i++) {
      t.DrawLatex( binmin-10, recobins[i]-2, to_string((int)recobins[i]).data() );
      if (i==0 || i==n-1) continue;
      l.DrawLine( binmin, recobins[i], binmax, recobins[i] );
    }
    t.SetNDC();
    t.SetTextSize(0.05);  t.DrawLatex(0.8, 0.96, "c/b");
  }

  h->GetXaxis()->SetTitle(xtitle);
  //h->GetXaxis()->SetNdivisions(9, 4, 0);
  h->GetXaxis()->SetNdivisions(5, 5, 0);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.1);

  //h->GetYaxis()->SetRangeUser(0.1, 10000);
  h->GetYaxis()->SetTitleSize(0.05);

  if ( plotname == "acceptance" ) {
    leg->SetTextSize(0.05);

    h->GetXaxis()->SetTickLength(0.03/t_scale);
    h->GetXaxis()->SetLabelSize(0);
    h->GetYaxis()->SetTitle("Acceptance");
    h->GetYaxis()->SetRangeUser(0, 1.2);
    h->GetYaxis()->SetTitleSize(0.05/t_scale);
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetLabelSize(0.04/t_scale);
  }

  map<TString, TGraphAsymmErrors*> m_graphs;
  map<TString, TH1D*> m_hists;

  //TGraphAsymmErrors* band = new TGraphAsymmErrors();
  //TGraph* nom = new TGraph();

  for (auto const& fpair : flavors) {
    TString flavor = fpair.first;

    if (plotname == "acceptance" || plotname == "convolute") {

      TH1D* h_allgen = (TH1D*) file->Get( flavor + "_genjet" );
      TH1D* h_gen    = (TH1D*) file->Get( flavor + "_genjet30" );
      cout << flavor << " " << h_gen->GetEntries() << " " << h_allgen->GetEntries() << endl;

      h_gen->Sumw2();
      h_allgen->Sumw2();

      m_graphs[flavor] = new TGraphAsymmErrors();
      m_graphs[flavor]->BayesDivide(h_gen, h_allgen);

      double *x = m_graphs[flavor]->GetX(), *y = m_graphs[flavor]->GetY(), *eyhigh = m_graphs[flavor]->GetEYhigh(), *eylow = m_graphs[flavor]->GetEYlow();

      m_hists[flavor] = new TH1D(flavor, flavor, nbins, 0, 5*nbins);
      for (int pt=0; pt<nbins; pt++) {
        int bin = m_hists[flavor]->FindBin( x[pt] );
        if (bin <= 1) continue;

        m_hists[flavor]->SetBinContent( bin, y[pt] );
        m_hists[flavor]->SetBinError( bin, std::max(eyhigh[pt], eylow[pt]) );
      }
  //    double error=0;
  //    cout << flavor << "\t" << m_hists[flavor]->IntegralAndError(1,nbins,error) << "\t" << error << endl;
      if (flavor == "incl" || plotname == "convolute") continue;
      m_hists[flavor]->SetLineColor(fpair.second);
      m_hists[flavor]->SetMarkerColor(fpair.second);
      leg->AddEntry(m_graphs[flavor], flavor, "LE");

      m_graphs[flavor]->SetLineColor(fpair.second);
      m_graphs[flavor]->Draw("same");
    }
    else {
      TH1D* h2 = (TH1D*) file->Get( flavor + "_" + plotname );

      TH1D* h1J = (TH1D*) file_1J->Get( flavor + "_" + plotname );
      h2->Add(h1J);
      TH1D* h2J = (TH1D*) file_2J->Get( flavor + "_" + plotname );
      h2->Add(h2J);

      h2 = (TH1D*) h2->Rebin(recobins.size()-1, flavor+"_rebin", &recobins[0]);
      m_hists[flavor] = h2;
    }

/*    else {
      if (flavor == "udsg") continue;
      TFile* f = TFile::Open( "mcfm_plots/zbjets" + flavor + ".root" );
      TH1D* hist = (TH1D*) f->Get( "id5" );

      int start = 7;
      for (int i=start,n=hist->GetNbinsX(); i<=n; i++) {

        if      (flavor == "")      nom->SetPoint( i-start, hist->GetBinCenter(i), hist->GetBinContent(i) );
        else if (flavor == "_low")  band->SetPointEYlow( i-start, -1 * (hist->GetBinContent(i)+hist->GetBinError(i)) );
        else if (flavor == "_high") {
          band->SetPoint( i-start, hist->GetBinCenter(i), 0 );
          band->SetPointEYhigh( i-start, hist->GetBinContent(i)+hist->GetBinError(i) );
        }
      }
    }*/
  }
/*
  double *x = nom->GetX(), *center = nom->GetY(), *high = band->GetEYhigh(), *low = band->GetEYlow();
  for (int i=0; i<nom->GetN(); i++) {

    double max = std::max( std::max(high[i], center[i]), -1*low[i]);
    double min = std::min( std::min(high[i], center[i]), -1*low[i]);
    double mid = 0;

    if      (center[i] != max && center[i] != min) mid = center[i];
    else if (high[i] != max && high[i] != min)     mid = high[i];
    else                                           mid = -1*low[i];

    nom->SetPoint( i, x[i], mid );
    band->SetPointEYhigh( i, max );
    band->SetPointEYlow( i, -1*min );
  }

  nom->SetLineStyle(7);
  nom->SetLineWidth(3);
  band->SetFillColor(kGreen+1);
  band->Draw("sameE3");
  nom->Draw("sameL");
*/
  if (plotname == "convolute") {
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetRangeUser(0, 200);
    h->GetYaxis()->SetTitle("d#sigma/dp_{T} [fb/5 GeV]");
    h->GetYaxis()->SetTitleOffset(1.2);

    map< TString, TH1D* > m_con, m_30;
    map< TString, TString > flavs = { {"b", "261"}, {"c", "262"}, {"incl", "41"} };
    for (auto const& fpair : flavs) {
      TH1D* h_pt5 = getpT5("mcfm/" + fpair.second + "_pt15.root");
      if      (fpair.first == "b") {
        TH1D* h2 = getpT5("mcfm/263_pt15.root");
        h_pt5->Add(h2);
        TH1D* h3 = getpT5("mcfm/266_pt15.root");
        h_pt5->Add(h3);
      }
      else if (fpair.first == "c") {
        TH1D* h2 = getpT5("mcfm/264_pt15.root");
        h_pt5->Add(h2);
        TH1D* h3 = getpT5("mcfm/267_pt15.root");
        h_pt5->Add(h3);
      }

      TH1D* h_con = (TH1D*) h_pt5->Clone("h_con");
      TH1D* h_acpt = (TH1D*) h_pt5->Clone("h_acpt");
      for (int i=1; i<=40; i++) { h_acpt->SetBinContent( i, m_hists[fpair.first]->GetBinContent(i) );  h_acpt->SetBinError( i, m_hists[fpair.first]->GetBinError(i) ); }
      h_con->Multiply(h_acpt);

      h_con = (TH1D*) h_con->Rebin(recobins.size()-1, "rebin_con", &recobins[0]);
      m_con[fpair.first] = h_con;

      double error=0;
      cout << fpair.first << "\t" << h_con->IntegralAndError(1,40,error,"width") << "\t" << error << endl;

      h->Draw();
      h->GetYaxis()->SetRangeUser(0, h_pt5->GetMaximum()*1.2);

      h_pt5->Draw("same");
      h_con->SetMarkerColor(kRed);
      h_con->SetLineColor(kRed);
      h_con->Draw("same");
      leg->SetHeader(fpair.first + " jets");
      leg->AddEntry(h_pt5, "mcfm", "LE");
      leg->AddEntry(h_con, "convoluted", "LE");
      leg->Draw();

      TLatex text;
      text.SetNDC();

      text.SetTextColor(kBlack);
      text.SetTextSize(0.05);
      text.SetTextFont(61);
      text.DrawLatex(0.18, 0.96, "CMS");
      c->Print(fpair.first + "_convolute.pdf");
      c->Clear("D");  leg->Clear();

      h_pt5 = getpT5("mcfm/" + fpair.second + "_pt30.root");
      if      (fpair.first == "b") {
        TH1D* h2 = getpT5("mcfm/263_pt30.root");
        h_pt5->Add(h2);
        TH1D* h3 = getpT5("mcfm/266_pt30.root");
        h_pt5->Add(h3);
      }
      else if (fpair.first == "c") {
        TH1D* h2 = getpT5("mcfm/264_pt30.root");
        h_pt5->Add(h2);
        TH1D* h3 = getpT5("mcfm/267_pt30.root");
        h_pt5->Add(h3);
      }

      h_pt5 = (TH1D*) h_pt5->Rebin(recobins.size()-1, "rebin_30", &recobins[0]);
      m_30[fpair.first] = h_pt5;

      error=0;
      cout << fpair.first << "\t" << h_pt5->IntegralAndError(1,40,error,"width") << "\t" << error << endl;
    }

    h->GetYaxis()->SetTitle("Ratio");
    vector< pair<TString, TString> > ratios = { {"c", "incl"}, {"b", "incl"}, {"c", "b"} };
    for (auto const& rpair : ratios) {

      cout << rpair.first + "/" + rpair.second << endl;
      int bin200 = m_con[rpair.first]->FindBin(200);
      double x=0, xerr=0, y=0, yerr=0;
      x = m_con[rpair.first]->IntegralAndError(1,bin200,xerr);
      y = m_con[rpair.second]->IntegralAndError(1,bin200,yerr);
      divide(x, xerr, y, yerr);

      TH1D* r_con = (TH1D*) m_con[rpair.first]->Clone(rpair.first+rpair.second);
      r_con->Divide( m_con[rpair.second] );
      TH1D* r_30 = (TH1D*) m_30[rpair.first]->Clone(rpair.first+rpair.second);
      r_30->Divide( m_30[rpair.second] );

      h->Draw();
      h->GetYaxis()->SetRangeUser(0, r_con->GetMaximum()*2);

      r_con->Draw("same");  r_con->SetMarkerColor(kRed);  r_con->SetLineColor(kRed);
      r_30->Draw("same");  r_30->SetMarkerColor(kBlack);  r_30->SetLineColor(kBlack);
      leg->SetHeader(rpair.first + "/" + rpair.second);
      leg->AddEntry(r_30, "mcfm", "LE");
      leg->AddEntry(r_con, "convoluted", "LE");
      leg->Draw();

      TLatex text;
      text.SetNDC();

      text.SetTextColor(kBlack);
      text.SetTextSize(0.05);
      text.SetTextFont(61);
      text.DrawLatex(0.18, 0.96, "CMS");
      c->Print(rpair.first + "_" + rpair.second + "_mcfm.pdf");
      c->Clear("D");  leg->Clear();
    }
    return;
  }
  else if ( plotname != "acceptance" && plotname != "jet_genjet") {
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->GetXaxis()->SetRangeUser(0, 200);

    vector< pair<TString, pair<TString, float> > > ratios = { {"c", {"incl", 0.2}}, {"b", {"incl", 0.1}}, {"c", {"b", 3}} };
//    vector< pair<TString, pair<TString, float> > > ratios = { {"c", {"udsg", 0.2}}, {"b", {"udsg", 0.1}}, {"c", {"b", 3}} };
    for (auto const& rpair : ratios) {
      TString num = rpair.first, denom = rpair.second.first;
      TString ytitle = num + "/" + denom;
      cout << ytitle << endl;

      TH1D* r = (TH1D*) m_hists[num]->Clone(num+denom);
      r->Divide( m_hists[denom] );

      TFile* file2 = TFile::Open( "duong.root" );
      TH1D* h_measured = (TH1D*) file2->Get("r" + num + (denom=="incl"?"j":denom) + "_jet_pt_Muon");
//      TH1D* h_measured = (TH1D*) file2->Get("r" + num + (denom=="udsg"?"j":denom) + "_jet_pt_Muon");
//      TH1D* h_measured = (TH1D*) file2->Get("r" + num + (denom=="incl"?"j":denom) + "_V_pt_Muon");
      TH1D* h_meas2 = new TH1D("h_meas2","h_meas2",recobins.size()-1,&recobins[0]);
      for (int i=1; i<recobins.size(); i++) { h_meas2->SetBinContent( i, h_measured->GetBinContent(i) );  h_meas2->SetBinError( i, h_measured->GetBinError(i) ); }

      TF1* f_const = new TF1("f_const", "[0]");
      h_meas2->Fit(f_const, "N");

      int bin200 = m_hists[num]->FindBin(200);
      double x=0, xerr=0, y=0, yerr=0;
      x = m_hists[num]->IntegralAndError(1,bin200,xerr);
      y = m_hists[denom]->IntegralAndError(1,bin200,yerr);
      divide(x, xerr, y, yerr);

      h->Draw();
      h->GetYaxis()->SetRangeUser(0, rpair.second.second);
      h->GetYaxis()->SetTitle(ytitle);
      r->Draw("same");
      h_meas2->SetMarkerColor(kRed);
      h_meas2->SetLineColor(kRed);
      h_meas2->Draw("same");

      leg->AddEntry(h_meas2, "Measured", "LE");
      leg->AddEntry(r, "amc@NLO", "LE");
      leg->Draw();

      TLatex text;
      text.SetNDC();

      text.SetTextColor(kBlack);
      text.SetTextSize(0.05);
      text.SetTextFont(61);
      text.DrawLatex(0.18, 0.96, "CMS");
      c->Print(num + "_" + denom + "_" + plotname + ".pdf");
      c->Clear("D");  leg->Clear();
    }
    return;
  }
  leg->Draw();

  TLatex text;
  text.SetNDC();

  text.SetTextColor(kBlack);
  text.SetTextSize(0.05);
  text.SetTextFont(61);
  text.DrawLatex(0.18, 0.96, "CMS");

  text.SetTextSize(0.035);
  text.SetTextFont(42);

  //text.DrawLatex( 0.65, 0.85, "Zb scale uncertainty" );

  if ( plotname == "acceptance" ) {
    bottom->cd();
    TH1D* bh = (TH1D*) h->Clone("bh");
    TLegend* bleg = new TLegend(.3,.55,.45,.85);
    bleg->SetTextSize(0.07);
    bleg->SetBorderSize(0);
    bleg->SetFillColor(0);
    bleg->SetFillStyle(0);

    bh->GetXaxis()->SetTickLength(0.03/b_scale);
    bh->GetXaxis()->SetLabelSize(0.04/b_scale);
    bh->GetXaxis()->SetTitle(xtitle);
    bh->GetXaxis()->SetTitleSize(0.05/b_scale);
    bh->GetXaxis()->SetTitleOffset(0.9);

    bh->GetYaxis()->SetRangeUser(0.9, 1.1);
    bh->GetYaxis()->SetNdivisions(5, 3, 0);
    bh->GetYaxis()->SetLabelSize(0.04/b_scale);
    bh->GetYaxis()->SetTitle("Ratio");
    bh->GetYaxis()->SetTitleSize(0.05/b_scale);
    bh->GetYaxis()->SetTitleOffset(0.4);

    bh->Draw();

    TH1D* b_c = (TH1D*) m_hists["b"]->Clone("b_c");
    b_c->Divide( m_hists["c"] );
    b_c->SetLineColor(kBlue);    b_c->Draw("same");

    m_hists["b"]->Divide( m_hists["incl"] );  m_hists["b"]->Draw("same");
    m_hists["c"]->Divide( m_hists["incl"] );  m_hists["c"]->Draw("same");

    bleg->AddEntry(m_hists["b"], "b / incl", "LE");
    bleg->AddEntry(m_hists["c"], "c / incl", "LE");
    bleg->AddEntry(b_c, "b / c", "LE");
    bleg->Draw();
/*
    TF1* f_const = new TF1("f_const", "[0]");
    b_c->Fit(f_const);
    m_hists["b"]->Fit(f_const);
    m_hists["c"]->Fit(f_const);
*/
  }

  c->Print(plotname + ".pdf");
}

void setStyle() {

//Style//

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(0.9);

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  //tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

  tdrStyle->SetPaperSize(20.,20.);

  tdrStyle->SetHatchesLineWidth(5);
  tdrStyle->SetHatchesSpacing(0.05);

  tdrStyle->SetOptStat(0);

  tdrStyle->cd();

//End Style//
}
