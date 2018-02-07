// Chad Harrington 1/25/2018 - Draw histogram

void setStyle();

//plot acceptance, genjet, recojet, jet0_genjet0

void drawHist( TString fname = "plots.root", TString plotname = "acceptance", TString gencut = "0" ) {
  if (plotname != "acceptance" && plotname != "recojet" && gencut != "0") gencut = "0";

  TFile* file = TFile::Open(fname);
  vector< pair<TString, EColor> > flavors = { {"b", kRed},  {"c", kGreen}, {"udsg", kBlack} };

  setStyle();
  TCanvas* c = new TCanvas("c", "c", 600, 600);

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

  TLegend* leg = new TLegend(.6,.75,.85,.9);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);

  int nbins = 500;
  TH1D* h = new TH1D("h", "h", nbins, 0, nbins);
  //TH2D* h = new TH2D("h","h",nbins,0,nbins,nbins,0,nbins);

  TString xtitle = "p_{T}^{cut} (GeV)";
  if      (plotname == "genjet")       xtitle = "Leading genjet p_{T} (GeV)";
  else if (plotname == "recojet")      xtitle = "Leading recojet p_{T} (GeV)";
  else if (plotname == "jet0_genjet0") { xtitle = "Leading genjet p_{T} (GeV)";  h->GetYaxis()->SetTitle("Leading recojet p_{T} (GeV)"); }

  h->GetXaxis()->SetTitle(xtitle);
  h->GetXaxis()->SetNdivisions(5, 5, 0);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(1.1);

  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.3);

  if ( plotname == "acceptance" ) {
    leg->SetTextSize(0.05);

    h->GetXaxis()->SetTickLength(0.03/t_scale);
    h->GetXaxis()->SetLabelSize(0);
    h->GetYaxis()->SetTitle("Acceptance (p_{T} > p_{T}^{cut})");
    h->GetYaxis()->SetTitleSize(0.05/t_scale);
    h->GetYaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetLabelSize(0.04/t_scale);
  }

  h->Draw();

  //TH2D* h2  = (TH2D*) file->Get(plotname);
  //h->GetYaxis()->SetNdivisions(5, 5, 0);
  //h2->Draw("same");

  float ymax=0;
  map<TString, TH1D*> m_flav;

  for (auto const& fpair : flavors) {
    TString flavor = fpair.first;

    if (plotname == "acceptance") {
      TH1D* h_jet =    (TH1D*) file->Get( flavor + "_jet0_" + gencut );
      TH1D* h_genjet = (TH1D*) file->Get( flavor + "_genjet0" );
      double ngenjet = h_genjet->Integral(1,nbins+1);

      m_flav[flavor] = new TH1D(flavor, flavor, 500, 0, 500);
      for (int i=1; i<=nbins; i++) m_flav[flavor]->SetBinContent(i, h_jet->Integral(i, nbins+1)/ngenjet);
      leg->AddEntry(m_flav[flavor], flavor, "L");
    }
    else {
      TString hname = plotname == "genjet" ? "_genjet0" : "_jet0_" + gencut;
      m_flav[flavor] = (TH1D*) file->Get( flavor + hname );
      m_flav[flavor]->Scale( 1 / m_flav[flavor]->Integral() );
      leg->AddEntry(m_flav[flavor], Form("%s  (Mean %.1f)", flavor.Data(), m_flav[flavor]->GetMean() ), "L");
    }

    float m = m_flav[flavor]->GetMaximum();
    if (m > ymax) ymax = m;

    m_flav[flavor]->SetLineColor(fpair.second);
    m_flav[flavor]->Draw("same");
  }
  h->GetYaxis()->SetRangeUser(0, ymax*1.2);
  leg->Draw();

  TLatex text;
  text.SetNDC();

  text.SetTextColor(kBlack);
  text.SetTextSize(0.05);
  text.SetTextFont(61);
  text.DrawLatex(0.18, 0.96, "CMS");

  text.SetTextSize(0.035);
  text.SetTextFont(42);
  if (gencut != "0") text.DrawLatex( 0.65, 0.65, Form("#bf{genjet p_{T} > %s GeV}", gencut.Data()) );

  if ( plotname == "acceptance" ) {
    bottom->cd();
    TH1D* bh = (TH1D*) h->Clone("bh");
    TLegend* bleg = new TLegend(.6,.45,.85,.75);
    bleg->SetTextSize(0.07);
    bleg->SetBorderSize(0);
    bleg->SetFillColor(0);
    bleg->SetFillStyle(0);

    bh->GetXaxis()->SetTickLength(0.03/b_scale);
    bh->GetXaxis()->SetLabelSize(0.04/b_scale);
    bh->GetXaxis()->SetTitle(xtitle);
    bh->GetXaxis()->SetTitleSize(0.05/b_scale);
    bh->GetXaxis()->SetTitleOffset(0.9);

    bh->GetYaxis()->SetRangeUser(0, 5);
    bh->GetYaxis()->SetNdivisions(5, 3, 0);
    bh->GetYaxis()->SetLabelSize(0.04/b_scale);
    bh->GetYaxis()->SetTitle("Ratio");
    bh->GetYaxis()->SetTitleSize(0.05/b_scale);
    bh->GetYaxis()->SetTitleOffset(0.4);

    bh->Draw();

    TH1D* b_udsg = (TH1D*) m_flav["b"]->Clone("b_udsg");
    TH1D* c_udsg = (TH1D*) m_flav["c"]->Clone("c_udsg");
    TH1D* b_c =    (TH1D*) m_flav["b"]->Clone("b_c");
    b_c->SetLineColor(kBlue);

    b_udsg->Divide(m_flav["udsg"]);  b_udsg->Draw("same");
    c_udsg->Divide(m_flav["udsg"]);  c_udsg->Draw("same");
    b_c->Divide(m_flav["c"]);        b_c->Draw("same");

    bleg->AddEntry(b_udsg, "b / udsg", "L");
    bleg->AddEntry(c_udsg, "c / udsg", "L");
    bleg->AddEntry(b_c, "b / c", "L");
    bleg->Draw();
  }

  c->Print(plotname + (gencut != "0" ? gencut: "") + ".pdf");
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
