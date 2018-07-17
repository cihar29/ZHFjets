vector<double> genbins =  {30, 40, 60, 110, 200};
vector<double> recobins = {30, 35, 40, 50, 60, 80, 110, 140, 200};

TUnfoldBinning* generatorBinning = new TUnfoldBinning("generator");
TUnfoldBinning* detectorBinning = new TUnfoldBinning("detector");

vector<TString> labels = {"30_35", "35_40", "40_50", "50_60", "60_80", "80_110", "110_140", "140_200"};
map<TString, double> m_label = { {"30_35", 32.5}, {"35_40", 37.5}, {"40_50", 45}, {"50_60", 55},
                                 {"60_80", 70}, {"80_110", 95}, {"110_140", 125}, {"140_200", 170} };

void fillData( TH1D*& h_data_tunfold, TH1D*& h_data, bool mc ) {

  for (int i=1, n=h_data->GetNbinsX(); i<=n; i++) {

    double pt = mc ? h_data->GetBinCenter(i) : m_label[h_data->GetXaxis()->GetBinLabel(i)];
    int bin = h_data_tunfold->FindBin( pt );

    h_data_tunfold->SetBinContent( bin, h_data_tunfold->GetBinContent(bin)+h_data->GetBinContent(i) );
    h_data_tunfold->SetBinError( bin, sqrt( h_data_tunfold->GetBinError(bin)*h_data_tunfold->GetBinError(bin)+h_data->GetBinError(i)*h_data->GetBinError(i) ) );
  }
}

void fillResponse( TH2D*& h_response_tunfold, TH2D*& h_response, TH1D* h_fake, TH1D* h_miss ) {

  for (int i=1, nx=h_response->GetXaxis()->GetNbins(); i<=nx; i++) {
    for (int j=1, ny=h_response->GetYaxis()->GetNbins(); j<=ny; j++) {

      int genBin = h_response_tunfold->GetXaxis()->FindBin( h_response->GetXaxis()->GetBinCenter(i) );
      int recBin = h_response_tunfold->GetYaxis()->FindBin( h_response->GetYaxis()->GetBinCenter(j) );

      h_response_tunfold->SetBinContent( genBin, recBin, h_response_tunfold->GetBinContent(genBin,recBin)+h_response->GetBinContent(i,j) );
      h_response_tunfold->SetBinError( genBin, recBin,
        sqrt( h_response_tunfold->GetBinError(genBin,recBin)*h_response_tunfold->GetBinError(genBin,recBin)+h_response->GetBinError(i,j)*h_response->GetBinError(i,j) )
                                     );
    }
  }

  // clear reco underflow for misses
  for (int j=0, ny=h_response->GetYaxis()->GetNbins()+1; j<=ny; j++) {
    h_response_tunfold->SetBinContent( j, 0, 0. );
    h_response_tunfold->SetBinError( j, 0, 0. );
  }

  // misses
  for (int i=1, n=h_miss->GetNbinsX(); i<=n; i++) {

    int genBin = h_response_tunfold->GetXaxis()->FindBin( h_miss->GetBinCenter(i) );
    int recBin = 0;

    h_response_tunfold->SetBinContent( genBin, recBin, h_response_tunfold->GetBinContent(genBin,recBin)+h_miss->GetBinContent(i) );
    h_response_tunfold->SetBinError( genBin, recBin,
      sqrt( h_response_tunfold->GetBinError(genBin,recBin)*h_response_tunfold->GetBinError(genBin,recBin)+h_miss->GetBinError(i)*h_miss->GetBinError(i) ) );
  }
/*
  // fakes
  for (int i=1, n=h_fake->GetNbinsX(); i<=n; i++) {

    int genBin = bkgBinning->GetGlobalBinNumber( h_fake->GetBinCenter(i) );
    int recBin = detectorDistribution->GetGlobalBinNumber( h_fake->GetBinCenter(i) );

    h_response_tunfold->SetBinContent( genBin, recBin, h_response_tunfold->GetBinContent(genBin,recBin)+h_fake->GetBinContent(i) );
    h_response_tunfold->SetBinError( genBin, recBin,
      sqrt( h_response_tunfold->GetBinError(genBin,recBin)*h_response_tunfold->GetBinError(genBin,recBin)+h_fake->GetBinError(i)*h_fake->GetBinError(i) ) );
  }
*/
}

TH1* getUnfold( TH2D*& h_response, TH1D*& h_data ) {

  TUnfoldDensity unfold( h_response, TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeCurvature, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeBinWidth,
                         generatorBinning, detectorBinning, 0, "*[B]" );

//  TH2D* inputEmatrix = detectorBinning->CreateErrorMatrixHistogram("input_covar", true);
//  for (int i=1, n=inputEmatrix->GetNbinsX(); i<=n; i++) {
//    double e = h_data->GetBinError(i);

//    inputEmatrix->SetBinContent(i, i, e*e);
//    if (e<=0.) inputEmatrix->SetBinContent(i,i+1,1.0);
//  }
//  unfold.SetInput( h_data, 0.0, 0.0, inputEmatrix );

  unfold.SetInput( h_data );
/*
  TH2 *histL = unfold.GetL("L");
  for (int j=1, ny=histL->GetNbinsY(); j<=ny; j++) {
    cout<<"L["<<unfold.GetLBinning()->GetBinName(j)<<"]";
    for (int i=1, nx=histL->GetNbinsX(); i<=nx; i++) {
      double c = histL->GetBinContent(i,j);
      if(c!=0.0) cout<<" ["<<i<<"]="<<c;
    }
    cout<<"\n";
  }
*/
  int nScan=100;
  TSpline* rhoLogTau=0, *logTauX=0, *logTauY=0, *logTauCurvature=0;
  TGraph* lCurve=0;

  int iBest = unfold.ScanTau( nScan, 0, 0, &rhoLogTau, TUnfoldDensity::kEScanTauRhoAvg, "generator", 0, &lCurve, &logTauX, &logTauY );
//  int iBest = unfold.ScanLcurve( nScan, 0, 0, &lCurve );

  cout << "chi2\t" << unfold.GetChi2A() << "\t+\t" << unfold.GetChi2L() << "\t/\t" << unfold.GetNdf() << endl;
/*
  TH2D* cov_matrix = new TH2D("cov_matrix","cov_matrix",genbins.size()-1,&genbins[0],genbins.size()-1,&genbins[0]);
  unfold.GetEmatrix(cov_matrix);
  //unfold.GetEmatrixSysUncorr(cov_matrix, 0, false);

  Double_t t[1],rho[1],x[1],y[1];
  rhoLogTau->GetKnot(iBest,t[0],rho[0]);
  lCurve->GetPoint(iBest,x[0],y[0]);
  TGraph *bestRhoLogTau=new TGraph(1,t,rho);
  TGraph *bestLCurve=new TGraph(1,x,y);
  Double_t *tAll=new Double_t[nScan],*rhoAll=new Double_t[nScan];
  for(Int_t i=0;i<nScan;i++) {
     rhoLogTau->GetKnot(i,tAll[i],rhoAll[i]);
  }
  TGraph *knots=new TGraph(nScan,tAll,rhoAll);

  rhoLogTau->Draw();

  TCanvas* c1 = new TCanvas("c1");
  logTauX->Draw();

  TCanvas* c2 = new TCanvas("c2");
  logTauY->Draw();

  TCanvas* c3 = new TCanvas("c3");
  lCurve->Draw();

  TCanvas* c4 = new TCanvas("c4", "bestRhoLogTau");
  bestRhoLogTau->SetMarkerStyle(20);
  bestRhoLogTau->Draw();

  TCanvas* c5 = new TCanvas("c5", "bestLCurve");
  bestLCurve->SetMarkerStyle(20);
  bestLCurve->Draw();

  TCanvas* c6 = new TCanvas("c6", "knots");
  knots->SetMarkerStyle(20);
  knots->Draw();

  TCanvas* c7 = new TCanvas("c7");
  cov_matrix->Draw("textcolz");
*/
  TString name = h_response->GetName();
  return unfold.GetOutput( TString("h_unfolded_signal_")+h_data->GetName() );
}

void unfold( TString ratio = "rcb", bool mc = true ) {
  TString num, denom;
  if      (ratio == "rcj") { num = "c"; denom = "incl"; }
  else if (ratio == "rbj") { num = "b"; denom = "incl"; }
  else                     { num = "c"; denom = "b"; }

  TFile* file_rsp = TFile::Open( mc ? "plots_leading_rsp70.root" : "plots_leading_total.root");
  TFile* file_data = TFile::Open( mc ? "plots_leading_data30.root" : "duong.root" );

  TH2D* h_response_num   = (TH2D*) file_rsp->Get(num + "_Ljet_Lgenjet_mtchd");
  TH2D* h_response_denom = (TH2D*) file_rsp->Get(denom + "_Ljet_Lgenjet_mtchd");

  TH1D* h_fake_num   = (TH1D*) file_rsp->Get(num + "_Ljet_unmtchd");
  TH1D* h_fake_denom = (TH1D*) file_rsp->Get(denom + "_Ljet_unmtchd");

  TH1D* h_miss_num   = (TH1D*) file_rsp->Get(num + "_Lgenjet_unmtchd");
  TH1D* h_miss_denom = (TH1D*) file_rsp->Get(denom + "_Lgenjet_unmtchd");

  TH1D* h_data_num=0, *h_data_denom=0;
  if (mc) {
    // validation samples
    h_data_num   = (TH1D*) file_data->Get(num + "_Ljet");
    h_data_denom = (TH1D*) file_data->Get(denom + "_Ljet");
  }
  else {
    // DUONG'S HISTS //
    h_data_num   = (TH1D*) file_data->Get("jet_pt_N" + num + "_Muon");
    h_data_denom = (TH1D*) file_data->Get("jet_pt_N" + (denom=="incl"?"inc":denom) + "_Muon");

    TH1D* h_data_num_stat   = (TH1D*) file_data->Get("jet_pt_N" + num + "_stat_Muon");
    TH1D* h_data_denom_stat = (TH1D*) file_data->Get("jet_pt_N" + denom + "_stat_Muon");

    for (int i=1, n=h_data_num->GetNbinsX(); i<=n; i++) {
      h_data_num->SetBinError( i, h_data_num_stat->GetBinContent(i) );
      if (denom != "incl")
        h_data_denom->SetBinError( i, h_data_denom_stat->GetBinContent(i) );
    }
  }

  // UNFOLD //

  // binning

  generatorBinning->AddAxis( "genpt", genbins.size()-1, &genbins[0], true, true );
  detectorBinning->AddAxis( "recopt", recobins.size()-1, &recobins[0], true, true );

  TH1D* h_data_num_tunfold = new TH1D("h_data_num_tunfold", "h_data_num_tunfold", recobins.size()-1, &recobins[0]);
  TH1D* h_data_denom_tunfold = new TH1D("h_data_denom_tunfold", "h_data_denom_tunfold", recobins.size()-1, &recobins[0]);
  TH2D* h_response_num_tunfold = new TH2D("h_response_num_tunfold", "h_response_num_tunfold", genbins.size()-1, &genbins[0], recobins.size()-1, &recobins[0]);
  TH2D* h_response_denom_tunfold = new TH2D("h_response_denom_tunfold", "h_response_denom_tunfold", genbins.size()-1, &genbins[0], recobins.size()-1, &recobins[0]);

  // fill hists

  fillData(h_data_num_tunfold, h_data_num, mc);
  fillData(h_data_denom_tunfold, h_data_denom, mc);
  fillResponse(h_response_num_tunfold, h_response_num, h_fake_num, h_miss_num);
  fillResponse(h_response_denom_tunfold, h_response_denom, h_fake_denom, h_miss_denom);

  // unfold

  TH1* h_unfolded_num   = getUnfold( h_response_num_tunfold, h_data_num_tunfold );
  TH1* h_unfolded_denom = getUnfold( h_response_denom_tunfold, h_data_denom_tunfold );

  // Evaluate systematic uncertainty (tested systDN and same result)
  TH1* h_unfolded_num_systUP = 0, *h_unfolded_denom_systUP = 0;
  if (!mc) {
    TH1D* h_data_num_systUP = (TH1D*) file_data->Get("jet_pt_N" + num + "_syst_Muon");
    h_data_num_systUP->Add(h_data_num);
    TH1D* h_data_num_systUP_tunfold = new TH1D("h_data_num_systUP_tunfold", "h_data_num_systUP_tunfold", recobins.size()-1, &recobins[0]);
    fillData(h_data_num_systUP_tunfold, h_data_num_systUP, mc);
    h_unfolded_num_systUP = getUnfold( h_response_num_tunfold, h_data_num_systUP_tunfold );

    TH1D* h_data_denom_systUP=0;
    if (denom != "incl") h_data_denom_systUP = (TH1D*) file_data->Get("jet_pt_N" + denom + "_syst_Muon");
    else {
      // make incl syst = incl stat
      h_data_denom_systUP = new TH1D("h_data_denom_systUP", "h_data_denom_systUP", 8, 0, 8);
      for (int i=1, nreco=recobins.size(); i<nreco; i++) {
        h_data_denom_systUP->SetBinContent( i, h_data_denom->GetBinError(i) );
        h_data_denom_systUP->SetBinError( i, sqrt( h_data_denom->GetBinError(i) ) );
        h_data_denom_systUP->GetXaxis()->SetBinLabel(i, labels[i-1]);
      }
    }

    h_data_denom_systUP->Add(h_data_denom);
    TH1D* h_data_denom_systUP_tunfold = new TH1D("h_data_denom_systUP_tunfold", "h_data_denom_systUP_tunfold", recobins.size()-1, &recobins[0]);
    fillData(h_data_denom_systUP_tunfold, h_data_denom_systUP, mc);
    h_unfolded_denom_systUP = getUnfold( h_response_denom_tunfold, h_data_denom_systUP_tunfold );
  }

  // OUTPUT //

  TH1D* h_gen_num_rebin1=0, *h_gen_denom_rebin1=0, *h_gen_num_rebin2=0, *h_gen_denom_rebin2=0,
       *h_gen_num_rebin3=0, *h_gen_denom_rebin3=0, *h_unfolded_ratio=0;
  if (mc) {
    cout << "Reco" << endl;
    for (int i=1, nreco=recobins.size(); i<nreco; i++) cout << h_data_num_tunfold->GetBinContent(i)/h_data_denom_tunfold->GetBinContent(i) << endl;
    cout << endl;

    // VALIDATION //
    TH1D* h_gen_num1   = (TH1D*) file_data->Get(num + "_Lgenjet_mtchd");
    TH1D* h_gen_denom1 = (TH1D*) file_data->Get(denom + "_Lgenjet_mtchd");
    h_gen_num_rebin1   = (TH1D*) h_gen_num1->Rebin(genbins.size()-1, "h_gen_num1", &genbins[0]);
    h_gen_denom_rebin1 = (TH1D*) h_gen_denom1->Rebin(genbins.size()-1, "h_gen_denom1", &genbins[0]);

//    cout << "Gen (matched to leading reco)" << endl;
//    for (int i=1, ngen=genbins.size(); i<ngen; i++) cout << h_gen_num_rebin1->GetBinContent(i)/h_gen_denom_rebin1->GetBinContent(i) << endl;
//    cout << endl;

    TH1D* h_gen_num2   = (TH1D*) file_data->Get(num + "_Lgenjet_unmtchd");
    TH1D* h_gen_denom2 = (TH1D*) file_data->Get(denom + "_Lgenjet_unmtchd");
    h_gen_num2->Add(h_gen_num1);
    h_gen_denom2->Add(h_gen_denom1);
    h_gen_num_rebin2   = (TH1D*) h_gen_num2->Rebin(genbins.size()-1, "h_gen_num2", &genbins[0]);
    h_gen_denom_rebin2 = (TH1D*) h_gen_denom2->Rebin(genbins.size()-1, "h_gen_denom2", &genbins[0]);

    vector<double> v_gen, v_unfolded;
    cout << "Gen (matched to leading reco + unmatched leading)" << endl;
    for (int i=1, ngen=genbins.size(); i<ngen; i++) {
      v_gen.push_back( h_gen_num_rebin2->GetBinContent(i)/h_gen_denom_rebin2->GetBinContent(i) );
      cout << v_gen.back() << endl;
    }
    cout << endl;

    TH1D* h_gen_num3   = (TH1D*) file_data->Get(num + "_Lgenjet");
    TH1D* h_gen_denom3 = (TH1D*) file_data->Get(denom + "_Lgenjet");
    h_gen_num_rebin3   = (TH1D*) h_gen_num3->Rebin(genbins.size()-1, "h_gen_num3", &genbins[0]);
    h_gen_denom_rebin3 = (TH1D*) h_gen_denom3->Rebin(genbins.size()-1, "h_gen_denom3", &genbins[0]);

//    cout << "Gen (leading gen)" << endl;
//    for (int i=1, ngen=genbins.size(); i<ngen; i++) cout << h_gen_num_rebin3->GetBinContent(i)/h_gen_denom_rebin3->GetBinContent(i) << endl;
//    cout << endl;

    cout << "Unfolded" << endl;
    for (int i=1, ngen=genbins.size(); i<ngen; i++) {
      v_unfolded.push_back( h_unfolded_num->GetBinContent(i)/h_unfolded_denom->GetBinContent(i) );
      cout << v_unfolded.back() << endl;
    }
    cout << endl;

    cout << "Percent Difference" << endl;
    for (int i=1, ngen=genbins.size(); i<ngen; i++) cout << (v_unfolded[i-1] - v_gen[i-1]) / v_gen[i-1] * 100 << endl;
    cout << endl;
  }
  else {
    map<TString, double> m_eff = { {"c", 0.094}, {"b", 0.54}, {"incl", 1} };
    TH1D* h_sf_num   = (TH1D*) file_data->Get("jet_pt_SF" + num + "_Muon");
    TH1D* h_sf_denom = 0;
    if (denom == "incl") {
      h_sf_denom = new TH1D("incl_sf", "incl_sf", recobins.size()-1,&recobins[0]);
      for (int i=1, nreco=recobins.size(); i<nreco; i++) h_sf_denom->SetBinContent(i, 1);
    }
    else h_sf_denom = (TH1D*) file_data->Get("jet_pt_SF" + denom + "_Muon");

    cout << "Reco" << endl;
    for (int i=1, nreco=recobins.size(); i<nreco; i++)
      cout << h_sf_num->GetBinContent(i)*h_data_num_tunfold->GetBinContent(i)*m_eff[denom]
             /h_sf_denom->GetBinContent(i)/h_data_denom_tunfold->GetBinContent(i)/m_eff[num] << endl;
    cout << endl;

    h_unfolded_ratio = (TH1D*) h_unfolded_num->Clone("h_unfolded_ratio");
    h_unfolded_ratio->Divide( h_unfolded_denom );
    h_unfolded_ratio->Scale( m_eff[denom]/m_eff[num] );

    cout << "Unfolded" << endl;
    for (int i=1, ngen=genbins.size(); i<ngen; i++) cout << h_unfolded_ratio->GetBinContent(i) << endl;
    cout << endl;

    h_unfolded_num_systUP->Divide( h_unfolded_denom_systUP );
    h_unfolded_num_systUP->Scale( m_eff[denom]/m_eff[num] );

    cout << "Unfolded Stat and Systematic" << endl;
    for (int i=1, ngen=genbins.size(); i<ngen; i++) {
      double err = h_unfolded_num_systUP->GetBinContent(i) - h_unfolded_ratio->GetBinContent(i);
      cout << h_unfolded_ratio->GetBinError(i) << "\t" << err << endl;
      h_unfolded_ratio->SetBinError( i, sqrt( h_unfolded_ratio->GetBinError(i)*h_unfolded_ratio->GetBinError(i) + err*err ) );
    }
  }

  // write to file
  TFile* outFile = new TFile(ratio+"_unfold.root","RECREATE");
  outFile->cd();

  if (mc) {
    h_gen_num_rebin1->Write();
    h_gen_denom_rebin1->Write();
    h_gen_num_rebin2->Write();
    h_gen_denom_rebin2->Write();
    h_gen_num_rebin3->Write();
    h_gen_denom_rebin3->Write();
  }
  else {
    h_unfolded_ratio->Write();
  }
  h_response_num_tunfold->Write();
  h_response_denom_tunfold->Write();
  h_unfolded_num->Write();
  h_unfolded_denom->Write();
  h_data_num_tunfold->Write();
  h_data_denom_tunfold->Write();

  outFile->Write();
  delete outFile;
  outFile = 0;
}
