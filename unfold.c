vector<double> genbins =  {30, 40, 60, 80, 110, 200};
vector<double> recobins = {30, 35, 40, 50, 60, 80, 110, 140, 200};

void rebin2D( TH2D*& h_response ) {

  int nx = h_response->GetXaxis()->GetNbins(), ny = h_response->GetYaxis()->GetNbins();
  //const double* xbins = h_response->GetXaxis()->GetXbins()->GetArray(), *ybins = h_response->GetYaxis()->GetXbins()->GetArray(); // this only works for variable binning
  vector<double> xbins(nx+1), ybins(ny+1);
  for (int i=1, nxbins=xbins.size(); i<=nxbins; i++) xbins[i-1] = h_response->GetXaxis()->GetBinLowEdge(i);
  for (int i=1, nybins=ybins.size(); i<=nybins; i++) ybins[i-1] = h_response->GetYaxis()->GetBinLowEdge(i);

  // PRINT MATRIX BEFORE REBINNING //
//  for (int i=1, nxbins=xbins.size(); i<nxbins; i++) {
//    for (int j=1, nybins=ybins.size(); j<nybins; j++)
//      cout << xbins[i-1] << "-" << xbins[i] << " " << ybins[j-1] << "-" << ybins[j] << "\t"
//           << h_response->GetBinContent(i,j) << "\t" << h_response->GetBinError(i,j) << endl;
//  }

  TH2D* h_response2 = new TH2D((TString)h_response->GetName()+"_rebin",(TString)h_response->GetName()+"_rebin",genbins.size()-1,&genbins[0],recobins.size()-1,&recobins[0]);
  for (int i=1, ngen=genbins.size(); i<ngen; i++) {
    for (int j=1, nreco=recobins.size(); j<nreco; j++) {

      double merged=0, merged_err2=0;
      for (int g=1, nxbins=xbins.size(); g<nxbins; g++) {
        if (xbins[g-1]<genbins[i-1] || xbins[g-1]>=genbins[i]) continue;

        for (int h=1, nybins=ybins.size(); h<nybins; h++) {
          if (ybins[h-1]<recobins[j-1] || ybins[h-1]>=recobins[j]) continue;
          merged      += h_response->GetBinContent(g, h);
          merged_err2 += h_response->GetBinError(g, h) * h_response->GetBinError(g, h);
        }
      }
      h_response2->SetBinContent(i, j, merged);          //merged<=0?0:merged);
      h_response2->SetBinError(i, j, sqrt(merged_err2)); //merged<=0?0:sqrt(merged_err2));
//      h_response2->SetBinContent(i, j, i==j?1:0.01);          //merged<=0?0:merged);
//      h_response2->SetBinError(i, j, 0.001); //merged<=0?0:sqrt(merged_err2));
    }
  }
  h_response = h_response2;
}

void unfold(TString ratio = "rcb") {

  TString num, denom;
  if      (ratio == "rcj") { num = "c"; denom = "incl"; }
  else if (ratio == "rbj") { num = "b"; denom = "incl"; }
  else                     { num = "c"; denom = "b"; }

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TFile* file1 = TFile::Open( "plots_wgt.root" );
  TFile* file2 = TFile::Open( "plots_wgt.root" );
//  TFile* file2 = TFile::Open( "duong.root" );

  TH2D* h_response_num   = (TH2D*) file1->Get(num + "_jet_genjet");
  TH2D* h_response_denom = (TH2D*) file1->Get(denom + "_jet_genjet");

  TH1D* h_data_num   = (TH1D*) file2->Get(num + "_jet30");    //c_jet30, c_Ljet
  h_data_num = (TH1D*) h_data_num->Rebin(recobins.size()-1, "h_data_num_rebin", &recobins[0]);
  TH1D* h_data_denom = (TH1D*) file2->Get(denom + "_jet30"); //c_jet30, c_Ljet
  h_data_denom = (TH1D*) h_data_denom->Rebin(recobins.size()-1, "h_data_denom_rebin", &recobins[0]);
  TH1D* h_data = (TH1D*) h_data_num->Clone("h_data");
  h_data->Divide(h_data_denom);

//  TH1D* h_data = (TH1D*) file2->Get(ratio + "_jet_pt_mm");

  // CONVERT DUONG'S HISTS //

//  TH1D* h_data2 = new TH1D("h_data_rebin","h_data_rebin",recobins.size()-1,&recobins[0]);
//  for (int i=1; i<recobins.size(); i++) { h_data2->SetBinContent( i, h_data->GetBinContent(i) );  h_data2->SetBinError( i, h_data->GetBinError(i) ); }
//  h_data = h_data2;

  // REBIN RESPONSE MATRIX //

  rebin2D(h_response_num);
  rebin2D(h_response_denom);
  TH2D* h_response = (TH2D*) h_response_num->Clone("h_response");
  h_response->Divide(h_response_denom);

  // PRINT MATRIX AFTER REBINNING //
//  for (int i=1, ngen=genbins.size(); i<ngen; i++) {
//    for (int j=1, nreco=recobins.size(); j<nreco; j++)
//      cout << genbins[i-1] << "-" << genbins[i] << " " << recobins[j-1] << "-" << recobins[j] << "\t"
//           << h_response->GetBinContent(i,j) << "\t" << h_response->GetBinError(i,j) << endl;
//  }

  // UNFOLD //

  TUnfoldBinning genBinning("signal"), recoBinning("background");
  genBinning.AddAxis( "genpt", genbins.size()-1, &genbins[0], false, false );
  recoBinning.AddAxis( "recopt", recobins.size()-1, &recobins[0], false, false );

  TUnfoldDensity unfold( h_response, TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeSize, TUnfold::kEConstraintArea, TUnfoldDensity::kDensityModeBinWidth,
                         &genBinning, &recoBinning, 0, "*[UOB]" );

  unfold.SetInput( h_data );
  //unfold.DoUnfold(0.00001);

  int nScan=100;
  TSpline* rhoLogTau=0, *logTauX=0, *logTauY=0;
  TGraph* lCurve=0;

  unfold.ScanTau( nScan, 0, 0, &rhoLogTau, TUnfoldDensity::kEScanTauRhoMax, "signal", "*[UOB]", &lCurve, &logTauX, &logTauY );

  TH1* h_unfolded = unfold.GetOutput("unfolded_signal");

  // VALIDATION //

  TH1D* h_gen_num   = (TH1D*) file2->Get(num + "_genjet30");    //c_genjet30, c_Lgenjet
  h_gen_num = (TH1D*) h_gen_num->Rebin(genbins.size()-1, "h_gen_num_rebin", &genbins[0]);
  TH1D* h_gen_denom = (TH1D*) file2->Get(denom + "_genjet30"); //c_genjet30, c_Lgenjet
  h_gen_denom = (TH1D*) h_gen_denom->Rebin(genbins.size()-1, "h_gen_denom_rebin", &genbins[0]);
  TH1D* h_gen = (TH1D*) h_gen_num->Clone("h_gen");
  h_gen->Divide(h_gen_denom);

  for (int i=1, ngen=genbins.size(); i<ngen; i++) cout << h_gen->GetBinContent(i) << "\t" << h_unfolded->GetBinContent(i) << endl;

  h_gen->GetYaxis()->SetRangeUser(0,3);
  h_gen->Draw();
  h_data->SetLineColor(kGreen);
  h_data->Draw("same");
  h_unfolded->SetLineColor(kRed);
  h_unfolded->Draw("same");
return;

  // OUTPUT //

  TFile* outFile = new TFile("unfold.root","RECREATE");
  outFile->cd();

  h_response->Write();
  h_response_num->Write();
  h_response_denom->Write();
  h_data->Write();
  h_unfolded->Write();

  outFile->Write();
  delete outFile;
  outFile = 0;
}
