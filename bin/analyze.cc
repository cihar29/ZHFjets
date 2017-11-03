//chad harrington - 11/2/2017

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <time.h>

using namespace std;

//void FillHist1D(const TString& histName, const Double_t& value, const double& weight);
//void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight);
void setPars(const string& parFile); 
//void setWeight(const string& parFile); 

map<TString, TH1*> m_Histos1D;
map<TString, TH2*> m_Histos2D;

TString inName, outName;
string channel;
double weight0, weight;

const int MAXJET = 50;
const int MAXLEP = 2;
//const float MUONMASS = 0.10566;
//const float ELEMASS = 0.;
//const float btagWP_L = 0.5426;
//const float btagWP_M = 0.8484;

int main(int argc, char* argv[]) {

  if (argc == 1) { cout << "Please provide a parameter file." << endl; return -1; }

  string parFile = argv[1];
  setPars(parFile);

  weight0 = 1.;
  //weight0 = -1;
  //if (argc == 3)    { string wFile = argv[2]; setWeight(wFile); }
  //if (weight0 == -1) { cout << "Weight set to 1" << endl; weight0 = 1.; }
  //else                cout << "Weight set to " << weight0 << endl;

  //Open Files//

  TFile* inFile = TFile::Open(inName);
  TTree* T = (TTree*) inFile->Get("tree");
  Long64_t nEntries = T->GetEntries();
  cout << endl << nEntries << " Events" << endl;
  cout << "Processing " + inName << endl;
  cout << "Channel: " + channel << endl;

  //Cuts//

  enum Cuts{
    channelcut, trigcut, selcut, zjetmasscut, njetcut, metcut, hfjetcut, numCuts
  };
  vector<pair<string, double> > v_cuts(numCuts);

  v_cuts[channelcut]=make_pair("Channel",0.); v_cuts[trigcut]=make_pair("Trigger",0.); v_cuts[selcut]=make_pair("Selection",0.);
  v_cuts[zjetmasscut]=make_pair("Z-Jet Mass",0.); v_cuts[njetcut]=make_pair("nJet",0.); v_cuts[metcut]=make_pair("MET",0.);
  v_cuts[hfjetcut]=make_pair("HF Jet",0.);

  //Histograms//
/*
  int nDirs = 2;
  for (int i=0; i<nDirs; i++) {
    TString hname = Form("%i_lep0pt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,200);
    hname = Form("%i_lep1pt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,200);
    hname = Form("%i_lep0eta",i);
    m_Histos1D[hname] = new TH1D(hname,hname,60,-3,3);
    hname = Form("%i_lep1eta",i);
    m_Histos1D[hname] = new TH1D(hname,hname,60,-3,3);
    hname = Form("%i_lep0phi",i);
    m_Histos1D[hname] = new TH1D(hname,hname,60,-3.2,3.2);
    hname = Form("%i_lep1phi",i);
    m_Histos1D[hname] = new TH1D(hname,hname,60,-3.2,3.2);
    hname = Form("%i_lep0iso03",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,1);
    hname = Form("%i_lep1iso03",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,1);
    hname = Form("%i_lep0iso04",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,1);
    hname = Form("%i_lep1iso04",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,1);

    hname = Form("%i_jet_pt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,200,0,200);
    hname = Form("%i_jet_eta",i);
    m_Histos1D[hname] = new TH1D(hname,hname,60,-3,3);
    hname = Form("%i_jet_phi",i);
    m_Histos1D[hname] = new TH1D(hname,hname,60,-3.2,3.2);
    hname = Form("%i_jet_csv",i);
    m_Histos1D[hname] = new TH1D(hname,hname,1000,0,1);
    hname = Form("%i_jet_vtxmass",i);
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,10);
    hname = Form("%i_njet",i);
    m_Histos1D[hname] = new TH1D(hname,hname,10,0,10);
    hname = Form("%i_njet_csvm_vtxmass",i);
    m_Histos1D[hname] = new TH1D(hname,hname,10,0,10);

    hname = Form("%i_V_pt",i);
    m_Histos1D[hname] = new TH1D(hname,hname,500,0,500);
    hname = Form("%i_V_eta",i);
    m_Histos1D[hname] = new TH1D(hname,hname,60,-3,3);
    hname = Form("%i_V_phi",i);
    m_Histos1D[hname] = new TH1D(hname,hname,60,-3.2,3.2);
    hname = Form("%i_V_mass",i);
    m_Histos1D[hname] = new TH1D(hname,hname,500,0,500);

//    hname = Form("%i_V_pt_jet_pt",i);
//    m_Histos2D[hname] = new TH2D(hname,hname,9,{30, 35, 40, 50, 60, 80, 110, 140, 200},9,{0, 20, 30, 40, 50, 70, 90, 120, 200});
  }
*/
  //Set Branches//

  int HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v, HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v,
      HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v, HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v,
      HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v, HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;

  T->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v", &HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v);
  T->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", &HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v);
  T->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v", &HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v);
  T->SetBranchAddress("HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", &HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v);
  T->SetBranchAddress("HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
  T->SetBranchAddress("HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", &HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v);

  int nLep = MAXLEP;
  float vLeptons_new_pt[nLep], vLeptons_new_eta[nLep], vLeptons_new_phi[nLep], vLeptons_new_mass[nLep],
        vLeptons_new_pfRelIso03[nLep], vLeptons_new_pfRelIso04[nLep], vLeptons_new_etaSc[nLep];
  float Vtype_new;

  T->SetBranchAddress("vLeptons_new_pt", vLeptons_new_pt);
  T->SetBranchAddress("vLeptons_new_eta", vLeptons_new_eta);
  T->SetBranchAddress("vLeptons_new_phi", vLeptons_new_phi);
  T->SetBranchAddress("vLeptons_new_mass", vLeptons_new_mass);
  T->SetBranchAddress("vLeptons_new_pfRelIso03", vLeptons_new_pfRelIso03);
  T->SetBranchAddress("vLeptons_new_pfRelIso04", vLeptons_new_pfRelIso04);
  T->SetBranchAddress("vLeptons_new_etaSc", vLeptons_new_etaSc);
  T->SetBranchAddress("Vtype_new", &Vtype_new);

  int nJet=MAXJET;
  int idxJet_passCSV_SVT[2];
  float Jet_pt[nJet], Jet_eta[nJet], Jet_gcc_weight[nJet], Jet_gbb_weight[nJet];

  T->SetBranchAddress("nJet", &nJet);
  T->SetBranchAddress("idxJet_passCSV_SVT", idxJet_passCSV_SVT);
  T->SetBranchAddress("Jet_pt", Jet_pt);
  T->SetBranchAddress("Jet_eta", Jet_eta);
  T->SetBranchAddress("Jet_gcc_weight", Jet_gcc_weight);
  T->SetBranchAddress("Jet_gbb_weight", Jet_gbb_weight);

  float met_pt;

  T->SetBranchAddress("met_pt", &met_pt);

  //Loop Over Entries//
  time_t start = time(NULL);

  for (Long64_t n=0; n<nEntries; n++) {
    T->GetEntry(n);

    TLorentzVector lep0, lep1;
    weight = weight0;

    if (channel == "mm") {
      if (Vtype_new == 0) {
        v_cuts[channelcut].second += weight;

        if ( !HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v && !HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v &&
             !HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v && !HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v ) continue;
        v_cuts[trigcut].second += weight;

        if (vLeptons_new_pt[0] < 25 || vLeptons_new_pt[1] < 25) continue;
        if (fabs(vLeptons_new_eta[0]) > 2.4 || fabs(vLeptons_new_eta[1]) > 2.4) continue;
        if (vLeptons_new_pfRelIso04[0] > 0.25 || vLeptons_new_pfRelIso04[1] > 0.25) continue;

        v_cuts[selcut].second += weight;
      }
      else continue;
    }
    else if (channel == "ee") {
      if (Vtype_new == 1) {
        v_cuts[channelcut].second += weight;

        if (!HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v && !HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v) continue;
        v_cuts[trigcut].second += weight;

        if (vLeptons_new_pt[0] < 25 || vLeptons_new_pt[1] < 25) continue;
        if (fabs(vLeptons_new_eta[0]) > 2.4 || fabs(vLeptons_new_eta[1]) > 2.4) continue;
        if (1.4442 < fabs(vLeptons_new_etaSc[0]) && fabs(vLeptons_new_etaSc[0]) < 1.5660) continue;
        if (1.4442 < fabs(vLeptons_new_etaSc[1]) && fabs(vLeptons_new_etaSc[1]) < 1.5660) continue;
        if (vLeptons_new_pfRelIso03[0] > 0.25 || vLeptons_new_pfRelIso03[1] > 0.25) continue;

        v_cuts[selcut].second += weight;
      }
      else continue;
    }
    lep0.SetPtEtaPhiM(vLeptons_new_pt[0],vLeptons_new_eta[0],vLeptons_new_phi[0],vLeptons_new_mass[0]);
    lep1.SetPtEtaPhiM(vLeptons_new_pt[1],vLeptons_new_eta[1],vLeptons_new_phi[1],vLeptons_new_mass[1]);

    double V_mass = (lep0+lep1).M();
    if (V_mass<70 || V_mass>110) continue;
    v_cuts[zjetmasscut].second += weight;

    int nGoodJet = 0;
    for (int i=0; i<nJet; i++) {
      if (Jet_pt[i]>30 && fabs(Jet_eta[i])<2.4) nGoodJet++;
    }
    if (nGoodJet < 1) continue;
    v_cuts[njetcut].second += weight;

    if (met_pt > 40) continue;
    v_cuts[metcut].second += weight;

    int idx = idxJet_passCSV_SVT[1];
    if (idx < 0) continue;
    if (Jet_pt[idx] < 30 || fabs(Jet_eta[idx]) > 2.4) continue;
    if (Jet_gcc_weight[idx] > 1.1 || Jet_gbb_weight[idx] > 1.1) continue;

    v_cuts[hfjetcut].second += weight;
  }
  cout << difftime(time(NULL), start) << " s" << endl;

  //TH1D* cuts = new TH1D("cuts","cuts",numCuts,-0.5,float(numCuts)-0.5);

  //Cutflow Table//
  cout<<"===================================================================================================\n";
  cout<<"                                Cut Flow Table: " + inName( inName.Last('/')+1, inName.Index('.')-inName.Last('/')-1 ) + "\n";
  cout<<"===================================================================================================\n";

  cout<<      "                          |||          Nevent          |||     Efficiency (Relative Efficiency)\n";

  for (int i=0; i<numCuts; i++) {
    if (i == 0)
      cout << Form("%-25s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[0].second) << endl;
    else
      cout << Form("%-25s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[i-1].second) << endl;

    //cuts->SetBinContent(i+1, v_cuts[i].second);
    //cuts->GetXaxis()->SetBinLabel(i+1,Form("%i",i+1));
  }
  cout << endl;

  //Write Histograms//
/*
  TFile* outFile = new TFile(outName,"RECREATE");
  outFile->cd();

  cuts->Write();

  for (int i=0; i<nDirs; i++) outFile->mkdir( Form("%i/", i) );

  for (map<TString, TH1*>::iterator hid = m_Histos1D.begin(); hid != m_Histos1D.end(); hid++) {
    outFile->cd();
    TString prefix = hid->first(0, 1);

    if ( prefix.IsDigit() ) outFile->cd(outName + ":/" + prefix);

    hid->second->Write();
  }
  for (map<TString, TH2*>::iterator hid = m_Histos2D.begin(); hid != m_Histos2D.end(); hid++) {
    outFile->cd();
    TString prefix = hid->first(0, 1);

    if ( prefix.IsDigit() ) outFile->cd(outName + ":/" + prefix);

    hid->second->Write();
  }

  outFile->Write();
  delete outFile;
  outFile = 0;
*/
}
/*
void FillHist1D(const TString& histName, const Double_t& value, const double& weight) {
  map<TString, TH1*>::iterator hid=m_Histos1D.find(histName);
  if (hid==m_Histos1D.end())
    cout << "%FillHist1D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value, weight);
}

void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight) {
  map<TString, TH2*>::iterator hid=m_Histos2D.find(histName);
  if (hid==m_Histos2D.end())
    cout << "%FillHist2D -- Could not find histogram with name: " << histName << endl;
  else
    hid->second->Fill(value1, value2, weight);
}

void FillHists(const TString& prefix, const int& nEle, const int& nGoodEle, const int& nMuon, const int& nGoodMuon, const int& nJet, const int& nGoodJet,
               const float& dphi_jet0met, const float& dphi_jet1met, const int& nPV) {

  FillHist1D(prefix+"nEleDiff", nEle-nGoodEle, weight);
  FillHist1D(prefix+"nMuonDiff", nMuon-nGoodMuon, weight);
}

void setWeight(const string& wFile) {

  ifstream file(wFile);
  string line;

  TString name = inName( inName.Last('/')+1, inName.Index('.')-inName.Last('/')-1 );
  while (getline(file, line)) {

    if (line.length() > 0) {
      while (line.at(0) == ' ') line.erase(0, 1);
    }

    int delim_pos = line.find(' ');
    if (delim_pos == -1) continue;

    TString dataset = line.substr(0, delim_pos).data();
    if ( name.EqualTo(dataset, TString::kIgnoreCase) ) {

      while (line.at(line.length()-1) == ' ') line.erase(line.length()-1, line.length());

      //weight is found in the last column
      delim_pos = line.length()-1;
      while (line.at(delim_pos) != ' ') delim_pos--;

      line.erase(0, delim_pos+1);
      weight0 = stod(line);
      break;
    }
  }
  file.close();
}
*/
void setPars(const string& parFile) {

  ifstream file(parFile);
  string line;

  while (getline(file, line)) {

    if (line.length() > 0){
      while (line.at(0) == ' ') line.erase(0, 1);
    }

    int delim_pos = line.find(' ');
    if (delim_pos == -1) continue;

    string var = line.substr(0, delim_pos);
    line.erase(0, delim_pos + 1);

    while (line.at(0) == ' ') line.erase(0, 1);
    while (line.at(line.length()-1) == ' ') line.erase(line.length()-1, line.length());

    if      (var == "inName")  inName = line.data();
    else if (var == "outName") outName = line.data();
    else if (var == "channel") channel = line;
  }
  file.close();
}
