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

bool isMC;
TString inName;
string channel;
double weight0, weight;

const int MAXJET = 50;
const int MAXLEP = 2;

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
  cout << nEntries << " Events" << endl;
  cout << "Processing " + inName << endl;
  cout << "Channel: " + channel << endl;

  //Cuts//

  enum Cuts{
    channelcut, trigcut, selcut, zjetmasscut, metcut, jetcut, bjet, cjet, ljet, genjetcut, bgenjet, cgenjet, lgenjet, numCuts
  };
  vector<pair<string, double> > v_cuts(numCuts);

  v_cuts[channelcut]=make_pair("Channel",0.); v_cuts[trigcut]=make_pair("Trigger",0.); v_cuts[selcut]=make_pair("Selection",0.);
  v_cuts[zjetmasscut]=make_pair("Z-Jet Mass",0.); v_cuts[jetcut]=make_pair("GoodJet",0.); v_cuts[metcut]=make_pair("MET",0.);
  v_cuts[bjet]=make_pair("bJet",0.); v_cuts[cjet]=make_pair("cJet",0.); v_cuts[ljet]=make_pair("lJet",0.); v_cuts[genjetcut]=make_pair("GoodGenJet",0.);
  v_cuts[bgenjet]=make_pair("bgenjet",0.); v_cuts[cgenjet]=make_pair("cgenjet",0.); v_cuts[lgenjet]=make_pair("lgenjet",0.);

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

  int nJet = MAXJET;
  int idxJet_passCSV_SVT[2], Jet_hadronFlavour[nJet];
  float Jet_pt[nJet], Jet_eta[nJet], Jet_phi[nJet], Jet_mass[nJet];

  T->SetBranchAddress("nJet", &nJet);
  T->SetBranchAddress("idxJet_passCSV_SVT", idxJet_passCSV_SVT);
  T->SetBranchAddress("Jet_pt", Jet_pt);
  T->SetBranchAddress("Jet_eta", Jet_eta);
  T->SetBranchAddress("Jet_phi", Jet_phi);
  T->SetBranchAddress("Jet_mass", Jet_mass);

  int nGenJet = MAXJET;
  int GenJet_numBHadrons[nGenJet], GenJet_numCHadrons[nGenJet]; //GenJet_pdgId[nGenJet];
  float GenJet_pt[nGenJet], GenJet_eta[nGenJet], GenJet_phi[nGenJet], GenJet_mass[nGenJet];

  if (isMC) {
    T->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour);

    T->SetBranchAddress("nGenJet", &nGenJet);
    T->SetBranchAddress("GenJet_pt", GenJet_pt);
    T->SetBranchAddress("GenJet_eta", GenJet_eta);
    T->SetBranchAddress("GenJet_phi", GenJet_phi);
    T->SetBranchAddress("GenJet_mass", GenJet_mass);
    //T->SetBranchAddress("GenJet_pdgId", GenJet_pdgId);
    T->SetBranchAddress("GenJet_numBHadrons", GenJet_numBHadrons);
    T->SetBranchAddress("GenJet_numCHadrons", GenJet_numCHadrons);
  }

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

    if (met_pt > 40) continue;
    v_cuts[metcut].second += weight;

    map<int, int> genjet_jet_idx;
    int hfidx = idxJet_passCSV_SVT[1], genhfidx = -1;

    for (int i=0; i<nJet; i++) {
      TLorentzVector jet;
      jet.SetPtEtaPhiM(Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i]);

      double rmin = 99.;
      int matched_idx = -1;

      for (int j=0; j<nGenJet; j++) {
        if (genjet_jet_idx.find(j) != genjet_jet_idx.end()) continue;  //this genjet is already matched to a jet

        TLorentzVector genJet;
        genJet.SetPtEtaPhiM(GenJet_pt[j], GenJet_eta[j], GenJet_phi[j], GenJet_mass[j]);

        double dR = jet.DeltaR(genJet);
        if ( dR<rmin && dR<0.2 ) {
          rmin = dR;
          matched_idx = j;
        }
      }
      if (matched_idx != -1) genjet_jet_idx[matched_idx] = i;
      if (hfidx == i)        genhfidx = matched_idx;
    }

    //loop over matched jets
    bool goodJet = false;
    for (int i=0; i<nJet; i++) {

      bool matched = false;
      for (auto const& m : genjet_jet_idx) {
        if (i==m.second) {matched = true; break;}
      }
      if ( matched && Jet_pt[i]>30 && fabs(Jet_eta[i])<2.4 ) {goodJet = true; break;}
    }
    if (goodJet) {
      v_cuts[jetcut].second += weight;

      if (hfidx >=0 ) {
        if (Jet_pt[hfidx]<30 || fabs(Jet_eta[hfidx])>2.4) continue;

        if      (Jet_hadronFlavour[hfidx] == 5) v_cuts[bjet].second += weight;
        else if (Jet_hadronFlavour[hfidx] == 4) v_cuts[cjet].second += weight;
        else                                    v_cuts[ljet].second += weight;
      }
    }

    //loop over matched genjets
    bool goodGenJet = false;
    for (int j=0; j<nGenJet; j++) {
      if (genjet_jet_idx.find(j) == genjet_jet_idx.end()) continue;

      if ( GenJet_pt[j]>30 && fabs(GenJet_eta[j])<2.4 ) {goodGenJet = true; break;}
    }
    if (goodGenJet) {
      v_cuts[genjetcut].second += weight;

      if (genhfidx >=0 ) {
        if (GenJet_pt[genhfidx]<30 || fabs(GenJet_eta[genhfidx])>2.4) continue;

        //if      (abs(GenJet_pdgId[genhfidx]) == 5) v_cuts[bgenjet].second += weight;
        //else if (abs(GenJet_pdgId[genhfidx]) == 4) v_cuts[cgenjet].second += weight;
        if      (GenJet_numBHadrons[genhfidx] >= 1) v_cuts[bgenjet].second += weight;
        else if (GenJet_numCHadrons[genhfidx] >= 1) v_cuts[cgenjet].second += weight;
        else                                        v_cuts[lgenjet].second += weight;
      }
    }
  }
  cout << difftime(time(NULL), start) << " s" << endl;

  //Cutflow Table//
  cout<<"===================================================================================================\n";
  cout<<"                                Cut Flow Table: " + inName( inName.Last('/')+1, inName.Index('.')-inName.Last('/')-1 ) + "\n";
  cout<<"===================================================================================================\n";

  cout<<      "                          |||          Nevent          |||     Efficiency (Relative Efficiency)\n";

  for (int i=0; i<numCuts; i++) {
    if (i == 0)
      cout << Form("%-25s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[0].second) << endl;

    else if (i>genjetcut)
      cout << Form("%-25s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[genjetcut].second) << endl;
    else if (i>jetcut)
      cout << Form("%-25s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[jetcut].second) << endl;
    else
      cout << Form("%-25s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[i-1].second) << endl;

   if (i==metcut || i==ljet)
      cout << "---------------------------------------------------------------------------------------------------" << endl;
  }
  cout << endl;
}
/*
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

    if (var == "isMC") {
      if (line == "true") isMC = true;
      else isMC = false;
    }
    else if (var == "inName")  inName = line.data();
    else if (var == "channel") channel = line;
  }
  file.close();
}
