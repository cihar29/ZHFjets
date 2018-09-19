//chad harrington - 11/2/2017

#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <time.h>

using namespace std;

void FillHist1D(const TString& histName, const Double_t& value, const double& weight);
void FillHist2D(const TString& histName, const Double_t& value1, const Double_t& value2, const double& weight);
void setPars(const string& parFile);
//void setWeight(const string& parFile);

map<TString, TH1*> m_Histos1D;
map<TString, TH2*> m_Histos2D;

bool isMC = false;
TString inName, outName;
string channel;
double xs, lumi, posWgt, negWgt;

const int MAXJET = 50;
const int MAXLEP = 4;

vector<double> genbins =  {30, 40, 50, 60, 80, 110, 140, 200};                                    //Jet pt
vector<double> recobins = {30, 35, 40, 45, 50, 55, 60, 70, 80, 95, 110, 125, 140, 170, 200, 250}; //Jet pt
//vector<double> recobins = {30, 35, 40, 45, 50, 55, 60, 70, 80, 95, 110, 125, 140, 170, 200}; //Jet pt

vector<double> genbins2 = {0, 20, 30, 40, 50, 70, 90, 120, 200};                                         //Vpt
vector<double> recobins2 = {0, 10, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 105, 120, 160, 200, 250}; //Vpt
//vector<double> recobins2 = {0, 10, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 105, 120, 160, 200}; //Vpt
float max_val = 225;

int main(int argc, char* argv[]) {

  if (argc == 1) { cout << "Please provide a parameter file." << endl; return -1; }

  string parFile = argv[1];
  setPars(parFile);

  double weight0 = 1., weight;
  if (xs != 0) weight0 = lumi * xs * 1000 / (posWgt-negWgt);

  //Open Files//

  TChain* T = new TChain("tree");
  if (inName == "madgraph") {
    TString path = "root://131.225.204.161:1094//store/user/yokugawa/ZPJ_datasets_2017-03_MC_SUMMER2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MC_SUMMER2016_PRIMARY_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v2/170224_012042/";
    for (int i=1; i<480; i++)
      T->Add(path + "0000/tree_" + to_string(i) + ".root");
    path = "root://131.225.204.161:1094//store/user/yokugawa/ZPJ_datasets_2017-03_MC_SUMMER2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MC_SUMMER2016_PRIMARY_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext2-v1/170224_012132/";
    for (int i=1; i<840; i++)
      T->Add(path + "0000/tree_" + to_string(i) + ".root");
  }
  else if (inName == "DY0J") {
    TString path = "root://131.225.204.161:1094//store/group/lpcphys/noreplica/duong/VHbbTuples/MC_V25/DY_amcnlo/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/MCSUMMER2016_TEST_DYToLL_0J_13TeV-amcatnloFXFX-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v1/180126_003912/";
    for (int i=1; i<665; i++)
      T->Add(path + "0000/tree_" + to_string(i) + ".root");
    path = "root://131.225.204.161:1094//store/group/lpcphys/noreplica/duong/VHbbTuples/MC_V25/DY_amcnlo/DYToLL_0J_13TeV-amcatnloFXFX-pythia8/MCSUMMER2016_TEST_DYToLL_0J_13TeV-amcatnloFXFX-Py8__RunIISummer16MAv2-PUMoriond17_backup_80r2as_2016_TrancheIV_v6-v1/180126_030824/";
    for (int i=1; i<605; i++)
      T->Add(path + "0000/tree_" + to_string(i) + ".root");
  }
  else if (inName == "DY1J") {
    TString path = "root://131.225.204.161:1094//store/group/lpcphys/noreplica/duong/VHbbTuples/MC_V25/DY_amcnlo/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/MCSUMMER2016_TEST_DYToLL_1J_13TeV-amcatnloFXFX-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v1/180126_030917/";
    for (int i=1; i<800; i++)
      T->Add(path + "0000/tree_" + to_string(i) + ".root");
    path = "root://131.225.204.161:1094//store/group/lpcphys/noreplica/duong/VHbbTuples/MC_V25/DY_amcnlo/DYToLL_1J_13TeV-amcatnloFXFX-pythia8/MCSUMMER2016_TEST_DYToLL_1J_13TeV-amcatnloFXFX-Py8__RunIISummer16MAv2-PUMoriond17_backup_80r2as_2016_TrancheIV_v6-v1/180126_031023/";
    for (int i=1; i<560; i++)
      T->Add(path + "0000/tree_" + to_string(i) + ".root");
  }
  else if (inName == "DY2J") {
    TString path = "root://131.225.204.161:1094//store/group/lpcphys/noreplica/duong/VHbbTuples/MC_V25/DY_amcnlo/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/MCSUMMER2016_TEST_DYToLL_2J_13TeV-amcatnloFXFX-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6-v2/180126_031129/";
    for (int i=1; i<755; i++)
      T->Add(path + "0000/tree_" + to_string(i) + ".root");
    path = "root://131.225.204.161:1094//store/group/lpcphys/noreplica/duong/VHbbTuples/MC_V25/DY_amcnlo/DYToLL_2J_13TeV-amcatnloFXFX-pythia8/MCSUMMER2016_TEST_DYToLL_2J_13TeV-amcatnloFXFX-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext1-v1/180126_031239/";
    for (int i=1; i<671; i++)
      T->Add(path + "0000/tree_" + to_string(i) + ".root");
  }
  else {
    for (int i=1; i<1000; i++)
      T->Add(inName + "0000/tree_" + to_string(i) + ".root");
    for (int i=1000; i<1455; i++)
      T->Add(inName + "0001/tree_" + to_string(i) + ".root");
  }

  Long64_t nEntries = T->GetEntries();
  cout << nEntries << " Events" << endl;
  cout << "Processing " + inName << endl;
  cout << "Channel: " + channel << endl;
  cout << "weight0: " << weight0  << endl;

  //Cuts//

  enum Cuts{
    channelcut, trigcut, selcut, zjetmasscut, metcut, jetcut, bjet, cjet, ljet, jetcut30, bjetcut30, cjetcut30, ljetcut30, numCuts
  };
  vector<pair<string, double> > v_cuts(numCuts);

  v_cuts[channelcut]=make_pair("Channel",0.);     v_cuts[trigcut]=make_pair("Trigger",0.); v_cuts[selcut]=make_pair("Selection",0.);
  v_cuts[zjetmasscut]=make_pair("Z-Jet Mass",0.); v_cuts[metcut]=make_pair("MET",0.);      v_cuts[jetcut]=make_pair("Matched Jet",0.);
  v_cuts[bjet]=make_pair("bJet",0.);              v_cuts[cjet]=make_pair("cJet",0.);       v_cuts[ljet]=make_pair("lJet",0.);
  v_cuts[jetcut30]=make_pair("Matched Jet, pt>30GeV",0.);
  v_cuts[bjetcut30]=make_pair("bJetcut30",0.); v_cuts[cjetcut30]=make_pair("cJetcut30",0.); v_cuts[ljetcut30]=make_pair("lJetcut30",0.);

  //Histograms//
  // Gen Tagged for acceptance study //

  TString hname = "b_genjet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_genjet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_genjet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_genjet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "b_jet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_jet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_jet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_jet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "b_genjet30";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_genjet30";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_genjet30";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_genjet30";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "b_jet30";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_jet30";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_jet30";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_jet30";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "b_jet_genjet";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
  hname = "c_jet_genjet";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
  hname = "udsg_jet_genjet";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
  hname = "incl_jet_genjet";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);

  // Reco tagged for unfolding //

  hname = "b_Ljet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_Ljet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_Ljet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_Ljet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "b_Ljet_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_Ljet_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_Ljet_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_Ljet_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "b_Lgenjet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_Lgenjet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_Lgenjet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_Lgenjet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "b_Lgenjet_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_Lgenjet_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_Lgenjet_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_Lgenjet_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "b_Lgenjet_unmtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_Lgenjet_unmtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_Lgenjet_unmtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_Lgenjet_unmtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

//  hname = "b_Ljet_Lgenjet_mtchd";
//  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
//  hname = "c_Ljet_Lgenjet_mtchd";
//  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
//  hname = "udsg_Ljet_Lgenjet_mtchd";
//  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
//  hname = "incl_Ljet_Lgenjet_mtchd";
//  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);

  hname = "b_Ljet_Lgenjet";
  m_Histos2D[hname] = new TH2D(hname,hname, genbins.size()-1, &genbins[0], recobins.size()-1, &recobins[0]);
  hname = "c_Ljet_Lgenjet";
  m_Histos2D[hname] = new TH2D(hname,hname, genbins.size()-1, &genbins[0], recobins.size()-1, &recobins[0]);
  hname = "udsg_Ljet_Lgenjet";
  m_Histos2D[hname] = new TH2D(hname,hname, genbins.size()-1, &genbins[0], recobins.size()-1, &recobins[0]);
  hname = "incl_Ljet_Lgenjet";
  m_Histos2D[hname] = new TH2D(hname,hname, genbins.size()-1, &genbins[0], recobins.size()-1, &recobins[0]);

  hname = "b_Vpt";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_Vpt";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_Vpt";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_Vpt";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "b_Vpt_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_Vpt_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_Vpt_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_Vpt_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "b_genVpt_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_genVpt_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_genVpt_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_genVpt_mtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "b_genVpt_unmtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "c_genVpt_unmtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "udsg_genVpt_unmtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_genVpt_unmtchd";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

//  hname = "b_Vpt_genVpt_mtchd";
//  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
//  hname = "c_Vpt_genVpt_mtchd";
//  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
//  hname = "udsg_Vpt_genVpt_mtchd";
//  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
//  hname = "incl_Vpt_genVpt_mtchd";
//  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);

  hname = "b_Vpt_genVpt";
  m_Histos2D[hname] = new TH2D(hname,hname, genbins2.size()-1, &genbins2[0], recobins2.size()-1, &recobins2[0]);
  hname = "c_Vpt_genVpt";
  m_Histos2D[hname] = new TH2D(hname,hname, genbins2.size()-1, &genbins2[0], recobins2.size()-1, &recobins2[0]);
  hname = "udsg_Vpt_genVpt";
  m_Histos2D[hname] = new TH2D(hname,hname, genbins2.size()-1, &genbins2[0], recobins2.size()-1, &recobins2[0]);
  hname = "incl_Vpt_genVpt";
  m_Histos2D[hname] = new TH2D(hname,hname, genbins2.size()-1, &genbins2[0], recobins2.size()-1, &recobins2[0]);

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
  float vLeptons_pt[nLep], vLeptons_eta[nLep], vLeptons_phi[nLep], vLeptons_mass[nLep],
        vLeptons_pfRelIso03[nLep], vLeptons_pfRelIso04[nLep], vLeptons_etaSc[nLep];
  float Vtype;

  T->SetBranchAddress("vLeptons_pt", vLeptons_pt);
  T->SetBranchAddress("vLeptons_eta", vLeptons_eta);
  T->SetBranchAddress("vLeptons_phi", vLeptons_phi);
  T->SetBranchAddress("vLeptons_mass", vLeptons_mass);
  T->SetBranchAddress("vLeptons_pfRelIso03", vLeptons_pfRelIso03);
  T->SetBranchAddress("vLeptons_pfRelIso04", vLeptons_pfRelIso04);
  T->SetBranchAddress("vLeptons_etaSc", vLeptons_etaSc);
  T->SetBranchAddress("Vtype", &Vtype);

  int nJet = MAXJET;
  int Jet_hadronFlavour[nJet];
  float Jet_pt[nJet], Jet_eta[nJet], Jet_phi[nJet], Jet_mass[nJet]; //, Jet_btagCSV[nJet], Jet_ctagVsL[nJet], Jet_ctagVsB[nJet];

  T->SetBranchAddress("nJet", &nJet);
  T->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);
  T->SetBranchAddress("Jet_pt", Jet_pt);
  T->SetBranchAddress("Jet_eta", Jet_eta);
  T->SetBranchAddress("Jet_phi", Jet_phi);
  T->SetBranchAddress("Jet_mass", Jet_mass);
  //T->SetBranchAddress("Jet_btagCSV", Jet_btagCSV);
  //T->SetBranchAddress("Jet_ctagVsL", Jet_ctagVsL);
  //T->SetBranchAddress("Jet_ctagVsB", Jet_ctagVsB);

  int nGenJet = MAXJET;
  int GenJet_numBHadrons[nGenJet], GenJet_numCHadrons[nGenJet]; //GenJet_pdgId[nGenJet];
  float GenJet_pt[nGenJet], GenJet_eta[nGenJet], GenJet_phi[nGenJet], GenJet_mass[nGenJet];

  int nGenLep = MAXLEP;
  float GenLep_pt[nGenLep], GenLep_eta[nGenLep], GenLep_phi[nGenLep], GenLep_mass[nGenLep];

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

    T->SetBranchAddress("nGenLep", &nGenLep);
    T->SetBranchAddress("GenLep_pt", GenLep_pt);
    T->SetBranchAddress("GenLep_eta", GenLep_eta);
    T->SetBranchAddress("GenLep_phi", GenLep_phi);
    T->SetBranchAddress("GenLep_mass", GenLep_mass);
  }

  float met_pt, puWeight, genWeight;

  T->SetBranchAddress("met_pt", &met_pt);
  T->SetBranchAddress("puWeight", &puWeight);
  T->SetBranchAddress("genWeight", &genWeight);

  //Loop Over Entries//
  time_t start = time(NULL);

  for (Long64_t n=0; n<nEntries; n++) {
    T->GetEntry(n);

    TLorentzVector lep0, lep1;
    weight = weight0*puWeight*TMath::Sign(1,genWeight);

    if (channel == "mm") {
      if (Vtype == 0) {
        v_cuts[channelcut].second += weight;

        if ( !HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v && !HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v &&
             !HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v && !HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v ) continue;
        v_cuts[trigcut].second += weight;

        if (vLeptons_pt[0] < 25 || vLeptons_pt[1] < 25) continue;
        if (fabs(vLeptons_eta[0]) > 2.4 || fabs(vLeptons_eta[1]) > 2.4) continue;
        if (vLeptons_pfRelIso04[0] > 0.25 || vLeptons_pfRelIso04[1] > 0.25) continue;

        v_cuts[selcut].second += weight;
      }
      else continue;
    }
    else if (channel == "ee") {
      if (Vtype == 1) {
        v_cuts[channelcut].second += weight;

        if (!HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v && !HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v) continue;
        v_cuts[trigcut].second += weight;

        if (vLeptons_pt[0] < 25 || vLeptons_pt[1] < 25) continue;
        if (fabs(vLeptons_eta[0]) > 2.4 || fabs(vLeptons_eta[1]) > 2.4) continue;
        if (1.4442 < fabs(vLeptons_etaSc[0]) && fabs(vLeptons_etaSc[0]) < 1.5660) continue;
        if (1.4442 < fabs(vLeptons_etaSc[1]) && fabs(vLeptons_etaSc[1]) < 1.5660) continue;
        if (vLeptons_pfRelIso03[0] > 0.25 || vLeptons_pfRelIso03[1] > 0.25) continue;

        v_cuts[selcut].second += weight;
      }
      else continue;
    }
    lep0.SetPtEtaPhiM(vLeptons_pt[0],vLeptons_eta[0],vLeptons_phi[0],vLeptons_mass[0]);
    lep1.SetPtEtaPhiM(vLeptons_pt[1],vLeptons_eta[1],vLeptons_phi[1],vLeptons_mass[1]);

    double V_mass = (lep0+lep1).M();
    if (V_mass<71 || V_mass>111) continue;
    v_cuts[zjetmasscut].second += weight;

    if (met_pt > 40) continue;
    v_cuts[metcut].second += weight;

    // Gen matching //
    map<int, int> matched0, matched30; //<reco idx, gen idx>
    for (int igen=0; igen<nGenJet; igen++) {

      TString flav = "udsg";
      if      (GenJet_numBHadrons[igen] >= 1) flav = "b";
      else if (GenJet_numCHadrons[igen] >= 1) flav = "c";

      TLorentzVector genJet;
      genJet.SetPtEtaPhiM(GenJet_pt[igen], GenJet_eta[igen], GenJet_phi[igen], GenJet_mass[igen]);

      // No cut on pt //
      double rmin = 99.;
      int idx = -1;
      for (int ireco=0; ireco<nJet; ireco++) {
        if (fabs(Jet_eta[ireco])>2.4) continue;
        if ( matched0.find(ireco) != matched0.end() ) continue;  // this jet is already matched to a genjet

        TLorentzVector jet;
        jet.SetPtEtaPhiM(Jet_pt[ireco], Jet_eta[ireco], Jet_phi[ireco], Jet_mass[ireco]);

        double dR = jet.DeltaR(genJet);
        if ( dR<rmin && dR<0.2 ) {
          rmin = dR;
          idx = ireco;
        }
      }
      if (idx != -1) {
        matched0[idx] = igen;

        FillHist1D(flav + "_genjet", GenJet_pt[igen], weight);
        FillHist1D(flav + "_jet", Jet_pt[idx], weight);
        FillHist1D("incl_genjet", GenJet_pt[igen], weight);
        FillHist1D("incl_jet", Jet_pt[idx], weight);
      }
/*
      // cut on pt<30 GeV //
      // for acceptance studies //
      rmin = 99.;
      idx = -1;
      for (int ireco=0; ireco<nJet; ireco++) {
        if (Jet_pt[ireco]<30 || fabs(Jet_eta[ireco])>2.4) continue;
        if ( matched30.find(ireco) != matched30.end() ) continue;  // this jet is already matched to a genjet

        TLorentzVector jet;
        jet.SetPtEtaPhiM(Jet_pt[ireco], Jet_eta[ireco], Jet_phi[ireco], Jet_mass[ireco]);

        double dR = jet.DeltaR(genJet);
        if ( dR<rmin && dR<0.2 ) {
          rmin = dR;
          idx = ireco;
        }
      }
      if (idx != -1) {
        matched30[idx] = igen;

        FillHist1D(flav + "_genjet30", GenJet_pt[igen], weight);
        FillHist1D(flav + "_jet30", Jet_pt[idx], weight);
        FillHist2D(flav + "_jet_genjet", GenJet_pt[igen], Jet_pt[idx], weight);

        FillHist1D("incl_genjet30", GenJet_pt[igen], weight);
        FillHist1D("incl_jet30", Jet_pt[idx], weight);
        FillHist2D("incl_jet_genjet", GenJet_pt[igen], Jet_pt[idx], weight);
      }*/
    }

    // Reco flavor //
    int idx = -1;
    double ptmax = -1;
    TString flav = "udsg";

    for (int ireco=0; ireco<nJet; ireco++) {
      if ( fabs(Jet_eta[ireco])>2.4 ) continue;

      if      ( Jet_hadronFlavour[ireco]==5 )                flav = "b";
      else if ( flav != "b" && Jet_hadronFlavour[ireco]==4 ) flav = "c";

      if ( Jet_pt[ireco]>ptmax ) {
        ptmax = Jet_pt[ireco];
        idx = ireco;
      }
    }

    // We have the leading reco jet within eta acceptance
    if (idx != -1) {
      //FillHist1D(flav + "_Ljet", Jet_pt[idx], weight);
      //FillHist1D("incl_Ljet", Jet_pt[idx], weight);
      FillHist1D(flav + "_Ljet", std::min(Jet_pt[idx],max_val), weight);
      FillHist1D("incl_Ljet", std::min(Jet_pt[idx],max_val), weight);

      if (nGenJet > 0) {
        FillHist1D(flav + "_Lgenjet", GenJet_pt[0], weight);
        FillHist1D("incl_Lgenjet", GenJet_pt[0], weight);
      }

      // This reco jet is matched to a gen jet, reco underflow reserved for unmatched jets
      if ( matched0.find(idx) != matched0.end() && Jet_pt[idx]>30 ) {
        int genidx = matched0[idx];

        FillHist1D(flav + "_Ljet_mtchd", std::min(Jet_pt[idx],max_val), weight);
        FillHist1D("incl_Ljet_mtchd", std::min(Jet_pt[idx],max_val), weight);

        FillHist1D(flav + "_Lgenjet_mtchd", GenJet_pt[genidx], weight);
        FillHist1D("incl_Lgenjet_mtchd", GenJet_pt[genidx], weight);

        FillHist2D(flav + "_Ljet_Lgenjet", GenJet_pt[genidx], std::min(Jet_pt[idx],max_val), weight);
        FillHist2D("incl_Ljet_Lgenjet", GenJet_pt[genidx], std::min(Jet_pt[idx],max_val), weight);
      }
      else if ( matched0.find(idx) == matched0.end() && GenJet_pt[0]>30 ) {

        FillHist1D(flav + "_Lgenjet_unmtchd", GenJet_pt[0], weight);
        FillHist1D("incl_Lgenjet_unmtchd", GenJet_pt[0], weight);

        FillHist2D(flav + "_Ljet_Lgenjet", GenJet_pt[0], 0, weight);
        FillHist2D("incl_Ljet_Lgenjet", GenJet_pt[0], 0, weight);
      }

      ////////////////////Vpt
      if (Jet_pt[idx] > 30) {
        float V_pt = (lep0+lep1).Pt();

        //FillHist1D(flav + "_Vpt", V_pt, weight);
        //FillHist1D("incl_Vpt", V_pt, weight);
        FillHist1D(flav + "_Vpt", std::min(V_pt,max_val), weight);
        FillHist1D("incl_Vpt", std::min(V_pt,max_val), weight);

        if ( nGenLep >= 2) {
          //Lep matching
          map<int, int> matched = { {0,-1}, {1,-1} }; //<reco idx, gen idx>
          //this is sloppy but w/e
          vector< TLorentzVector > v_gLep;
          for (int i=0; i<nGenLep; i++) {

            TLorentzVector gLep;
            gLep.SetPtEtaPhiM( GenLep_pt[i], GenLep_eta[i], GenLep_phi[i], GenLep_mass[i] );
            v_gLep.push_back( gLep );

            double rmin = 99.;
            int idx = -1;

            double dR = lep0.DeltaR( gLep );
            if ( dR<rmin && dR<0.2 && matched[0]==-1) {
              rmin = dR;
              idx = 0;
            }

            dR = lep1.DeltaR( gLep );
            if ( dR<rmin && dR<0.2 && matched[1]==-1 ) {
              rmin = dR;
              idx = 1;
            }

            if (idx != -1) matched[idx] = i;
          }

          //both reco leps are matched
          if ( matched[0]!=-1 && matched[1]!=-1 ) {

            float genVpt = (v_gLep[ matched[0] ] + v_gLep[ matched[1] ]).Pt();

            FillHist1D(flav + "_Vpt_mtchd", std::min(V_pt,max_val), weight);
            FillHist1D("incl_Vpt_mtchd", std::min(V_pt,max_val), weight);

            FillHist1D(flav + "_genVpt_mtchd", genVpt, weight);
            FillHist1D("incl_genVpt_mtchd", genVpt, weight);

            FillHist2D(flav + "_Vpt_genVpt", genVpt, std::min(V_pt,max_val), weight);
            FillHist2D("incl_Vpt_genVpt", genVpt, std::min(V_pt,max_val), weight);
          }
          else {

            float genVpt = (v_gLep[0] + v_gLep[1]).Pt();

            FillHist1D(flav + "_genVpt_unmtchd", genVpt, weight);
            FillHist1D("incl_genVpt_unmtchd", genVpt, weight);

            FillHist2D(flav + "_Vpt_genVpt", genVpt, 0, weight);
            FillHist2D("incl_Vpt_genVpt", genVpt, 0, weight);
          }
        }
      }
      ////////////////////End Vpt
    }

  }
  cout << difftime(time(NULL), start) << " s" << endl;

  v_cuts[jetcut].second = m_Histos1D["incl_jet"]->Integral();
  v_cuts[bjet].second = m_Histos1D["b_jet"]->Integral();
  v_cuts[cjet].second = m_Histos1D["c_jet"]->Integral();
  v_cuts[ljet].second = m_Histos1D["udsg_jet"]->Integral();

  v_cuts[jetcut30].second = m_Histos1D["incl_jet30"]->Integral();
  v_cuts[bjetcut30].second = m_Histos1D["b_jet30"]->Integral();
  v_cuts[cjetcut30].second = m_Histos1D["c_jet30"]->Integral();
  v_cuts[ljetcut30].second = m_Histos1D["udsg_jet30"]->Integral();

  //Cutflow Table//
  cout<<"===================================================================================================\n";
  cout<<"                                Cut Flow Table: " + inName( inName.Last('/')+1, inName.Index('.')-inName.Last('/')-1 ) + "\n";
  cout<<"===================================================================================================\n";

  cout<<      "                          |||          Nevent          |||     Efficiency (Relative Efficiency)\n";

  for (int i=0; i<numCuts; i++) {
    if (i == 0)
      cout << Form("%-25s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[0].second) << endl;

    else if (i==jetcut30)
      cout << Form("%-25s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[metcut].second) << endl;
    else if (i>jetcut30)
      cout << Form("%-25s |||       %12.1f       |||       %1.6f (%1.4f)",
                   v_cuts[i].first.data(), v_cuts[i].second, v_cuts[i].second/v_cuts[0].second, v_cuts[i].second/v_cuts[jetcut30].second) << endl;
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

  TFile* outFile = new TFile(outName,"RECREATE");
  outFile->cd();

  for (map<TString, TH1*>::iterator hid = m_Histos1D.begin(); hid != m_Histos1D.end(); hid++) hid->second->Write();
  for (map<TString, TH2*>::iterator hid = m_Histos2D.begin(); hid != m_Histos2D.end(); hid++) hid->second->Write();

  delete outFile;
  outFile = 0;
}

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

    if      (var == "isMC")   { if (line == "true") isMC = true; }
    else if (var == "inName")   inName = line.data();
    else if (var == "outName")  outName = line.data();
    else if (var == "channel")  channel = line;
    else if (var == "xs")       xs = stod( line );
    else if (var == "lumi")     lumi = stod( line );
    else if (var == "posWgt")   posWgt = stod( line );
    else if (var == "negWgt")   negWgt = stod( line );
  }
  file.close();
}
