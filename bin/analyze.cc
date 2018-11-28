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
void FillHist2D(const TString& histName, const TString& value1, const TString& value2, const double& weight);
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
  else if (inName == "amcnlo") {
//    TString path = "root://131.225.204.161:1094//store/user/yokugawa/ZPJ_datasets_2017-03_MC_SUMMER2016/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MC_SUMMER2016_PRIMARY_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext2-v1/170224_012232/";
    TString path = "root://131.225.204.161:1094//store/user/duong/noreplica/VHbbTuples/MC_V25/DY_amcnlo/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MCSUMMER2016_TEST_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-Py8__RunIISummer16MAv2-PUMoriond17_80r2as_2016_TrancheIV_v6_ext2-v1/181125_101328/";
    for (int i=1; i<1000; i++)
      T->Add(path + "0000/tree_" + to_string(i) + ".root");
    for (int i=1000; i<1455; i++)
      T->Add(path + "0001/tree_" + to_string(i) + ".root");
  }
  else T->Add( inName );

  Long64_t nEntries = T->GetEntries();
  cout << nEntries << " Events" << endl;
  cout << "Processing " + inName << endl;
  cout << "Channel: " + channel << endl;
  cout << "weight0: " << weight0  << endl;

  //Histograms//
  vector<TString> flabels = { "b", "c", "udsg", "x" };
  int nFlav = flabels.size();

  TString hname = "flavor_jet_genjet";
  m_Histos2D[hname] = new TH2D(hname,hname,nFlav,0,nFlav,nFlav,0,nFlav);

  hname = "flavor_jet_partjet";
  m_Histos2D[hname] = new TH2D(hname,hname,nFlav,0,nFlav,nFlav,0,nFlav);

  hname = "flavor_jet_mcfmjet";
  m_Histos2D[hname] = new TH2D(hname,hname,nFlav,0,nFlav,nFlav,0,nFlav);

  for ( int i=0; i<nFlav; i++) {
    m_Histos2D["flavor_jet_genjet"]->GetXaxis()->SetBinLabel( i+1, flabels[i] );
    m_Histos2D["flavor_jet_genjet"]->GetYaxis()->SetBinLabel( i+1, flabels[i] );

    m_Histos2D["flavor_jet_partjet"]->GetXaxis()->SetBinLabel( i+1, flabels[i] );
    m_Histos2D["flavor_jet_partjet"]->GetYaxis()->SetBinLabel( i+1, flabels[i] );

    m_Histos2D["flavor_jet_mcfmjet"]->GetXaxis()->SetBinLabel( i+1, flabels[i] );
    m_Histos2D["flavor_jet_mcfmjet"]->GetYaxis()->SetBinLabel( i+1, flabels[i] );

    if ( flabels[i] == "x" ) continue;

    hname = flabels[i] + "_recojet";
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
    hname = flabels[i] + "_vpt";
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

    hname = flabels[i] + "_genjet";
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
    hname = flabels[i] + "_genvpt";
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

    hname = flabels[i] + "_partjet";
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
    hname = flabels[i] + "_partvpt";
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
    hname = flabels[i] + "_mcfmjet";
    m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

    hname = flabels[i] + "_recojet_genjet_mtchd";
    m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
    hname = flabels[i] + "_recojet_partjet_mtchd";
    m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
    hname = flabels[i] + "_recojet_mcfmjet_mtchd";
    m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);

    hname = flabels[i] + "_genjet_partjet_mtchd";
    m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
    hname = flabels[i] + "_genjet_mcfmjet_mtchd";
    m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);

    hname = flabels[i] + "_partjet_mcfmjet_mtchd";
    m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);

    hname = flabels[i] + "_vpt_genvpt_mtchd";
    m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
    hname = flabels[i] + "_vpt_partvpt_mtchd";
    m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
    hname = flabels[i] + "_genvpt_partvpt_mtchd";
    m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
  }

  hname = "incl_recojet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_vpt";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "incl_genjet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_genvpt";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "incl_partjet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_partvpt";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);
  hname = "incl_mcfmjet";
  m_Histos1D[hname] = new TH1D(hname,hname,100,0,500);

  hname = "incl_recojet_genjet_mtchd";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
  hname = "incl_recojet_partjet_mtchd";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
  hname = "incl_recojet_mcfmjet_mtchd";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);

  hname = "incl_genjet_partjet_mtchd";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
  hname = "incl_genjet_mcfmjet_mtchd";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);

  hname = "incl_partjet_mcfmjet_mtchd";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);

  hname = "incl_vpt_genvpt_mtchd";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
  hname = "incl_vpt_partvpt_mtchd";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);
  hname = "incl_genvpt_partvpt_mtchd";
  m_Histos2D[hname] = new TH2D(hname,hname,100,0,500,100,0,500);

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

  ULong64_t evt;
  T->SetBranchAddress("evt", &evt);

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
  float Jet_pt[nJet], Jet_eta[nJet], Jet_phi[nJet], Jet_mass[nJet];

  T->SetBranchAddress("nJet", &nJet);
  T->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);
  T->SetBranchAddress("Jet_pt", Jet_pt);
  T->SetBranchAddress("Jet_eta", Jet_eta);
  T->SetBranchAddress("Jet_phi", Jet_phi);
  T->SetBranchAddress("Jet_mass", Jet_mass);

  int nGenJet = MAXJET;
//  int GenJet_numBHadrons[nGenJet], GenJet_numCHadrons[nGenJet];
  int GenJet_partonFlavour[nGenJet]; //, GenJet_hadronFlavour[nGenJet];
  float GenJet_pt[nGenJet], GenJet_eta[nGenJet], GenJet_phi[nGenJet], GenJet_mass[nGenJet];

  int nGenLep = MAXLEP;
  int GenLep_sourceId[nGenLep];
  float GenLep_pt[nGenLep], GenLep_eta[nGenLep], GenLep_phi[nGenLep], GenLep_mass[nGenLep];

  if (isMC) {
    T->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour);

    T->SetBranchAddress("nGenJet", &nGenJet);
    T->SetBranchAddress("GenJet_pt", GenJet_pt);
    T->SetBranchAddress("GenJet_eta", GenJet_eta);
    T->SetBranchAddress("GenJet_phi", GenJet_phi);
    T->SetBranchAddress("GenJet_mass", GenJet_mass);
    //T->SetBranchAddress("GenJet_numBHadrons", GenJet_numBHadrons);
    //T->SetBranchAddress("GenJet_numCHadrons", GenJet_numCHadrons);
    T->SetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour);
//    T->SetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour);

    T->SetBranchAddress("nGenLep", &nGenLep);
    T->SetBranchAddress("GenLep_pt", GenLep_pt);
    T->SetBranchAddress("GenLep_eta", GenLep_eta);
    T->SetBranchAddress("GenLep_phi", GenLep_phi);
    T->SetBranchAddress("GenLep_mass", GenLep_mass);
    T->SetBranchAddress("GenLep_sourceId", GenLep_sourceId);
  }

  float met_pt, puWeight, genWeight;

  T->SetBranchAddress("met_pt", &met_pt);
  T->SetBranchAddress("puWeight", &puWeight);
  T->SetBranchAddress("genWeight", &genWeight);

  // AOD tree with parton jets //

  TFile* maskFile = TFile::Open("root://131.225.204.161:1094//store/user/cihar29/dy_mcfmjets/dy_mcfmjets.root");
  TTree* maskT = (TTree*) maskFile->Get("T");
  maskT->BuildIndex("event");

  int nPartJet=MAXJET;
  float PartJet_pt[nPartJet], PartJet_eta[nPartJet], PartJet_phi[nPartJet], PartJet_mass[nPartJet];
  int PartJet_flav[nPartJet];

  maskT->SetBranchAddress("npGenJet", &nPartJet);
  maskT->SetBranchAddress("pGenJet_pt", PartJet_pt);
  maskT->SetBranchAddress("pGenJet_eta", PartJet_eta);
  maskT->SetBranchAddress("pGenJet_phi", PartJet_phi);
  maskT->SetBranchAddress("pGenJet_mass", PartJet_mass);
  maskT->SetBranchAddress("pGenJet_flav", PartJet_flav);

  int nmcfmJet=MAXJET;
  float mcfmJet_pt[nmcfmJet], mcfmJet_eta[nmcfmJet], mcfmJet_phi[nmcfmJet], mcfmJet_mass[nmcfmJet];
  int mcfmJet_flav[nmcfmJet];

  maskT->SetBranchAddress("nstat23Jet", &nmcfmJet);
  maskT->SetBranchAddress("stat23Jet_pt", mcfmJet_pt);
  maskT->SetBranchAddress("stat23Jet_eta", mcfmJet_eta);
  maskT->SetBranchAddress("stat23Jet_phi", mcfmJet_phi);
  maskT->SetBranchAddress("stat23Jet_mass", mcfmJet_mass);
  maskT->SetBranchAddress("stat23Jet_flav", mcfmJet_flav);

  int nGen=MAXLEP;
  float gen_pt[nGen];

  maskT->SetBranchAddress("nGen", &nGen);
  maskT->SetBranchAddress("gen_pt", &gen_pt);

  //Loop Over Entries//
  time_t start = time(NULL);

  for (Long64_t n=0; n<nEntries; n++) {
    T->GetEntry(n);

    if (maskT->GetEntryWithIndex(evt) == -1) continue;

    TLorentzVector lep0, lep1;
    weight = weight0*puWeight*TMath::Sign(1,genWeight);

    if (channel == "mm") {
      if (Vtype == 0) {

        if ( !HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v && !HLT_BIT_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v &&
             !HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v && !HLT_BIT_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v ) continue;

        if (vLeptons_pt[0] < 25 || vLeptons_pt[1] < 25) continue;
        if (fabs(vLeptons_eta[0]) > 2.4 || fabs(vLeptons_eta[1]) > 2.4) continue;
        if (vLeptons_pfRelIso04[0] > 0.25 || vLeptons_pfRelIso04[1] > 0.25) continue;
      }
      else continue;
    }
    else if (channel == "ee") {
      if (Vtype == 1) {

        if (!HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v && !HLT_BIT_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v) continue;

        if (vLeptons_pt[0] < 25 || vLeptons_pt[1] < 25) continue;
        if (fabs(vLeptons_eta[0]) > 2.4 || fabs(vLeptons_eta[1]) > 2.4) continue;
        if (1.4442 < fabs(vLeptons_etaSc[0]) && fabs(vLeptons_etaSc[0]) < 1.5660) continue;
        if (1.4442 < fabs(vLeptons_etaSc[1]) && fabs(vLeptons_etaSc[1]) < 1.5660) continue;
        if (vLeptons_pfRelIso03[0] > 0.25 || vLeptons_pfRelIso03[1] > 0.25) continue;
      }
      else continue;
    }
    lep0.SetPtEtaPhiM(vLeptons_pt[0],vLeptons_eta[0],vLeptons_phi[0],vLeptons_mass[0]);
    lep1.SetPtEtaPhiM(vLeptons_pt[1],vLeptons_eta[1],vLeptons_phi[1],vLeptons_mass[1]);

    double V_mass = (lep0+lep1).M();
    if (V_mass<71 || V_mass>111) continue;

    if (met_pt > 40) continue;

    // gen leptons
    vector<TLorentzVector> gleps;
    for ( int i=0; i<nGenLep; i++ ) {
      if ( GenLep_sourceId[i] == 23 ) {
        gleps.push_back( TLorentzVector() );
        gleps.back().SetPtEtaPhiM( GenLep_pt[i],GenLep_eta[i],GenLep_phi[i],GenLep_mass[i] );
      }
    }
    if ( gleps.size() < 2 ) continue;
    float genvpt = (gleps[0]+gleps[1]).Pt();
    float vpt = (lep0+lep1).Pt();

    // fill jet histos //
    // reco jets
    int ridx = -1;
    TString rflav = "x";

    for (int ireco=0; ireco<nJet; ireco++) {
      if ( fabs(Jet_eta[ireco])>2.4 || Jet_pt[ireco]<30 ) continue;

      if      ( rflav!="b" &&               Jet_hadronFlavour[ireco]==5 ) { rflav="b"; ridx=ireco; }
      else if ( rflav!="b" && rflav!="c" && Jet_hadronFlavour[ireco]==4 ) { rflav="c"; ridx=ireco; }
      else if ( rflav=="x" )                                              { rflav="udsg"; ridx=ireco; }
    }

    TLorentzVector recojet;
    if (ridx != -1) {
      FillHist1D(rflav + "_recojet", Jet_pt[ridx], weight);
      FillHist1D("incl_recojet", Jet_pt[ridx], weight);

      FillHist1D(rflav + "_vpt", vpt, weight);
      FillHist1D("incl_vpt", vpt, weight);

      recojet.SetPtEtaPhiM( Jet_pt[ridx], Jet_eta[ridx], Jet_phi[ridx], Jet_mass[ridx] );
    }

    // gen jets
    int gidx = -1;
    TString gflav = "x";

    for (int igen=0; igen<nGenJet; igen++) {
      if ( fabs(GenJet_eta[igen])>2.4 || GenJet_pt[igen]<30 ) continue;

      TLorentzVector genjet;
      genjet.SetPtEtaPhiM( GenJet_pt[igen], GenJet_eta[igen], GenJet_phi[igen], GenJet_mass[igen] );

      // skip genjets which are within 0.2 of genlepton
      if ( gleps[0].DeltaR( genjet ) < 0.2 || gleps[1].DeltaR( genjet ) < 0.2 ) continue;

      if      ( gflav!="b" &&               GenJet_partonFlavour[igen]==5 ) { gflav="b"; gidx=igen; }
      else if ( gflav!="b" && gflav!="c" && GenJet_partonFlavour[igen]==4 ) { gflav="c"; gidx=igen; }
      else if ( gflav=="x" )                                                { gflav="udsg"; gidx=igen; }
    }
    FillHist2D("flavor_jet_genjet", gflav, rflav, weight);

    TLorentzVector genjet;
    if (gidx != -1) {
      FillHist1D(gflav + "_genjet", GenJet_pt[gidx], weight);
      FillHist1D("incl_genjet", GenJet_pt[gidx], weight);

      FillHist1D(gflav + "_genvpt", genvpt, weight);
      FillHist1D("incl_genvpt", genvpt, weight);

      genjet.SetPtEtaPhiM( GenJet_pt[gidx], GenJet_eta[gidx], GenJet_phi[gidx], GenJet_mass[gidx] );

      if (gflav == rflav) {
        FillHist2D(gflav + "_vpt_genvpt_mtchd", genvpt, vpt, weight);
        FillHist2D("incl_vpt_genvpt_mtchd", genvpt, vpt, weight);

        if (genjet.DeltaR( recojet ) < 0.2) {
          FillHist2D(gflav + "_recojet_genjet_mtchd", GenJet_pt[gidx], Jet_pt[ridx], weight);
          FillHist2D("incl_recojet_genjet_mtchd", GenJet_pt[gidx], Jet_pt[ridx], weight);
        }
      }
    }

    // parton jets
    int pidx = -1;
    TString pflav = "x";

    for (int ipart=0; ipart<nPartJet; ipart++) {
      if ( fabs(PartJet_eta[ipart])>2.4 || PartJet_pt[ipart]<30 ) continue;

      TLorentzVector partjet;
      partjet.SetPtEtaPhiM( PartJet_pt[ipart], PartJet_eta[ipart], PartJet_phi[ipart], PartJet_mass[ipart] );

      // skip genjets which are within 0.2 of genlepton
      if ( gleps[0].DeltaR( partjet ) < 0.2 || gleps[1].DeltaR( partjet ) < 0.2 ) continue;

      if      ( pflav!="b" &&               PartJet_flav[ipart]==5 ) { pflav="b"; pidx=ipart; }
      else if ( pflav!="b" && pflav!="c" && PartJet_flav[ipart]==4 ) { pflav="c"; pidx=ipart; }
      else if ( pflav=="x" )                                         { pflav="udsg"; pidx=ipart; }
    }
    FillHist2D("flavor_jet_partjet", pflav, rflav, weight);

    TLorentzVector partjet;
    if (pidx != -1) {
      FillHist1D(pflav + "_partjet", PartJet_pt[pidx], weight);
      FillHist1D("incl_partjet", PartJet_pt[pidx], weight);

      FillHist1D(pflav + "_partvpt", gen_pt[0], weight);
      FillHist1D("incl_partvpt", gen_pt[0], weight);

      partjet.SetPtEtaPhiM( PartJet_pt[pidx], PartJet_eta[pidx], PartJet_phi[pidx], PartJet_mass[pidx] );

      if (pflav == rflav) {
        FillHist2D(pflav + "_vpt_partvpt_mtchd", gen_pt[0], vpt, weight);
        FillHist2D("incl_vpt_partvpt_mtchd", gen_pt[0], vpt, weight);

        if (partjet.DeltaR( recojet ) < 0.2) {
          FillHist2D(pflav + "_recojet_partjet_mtchd", PartJet_pt[pidx], Jet_pt[ridx], weight);
          FillHist2D("incl_recojet_partjet_mtchd", PartJet_pt[pidx], Jet_pt[ridx], weight);
        }
      }
      if (pflav == gflav) {
        FillHist2D(pflav + "_genvpt_partvpt_mtchd", gen_pt[0], genvpt, weight);
        FillHist2D("incl_genvpt_partvpt_mtchd", gen_pt[0], genvpt, weight);

        if (partjet.DeltaR( genjet ) < 0.2) {
          FillHist2D(pflav + "_genjet_partjet_mtchd", PartJet_pt[pidx], GenJet_pt[gidx], weight);
          FillHist2D("incl_genjet_partjet_mtchd", PartJet_pt[pidx], GenJet_pt[gidx], weight);
        }
      }
    }

    // status 23 (mcfm) jets
    int midx = -1;
    TString mflav = "x";

    for (int imcfm=0; imcfm<nmcfmJet; imcfm++) {
      if ( fabs(mcfmJet_eta[imcfm])>2.4 || mcfmJet_pt[imcfm]<30 ) continue;

      TLorentzVector mcfmjet;
      mcfmjet.SetPtEtaPhiM( mcfmJet_pt[imcfm], mcfmJet_eta[imcfm], mcfmJet_phi[imcfm], mcfmJet_mass[imcfm] );

      // skip genjets which are within 0.2 of genlepton
      if ( gleps[0].DeltaR( mcfmjet ) < 0.2 || gleps[1].DeltaR( mcfmjet ) < 0.2 ) continue;

      if      ( mflav!="b" &&               mcfmJet_flav[imcfm]==5 ) { mflav="b"; midx=imcfm; }
      else if ( mflav!="b" && mflav!="c" && mcfmJet_flav[imcfm]==4 ) { mflav="c"; midx=imcfm; }
      else if ( mflav=="x" )                                         { mflav="udsg"; midx=imcfm; }
    }
    FillHist2D("flavor_jet_mcfmjet", mflav, rflav, weight);

    TLorentzVector mcfmjet;
    if (midx != -1) {
      FillHist1D(mflav + "_mcfmjet", mcfmJet_pt[midx], weight);
      FillHist1D("incl_mcfmjet", mcfmJet_pt[midx], weight);

      mcfmjet.SetPtEtaPhiM( mcfmJet_pt[midx], mcfmJet_eta[midx], mcfmJet_phi[midx], mcfmJet_mass[midx] );

      if (mflav == rflav && mcfmjet.DeltaR( recojet ) < 0.2) {
        FillHist2D(mflav + "_recojet_mcfmjet_mtchd", mcfmJet_pt[midx], Jet_pt[ridx], weight);
        FillHist2D("incl_recojet_mcfmjet_mtchd", mcfmJet_pt[midx], Jet_pt[ridx], weight);
      }
      if (mflav == gflav && mcfmjet.DeltaR( genjet ) < 0.2) {
        FillHist2D(mflav + "_genjet_mcfmjet_mtchd", mcfmJet_pt[midx], GenJet_pt[gidx], weight);
        FillHist2D("incl_genjet_mcfmjet_mtchd", mcfmJet_pt[midx], GenJet_pt[gidx], weight);
      }
      if (mflav == pflav && mcfmjet.DeltaR( partjet ) < 0.2) {
        FillHist2D(mflav + "_partjet_mcfmjet_mtchd", mcfmJet_pt[midx], PartJet_pt[pidx], weight);
        FillHist2D("incl_partjet_mcfmjet_mtchd", mcfmJet_pt[midx], PartJet_pt[pidx], weight);
      }
    }

  }
  cout << difftime(time(NULL), start) << " s" << endl;

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

void FillHist2D(const TString& histName, const TString& value1, const TString& value2, const double& weight) {
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
