#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <exception>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TGraphAsymmErrors.h"
#include <stdlib.h>
#include "DesyTauAnalyses/NTupleMaker/interface/json.h"
#include "DesyTauAnalyses/NTupleMaker/interface/PileUp.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "DesyTauAnalyses/NTupleMaker/interface/Jets.h"
#include "DesyTauAnalyses/NTupleMaker/interface/AnalysisMacro.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "HTT-utilities/RecoilCorrections_KIT/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections_KIT/interface/MEtSys.h"
#include "RooWorkspace.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"


int main(int argc, char * argv[]) {

  // **** configuration
  Config cfg(argv[1]);
  string Channel="eltau";

  // kinematic cuts on electrons
  const bool isData = cfg.get<bool>("IsData");
  cout<<"Start Here"<<endl;
  string TauId = cfg.get<string>("TauId");

/*
   bool ApplyTauEnergyScaleUnc     = cfg.get<bool>("ApplyTauEnergyScaleUnc");
   double TauEnergyScaleUnc   = cfg.get<double>("TauEnergyScaleUnc");
   bool ApplyTauCorrectionUncSignPositive     = cfg.get<bool>("ApplyTauCorrectionUncSignPositive");

   bool ApplyMuEnergyScaleUnc     = cfg.get<bool>("ApplyMuEnergyScaleUnc");
   double MuEnergyScaleUnc   = cfg.get<double>("MuEnergyScaleUnc");
   bool ApplyMuonCorrectionUncSignPositive     = cfg.get<bool>("ApplyMuonCorrectionUncSignPositive");

   bool ApplyElEnergyScaleUnc     = cfg.get<bool>("ApplyElEnergyScaleUnc");
   double ElEnergyScaleUncBarrel   = cfg.get<double>("ElEnergyScaleUncBarrel");
   double ElEnergyScaleUncEndcaps   = cfg.get<double>("ElEnergyScaleUncEndcaps");
   bool ApplyElectronCorrectionUncSignPositive     = cfg.get<bool>("ApplyElectronCorrectionUncSignPositive");

   bool ApplyJetEnergyCorrectionUnc     = cfg.get<bool>("ApplyJetEnergyCorrectionUnc");
   bool ApplyJetEnergyCorrectionUncSignPositive     = cfg.get<bool>("ApplyJetEnergyCorrectionUncSignPositive");
*/

  cout<<"cont "<<endl;
  bool ApplyTauEnergyScaleUnc = false;
  bool ApplyTauCorrectionUncSignPositive = false;
  bool ApplyElEnergyScaleUnc = false;
  bool ApplyElectronCorrectionUncSignPositive = false;
  bool ApplyMuEnergyScaleUnc = false;
  bool ApplyMuonCorrectionUncSignPositive = false;
  bool ApplyJetEnergyCorrectionUnc = false;
  bool ApplyJetEnergyCorrectionUncSignPositive = false;



  double TauEnergyScaleUnc   = cfg.get<double>("TauEnergyScaleUnc");
  double MuEnergyScaleUnc   = cfg.get<double>("MuEnergyScaleUnc");
  double ElEnergyScaleUncBarrel   = cfg.get<double>("ElEnergyScaleUncBarrel");
  double ElEnergyScaleUncEndcaps   = cfg.get<double>("ElEnergyScaleUncEndcaps");
// string Systematic=argv[5];
 string BTag_ = "central";
 string BtagCVS = "DeepCSV_2017data.csv" ;
/*
 if (Systematic=="1" || Systematic=="" || isData) Systematic = "Nominal";
  if (string::npos != Systematic.find("TauEnUp")){ ApplyTauEnergyScaleUnc = true; ApplyTauCorrectionUncSignPositive = true;}
	if (string::npos != Systematic.find("TauEnDown")){ ApplyTauEnergyScaleUnc = true; ApplyTauCorrectionUncSignPositive = false;}

	if (string::npos != Systematic.find("ElEnUp")){ ApplyElEnergyScaleUnc = true; ApplyElectronCorrectionUncSignPositive = true;}
	if (string::npos != Systematic.find("ElEnDown")){ ApplyElEnergyScaleUnc = true; ApplyElectronCorrectionUncSignPositive = false;}

	if (string::npos != Systematic.find("MuEnUp")){ ApplyMuEnergyScaleUnc = true; ApplyMuonCorrectionUncSignPositive = true;}
	if (string::npos != Systematic.find("MuEnDown")){ ApplyMuEnergyScaleUnc = true; ApplyMuonCorrectionUncSignPositive = false;}

	if (string::npos != Systematic.find("JetEnUp")){ ApplyJetEnergyCorrectionUnc = true; ApplyJetEnergyCorrectionUncSignPositive = true;}
	if (string::npos != Systematic.find("JetEnDown")){ ApplyJetEnergyCorrectionUnc = true; ApplyJetEnergyCorrectionUncSignPositive = false;}

	if (string::npos != Systematic.find("BTagUp")){ BTag_ = "up";}
	if (string::npos != Systematic.find("BTagDown")){ BTag_ = "down";}
*/

  // kinematics electrons
  const float  ptElectronCut       = cfg.get<float>("ptElectronCuteltau");
  const double etaElectronCut     = cfg.get<double>("etaElectronCuteltau");
  const double dxyElectronCut     = cfg.get<double>("dxyElectronCuteltau");
  const double dzElectronCut      = cfg.get<double>("dzElectronCuteltau");
  const double isoElectronLowCut  = cfg.get<double>("isoElectronLowCuteltau");
  const double isoElectronHighCut  = cfg.get<double>("isoElectronHighCuteltau");
  const bool applyElectronId     = cfg.get<bool>("ApplyElectronId");

  //Tau
  cout<<" Tau :"<<endl;
  const double ptTauCut = cfg.get<double>("ptTauCut");
  const double etaTauCut = cfg.get<double>("etaTauCut");
  const double decayModeFinding    = cfg.get<double>("decayModeFinding");
  const double  againstElectronVLooseMVA6  = cfg.get<double>("againstElectronVLooseMVA6");
  const double  againstMuonTight3  = cfg.get<double>("againstMuonTight3");
  const double  vertexz =  cfg.get<double>("vertexz");
  const double  byCombinedIsolationDeltaBetaCorrRaw3Hits = cfg.get<double>("byCombinedIsolationDeltaBetaCorrRaw3Hits");

  // veto muons
  cout<<" Veto Muons :"<<endl;
  const float ptVetoMuonCut   = cfg.get<float>("ptVetoMuonCut");
  const float etaVetoMuonCut  = cfg.get<float>("etaVetoMuonCut");
  const float dxyVetoMuonCut   = cfg.get<float>("dxyVetoMuonCut");
  const float dzVetoMuonCut   = cfg.get<float>("dzVetoMuonCut");
  const float isoVetoMuonCut   = cfg.get<float>("isoVetoMuonCut");
  const bool applyVetoMuonId     = cfg.get<bool>("ApplyVetoMuonId");

  // dilepton veto
  cout<<" DiEl Veto :"<<endl;
  const float ptDilepElectronCut = cfg.get<float>("ptDilepElectronCuteltau");
  const float etaDilepElectronCut = cfg.get<float>("etaDilepElectronCuteltau");
  const float dxyDilepElectronCut = cfg.get<float>("dxyDilepElectronCuteltau");
  const float dzDilepElectronCut = cfg.get<float>("dzDilepElectronCuteltau");
  const float isoDilepElectronCut = cfg.get<float>("isoDilepElectronCuteltau");
  const float dRDilepVetoCut = cfg.get<float>("dRDilepVetoCuteltau");

  //Veto Electrons
  cout<<" El Veto :"<<endl;
  const float ptVetoElectronCut   = cfg.get<float>("ptVetoElectronCut");
  const float etaVetoElectronCut  = cfg.get<float>("etaVetoElectronCut");
  const float dxyVetoElectronCut  = cfg.get<float>("dxyVetoElectronCut");
  const float dzVetoElectronCut   = cfg.get<float>("dzVetoElectronCut");
  const float isoVetoElectronCut  = cfg.get<float>("isoVetoElectronCut");
  const bool applyVetoElectronId  = cfg.get<bool>("ApplyVetoElectronId");

  //Vertex Cuts
  cout<<" vertex cut :"<<endl;
  const double ndofVertexCut  = cfg.get<double>("NdofVertexCut");
  const double zVertexCut     = cfg.get<double>("ZVertexCut");
  const double dVertexCut     = cfg.get<double>("DVertexCut");
  const string dataBaseDir = cfg.get<string>("DataBaseDir");



  const string Region  = cfg.get<string>("Region");
  const string Sign  = cfg.get<string>("Sign");

  const double leadchargedhadrcand_dz = cfg.get<double>("leadchargedhadrcand_dz");
  const double leadchargedhadrcand_dxy = cfg.get<double>("leadchargedhadrcand_dxy");

  //Kinematic cuts on Jets
  const double etaJetCut   = cfg.get<double>("etaJetCut");
  const double ptJetCut   = cfg.get<double>("ptJetCut");


  // topological cuts
  cout<<" topological cuts"<<endl;
  const double dRleptonsCuteltau   = cfg.get<double>("dRleptonsCuteltau");
  const double deltaRTrigMatch = cfg.get<double>("DRTrigMatch");
  const bool isIsoR03 = cfg.get<bool>("IsIsoR03");
  string TrigLeg = cfg.get<string>("SingleElectronFilterName");
  const float SingleElectronTriggerPtCut = cfg.get<float>("SingleElectronTriggerPtCut");


  //Vertex distributions filenames and histname
  //Attach the year;

  string cmsswBase = (getenv ("CMSSW_BASE"));

  const string _year   = cfg.get<string>("year");
  const string jsonFile = cfg.get<string>("jsonFile");
  cout<<" Start Json Parcing"<<endl; 
  string fullPathToJsonFile = cmsswBase + "/src/DesyTauAnalyses/NTupleMaker/test/json/" + jsonFile;

  const string ElectronIdIsoEffFile = cfg.get<string>("ElectronIdIsoEffFile"+_year);
  const string ElectronTriggerEffFile = cfg.get<string>("ElectrontrigEffFile"+_year);
  

  const string recoilFileName   = cfg.get<string>("RecoilFileName"+_year);
  TString RecoilFileName(recoilFileName);

  const string recoilFileNamePuppi   = cfg.get<string>("RecoilFileNamePuppi"+_year);
  TString RecoilFileNamePuppi(recoilFileNamePuppi);

  const string metSysFileName   = cfg.get<string>("MetSysFileName"+_year);
  TString MetSysFileName(metSysFileName);

  const string metSysFileNamePuppi   = cfg.get<string>("MetSysFileNamePuppi"+_year);
  TString MetSysFileNamePuppi(metSysFileNamePuppi);

  const string zMassPtWeightsFileName   = cfg.get<string>("ZMassPtWeightsFileName"+_year);
  TString ZMassPtWeightsFileName(zMassPtWeightsFileName);

  const string zMassPtWeightsHistName   = cfg.get<string>("ZMassPtWeightsHistName"+_year);
  TString ZMassPtWeightsHistName(zMassPtWeightsHistName);

  //Definition of the Corrector
  kit::MEtSys metSys(MetSysFileName);
  kit::RecoilCorrector recoilMetCorrector(RecoilFileName);
  kit::MEtSys metSysPuppi(MetSysFileNamePuppi);
  kit::RecoilCorrector recoilMetCorrectorPuppi(RecoilFileNamePuppi);


  //  const string TauFakeRateFile = cfg.get<string>("TauFakeRateEff");
  // Run-lumi selector
  std::vector<Period> periods;

  if (isData) { // read the good runs
    std::fstream inputFileStream(fullPathToJsonFile.c_str(), std::ios::in);
    if (inputFileStream.fail() ) {
      std::cout << "Error: cannot find json file " << fullPathToJsonFile << std::endl;
      std::cout << "please check" << std::endl;
      std::cout << "quitting program" << std::endl;
      exit(-1);
    }

    for(std::string s; std::getline(inputFileStream, s); ) {
      periods.push_back(Period());
      std::stringstream ss(s);
      ss >> periods.back();
    }
  }
  TString MainTrigger(TrigLeg);

   const double bTag   = cfg.get<double>("bTag");

   xs=1;fact=1;fact2=1;
   unsigned int RunMin = 9999999;
   unsigned int RunMax = 0;

   ifstream ifs("xsecs");
   string line;
   XSec=1.;
   xsecs=XSec;
   std::vector<unsigned int> allRuns; allRuns.clear();

  bool doThirdLeptVeto=true;
  bool doMuVeto=true;
  bool SUSY = false;
  float SusyMotherMassF;
  float SusyLSPMassF;
  char ff[200];
  sprintf(ff,"%s/%s",argv[3],argv[2]);

  std::string rootFileName(argv[2]);
  std::string NrootFile(argv[4]);
  std::ifstream fileList(ff);
  std::ifstream fileList0(ff);
  std::string ntupleName("makeroottree/AC1B");
  std::string initNtupleName("initroottree/AC1B");

  int nTotalFiles = 0;

  string SaveDir=argv[3];

  TString TStrName(rootFileName+"_"+Region+"_"+Sign);
  datasetName = rootFileName.c_str();

  if (string::npos != datasetName.find("SMS-") || string::npos != datasetName.find("stau") || string::npos != datasetName.find("C1") || string::npos != datasetName.find("Chi") || string::npos != datasetName.find("Stop")) SUSY = true;

  //PU reweighting
  PileUp * PUofficial = new PileUp();

  TString Pileup_file_name_data = cfg.get<string>("Pileup_file_name_data"+_year);
  TString Pileup_file_name_mc = cfg.get<string>("Pileup_file_name_mc"+_year);

  TFile * filePUdistribution_data = new TFile(TString(cmsswBase) +
     "/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+Pileup_file_name_data,"read");
TFile * filePUdistribution_MC = new TFile (TString(cmsswBase) +
    "/src/DesyTauAnalyses/NTupleMaker/data/PileUpDistrib/"+Pileup_file_name_mc, "read");

  TH1D * PU_data = (TH1D *)filePUdistribution_data->Get("pileup");
  TH1D * PU_mc = (TH1D *)filePUdistribution_MC->Get("pileup");
  PUofficial->set_h_data(PU_data);
  PUofficial->set_h_MC(PU_mc);

  //BTags:==============================================================

  BTagCalibration calib("csvv2", cmsswBase+"/src/DesyTauAnalyses/NTupleMaker/data/"+BtagCVS);
  BTagCalibrationReader reader_B(BTagEntry::OP_MEDIUM,BTag_);
  BTagCalibrationReader reader_C(BTagEntry::OP_MEDIUM,BTag_);
  BTagCalibrationReader reader_Light(BTagEntry::OP_MEDIUM,BTag_);

  if (!SUSY){reader_B.load(calib,BTagEntry::FLAV_B,"comb");
  reader_C.load(calib,BTagEntry::FLAV_C,"comb");
  reader_Light.load(calib,BTagEntry::FLAV_UDSG,"incl");}

  if (SUSY){reader_B.load(calib,BTagEntry::FLAV_B,"fastsim");
  reader_C.load(calib,BTagEntry::FLAV_C,"fastsim");
  reader_Light.load(calib,BTagEntry::FLAV_UDSG,"fastsim");}


  float etaBTAG[2] = {0.5,2.1};
  float ptBTAG[5] = {25.,35.,50.,100.,200.};

  for (int iEta=0; iEta<2; ++iEta) {
    for (int iPt=0; iPt<5; ++iPt) {
      float sfB = reader_B.eval_auto_bounds(BTag_,BTagEntry::FLAV_B, etaBTAG[iEta], ptBTAG[iPt]);
      float sfC = reader_C.eval_auto_bounds(BTag_,BTagEntry::FLAV_C, etaBTAG[iEta], ptBTAG[iPt]);
      float sfLight = reader_Light.eval_auto_bounds(BTag_,BTagEntry::FLAV_UDSG, etaBTAG[iEta], ptBTAG[iPt]);
      printf("pT = %3.0f   eta = %3.1f  ->  SFb = %5.3f   SFc = %5.3f   SFl = %5.3f\n",ptBTAG[iPt],etaBTAG[iEta],sfB,sfC,sfLight);
    }
  }

   //TODO add the source of the efficiencies;
  TFile * fileTagging = new TFile(TString(cmsswBase)+TString("/src/DesyTauAnalyses/NTupleMaker/data/tagging_efficiencies_Moriond2017.root"));
  TH1F * tagEff_B = (TH1F*)fileTagging->Get("btag_eff_b");
  TH1F * tagEff_C = (TH1F*)fileTagging->Get("btag_eff_c");
  TH1F * tagEff_Light = (TH1F*)fileTagging->Get("btag_eff_oth");
  TRandom3 rand;

  float MaxBJetPt = 1000.;
  float MaxLJetPt = 1000.;
  float MinLJetPt = 20.;
  float MinBJetPt = 20.;

  TFile * fileZMassPtWeights = new TFile(TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/zptm_weights_2018_kit.root");
  if (fileZMassPtWeights->IsZombie()) {
    exit(-1);
  }
  TH2D * histZMassPtWeights = (TH2D*)fileZMassPtWeights->Get("zptmass_histo");
  if (histZMassPtWeights==NULL) {
  std::cout << " ZMassPT Weights histogram cannot found in file " << TString(cmsswBase) << "/src/DesyTauAnalyses/NTupleMaker/data/" << std::endl;
  exit(-1);
  }

    ScaleFactor * SF_elIdIso = new ScaleFactor();
   SF_elIdIso->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronIdIsoEffFile));

  ScaleFactor * SF_electronTrigger = new ScaleFactor();
  SF_electronTrigger->init_ScaleFactor(TString(cmsswBase)+"/src/"+TString(ElectronTriggerEffFile));

   TString correctionsWorkspaceFileName = TString(cmsswBase)+"/src/DesyTauAnalyses/NTupleMaker/data/htt_scalefactors_legacy_2018.root";
   TFile * correctionWorkSpaceFile = new TFile(correctionsWorkspaceFileName);
   RooWorkspace *correctionWS = (RooWorkspace*)correctionWorkSpaceFile->Get("w");

  TFile * file;
  if (isData) file = new TFile(SaveDir+"/"+TStrName+TString("_DataDriven.root"),"update");
  else file = new TFile(SaveDir+"/"+TStrName+TString(".root"),"update");

  file->mkdir(Channel.c_str());
  file->cd(Channel.c_str());

  int nFiles = 0;
  int nEvents = 0;
  int selEvents = 0;
  bool lumi=false;


  std::string dummy;
  while (fileList0 >> dummy) nTotalFiles++;

  SetupTree();

  if (argv[4] != NULL  && atoi(argv[4])< nTotalFiles) nTotalFiles=atoi(argv[4]);


  cout<<" there are in total "<<nTotalFiles<<" files for this dataset "<<endl;
  for (int iF=0; iF<nTotalFiles; ++iF) {

    std::string filen;
    fileList >> filen;

    TFile * file_ = TFile::Open(TString(filen));


	bool WithInit = true;
	if (SUSY) WithInit=false;

	if (WithInit) cout << "With initroottree"<<endl;
	if (!WithInit) cout << "Without initroottree"<<endl;


    TTree * _inittree = NULL;
    if (!WithInit)  _inittree = (TTree*)file_->Get(TString(ntupleName));
    if (WithInit)  _inittree = (TTree*)file_->Get(TString(initNtupleName));
    if (_inittree==NULL) continue;

    Float_t genweight;
    if (!isData)_inittree->SetBranchAddress("genweight",&genweight);
    Long64_t numberOfEntriesInitTree = _inittree->GetEntries();

    for (Long64_t iEntry=0; iEntry<numberOfEntriesInitTree; ++iEntry) {
        _inittree->GetEntry(iEntry);
        if (isData)histWeightsH->Fill(0.,1.);
    }


    TTree * _tree = NULL;
    _tree = (TTree*)file_->Get(TString(ntupleName));

    if (_tree==NULL) continue;
    Long64_t numberOfEntries = _tree->GetEntries();
    AC1B analysisTree(_tree);
    TLorentzVector BlobA; BlobA.SetXYZT(0,0,0,0);
    TLorentzVector BlobB; BlobB.SetXYZT(0,0,0,0);
    TLorentzVector PairLV; PairLV.SetXYZT(0,0,0,0);

    if (!isData && !WithInit){
		for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) {
			analysisTree.GetEntry(iEntry);
			histWeightsH->Fill(0.,analysisTree.genweight);

	    if (SUSY){
            BlobA.SetXYZT(0,0,0,0);
            BlobB.SetXYZT(0,0,0,0);
            PairLV.SetXYZT(0,0,0,0);
	    //Vector clear	
	    for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {

	        TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
					      analysisTree.genparticles_py[igen],
					      analysisTree.genparticles_pz[igen],
					      analysisTree.genparticles_e[igen]);

	  if (string::npos != datasetName.find("C1N2") || string::npos != datasetName.find("Chi")) {
			  if (abs(analysisTree.genparticles_pdgid[igen])==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==1000023 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }

	  if (string::npos != datasetName.find("C1C1")) {
			  if (analysisTree.genparticles_pdgid[igen]==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==-1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }

	  if (string::npos != datasetName.find("left") || string::npos != datasetName.find("max") ) {
			  if (analysisTree.genparticles_pdgid[igen]==1000015 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==-1000015 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }
	  if (string::npos != datasetName.find("right")) {
			  if (analysisTree.genparticles_pdgid[igen]==2000015 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==-2000015 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }

	}

	if (BlobA.M()>0 && BlobB.M()>0) PairLV = BlobA+BlobB;
	if (PairLV.Pt()>0) histPt->Fill(PairLV.Pt());

			}

		}
	}



  	float genweights=1.;

    if(!isData && WithInit)
      {
	TTree *genweightsTree = (TTree*)file_->Get("initroottree/AC1B");
	genweightsTree->SetBranchAddress("genweight",&genweights);


	Long64_t numberOfEntriesInit = genweightsTree->GetEntries();
	for (Long64_t iEntryInit=0; iEntryInit<numberOfEntriesInit; ++iEntryInit) {

	        genweightsTree->GetEntry(iEntryInit);
	        histWeightsH->Fill(0.,genweights);
            }
    }


    for (Long64_t iEntry=0; iEntry<numberOfEntries; ++iEntry) {
      Float_t weight = 1.;
      Float_t puweight = 1.;
      Float_t topptweight = 1.;
      analysisTree.GetEntry(iEntry);
      nEvents++;

      if (nEvents%100==0) cout << "      processed " << nEvents << " events" << endl;

      if (fabs(analysisTree.primvertex_z)>zVertexCut) continue;
      if (analysisTree.primvertex_ndof<ndofVertexCut) continue;
      double dVertex = (analysisTree.primvertex_x*analysisTree.primvertex_x+
			            analysisTree.primvertex_y*analysisTree.primvertex_y);
      if (dVertex>dVertexCut) continue;
      if (analysisTree.primvertex_count<2) continue;

      bool lumi=false;
      bool CutBasedTauId = false;

      if (ApplyTauEnergyScaleUnc && !isData)
		{
			double ApplyTauCorrectionUncSign=1;
			if (!ApplyTauCorrectionUncSignPositive) ApplyTauCorrectionUncSign = -1;

		    for (unsigned int it = 0; it<analysisTree.tau_count; ++it)
			{

                analysisTree.tau_pt[it] *= (1 +ApplyTauCorrectionUncSign*TauEnergyScaleUnc);
                analysisTree.tau_px[it] *= (1 +ApplyTauCorrectionUncSign*TauEnergyScaleUnc);
                analysisTree.tau_py[it] *= (1 +ApplyTauCorrectionUncSign*TauEnergyScaleUnc);
                analysisTree.tau_pz[it] *= (1 +ApplyTauCorrectionUncSign*TauEnergyScaleUnc);
                analysisTree.tau_e[it] *= (1 +ApplyTauCorrectionUncSign*TauEnergyScaleUnc);

			//analysisTree.pfmet_ex = analysisTree.pfmet_ex+((analysisTree.tau_px[it]/TauEnergyScaleUnc)-analysisTree.tau_px[it]);
			//analysisTree.pfmet_ey = analysisTree.pfmet_ey+((analysisTree.tau_py[it]/TauEnergyScaleUnc)-analysisTree.tau_py[it]);
			}

		}

	if (ApplyJetEnergyCorrectionUnc && !isData)
		{
		double ApplyJetEnergyCorrectionUncSign=1;
		for (unsigned int it = 0; it<analysisTree.pfjet_count; ++it)
			{
			if (!ApplyJetEnergyCorrectionUncSignPositive) ApplyJetEnergyCorrectionUncSign = -1;
			analysisTree.pfjet_pt[it] *=(1 + ApplyJetEnergyCorrectionUncSign*analysisTree.pfjet_jecUncertainty[it]);
			analysisTree.pfjet_px[it] *=(1 + ApplyJetEnergyCorrectionUncSign*analysisTree.pfjet_jecUncertainty[it]);
			analysisTree.pfjet_py[it] *=(1 + ApplyJetEnergyCorrectionUncSign*analysisTree.pfjet_jecUncertainty[it]);
			analysisTree.pfjet_pz[it] *=(1 + ApplyJetEnergyCorrectionUncSign*analysisTree.pfjet_jecUncertainty[it]);
			analysisTree.pfjet_e[it] *=(1 + ApplyJetEnergyCorrectionUncSign*analysisTree.pfjet_jecUncertainty[it]);
			}

		}

	if (ApplyElEnergyScaleUnc && !isData)
		{
		double ElEnergyScaleUnc=1;
			double ApplyElectronCorrectionUncSign=1;
			if (!ApplyElectronCorrectionUncSignPositive) ApplyElectronCorrectionUncSign = -1;

		for (unsigned int it = 0; it<analysisTree.electron_count; ++it)
			{

			if (analysisTree.electron_eta[it] < 1.48)  ElEnergyScaleUnc = ElEnergyScaleUncBarrel;
			if (analysisTree.electron_eta[it] > 1.48)  ElEnergyScaleUnc = ElEnergyScaleUncEndcaps;

			analysisTree.electron_pt[it] *=(1 + ApplyElectronCorrectionUncSign*ElEnergyScaleUnc);
			analysisTree.electron_px[it] *=(1 + ApplyElectronCorrectionUncSign*ElEnergyScaleUnc);
			analysisTree.electron_py[it] *=(1 + ApplyElectronCorrectionUncSign*ElEnergyScaleUnc);
			analysisTree.electron_pz[it] *=(1 + ApplyElectronCorrectionUncSign*ElEnergyScaleUnc);

			}

		}

	if (ApplyMuEnergyScaleUnc && !isData)
		{
			double ApplyMuonCorrectionUncSign=1;
			if (!ApplyMuonCorrectionUncSignPositive) ApplyMuonCorrectionUncSign = -1;

		for (unsigned int it = 0; it<analysisTree.muon_count; ++it)
			{

			analysisTree.muon_pt[it] *= (1 + ApplyMuonCorrectionUncSign*MuEnergyScaleUnc);
			analysisTree.muon_px[it] *= (1 + ApplyMuonCorrectionUncSign*MuEnergyScaleUnc);
			analysisTree.muon_py[it] *= (1 + ApplyMuonCorrectionUncSign*MuEnergyScaleUnc);
			analysisTree.muon_pz[it] *= (1 + ApplyMuonCorrectionUncSign*MuEnergyScaleUnc);
			}

		}

      //Systematic study end

      float topPt = 0;
      float antitopPt = 0;
      LSF_weight = 1.;
      TFR_weight = 1.;
      top_weight = 1.;
      all_weight = 1.;
      pu_weight = 1.;
      gen_weight = 1.;
      trig_weight = 1.;
      //needed for Recoil
      bool isW = false;
      bool isDY = false;
      bool isZTT = false;
      bool isZMM = false;
      bool isZEE = false;
      bool isTOP = false;
      if (!isData &&  string::npos != filen.find("JetsToLNu") ) isW=true;
      if (!isData &&  string::npos != filen.find("TTWJetsToLNu") ) isW=false;
      if (!isData &&  string::npos != filen.find("JetsToLL_M") )  isDY=true;
      if (!isData &&  string::npos != filen.find("TTTo") ) isTOP=true;

      float nuPx = 0;
      float nuPy = 0;
      float nuPz = 0;
      float nuPt = 0;
      float nuPhi = 0;

      float nuPx_msv = 0;
      float nuPy_msv = 0;
      float nuPz_msv = 0;
      float nuPt_msv = 0;
      float nuPhi_msv = 0;

      float lepPx = 0;
      float lepPy = 0;
      float lepPz = 0;
      float bosonPx = 0;
      float bosonPy = 0;
      float bosonPz = 0;
      float bosonPt = 0;
      float bosonEta = 0;
      float bosonMass = -1;

      bool isZfound = false;
      bool isWfound = false;
      bool isHfound = false;
      bool isGSfound = false;

      std::vector<TLorentzVector> promptTausFirstCopy; promptTausFirstCopy.clear();
      std::vector<TLorentzVector> promptTausLastCopy;  promptTausLastCopy.clear();
      std::vector<TLorentzVector> promptElectrons; promptElectrons.clear();
      std::vector<TLorentzVector> promptMuons; promptMuons.clear();
      std::vector<TLorentzVector> promptNeutrinos; promptNeutrinos.clear();
      std::vector<TLorentzVector> tauNeutrinos; tauNeutrinos.clear();

      TLorentzVector promptTausLV; promptTausLV.SetXYZT(0.001,0.001,0,0);
      TLorentzVector promptVisTausLV; promptVisTausLV.SetXYZT(0.001,0.001,0,0);
      TLorentzVector zBosonLV; zBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector wBosonLV; wBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector hBosonLV; hBosonLV.SetXYZT(0,0,0,0);
      TLorentzVector promptElectronsLV; promptElectronsLV.SetXYZT(0.001,0.001,0,0);
      TLorentzVector promptMuonsLV; promptMuonsLV.SetXYZT(0.001,0.001,0,0);
      TLorentzVector promptNeutrinosLV;  promptNeutrinosLV.SetXYZT(0,0,0,0);
      TLorentzVector tauNeutrinosLV;  tauNeutrinosLV.SetXYZT(0,0,0,0);
      TLorentzVector wDecayProductsLV; wDecayProductsLV.SetXYZT(0,0,0,0);
      TLorentzVector fullVLV; fullVLV.SetXYZT(0,0,0,0);
      TLorentzVector visVLV; visVLV.SetXYZT(0,0,0,0);

    if (!isData) {

	  for (unsigned int igentau=0; igentau < analysisTree.gentau_count; ++igentau) {

	  TLorentzVector tauLV; tauLV.SetXYZT(analysisTree.gentau_px[igentau],
					      analysisTree.gentau_py[igentau],
					      analysisTree.gentau_pz[igentau],
					      analysisTree.gentau_e[igentau]);

	  TLorentzVector tauVisLV; tauVisLV.SetXYZT(analysisTree.gentau_visible_px[igentau],
						    analysisTree.gentau_visible_py[igentau],
						    analysisTree.gentau_visible_pz[igentau],
						    analysisTree.gentau_visible_e[igentau]);
	  if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isFirstCopy[igentau]) {
	    promptTausFirstCopy.push_back(tauLV);
	    promptTausLV += tauLV;
	    wDecayProductsLV += tauLV;
	  }
	  if (analysisTree.gentau_isPrompt[igentau]&&analysisTree.gentau_isLastCopy[igentau]) {
	    promptTausLastCopy.push_back(tauVisLV);
	    promptVisTausLV += tauVisLV;
	  }

	}
    BlobA.SetXYZT(0,0,0,0);
    BlobB.SetXYZT(0,0,0,0);
    PairLV.SetXYZT(0,0,0,0);

	for (unsigned int igen=0; igen < analysisTree.genparticles_count; ++igen) {


	  TLorentzVector genLV; genLV.SetXYZT(analysisTree.genparticles_px[igen],
					      analysisTree.genparticles_py[igen],
					      analysisTree.genparticles_pz[igen],
					      analysisTree.genparticles_e[igen]);

	if (SUSY){

	  if (string::npos != datasetName.find("C1N2") || string::npos != datasetName.find("Chi")) {
			  if (abs(analysisTree.genparticles_pdgid[igen])==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==1000023 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }

	  if (string::npos != datasetName.find("C1C1")) {
			  if (analysisTree.genparticles_pdgid[igen]==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==-1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }

	  if (string::npos != datasetName.find("left") || string::npos != datasetName.find("max") ) {
			  if (analysisTree.genparticles_pdgid[igen]==1000015 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==-1000015 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }
	  if (string::npos != datasetName.find("right")) {
			  if (analysisTree.genparticles_pdgid[igen]==2000015 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (analysisTree.genparticles_pdgid[igen]==-2000015 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }
	  if (string::npos != datasetName.find("Stop")) {
			  if (abs(analysisTree.genparticles_pdgid[igen])==1000006 && abs(analysisTree.genparticles_status[igen])==62) BlobA = genLV;
	  		  if (abs(analysisTree.genparticles_pdgid[igen])==1000024 && abs(analysisTree.genparticles_status[igen])==62) BlobB = genLV;
			  }

	if (BlobA.M()>0 && BlobB.M()>0) PairLV = BlobA+BlobB;


	PtSystem = PairLV.Pt();
	}

	  if (analysisTree.genparticles_pdgid[igen]==22 && analysisTree.genparticles_status[igen]==44)
	    isGSfound = true;

	  if (analysisTree.genparticles_pdgid[igen]==23) { 
	    isZfound = true;
	    zBosonLV = genLV;
	  }
	  if (analysisTree.genparticles_pdgid[igen]==25||
	      analysisTree.genparticles_pdgid[igen]==35||
	      analysisTree.genparticles_pdgid[igen]==36) { 
	    isHfound = true;
	    hBosonLV = genLV;
	  }
	  if (abs(analysisTree.genparticles_pdgid[igen])==24) { 
	    isWfound = true;
	    wBosonLV = genLV;
	  }

	  if (fabs(analysisTree.genparticles_pdgid[igen])==11) { 
	    if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
	      promptElectrons.push_back(genLV);
	      promptElectronsLV += genLV;
	      wDecayProductsLV += genLV;
	    }
	  }
	  
	  if (fabs(analysisTree.genparticles_pdgid[igen])==13) { 
	    if (analysisTree.genparticles_fromHardProcess[igen]&&analysisTree.genparticles_status[igen]==1) {
	      promptMuons.push_back(genLV);
	      promptMuonsLV += genLV;
	      wDecayProductsLV += genLV;
	    }
	  }

	  if (fabs(analysisTree.genparticles_pdgid[igen])==12||
	      fabs(analysisTree.genparticles_pdgid[igen])==14||
	      fabs(analysisTree.genparticles_pdgid[igen])==16)  {
	    if ((analysisTree.genparticles_fromHardProcess[igen]||analysisTree.genparticles_isPrompt[igen])&&
		!analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
		analysisTree.genparticles_status[igen]==1) {
	      promptNeutrinos.push_back(genLV);
	      promptNeutrinosLV += genLV;
	      wDecayProductsLV += genLV;
	    }
	    if (analysisTree.genparticles_isDirectHardProcessTauDecayProduct[igen]&&
		analysisTree.genparticles_status[igen]==1) {
	      tauNeutrinos.push_back(genLV);
	      tauNeutrinosLV += genLV;
	    }
	  }

	}

	isDYTT=false;
	isDYLL=false;
	isDYEE=false;
	isDYMM=false;
    //cout<<"Matching ISR Jets"<<endl;

	if (isDY) {

	  if (promptTausFirstCopy.size()==2) {
	    isZTT = true; isZMM = false; isZEE = false;
	    isDYTT=true;
	    bosonPx = promptTausLV.Px(); bosonPy = promptTausLV.Py(); bosonPz = promptTausLV.Pz();
	    bosonMass = promptTausLV.M();
	    bosonEta  = promptTausLV.Eta();
	    lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
	    //mtBoson_gen = mT(promptTausFirstCopy[0],promptTausFirstCopy[1]);
	  }
	  else if (promptMuons.size()==2) {
	    isZTT = false; isZMM = true; isZEE = false;
	    isDYMM=true;
	    bosonPx = promptMuonsLV.Px(); bosonPy = promptMuonsLV.Py(); bosonPz = promptMuonsLV.Pz();
	    bosonMass = promptMuonsLV.M();
	    bosonEta = promptMuonsLV.Eta();
	    lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
	    //mtBoson_gen = mT(promptMuons[0],promptMuons[1]);
	  }
	  else {
	    isZTT = false; isZMM = false; isZEE = true;
	    isDYEE=true;
	    bosonPx = promptElectronsLV.Px(); bosonPy = promptElectronsLV.Py(); bosonPz = promptElectronsLV.Pz();
	    bosonMass = promptElectronsLV.M();
	    bosonEta = promptElectronsLV.Eta();
	    lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
	    //if (promptElectrons.size()==2)
	    //  mtBoson_gen = mT(promptElectrons[0],promptElectrons[1]);
	  }
	  nuPx = tauNeutrinosLV.Px(); nuPy = tauNeutrinosLV.Py(); nuPz = tauNeutrinosLV.Pz();
	}

	else if (isW) {
	  bosonPx = wDecayProductsLV.Px(); bosonPy = wDecayProductsLV.Py(); bosonPz = wDecayProductsLV.Pz();
	  bosonMass = wDecayProductsLV.M();
	  if (promptTausLastCopy.size()==1) {
	    lepPx = promptVisTausLV.Px(); lepPy = promptVisTausLV.Py(); lepPz = promptVisTausLV.Pz();
	  }
	  else if (promptMuons.size()==1) {
	    lepPx = promptMuonsLV.Px(); lepPy = promptMuonsLV.Py(); lepPz = promptMuonsLV.Pz();
	  }
	  else {
	    lepPx = promptElectronsLV.Px(); lepPy = promptElectronsLV.Py(); lepPz = promptElectronsLV.Pz();
	  }
	  nuPx = promptNeutrinosLV.Px(); nuPy = promptNeutrinosLV.Py(); nuPz = promptNeutrinosLV.Pz();
	}
	else {
	  TLorentzVector bosonLV = promptTausLV + promptMuonsLV + promptElectronsLV + promptNeutrinosLV;
	  bosonPx = bosonLV.Px(); bosonPy = bosonLV.Py(); bosonPz = bosonLV.Pz();
	  TLorentzVector lepLV = promptVisTausLV + promptMuonsLV + promptElectronsLV;
	  lepPx = lepLV.Px(); lepPy = lepLV.Py(); lepPz = lepLV.Pz();
	  nuPx = promptNeutrinosLV.Px(); nuPy = promptNeutrinosLV.Py(); nuPz = promptNeutrinosLV.Pz();
	}

	nuPt = TMath::Sqrt(nuPx*nuPx+nuPy*nuPy);
	nuPhi = TMath::ATan2(nuPy,nuPx);

	bosonPt = TMath::Sqrt(bosonPx*bosonPx+bosonPy*bosonPy);

      }

	if (isDY) {
	  // applying Z pt mass weights
	  zptmassweight = 1;
	  if (bosonMass>50.0) {
	    float bosonMassX = bosonMass;
	    float bosonPtX = bosonPt;
	    if (bosonMassX>1000.) bosonMassX = 1000.;
	    if (bosonPtX<1.)      bosonPtX = 1.;
	    if (bosonPtX>1000.)   bosonPtX = 1000.;
	    zptmassweight = histZMassPtWeights->GetBinContent(histZMassPtWeights->GetXaxis()->FindBin(bosonMassX),
						histZMassPtWeights->GetYaxis()->FindBin(bosonPtX));
	  }
	}

    if (!isData && ( string::npos != filen.find("TTTo")  || string::npos != filen.find("TT_TuneCUETP8M2T4_13TeV-powheg-pythia8") || string::npos != filen.find("TT_TuneCUETP8M1_13TeV-powheg-pythia8")) )
	{
	  for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {
	    // 		cout<< "  info = " <<  int(analysisTree.genparticles_count) <<"  "<<int(analysisTree.genparticles_pdgid[igen])<<endl;

	    if (analysisTree.genparticles_pdgid[igen]==6)
	      topPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				  analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);
	    if (analysisTree.genparticles_pdgid[igen]==-6)
	      antitopPt = TMath::Sqrt(analysisTree.genparticles_px[igen]*analysisTree.genparticles_px[igen]+
				      analysisTree.genparticles_py[igen]*analysisTree.genparticles_py[igen]);


	  }

	  if (topPt>0.&&antitopPt>0.) {
	    topptweight = topPtWeight(topPt,antitopPt);
	    weight *= topptweight;
	    top_weight = topptweight;
	     // cout<<"  "<<topPt<<"  "<<antitopPt<<"  "<<topptweight<<endl;
	  }


      histTopPt->Fill(0.,topptweight*analysisTree.genweight);
      histTopPtSq->Fill(0.,topptweight*topptweight*analysisTree.genweight);

	}
	  if (!isData ) {
	    weight *= analysisTree.genweight;
	    gen_weight *=analysisTree.genweight;

	    double cLower, cUpper;
	    vector <double> ScalesV; ScalesV.clear();
		//	ScalesV.push_back(wScale0);
	    ScalesV.push_back(analysisTree.weightScale1);
	    ScalesV.push_back(analysisTree.weightScale2);
	    ScalesV.push_back(analysisTree.weightScale3);
	    ScalesV.push_back(analysisTree.weightScale4);
	    ScalesV.push_back(analysisTree.weightScale5);
	    ScalesV.push_back(analysisTree.weightScale6);
	    ScalesV.push_back(analysisTree.weightScale7);
	    ScalesV.push_back(analysisTree.weightScale8);

	    cLower = *min_element(ScalesV.begin(), ScalesV.end());
		cUpper = *max_element(ScalesV.begin(), ScalesV.end());
	    histWeightsScalesUp->Fill(0.,analysisTree.genweight*cUpper);
	    histWeightsScalesDown->Fill(0.,analysisTree.genweight*cLower);
	    histWeightsPDFUp->Fill(0.,analysisTree.genweight*analysisTree.weightPDFup);
	    histWeightsPDFDown->Fill(0.,analysisTree.genweight*analysisTree.weightPDFdown);

		lumi=true;
		wScale0 = analysisTree.weightScale0;
		wScale1 = analysisTree.weightScale1;
		wScale2 = analysisTree.weightScale2;
		wScale3 = analysisTree.weightScale3;
		wScale4 = analysisTree.weightScale4;
		wScale5 = analysisTree.weightScale5;
		wScale6 = analysisTree.weightScale6;
		wScale7 = analysisTree.weightScale7;
		wScale8 = analysisTree.weightScale8;
		wPDFUp = analysisTree.weightPDFup;
		wPDFDown = analysisTree.weightPDFdown;

	  }

    if (isData)  {
		XSec = 1.;
		histRuns->Fill(analysisTree.event_run);
		int n=analysisTree.event_run;
		int lum = analysisTree.event_luminosityblock;

		std::string num = std::to_string(n);
		std::string lnum = std::to_string(lum);
		for(const auto& a : periods)
	  {

	    if ( num.c_str() ==  a.name ) {

	      for(auto b = a.ranges.begin(); b != std::prev(a.ranges.end()); ++b) {
		  cout<<b->lower<<"  "<<b->bigger<<endl;
		if (lum  >= b->lower && lum <= b->bigger ) lumi = true;
	      }
	      auto last = std::prev(a.ranges.end());
	          std::cout <<"min "<< last->lower << "- max last " << last->bigger << std::endl;
	      if (  (lum >=last->lower && lum <= last->bigger )) lumi=true;


	    }

	  }
	    if (!lumi) continue;

    }
    if (analysisTree.event_run<RunMin) RunMin = analysisTree.event_run;
    if (analysisTree.event_run>RunMax) RunMax = analysisTree.event_run;
    bool isNewRun = true;
    if (allRuns.size()>0) {
	  for (unsigned int iR=0; iR<allRuns.size(); ++iR) {
	    if (analysisTree.event_run==allRuns.at(iR)) {
	    isNewRun = false;
	    break;
	    }
	  }
      }

      if (isNewRun) allRuns.push_back(analysisTree.event_run);
      if (!lumi) continue;
	  std::vector<TString> metFlags; metFlags.clear();

        metFlags.push_back("Flag_goodVertices");
	  metFlags.push_back("Flag_globalSuperTightHalo2016Filter");
	  metFlags.push_back("Flag_HBHENoiseFilter");
	  metFlags.push_back("Flag_HBHENoiseIsoFilter");
	  metFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
	  metFlags.push_back("Flag_BadPFMuonFilter");
	  if (isData){
	  metFlags.push_back("Flag_eeBadScFilter");
	  metFlags.push_back("ecalBadCalibReducedMINIAODFilter");
	  }

	 bool METflag = metFiltersPasses2(analysisTree, metFlags);
	 met_flag = METflag;

      //PU weights:
      if (!isData) {
            puweight = float(PUofficial->get_PUweight(double(analysisTree.numtruepileupinteractions)));
            weight *=puweight;
            pu_weight = puweight;
      }

      bool trigAccept = false;
      unsigned int nMainTrigger = 18;
      bool isMainTrigger = false;

      if (!isData){
        int nfilters = analysisTree.run_hltnames->size();//65
        for (int i=0; i<18; ++i) {
        try{
            TString HLTFilter(analysisTree.run_hltfilters->at(i));
            if (HLTFilter==MainTrigger) {
              nMainTrigger = i;
              isMainTrigger = true;
	      break;
	 }
        }
        catch (exception& e)
          {
            cout << "Standard exception: " << e.what() << endl;
          }

      }
      analysisTree.run_hltfilters->clear();
	}
    if (isData) isMainTrigger = true;

    if (!isMainTrigger) {
            return(-1);
    }
    unsigned int nBTagDiscriminant = -1;
    for (unsigned int iBTag=0; iBTag < analysisTree.run_btagdiscriminators->size(); ++iBTag) {
       TString discr(analysisTree.run_btagdiscriminators->at(iBTag));
       if (discr.Contains("pfCombinedInclusiveSecondaryVertexV2BJetTags"))
          nBTagDiscriminant = iBTag;
    }

    vector<int> electrons; electrons.clear();
      for (unsigned int ie = 0; ie<analysisTree.electron_count; ++ie) {
	if (analysisTree.electron_pt[ie]<SingleElectronTriggerPtCut) continue;
	if (fabs(analysisTree.electron_eta[ie])>etaElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[ie])>dxyElectronCut) continue;
	if (fabs(analysisTree.electron_dz[ie])>dzElectronCut) continue;
	bool electronMvaId = analysisTree.electron_mva_wp90_noIso_Fall17_v1[ie];
	if (applyElectronId && !electronMvaId) continue;
	if (applyElectronId && !analysisTree.electron_pass_conversion[ie]) continue;
	if (applyElectronId && analysisTree.electron_nmissinginnerhits[ie]>1) continue;
	if (fabs(analysisTree.electron_charge[ie]) !=1) continue;
	 electrons.push_back((int)ie);

      }
      if (electrons.size()==0 ) continue;

      vector<int> taus; taus.clear();
      for (unsigned int it = 0; it<analysisTree.tau_count; ++it) {

	if (analysisTree.tau_pt[it] < ptTauCut ) continue; 
	if (fabs(analysisTree.tau_eta[it])> etaTauCut) continue;
	if (analysisTree.tau_decayModeFinding[it]<decayModeFinding) continue;
	if ( fabs(analysisTree.tau_leadchargedhadrcand_dz[it])> leadchargedhadrcand_dz) continue;
        if ( fabs(analysisTree.tau_charge[it]) != 1 ) continue;
	  taus.push_back((int)it);

	}

      if (taus.size()==0) continue;

      cout<<"  electrons  "<<electrons.size()<< "  "<<taus.size()<<endl;

      int tau_index = -1;
      int el_index = -1;
      int mu_index = -1;

      float isoElecMin  = 1e+10;
      float isoTauMin = -1; 
      float isoTau =1.;
      if (CutBasedTauId) isoTauMin = -1;
      if (!CutBasedTauId) isoTauMin = -10;
      float ptEl = 0;
      float ptTau = 0;
      for (unsigned int im=0; im<electrons.size(); ++im) {
	bool isLegMatch = false;
	unsigned int eIndex  = electrons.at(im);
	float neutralHadIsoElec = analysisTree.electron_neutralHadIso[eIndex];
	float photonIsoElec = analysisTree.electron_photonIso[eIndex];
	float chargedHadIsoElec = analysisTree.electron_chargedHadIso[eIndex];
	float puIsoElec = analysisTree.electron_puIso[eIndex];
	if (isIsoR03) {
	  neutralHadIsoElec = analysisTree.electron_r03_sumNeutralHadronEt[eIndex];
	  photonIsoElec = analysisTree.electron_r03_sumPhotonEt[eIndex];
	  chargedHadIsoElec = analysisTree.electron_r03_sumChargedHadronPt[eIndex];
	  puIsoElec = analysisTree.electron_r03_sumPUPt[eIndex];
	}
	double neutralIsoElecN = neutralHadIsoElec + photonIsoElec - 0.5*puIsoElec;
	double neutralIsoElec = max(double(0),neutralIsoElecN); 
	float absIsoElec = chargedHadIsoElec + neutralIsoElec;
	float relIsoElec = absIsoElec/float(analysisTree.electron_pt[eIndex]);


	if (!SUSY)
	{ for (unsigned int iT=0; iT<analysisTree.trigobject_count; ++iT) {
	  if (analysisTree.trigobject_filters[iT][nMainTrigger]
	      &&analysisTree.electron_pt[eIndex]>SingleElectronTriggerPtCut &&
	      analysisTree.trigobject_pt[iT]>SingleElectronTriggerPtCut) { // IsoElec Leg
	    float dRtrig = deltaR(analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex],
				  analysisTree.trigobject_eta[iT],analysisTree.trigobject_phi[iT]);
	    if (dRtrig<deltaRTrigMatch) 
	      isLegMatch = true;
	    
	  }
	 }
	}
      cout<<"isLegMatch: "<<isLegMatch<<endl;
      if (SUSY && analysisTree.electron_pt[eIndex]>SingleElectronTriggerPtCut) isLegMatch = true;

      if (!isLegMatch) continue;

	for (unsigned int it=0; it<taus.size(); ++it) {

	  unsigned int tIndex = taus.at(it);

	  float dR = deltaR(analysisTree.tau_eta[tIndex],analysisTree.tau_phi[tIndex],
			    analysisTree.electron_eta[eIndex],analysisTree.electron_phi[eIndex]);

	  if (dR<dRleptonsCuteltau) continue;

	

	if (CutBasedTauId){isoTau = analysisTree.tau_chargedIsoPtSum[tIndex];

	  if ( (int) eIndex!= (int)el_index) {
	    if (relIsoElec==isoElecMin) {
	      if (analysisTree.electron_pt[eIndex]>ptEl) {
		isoElecMin  = relIsoElec;
		ptEl = analysisTree.electron_pt[eIndex];
		el_index = (int)eIndex;
		isoTauMin = isoTau;
		ptTau = analysisTree.tau_pt[tIndex];
		tau_index = (int)tIndex;
	      }
	    }
	    else if (relIsoElec<isoElecMin) {
	      isoElecMin  = relIsoElec;
	      ptEl = analysisTree.electron_pt[eIndex];
	      el_index = (int)eIndex;
	      isoTauMin = isoTau;
	      ptTau = analysisTree.tau_pt[tIndex];
	      tau_index = (int)tIndex;
	    }
	  }
	  else {
	    if (isoTau==isoTauMin) {
	      if (analysisTree.tau_pt[tIndex]>ptTau) {
		ptTau = analysisTree.tau_pt[tIndex];
		isoTauMin = isoTau;
		tau_index = (int)tIndex;
	      }
	    }
	    else if (isoTau>isoTauMin) {
	      ptTau = analysisTree.tau_pt[tIndex];
	      isoTauMin = isoTau;
	      tau_index = (int)tIndex;
	    }
	  }
	  
	}

	if (!CutBasedTauId){
	isoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tIndex];
        cout<<"isoTau : "<<isoTau<<endl;
	  if ( (int)eIndex!= (int)el_index) {
	    if (relIsoElec==isoElecMin) {
	      if (analysisTree.electron_pt[eIndex]>ptEl) {
		isoElecMin  = relIsoElec;
		ptEl = analysisTree.electron_pt[eIndex];
		el_index = (int)eIndex;
		isoTauMin = isoTau;
		ptTau = analysisTree.tau_pt[tIndex];
		tau_index = (int)tIndex;
	      }
	    }
	    else if (relIsoElec<isoElecMin) {
	      isoElecMin  = relIsoElec;
	      ptEl = analysisTree.electron_pt[eIndex];
	      el_index = (int)eIndex;
	      isoTauMin = isoTau;
	      ptTau = analysisTree.tau_pt[tIndex];
	      tau_index = (int)tIndex;
	    }
          }
          else {
            if (isoTau==isoTauMin) {
              if (analysisTree.tau_pt[tIndex]>ptTau) {
                ptTau = analysisTree.tau_pt[tIndex];
                isoTauMin = isoTau;
                tau_index = (int)tIndex;
              }
            }
            else if (isoTau>isoTauMin) {
              ptTau = analysisTree.tau_pt[tIndex];
              isoTauMin = isoTau;
              tau_index = (int)tIndex;
            }
          }

        }

      }
      }

      bool isTauId = false;

//      if ( analysisTree.tau_againstElectronVLooseMVA6[tau_index]>0.5 &&   analysisTree.tau_againstMuonTight3[tau_index]>0.5) isTauId = true;
        cout<<"Check index"<<endl;
	if (tau_index<0)continue;
       bool tauPass =false;
       cout<<"Start Tau Id"<<endl;
	
	if (!CutBasedTauId){
			if(string::npos != TauId.find("MVA")){

				tauPass=analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tau_index] > 0.5;
				TauDecayMode= analysisTree.tau_decayMode[tau_index];
				ta_IsoFlag=analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tau_index];
				ta_IsoFlagVTight[0]=analysisTree.tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[tau_index];
				ta_IsoFlagTight[0]=analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tau_index];
				ta_IsoFlagMedium[0]=analysisTree.tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[tau_index];
				ta_IsoFlagLoose[0]=analysisTree.tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[tau_index];
				ta_IsoFlagVLoose[0]=analysisTree.tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[tau_index];
				isoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tau_index];
			}
			if(string::npos != TauId.find("MVAv2")){

				tauPass=analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tau_index] > 0.5;
				TauDecayMode= analysisTree.tau_decayMode[tau_index];
				ta_IsoFlag=analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tau_index];
				ta_IsoFlagVTight[0]=analysisTree.tau_byVTightIsolationMVArun2017v2DBoldDMwLT2017[tau_index];
				ta_IsoFlagTight[0]=analysisTree.tau_byTightIsolationMVArun2017v2DBoldDMwLT2017[tau_index];
				ta_IsoFlagMedium[0]=analysisTree.tau_byMediumIsolationMVArun2017v2DBoldDMwLT2017[tau_index];
				ta_IsoFlagLoose[0]=analysisTree.tau_byLooseIsolationMVArun2017v2DBoldDMwLT2017[tau_index];
				ta_IsoFlagVLoose[0]=analysisTree.tau_byVLooseIsolationMVArun2017v2DBoldDMwLT2017[tau_index];
				isoTau = analysisTree.tau_byIsolationMVArun2v1DBoldDMwLTraw[tau_index];
			}
			if(string::npos != TauId.find("DeepTauId")){

				tauPass=analysisTree.tau_byTightDeepTau2017v2p1VSjet[tau_index] > 0.5;
				TauDecayMode= analysisTree.tau_decayMode[tau_index];
				ta_IsoFlag=analysisTree.tau_byTightDeepTau2017v2p1VSjet[tau_index];
				ta_IsoFlagVTight[0]=analysisTree.tau_byVTightDeepTau2017v2p1VSjet[tau_index];
				ta_IsoFlagTight[0]=analysisTree.tau_byTightDeepTau2017v2p1VSjet[tau_index];
				ta_IsoFlagMedium[0]=analysisTree.tau_byMediumDeepTau2017v2p1VSjet[tau_index];
				ta_IsoFlagLoose[0]=analysisTree.tau_byLooseDeepTau2017v2p1VSjet[tau_index];
				ta_IsoFlagVLoose[0]=analysisTree.tau_byVLooseDeepTau2017v2p1VSjet[tau_index];
				isoTau = analysisTree.tau_byDeepTau2017v2p1VSjetraw[tau_index];


		}


	 }
	
	if (CutBasedTauId){
		tauPass=
	         analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau_index] > 0.5;

          isoTau = analysisTree.tau_byCombinedIsolationDeltaBetaCorrRaw3Hits[tau_index];
          ta_IsoFlag=analysisTree.tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[tau_index];
	}

	ta_relIso[0]= analysisTree.tau_byIsolationMVArun2017v2DBoldDMwLTraw2017[tau_index];
	ta_relIso[0]= isoTau;//isoTauMin;
//        el_relIso[0] = isoElecMin;
	//el_relIso[0]= 1.0;//isoMuMin;

		double q = analysisTree.tau_charge[tau_index] * analysisTree.electron_charge[el_index];
		event_sign  = q;



      cout<<"All right with Taus"<<endl;

      bool          dilepton_veto=false;
      bool          extraelec_veto=false;
      bool          extramuon_veto=false;

      event_secondLeptonVeto = false;
      event_thirdLeptonVeto = false;

      // looking for extra muon
      bool foundExtraMuon = false;
      for (unsigned int ie = 0; ie<analysisTree.muon_count; ++ie) {
	if (isData && analysisTree.muon_isDuplicate[ie]) continue;
	if (isData && analysisTree.muon_isBad[ie]) continue;
	if (analysisTree.muon_pt[ie]<ptVetoMuonCut) continue;
	if (fabs(analysisTree.muon_eta[ie])>etaVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dxy[ie])>dxyVetoMuonCut) continue;
	if (fabs(analysisTree.muon_dz[ie])>dzVetoMuonCut) continue;
	if (applyVetoMuonId && !analysisTree.muon_isMedium[ie]) continue;

	float neutralHadIsoMu = analysisTree.muon_neutralHadIso[ie];
        float photonIsoMu = analysisTree.muon_photonIso[ie];
        float chargedHadIsoMu = analysisTree.muon_chargedHadIso[ie];
        float puIsoMu = analysisTree.muon_puIso[ie];
        if (isIsoR03) {
          neutralHadIsoMu = analysisTree.muon_r04_sumNeutralHadronEt[ie];
          photonIsoMu = analysisTree.muon_r04_sumPhotonEt[ie];
          chargedHadIsoMu = analysisTree.muon_r04_sumChargedHadronPt[ie];
          puIsoMu = analysisTree.muon_r04_sumPUPt[ie];
        }
        double neutralIsoMuN = neutralHadIsoMu + photonIsoMu - 0.5*puIsoMu;
        double neutralIsoMu = max(double(0),neutralIsoMuN);
        float absIsoMu = chargedHadIsoMu + neutralIsoMu;
        float relIsoMu = absIsoMu/analysisTree.muon_pt[ie];
	if (relIsoMu>isoVetoMuonCut) continue;
	foundExtraMuon = true;
      }

      // looking for extra electron's (dielectron veto)
      bool foundExtraElectron = false;
      vector<unsigned int> e_dielectrons; e_dielectrons.clear(); 
      for (unsigned int im = 0; im<analysisTree.electron_count; ++im) {
      	   if ((int)im==(int)el_index) continue;

	float neutralHadIsoElec = analysisTree.electron_neutralHadIso[im];
	float photonIsoElec = analysisTree.electron_photonIso[im];
	float chargedHadIsoElec = analysisTree.electron_chargedHadIso[im];
	float puIsoElec = analysisTree.electron_puIso[im];
	if (isIsoR03) {
	  neutralHadIsoElec = analysisTree.electron_r03_sumNeutralHadronEt[im];
	  photonIsoElec = analysisTree.electron_r03_sumPhotonEt[im];
	  chargedHadIsoElec = analysisTree.electron_r03_sumChargedHadronPt[im];
	  puIsoElec = analysisTree.electron_r03_sumPUPt[im];
	}
	double neutralIsoElecN = neutralHadIsoElec + photonIsoElec - 0.5*puIsoElec;
	double neutralIsoElec = max(double(0),neutralIsoElecN); 
	float absIsoElec = chargedHadIsoElec + neutralIsoElec;
	float relIsoElec = absIsoElec/analysisTree.electron_pt[im];

	if (analysisTree.electron_pt[im]>ptDilepElectronCut&&
	    fabs(analysisTree.electron_eta[im])<etaDilepElectronCut&&
	    fabs(analysisTree.electron_dxy[im])<dxyDilepElectronCut&&
	    fabs(analysisTree.electron_dz[im])<dzDilepElectronCut&&
	    analysisTree.electron_cutId_veto_Summer16[im]&&
	    relIsoElec<isoDilepElectronCut && 
	    fabs(analysisTree.electron_charge[im]) ==1)
	{
	    
	float dRelectrons = deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
			   analysisTree.electron_eta[im],analysisTree.electron_phi[im]);

	    if (dRelectrons>dRDilepVetoCut && (analysisTree.electron_charge[el_index]*analysisTree.electron_charge[im]<0.)) 
	      dilepton_veto = true;

	}
     e_dielectrons.push_back(im);

	if (analysisTree.electron_pt[im]<ptVetoElectronCut) continue;
	if (fabs(analysisTree.electron_eta[im])>etaVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dxy[im])>dxyVetoElectronCut) continue;
	if (fabs(analysisTree.electron_dz[im])>dzVetoElectronCut) continue;
	bool electronMvaId = analysisTree.electron_mva_wp90_noIso_Fall17_v1[im];
	if (applyVetoElectronId && !electronMvaId) continue;
	if (applyVetoElectronId && !analysisTree.electron_pass_conversion[im]) continue;
	if (applyVetoElectronId && analysisTree.electron_nmissinginnerhits[im]>1) continue;
	if (relIsoElec>isoVetoElectronCut) continue;
	foundExtraElectron = true;
      }

      extraelec_veto = foundExtraElectron;
      extramuon_veto = foundExtraMuon;

//      if(extraelec_veto || extramuon_veto)   event_thirdLeptonVeto = true;
/*
      dilepton_veto = false;
      if (e_dielectrons.size()>1) {
	for (unsigned int i1=0; i1<e_dielectrons.size()-1; ++i1) {
	  unsigned int indx1 = e_dielectrons[i1];
	  for (unsigned int i2=i1+1; i2<e_dielectrons.size(); ++i2 ) {
	    unsigned int indx2 = e_dielectrons[i2];
	    float dRelectrons = deltaR(analysisTree.electron_eta[indx1],analysisTree.electron_phi[indx1],
				   analysisTree.electron_eta[indx2],analysisTree.electron_phi[indx2]);
	    if (dRelectrons>dRDilepVetoCut && (analysisTree.electron_charge[indx1]*analysisTree.electron_charge[indx2]<0)) 
	      dilepton_veto = true;
 	  }
	}
      }
*/
   	event_secondLeptonVeto = dilepton_veto;
//       if (dilepton_veto)  continue;

      if(extraelec_veto || extramuon_veto)   event_thirdLeptonVeto = true;


      /////////////////Make sure that the selected electron fired the trigger and it is within a Dr < 0.5


      /////////////////////////////////
      cout<<" Apply trigger SF"<<endl;
      // ////////////////////////////
      cout<<"el_index : "<<el_index<<endl;
      double ptEl1 = analysisTree.electron_pt[el_index];
      double etaEl1 = analysisTree.electron_eta[el_index];
      float trigweight=1.;
      float EffFromData = 1.;
      float EffFromMC = 1.;


      if (!isData) {
      EffFromData = (float)SF_electronTrigger->get_EfficiencyData(double(ptEl1),double(etaEl1));
      EffFromMC = (float)SF_electronTrigger->get_EfficiencyMC(double(ptEl1),double(etaEl1));

	if (EffFromMC>1e-6)trigweight = EffFromData / EffFromMC;
	if (SUSY)  trigweight = EffFromData ;

	weight *= trigweight;
	trig_weight = trigweight;
    double IdIsoSF_el = 1;
      
      
     cout<<"Corrections"<<endl; 
      correctionWS->var("e_pt")->setVal(ptEl1);
      correctionWS->var("e_eta")->setVal(etaEl1);
	  correctionWS->var("e_iso")->setVal(el_relIso[0]);
			
	IdIsoSF_el = correctionWS->function("e_idiso_binned_ratio")->getVal();
      
    IdIsoSF_el = SF_elIdIso->get_ScaleFactor(ptEl1, etaEl1);

	LSF_weight = IdIsoSF_el;
	weight *= LSF_weight;
      }


      if (!isData){
	genTauMatched = false;
	genLeptonMatched = false;
	genLeptonMatchedPromptEl = false;
	genLeptonMatchedPromptMu = false;
	genLeptonMatchedPromptTau = false;
	genElMatchedToTauDecay = false;
	genMuMatchedToTauDecay = false;
	genTauMatchedToTauDecay = false;
	genElMatchedHadrDecay = false;
	genMuMatchedHadrDecay = false;
	genTauMatchedHadrDecay = false;
	genLeptonMatchedHFQ = false;
	genLeptonMatchedLFQ = false;
	genLeptonMatchedGluon =false;
	matchedTauToPromptEl = false;
	matchedTauToPromptMu = false;
	matchedTauToTauDecEl =false;
	matchedTauToTauDecMu =false;
	matchedTauToElHadronDec = false;
	matchedTauToMuHadronDec = false;
	matchedTauToTauHadronDec = false;
	matchedTauToGluon = false;
	matchedTauToHFQ = false;
	matchedTauToLFQ = false;
	genTauDecayMode1=-1;
	genTauDecayMode2=-1;
	TLorentzVector genTauV;
	TLorentzVector genPartV;

	bool FoundFirstMatchedTau = false;

	for (unsigned int gt = 0 ; gt < analysisTree.gentau_count; ++gt){

	  genTauV.SetXYZT(analysisTree.gentau_px[gt], analysisTree.gentau_py[gt], analysisTree.gentau_pz[gt], analysisTree.gentau_e[gt]);

	  double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			    genTauV.Eta(), genTauV.Phi());

	  if (Drr < 0.5 && ( analysisTree.gentau_isPrompt[gt] > 0.5  ) && genTauV.Pt() > 15. ) genTauMatched = true;
	  if (genTauMatched && !FoundFirstMatchedTau) {genTauDecayMode1 = analysisTree.gentau_decayMode[gt]; FoundFirstMatchedTau = true;}
	  if (genTauMatched && FoundFirstMatchedTau)   genTauDecayMode2 = analysisTree.gentau_decayMode[gt];

	}

	  for (unsigned int igen=0; igen<analysisTree.genparticles_count; ++igen) {

      		  if ( (abs(analysisTree.genparticles_pdgid[igen])==11 || abs(analysisTree.genparticles_pdgid[igen])==13 || abs(analysisTree.genparticles_pdgid[igen])==15 || abs(analysisTree.genparticles_pdgid[igen])<6 || abs(analysisTree.genparticles_pdgid[igen])==21)){

	  genPartV.SetXYZT(analysisTree.genparticles_px[igen], analysisTree.genparticles_py[igen], analysisTree.genparticles_pz[igen], analysisTree.genparticles_e[igen]);

	  double Drl=deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
			  genPartV.Eta(),genPartV.Phi());

	double DrTauLepton=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
			  genPartV.Eta(),genPartV.Phi());

		if (Drl < 0.5 && genPartV.Pt() > 8){
        //if ( abs(analysisTree.genparticles_pdgid[igen])==13)
        //cout<<analysisTree.genparticles_pdgid[igen]<<" isDirectHadronDecayProduct "<<analysisTree.genparticles_isDirectHadronDecayProduct[igen]<<" isTauDecay "<<analysisTree.genparticles_isTauDecayProduct[igen]<<" isDirectTau  "<<analysisTree.genparticles_isDirectTauDecayProduct[igen]<<"  isPrompt  "<<analysisTree.genparticles_isPrompt[igen]<<"  isDecayedLeptonHadron "<<analysisTree.genparticles_isDecayedLeptonHadron[igen]<<" is there gentau matched ?  "<<genTauMatched<<endl;
        genLeptonMatched = true;

		if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] > 0.5) genLeptonMatchedPromptEl = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] > 0.5) genLeptonMatchedPromptMu = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==15 && analysisTree.genparticles_isPrompt[igen] > 0.5) genLeptonMatchedPromptTau = true;

		if (abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] > 0.5 ) genElMatchedToTauDecay = true;
		if (abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] > 0.5 ) genMuMatchedToTauDecay = true;
		if (abs(analysisTree.genparticles_pdgid[igen])==15 && analysisTree.genparticles_isDirectPromptTauDecayProduct[igen] > 0.5 ) genTauMatchedToTauDecay = true;

		if (abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5 ) genElMatchedHadrDecay = true;
		if (abs(analysisTree.genparticles_pdgid[igen])==13 && (analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5 ||  analysisTree.genparticles_isDecayedLeptonHadron[igen] >0.5 )) genMuMatchedHadrDecay = true;
		if (abs(analysisTree.genparticles_pdgid[igen])==15 && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5 ) genTauMatchedHadrDecay = true;

		if ( (abs(analysisTree.genparticles_pdgid[igen])==4 || abs(analysisTree.genparticles_pdgid[igen])==5) && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5) genLeptonMatchedHFQ = true;

		if ( (abs(analysisTree.genparticles_pdgid[igen])<4) && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5) genLeptonMatchedLFQ = true;
		if ( (abs(analysisTree.genparticles_pdgid[igen])==21 ) && analysisTree.genparticles_isDirectHadronDecayProduct[igen] > 0.5) genLeptonMatchedGluon = true;

		}

		if (DrTauLepton < 0.5 && genPartV.Pt() > 8. ) 
			
		{

		if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.genparticles_isPrompt[igen] > 0.5 ) matchedTauToPromptEl = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.genparticles_isPrompt[igen] > 0.5 ) matchedTauToPromptMu = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.gentau_isDirectPromptTauDecayProduct[igen] > 0.5 ) matchedTauToTauDecEl = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.gentau_isDirectPromptTauDecayProduct[igen] > 0.5 ) matchedTauToTauDecMu = true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==11 && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToElHadronDec= true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==13 && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToMuHadronDec= true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==15 && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToTauHadronDec= true;
		if ( abs(analysisTree.genparticles_pdgid[igen])==21 && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToGluon= true;
		if ( (abs(analysisTree.genparticles_pdgid[igen])==4  || abs(analysisTree.genparticles_pdgid[igen])==5 ) && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToHFQ= true;
		if ( (abs(analysisTree.genparticles_pdgid[igen])<4 ) && analysisTree.gentau_isDirectHadronDecayProduct[igen] > 0.5 ) matchedTauToLFQ= true;

		}

		}//loop for gen 11 13 15
      
	  }
      
      
      }//!isData



      electron_index = (int)el_index;
      taus_index = (int)tau_index;

      el_count=(int)electrons.size();


      for (unsigned int ie=0;ie<electrons.size(); ++ie){
	  unsigned int eIndex = electrons[ie];
	el_px[ie]=analysisTree.electron_px[eIndex];
	el_py[ie]=analysisTree.electron_py[eIndex];
	el_pz[ie]=analysisTree.electron_pz[eIndex];
	el_eta[ie]=analysisTree.electron_eta[eIndex];
	el_pt[ie]=analysisTree.electron_pt[eIndex];
	el_phi[ie]=analysisTree.electron_phi[eIndex];
	el_charge[ie]=analysisTree.electron_charge[eIndex];
	el_dxy[ie]=analysisTree.electron_dxy[eIndex];
	el_dz[ie]=analysisTree.electron_dz[eIndex];
	el_dxyerr[ie]=analysisTree.electron_dxyerr[eIndex];
	el_dzerr[ie]=analysisTree.electron_dzerr[eIndex];

        el_neutralHadIso[ie] = analysisTree.electron_r03_sumNeutralHadronEt[eIndex];
        el_photonIso[ie] = analysisTree.electron_r03_sumPhotonEt[eIndex];
        el_chargedHadIso[ie] = analysisTree.electron_r03_sumChargedHadronPt[eIndex];
        el_puIso[ie] = analysisTree.electron_r03_sumPUPt[eIndex];
 
        double neutralIso = el_neutralHadIso[eIndex] + el_photonIso[eIndex] - 0.5*el_puIso[eIndex];
        neutralIso = max(double(0),neutralIso);
        el_neutralIso[ie] = neutralIso ;
        el_absIsoEl[ie] = el_chargedHadIso[eIndex] + el_neutralIso[eIndex];
	if ( (int)eIndex>-1)el_relIsoEl[ie]  = el_absIsoEl[eIndex]/el_pt[eIndex] ;
}

      ta_count=(int)taus.size();
      for (unsigned int it=0;it<taus.size(); ++it){
            unsigned int itt = taus[it];
            ta_px[it]=analysisTree.tau_px[itt];
            ta_py[it]=analysisTree.tau_py[itt];
            ta_pz[it]=analysisTree.tau_pz[itt];
            ta_eta[it]=analysisTree.tau_eta[itt];
            ta_pt[it]=analysisTree.tau_pt[itt];
            ta_phi[it]=analysisTree.tau_phi[itt];
            ta_charge[it]=analysisTree.tau_charge[itt];
            ta_dxy[it]=analysisTree.tau_dxy[itt];
            ta_dz[it]=analysisTree.tau_dz[itt];
            //
            ta_puCorrPtSum[it] = analysisTree.tau_puCorrPtSum[itt];
            ta_chargedIsoPtSum[it] = analysisTree.tau_chargedIsoPtSum[itt];
            ta_neutralIsoPtSum[it] = analysisTree.tau_neutralIsoPtSum[itt];
      }

      jet_count=(int)analysisTree.pfjet_count;
      for (unsigned int jj=0;jj<analysisTree.pfjet_count; ++jj){
        jet_e[jj] = analysisTree.pfjet_e[jj];
        jet_px[jj] = analysisTree.pfjet_px[jj];
        jet_py[jj] = analysisTree.pfjet_py[jj];
        jet_pz[jj] = analysisTree.pfjet_pz[jj];
        jet_pt[jj] = analysisTree.pfjet_pt[jj];
        jet_eta[jj] = analysisTree.pfjet_eta[jj];
        jet_phi[jj] = analysisTree.pfjet_phi[jj];
        jet_flavour[jj] = analysisTree.pfjet_flavour[jj];
        jet_btag[jj] = analysisTree.pfjet_btag[jj][0];
      }

      float jetEta = 2.4;
      float DRmax = 0.5;

      float bJetEtaCut = jetEta;

      vector<unsigned int> jets; jets.clear();
      vector<unsigned int> jetspt20; jetspt20.clear();
      vector<unsigned int> bjets; bjets.clear();
      vector<unsigned int> bjetsTight; bjetsTight.clear();
      vector<unsigned int> bjets_nocleaned; bjets_nocleaned.clear();

        int indexLeadingJet = -1;
        float ptLeadingJet = -1;

        int indexSubLeadingJet = -1;
        float ptSubLeadingJet = -1;

        int indexLeadingBJet = -1;

        int counter_cleaned_jets = 0;


      for (unsigned int jet=0; jet<analysisTree.pfjet_count; ++jet) {

	double drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
						  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);


		if (drr < 0.5 && !isData) 
		{

	     if (analysisTree.pfjet_flavour[jet] == 21) matchedTauToGluon = true;
	     if (abs(analysisTree.pfjet_flavour[jet]) == 4 || abs(analysisTree.pfjet_flavour[jet]) == 5) matchedTauToHFQ = true;
	     if (abs(analysisTree.pfjet_flavour[jet]) < 4) matchedTauToLFQ = true;
	     if (abs(analysisTree.pfjet_flavour[jet]) == 1 ) matchedTauToDownQ = true;
             if (abs(analysisTree.pfjet_flavour[jet]) == 2 ) matchedTauToUpQ = true;
             if (abs(analysisTree.pfjet_flavour[jet]) == 3 ) matchedTauToStrangeQ = true;
             if (abs(analysisTree.pfjet_flavour[jet]) == 4 ) matchedTauToCharmQ = true;
             if (abs(analysisTree.pfjet_flavour[jet]) == 5 ) matchedTauToBottomQ = true;

			}

	  double drl=deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
						  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);


		if (drl < 0.5 && !isData) 
		{

	     if (analysisTree.pfjet_flavour[jet] == 21) genLeptonMatchedGluon = true;
	     if (abs(analysisTree.pfjet_flavour[jet]) == 4 || abs(analysisTree.pfjet_flavour[jet]) == 5) genLeptonMatchedHFQ = true;
	     if (abs(analysisTree.pfjet_flavour[jet]) < 4) genLeptonMatchedLFQ = true;

			}

	if (fabs(analysisTree.pfjet_pt[jet])<ptJetCut) continue;
        float absJetEta = fabs(analysisTree.pfjet_eta[jet]);
	if (absJetEta > etaJetCut) continue;

	float jetPt = analysisTree.pfjet_pt[jet];

	bool isPFJetId = false ;
	isPFJetId =tightJetiD_2017(analysisTree,jet);

	if (!isPFJetId) continue;
	bool cleanedJet = true;
	bool cleanedJetTau = true;
	bool cleanedJetMu = true;

       double Dr=deltaR(analysisTree.electron_eta[el_index],analysisTree.electron_phi[el_index],
			 analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);
	if (  Dr  < DRmax) {
		cleanedJet=false;
		cleanedJetMu=false;


	double Drr=deltaR(analysisTree.tau_eta[tau_index],analysisTree.tau_phi[tau_index],
						  analysisTree.pfjet_eta[jet],analysisTree.pfjet_phi[jet]);

	if ( Drr < DRmax) {

		cleanedJet=false;
		cleanedJetTau=false;
	}
	if (!cleanedJet) continue;

	if (absJetEta<bJetEtaCut) { // jet within b-tagging acceptance

      	bool btaggedTight  (analysisTree.pfjet_btag[jet][0]  > 0.8838) ;
      	bool btagged  (analysisTree.pfjet_btag[jet][0]  > bTag) ;

	  if (!isData) {
	    int flavor = abs(analysisTree.pfjet_flavour[jet]);

	    double jet_scalefactor = 1;
	    double JetPtForBTag = jetPt;
	    double tageff = 1;

	    if (flavor==5) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader_B.eval_auto_bounds(BTag_,BTagEntry::FLAV_B, absJetEta, JetPtForBTag);
	      tageff = tagEff_B->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else if (flavor==4) {
	      if (JetPtForBTag>MaxBJetPt) JetPtForBTag = MaxBJetPt - 0.1;
	      if (JetPtForBTag<MinBJetPt) JetPtForBTag = MinBJetPt + 0.1;
	      jet_scalefactor = reader_C.eval_auto_bounds(BTag_,BTagEntry::FLAV_C, absJetEta, JetPtForBTag);
	      tageff = tagEff_C->Interpolate(JetPtForBTag,absJetEta);
	    }
	    else {
	      if (JetPtForBTag>MaxLJetPt) JetPtForBTag = MaxLJetPt - 0.1;
	      if (JetPtForBTag<MinLJetPt) JetPtForBTag = MinLJetPt + 0.1;
	      jet_scalefactor = reader_Light.eval_auto_bounds(BTag_,BTagEntry::FLAV_UDSG, absJetEta, JetPtForBTag);
	      tageff = tagEff_Light->Interpolate(JetPtForBTag,absJetEta);
	    }

	    if (tageff<1e-5)      tageff = 1e-5;
	    if (tageff>0.99999)   tageff = 0.99999;
	    rand.SetSeed((int)((absJetEta+5)*100000));
	    double rannum = rand.Rndm();

	    if (jet_scalefactor<1 && btagged) { // downgrade
	      double fraction = 1-jet_scalefactor;
	      if (rannum<fraction) {
		btagged = false;
		//		std::cout << "downgrading " << std::endl;
	      }
	    }
	    if (jet_scalefactor>1 && !btagged) { // upgrade
	      double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
	      if (rannum<fraction) {
		btagged = true;
		//		std::cout << "upgrading " << std::endl;
	      }
	    }

	    if (jet_scalefactor<1 && btaggedTight) { // downgrade
	      double fraction = 1-jet_scalefactor;
	      if (rannum<fraction) {
		btaggedTight = false;
		//		std::cout << "downgrading " << std::endl;
	      }
	    }
	    if (jet_scalefactor>1 && !btaggedTight) { // upgrade
	      double fraction = (jet_scalefactor-1.0)/(1.0/tageff-1.0);
	      if (rannum<fraction) {
		btaggedTight = true;
		//		std::cout << "upgrading " << std::endl;
	      }
	    }


	  } //is Data

	  if (btagged && cleanedJet) 	   bjets.push_back(jet);
	  if (btaggedTight && cleanedJet) 	   bjetsTight.push_back(jet);

	  if (btaggedTight && cleanedJet) bjets_cleanedTight[counter_cleaned_jets]=(int)jet;
	  if (btagged && cleanedJet) bjets_cleaned[counter_cleaned_jets]=(int)jet;

	}

	if (cleanedJet){

	jets.push_back((int)jet);
	jets_cleaned[counter_cleaned_jets]=(int)jet;
	jet_jecUn[counter_cleaned_jets] = analysisTree.pfjet_jecUncertainty[jet];
	counter_cleaned_jets++;
	}

      }///loop in all jets
    }
      njets = jets.size();
      jet_count = jets.size();
      //njetspt20 = jetspt20.size();
      nbtag = bjets.size();
      nbtagTight = bjetsTight.size();
      //nbtag_nocleaned = bjets_nocleaned.size();

      npv =  analysisTree.primvertex_count;
      npu = analysisTree.numtruepileupinteractions;
      if(SUSY){	
      SusyMother = SusyMotherMassF;

      SusyLSP = SusyLSPMassF;
      }

      //cout<<"Recoil corrections"<<endl;
      int njetsforrecoil = njets;
      if (isW) njetsforrecoil = njets + 1;

    ////while using old MC ntuples, need to use proper MET collection
      float pfmet_corr_x = 1.;
      float pfmet_corr_y = 1.;
      float met_x = 1.;
      float met_y = 1.;
      float puppimet_ex = 1.;
      float puppimet_ey = 1.;

      pfmet_corr_x = analysisTree.pfmetcorr_ex;
      pfmet_corr_y = analysisTree.pfmetcorr_ey;

      puppimet_ex = analysisTree.pfmetcorr_ex;
      puppimet_ey = analysisTree.pfmetcorr_ex;

      met_x = analysisTree.pfmetcorr_ex;
      met_y = analysisTree.pfmetcorr_ey;

      if ((isW || isDY) && !isData) {
            recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,pfmet_corr_x,pfmet_corr_y);
            float met_corr_x=1.;
            float met_corr_y=1.;
            met_x= analysisTree.pfmetcorr_ex_JetEnUp;
            met_y= analysisTree.pfmetcorr_ey_JetEnUp;
            met_corr_x= analysisTree.pfmetcorr_ex_JetEnUp;
            met_corr_y= analysisTree.pfmetcorr_ey_JetEnUp;

	        recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,met_corr_x,met_corr_y);
            met_ex_JetEnUp_recoil = met_corr_x;
            met_ey_JetEnUp_recoil = met_corr_y;


            met_x= analysisTree.pfmetcorr_ex_JetEnDown;
            met_y= analysisTree.pfmetcorr_ey_JetEnDown;
            met_corr_x= analysisTree.pfmetcorr_ex_JetEnDown;
            met_corr_y= analysisTree.pfmetcorr_ey_JetEnDown;
            recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,met_corr_x,met_corr_y);
            met_ex_JetEnDown_recoil = met_corr_x;
            met_ey_JetEnDown_recoil = met_corr_y;


            met_x=analysisTree.pfmetcorr_ex_UnclusteredEnUp;
            met_y=analysisTree.pfmetcorr_ey_UnclusteredEnUp;
            met_corr_x=analysisTree.pfmetcorr_ex_UnclusteredEnUp;
            met_corr_y=analysisTree.pfmetcorr_ey_UnclusteredEnUp;
            recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,met_corr_x,met_corr_y);
            met_ex_UnclusteredEnUp_recoil = met_corr_x;
            met_ey_UnclusteredEnUp_recoil = met_corr_y;


            met_x=analysisTree.pfmetcorr_ex_UnclusteredEnDown;
            met_y=analysisTree.pfmetcorr_ey_UnclusteredEnDown;
            met_corr_x=analysisTree.pfmetcorr_ex_UnclusteredEnDown;
            met_corr_y=analysisTree.pfmetcorr_ey_UnclusteredEnDown;
            recoilMetCorrector.CorrectByMeanResolution(met_x,met_y,bosonPx,bosonPy,lepPx,lepPy,njetsforrecoil,met_corr_x,met_corr_y);
            met_ex_UnclusteredEnDown_recoil = met_corr_x;
            met_ey_UnclusteredEnDown_recoil = met_corr_y;


            met_x = pfmet_corr_x;
            met_y = pfmet_corr_y;

      // MEt related systematic uncertainties

      int bkgdType = 0;
      if (isDY||isW)
	bkgdType = kit::MEtSys::ProcessType::BOSON;
      else if (isTOP)
	bkgdType = kit::MEtSys::ProcessType::TOP;
      else
	bkgdType = kit::MEtSys::ProcessType::EWK;

      float met_scaleUp_x   = met_x;
      float met_scaleUp_y   = met_y;
      float met_scaleDown_x = met_x;
      float met_scaleDown_y = met_y;
      float met_resoUp_x    = met_x;
      float met_resoUp_y    = met_y;
      float met_resoDown_x  = met_x;
      float met_resoDown_y  = met_y;

      met_scaleUp = TMath::Sqrt(met_scaleUp_x*met_scaleUp_x+
				   met_scaleUp_y*met_scaleUp_y);
      metphi_scaleUp = TMath::ATan2(met_scaleUp_y,met_scaleUp_x);
      
      met_scaleDown = TMath::Sqrt(met_scaleDown_x*met_scaleDown_x+
				     met_scaleDown_y*met_scaleDown_y);
      metphi_scaleDown = TMath::ATan2(met_scaleDown_y,met_scaleDown_x);
      
      met_resoUp = TMath::Sqrt(met_resoUp_x*met_resoUp_x+
				  met_resoUp_y*met_resoUp_y);
      metphi_resoUp = TMath::ATan2(met_resoUp_y,met_resoUp_x);
      
      met_resoDown = TMath::Sqrt(met_resoDown_x*met_resoDown_x+
				    met_resoDown_y*met_resoDown_y);
      metphi_resoDown = TMath::ATan2(met_resoDown_y,met_resoDown_x);

      met_ex_recoil = pfmet_corr_x;
      met_ey_recoil = pfmet_corr_y;
      //revert back to uncorrected met



      }//if isW, isDY !isData

      met_ex = analysisTree.pfmetcorr_ex;
      met_ey = analysisTree.pfmetcorr_ey;
      met_ez = analysisTree.pfmetcorr_ez;
      met_pt = TMath::Sqrt(met_ex*met_ex + met_ey*met_ey);
      met_phi = TMath::ATan2(met_y,met_x);

      //PF Corrected MET:
      met_ex_JetEnUp = analysisTree.pfmetcorr_ex_JetEnUp;
      met_ey_JetEnUp = analysisTree.pfmetcorr_ey_JetEnUp;

      met_ex_JetEnDown = analysisTree.pfmetcorr_ex_JetEnDown;
      met_ey_JetEnDown = analysisTree.pfmetcorr_ey_JetEnDown;

      met_ex_UnclusteredEnUp = analysisTree.pfmetcorr_ex_UnclusteredEnUp;
      met_ey_UnclusteredEnUp = analysisTree.pfmetcorr_ey_UnclusteredEnUp;

      met_ex_UnclusteredEnDown = analysisTree.pfmetcorr_ex_UnclusteredEnDown;
      met_ey_UnclusteredEnDown = analysisTree.pfmetcorr_ey_UnclusteredEnDown;


      met_ex_JetResUp = analysisTree.pfmetcorr_ex_JetResUp;
      met_ey_JetResUp  = analysisTree.pfmetcorr_ey_JetResUp;


      met_ex_JetResDown = analysisTree.pfmetcorr_ex_JetResDown;
      met_ey_JetResDown = analysisTree.pfmetcorr_ey_JetResDown;

      //Puppi MET Unc.

      puppimet_pt = TMath::Sqrt(puppimet_ex*puppimet_ex+puppimet_ey*puppimet_ey);
      puppimet_phi = TMath::ATan2(puppimet_ey,puppimet_ex);

      puppimet_ex_JetEnUp = analysisTree.puppimet_ex_JetEnUp;
      puppimet_ey_JetEnUp = analysisTree.puppimet_ey_JetEnUp;

      puppimet_ex_JetEnDown = analysisTree.puppimet_ex_JetEnDown;
      puppimet_ey_JetEnDown = analysisTree.puppimet_ey_JetEnDown;

      puppimet_ex_UnclusteredEnUp = analysisTree.puppimet_ex_UnclusteredEnUp;
      puppimet_ey_UnclusteredEnUp = analysisTree.puppimet_ey_UnclusteredEnUp;

      puppimet_ex_UnclusteredEnDown = analysisTree.puppimet_ex_UnclusteredEnDown;
      puppimet_ey_UnclusteredEnDown = analysisTree.puppimet_ey_UnclusteredEnDown;

      puppimet_ex_JetResUp = analysisTree.puppimet_ex_JetResUp;
      puppimet_ey_JetResUp = analysisTree.puppimet_ey_JetResUp;

      puppimet_ex_JetResDown = analysisTree.puppimet_ex_JetResDown;
      puppimet_ey_JetResDown = analysisTree.puppimet_ey_JetResDown;

      //Puppi PT and Phi calculation:
      puppimet_pt_JetEnUp = TMath::Sqrt(puppimet_ex_JetEnUp*puppimet_ex_JetEnUp+
                                     puppimet_ey_JetEnUp*puppimet_ey_JetEnUp);

      puppimet_phi_JetEnUp = TMath::ATan2(puppimet_ex_JetEnUp, puppimet_ey_JetEnUp);

      puppimet_pt_JetEnDown = TMath::Sqrt(puppimet_ex_JetEnDown*puppimet_ex_JetEnDown+
                                        puppimet_ey_JetEnDown*puppimet_ey_JetEnDown);

      puppimet_pt_JetEnDown = TMath::ATan2(puppimet_ex_JetEnDown, puppimet_ey_JetEnDown);

      puppimet_pt_UnclusteredEnUp = TMath::Sqrt(puppimet_ex_UnclusteredEnUp*puppimet_ex_UnclusteredEnUp+
                                            puppimet_ey_UnclusteredEnUp*puppimet_ey_UnclusteredEnUp);

      puppimet_phi_UnclusteredEnUp = TMath::ATan2(puppimet_ex_UnclusteredEnUp, puppimet_ey_UnclusteredEnUp);

      puppimet_pt_UnclusteredEnDown = TMath::Sqrt(puppimet_ex_UnclusteredEnDown*puppimet_ex_UnclusteredEnDown+
                                               puppimet_ey_UnclusteredEnDown*puppimet_ey_UnclusteredEnDown);

      puppimet_phi_UnclusteredEnDown = TMath::ATan2(puppimet_ex_UnclusteredEnDown, puppimet_ey_UnclusteredEnDown);

      puppimet_pt_JetResUp = TMath::Sqrt(puppimet_ex_JetResUp*puppimet_ex_JetResUp +
                                     puppimet_ey_JetResUp*puppimet_ey_JetResUp);
      puppimet_phi_JetResUp = TMath::ATan2(puppimet_ex_JetResUp, puppimet_ey_JetResUp);

      puppimet_pt_JetResDown = TMath::Sqrt(puppimet_ex_JetResDown*puppimet_ex_JetResDown +
                                            puppimet_ey_JetResDown*puppimet_ey_JetResDown);

      puppimet_phi_JetResDown = TMath::ATan2(puppimet_ex_JetResDown, puppimet_ey_JetResDown);


      //TODO Remove it
      float genmet_ex = 0;// analysisTree.genmet_ex;
      float genmet_ey = 0;//analysisTree.genmet_ey;

      genmet = TMath::Sqrt(genmet_ex*genmet_ex + genmet_ey*genmet_ey);
      genmetphi = TMath::ATan2(genmet_ey,genmet_ex);

      if (!isData) npartons = analysisTree.genparticles_noutgoing;

      all_weight = weight;
      event_run = analysisTree.event_run;
      event_lumi = analysisTree.event_luminosityblock;
      NuPx = nuPx;
      NuPy = nuPy;
      NuPz = nuPz;
      NuPt = nuPt;
      NuPhi = nuPhi;
      cout<<"Fill"<<endl;
      //genmet_Ex = analysisTree.genmet_ex;
      //genmet_Ey = analysisTree.genmet_ey;
      //genHT = analysisTree.genparticles_lheHt;
      T->Fill();
	  selEvents++;

    } // end of file processing (loop over events in one file)


    nFiles++;
    delete _tree;
    file_->Close();
    delete file_;
  }


  //cout<<" Run range from -----> "<<RunMin<<" to  "<<RunMax<<endl;


  std::cout << std::endl;
  int allEvents = (int)inputEventsH->GetEntries();
  std::cout << "Total number of input events    = " << allEvents << std::endl;
  std::cout << "Total number of events in Tree  = " << nEvents << std::endl;
  std::cout << "Total number of selected events = " << selEvents << std::endl;
  std::cout << std::endl;


  file->cd(Channel.c_str());
  hxsec->Fill(XSec);
  hxsec->Write();
  inputEventsH->Write();
  histWeightsH->Write();
  histWeightsScalesUp->Write();
  histWeightsScalesDown->Write();
  histWeightsPDFUp->Write();
  histWeightsPDFDown->Write();
  histTopPt->Write();
  histTopPtSq->Write();
  histRuns->Write();
  CutFlowUnW->Write();
  histPt->Write();
  file->Write();
  file->Close();

  cout<<"done"<<endl;
  delete file;

}

