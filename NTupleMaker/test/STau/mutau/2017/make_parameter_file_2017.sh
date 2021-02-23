#!/bin/bash

CHANNEL=$1
echo "CONFIGFILE,FILELIST" > parameters.txt

if [[ $CHANNEL == "mt" ]]; then
  # data
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleMuon_Run2017B 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleMuon_Run2017C 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleMuon_Run2017D 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleMuon_Run2017E 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleMuon_Run2017F 4

  # Embedded
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedMuTau_Run2017B 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedMuTau_Run2017C 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedMuTau_Run2017D 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedMuTau_Run2017E 4
  ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedMuTau_Run2017F 4
else
  if [[ $CHANNEL == "et" ]]; then
    # data
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleElectron_Run2017B 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleElectron_Run2017C 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleElectron_Run2017D 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleElectron_Run2017E 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_data.conf SingleElectron_Run2017F 4

    # Embedded
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElTau_Run2017B 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElTau_Run2017C 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElTau_Run2017D 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElTau_Run2017E 4
    ./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_embedded.conf EmbeddedElTau_Run2017F 4
  else
    echo "channel is not mt or et - exiting make_parameter_file_2017.sh"
    exit
  fi
fi

# SUSY Signals
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_ctau-0p01to10_mLSP-50to100_TuneCP5_13TeV 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_ctau-0p01to10_mStau-250to500_TuneCP5_13TeV 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_ctau-0p01to10_mStau-90_TuneCP5_13TeV 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_ctau-0p01to10_TuneCP5_13TeV-madgraphMLM 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau-Ewkino_lefthanded_dM-10to40_genHT-80_genMET 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau-Ewkino_lefthanded_dM-50_genHT-80_genMET 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_lefthanded_dM-10to50_genHT-80_genMET 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_lefthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_lefthanded_mStau-275to500_TuneCP5_13TeV-madgraphMLM 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_lefthanded_mStau-90_mLSP-20_TuneCP5_13TeV 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_lefthanded_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIFall17MiniAODv2 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_righthanded_mLSP-30to200_TuneCP5_13TeV-madgraphMLM 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_righthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_righthanded_mStau-275to500_TuneCP5_13TeV-madgraphMLM 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_righthanded_mStau-90_mLSP-20_TuneCP5_13TeV 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf SMS-TStauStau_righthanded_TuneCP5_13TeV-madgraphMLM-pythia8 4

# Signals
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf GluGluHToTauTauUncorrDecays_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf VBFHToTauTauUncorrDecays_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WminusHToTauTauUncorrDecays_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WplusHToTauTauUncorrDecays_M125 4
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ZHToTauTauUncorrDecays_M125 4

# DY
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_DYJetsToLL_M-10to50_13TeV-12Apr2018.conf DYJetsToLL_M-10to50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DYJetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf DY1JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf DY2JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY3JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.conf DY4JetsToLL_M-50 10

# W+jets
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WJetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf W1JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf W2JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf W3JetsToLNu 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.conf W4JetsToLNu 20

# EWK
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf EWKWPlus2Jets_WToLNu_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf EWKWMinus2Jets_WToLNu_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf EWKZ2Jets_ZToLL_M-50 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf EWKZ2Jets_ZToNuNu 10

# VV
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_WW_TuneCP5_13TeV-pythia8.conf WW 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_WZ_TuneCP5_13TeV-pythia8.conf WZ 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_ZZ_TuneCP5_13TeV-pythia8.conf ZZ 20

# Exclusive VV
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf VVTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf WZTo3LNu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ZZTo2L2Q 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ZZTo4L 10

# H->WW
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf GluGluHToWWTo2L2Nu_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf VBFHToWWTo2L2Nu_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf GluGluZH_HToWW_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf HWminusJ_HToWW_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf HWplusJ_HToWW_M125 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf HZJ_HToWW_M125 10

# TT
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8.conf TTTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf TTToHadronic 10
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC_TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8.conf TTToSemiLeptonic 10

# Single Top
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ST_t-channel_antitop_4f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ST_t-channel_top_4f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ST_tW_antitop_5f 20
./split_filelist.sh analysisMacroSynch_${CHANNEL}_17_MC.conf ST_tW_top_5f 20
