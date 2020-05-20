#!/bin/bash

echo "CONFIGFILE,FILELIST" > parameters.txt

## Test for 2018
#data
./split_filelist.sh analysisMacroSynch_mt_18_data.conf EGamma_Run2018A 10
./split_filelist.sh analysisMacroSynch_mt_18_data.conf EGamma_Run2018B 10
./split_filelist.sh analysisMacroSynch_mt_18_data.conf EGamma_Run2018C 10
./split_filelist.sh analysisMacroSynch_mt_18_data.conf EGamma_Run2018D 10

# Signals
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf GluGluHToTauTau_M125 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf VBFHToTauTau_M125 10
##  commented since now switching to TauSpinner weights
# ./split_filelist.sh analysisMacroSynch_mt_SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8.conf SUSYGluGluToHToTauTau 4

# Tau Spinner
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf GluGluHToTauTauUncorrDecays_M125 4
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf VBFHToTauTauUncorrDecays_M125 4

# DY
# TODO: add low mass DY (need to switch to the other PU file)
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf DYJetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf DY1JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf DY2JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf DY3JetsToLL_M-50 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf DY4JetsToLL_M-50 10

# Embedded
./split_filelist.sh analysisMacroSynch_mt_18_embedded.conf EmbeddedElTau_Run2018A 10
./split_filelist.sh analysisMacroSynch_mt_18_embedded.conf EmbeddedElTau_Run2018B 10
./split_filelist.sh analysisMacroSynch_mt_18_embedded.conf EmbeddedElTau_Run2018C 10
./split_filelist.sh analysisMacroSynch_mt_18_embedded.conf EmbeddedElTau_Run2018D 10

# W+jets
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf WJetsToLNu 20
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf W1JetsToLNu 20
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf W2JetsToLNu 20
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf W3JetsToLNu 20
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf W4JetsToLNu 20

# VV
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf WW 20
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf WZ 20
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf ZZ 20

# TT
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf TTTo2L2Nu 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf TTToHadronic 10
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf TTToSemiLeptonic 10

# Single Top
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf ST_t-channel_antitop_4f 20
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf ST_t-channel_top_4f 20
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf ST_tW_antitop_5f 20
./split_filelist.sh analysisMacroSynch_mt_18_MC.conf ST_tW_top_5f 20
