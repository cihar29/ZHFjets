#!/bin/bash

# execute as ./analyze_all.sh

channels=( "mm" "ee" )

dir="/uscmst1b_scratch/lpc1/lpctrig/duong/Output_ZplusC_V25/Zll_V25/syst_fromEOS/"
files=(
  "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"
  "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
  "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8"
  "WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8"
  "WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
  "WW_TuneCUETP8M1_13TeV-pythia8"
  "WZ_TuneCUETP8M1_13TeV-pythia8"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016B-23Sep2016-v3"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016C-23Sep2016-v1"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016D-23Sep2016-v1"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016E-23Sep2016-v1"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016F-23Sep2016-v1"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016G-23Sep2016-v1"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016H-PromptReco-v1"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016H-PromptReco-v2"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleEG__Run2016H-PromptReco-v3"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016B-23Sep2016-v3"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016C-23Sep2016-v1"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016D-23Sep2016-v1"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016E-23Sep2016-v1"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016F-23Sep2016-v1"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016G-23Sep2016-v1"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016H-PromptReco-v1"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016H-PromptReco-v2"
  "ZC2017_03_RRSEP2016_DOUBLELEP_DoubleMuon__Run2016H-PromptReco-v3"
  "ZZ_TuneCUETP8M1_13TeV-pythia8"
)

start=$(date +%s.%N)

mkdir logs
for chan in "${channels[@]}" ; do
  mkdir $chan
  log="logs/log_${chan}.txt"
  echo "" > $log

  for file in "${files[@]}" ; do

    mc="true"
    if [[ $file == *"DOUBLELEP"* ]] ; then
      mc="false"
    fi

    lines=( "isMC           ${mc}"
            "inName         ${dir}${file}.root"
            "outName        ${chan}/${file}_${chan}.root"
            "channel        ${chan}"
            )
    out=""

    for line in "${lines[@]}" ; do
      out="$out$line\n"
    done

    echo -e "$out" | column -t > pars.txt

    echo "Processing" $file
    analyze pars.txt >> $log

  done

done

end=$(date +%s.%N)
echo "$(echo "$end - $start" | bc) seconds"
