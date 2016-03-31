#!/bin/bash
cd /afs/cern.ch/work/z/zhicaiz/private/regression/CMSSW_7_4_12/src/HiggsAnalysis/GBRLikelihood/
ulimit -c 0
eval `scramv1 runtime -sh`
date 
root -q -l -b eregtraining_13TeV_Pi0.C+\(false,false,false\) 
root -q -l -b eregtesting_13TeV_Pi0.C+\(false,false,false\) 
date 
