date
mkdir -p log
root -q -l -b eregtraining_13TeV_Pi0_bothGamma.C+(false,false) > log/reg_log_EE.log
root -q -l -b eregtraining_13TeV_Pi0_bothGamma.C+(true,false) > log/reg_log_EB.log 
root -q -l -b eregtesting_13TeV_Pi0_bothGamma.C+(false,false) > log/test_log_EE.log
root -q -l -b eregtesting_13TeV_Pi0_bothGamma.C+(true,false) > log/test_log_EB.log
date
