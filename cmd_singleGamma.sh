date
mkdir -p log
root -q -l -b eregtraining_13TeV_Pi0_singleGamma.C+(false,false,false) > log/reg_log_gamma1_EE.log
root -q -l -b eregtraining_13TeV_Pi0_singleGamma.C+(true,false,false) > log/reg_log_gamma1_EB.log 
root -q -l -b eregtraining_13TeV_Pi0_singleGamma.C+(false,false,true) > log/reg_log_gamma2_EE.log
root -q -l -b eregtraining_13TeV_Pi0_singleGamma.C+(true,false,true) > log/reg_log_gamma2_EB.log
root -q -l -b eregtesting_13TeV_Pi0_singleGamma.C+(false,false,false) > log/test_log_gamma1_EE.log
root -q -l -b eregtesting_13TeV_Pi0_singleGamma.C+(true,false,false) > log/test_log_gamma1_EB.log
root -q -l -b eregtesting_13TeV_Pi0_singleGamma.C+(false,false,true) > log/test_log_gamma2_EE.log
root -q -l -b eregtesting_13TeV_Pi0_singleGamma.C+(true,false,true) > log/test_log_gamma2_EB.log
date
