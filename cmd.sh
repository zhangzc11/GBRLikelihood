date
root -q -l -b eregtraining_13TeV_Pi0.C+(false,false,false) > reg_log_gamma1_EE.log
root -q -l -b eregtraining_13TeV_Pi0.C+(true,false,false) > reg_log_gamma1_EB.log 
root -q -l -b eregtraining_13TeV_Pi0.C+(false,false,true) > reg_log_gamma2_EE.log
root -q -l -b eregtraining_13TeV_Pi0.C+(true,false,true) > reg_log_gamma2_EB.log
root -q -l -b eregtesting_13TeV_Pi0.C+(false,false,false) > test_log_gamma1_EE.log
root -q -l -b eregtesting_13TeV_Pi0.C+(true,false,false) > test_log_gamma1_EB.log
root -q -l -b eregtesting_13TeV_Pi0.C+(false,false,true) > test_log_gamma2_EE.log
root -q -l -b eregtesting_13TeV_Pi0.C+(true,false,true) > test_log_gamma2_EB.log
date
