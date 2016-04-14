date
mkdir -p log
root -q -l -b eregtesting_13TeV_Pi0.C+(false,false,0) > log/test_log_EE_both.log
root -q -l -b eregtesting_13TeV_Pi0.C+(true,false,0) > log/test_log_EB_both.log
root -q -l -b eregtesting_13TeV_Pi0.C+(false,false,1) > log/test_log_EE_both.log
root -q -l -b eregtesting_13TeV_Pi0.C+(true,false,1) > log/test_log_EB_both.log
root -q -l -b eregtesting_13TeV_Pi0.C+(false,false,2) > log/test_log_EE_both.log
root -q -l -b eregtesting_13TeV_Pi0.C+(true,false,2) > log/test_log_EB_both.log
date
