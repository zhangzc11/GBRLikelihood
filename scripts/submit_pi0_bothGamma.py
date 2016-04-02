#this is the macro to submit the regression for pi0->gammagamma
#Zhicai Zhang, zzhang2@caltech.edu
#Mar. 30, 2016
import subprocess, time, sys, os, shlex

queue = "cmscaf1nd"

if __name__ == "__main__":
	EEorEB = sys.argv[1]
	
	dobarrel = "false"
	if EEorEB == "EB":
		dobarrel = "true"
	pwd = os.getcwd()
	work_directory = pwd.replace("scripts","")
	#os.system("rm -rf "+pwd+"/submit")
	os.system("mkdir -p "+pwd+"/submit")

	env_script_n = pwd + "/submit/"+"_"+EEorEB+".sh"
	env_script_f = open(env_script_n, 'w')
	env_script_f.write("#!/bin/bash\n")
 	env_script_f.write("cd " + work_directory + "\n")
       	env_script_f.write("ulimit -c 0\n")
       	env_script_f.write("eval `scramv1 runtime -sh`\n")
       	#env_script_f.write("cmsenv \n")
       	env_script_f.write("date \n")
       	env_script_f.write("root -q -l -b eregtraining_13TeV_Pi0_bothGamma.C+\\("+dobarrel+",false\\) \n")
       	env_script_f.write("root -q -l -b eregtesting_13TeV_Pi0_bothGamma.C+\\("+dobarrel+",false\\) \n")
       	env_script_f.write("date \n")

	changePermission = subprocess.Popen(['chmod 777 ' + env_script_n], stdout=subprocess.PIPE, shell=True);
	debugout = changePermission.communicate()
	submit_s = 'bsub -q '+queue+' -o ' + pwd + "/submit/"+"_"+EEorEB+".log" + ' "source ' + env_script_n + '"'
	print "[submit_pi0_regression]  '-- " + submit_s
	submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
	output = (submitJobs.communicate()[0]).splitlines()
	print "[submit_pi0_regression]  '-- " + output[0]


